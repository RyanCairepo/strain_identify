#include <cmath>
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include "omp.h"
#include <ctime>
#include <chrono>
using namespace std;

// Usage:
// ./smallRNA_propor -k 5 -f /home/xuanzhan/Data/miRNA/simu/data/m_d/m_d1/rawread.txt
// awk '{print $2}' rawread.txt | sort | uniq -c |sort -nk1 -r > freq_uniq.txt

const int MAX_CHAR_NUM = 1<<20;

int freq[MAX_CHAR_NUM] = {0};

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule
const char* inbase = "ATCG";

int K_value;
unsigned int t_FN = 4;
std::string f_FN;
std::ostream *v1logger;

//store info from fastq file,includes ID,raw read, phreh quality
vector<string> F_id;
vector<string> F_read;
vector<string> F_quality;
vector<int> F_flag;
vector<string> F_id2;
//output files
ofstream outfile("correct_read.fastq");
ofstream outlistfile("ID_read_quality_cor.txt");
ofstream correctedfile("changed_list.txt");//[ID  Changed_read] as input file for check.c

hash<string> str_hash;
int size_str_hash = MAX_CHAR_NUM; //2^30

struct ReadNode
{
    string read;
    int freq;
    vector<int> group;
	string quality;
};

struct CandNode
{
    string read;
    int freq;
    string quality;
};

//store info unique read frequency file for hash check
vector<ReadNode> *readfreq; //2^30

//store candidates seq
vector<CandNode> candidates;
//store CandNode for Top rank
vector<CandNode> candfreq;

//show help
inline void displayHelp(const char* prog) {
	printf("smallRNAEC v1.1, by XuanZhang, March 2021.\n");
	printf("Usage: %s -k <threshold_value> -f <fastq_fileName>[option parameters]\n", prog);
	printf("\t options:\n \t\t -t <threads>\n\n");
	// printf("-----------\n");
	printf("\t\t -f is the raw read info file from fastq\n");
	printf("\t\t -k is the threshold value, for erronous reads selected\n");
	printf("\t\t -t is the number of threads\n");

	printf("Example:\n\t\t");
	printf("./smallRNA_mix_fq -k 5 -f ID_read_quality.txt\n\n");
}

//show current params
inline void displayParams() {
	printf("threshold_value is k = %d\n", K_value);
	printf("the number of threads is t = %d\n", t_FN);
	printf("raw read info file is: %s\n", f_FN.c_str());
}

//get params
inline void getPars(int argc, char* argv[]) {
	v1logger = &std::cout;
	//bool is1 = false, is2 = false, is3 = false;
	//bool iskmer = false; //four
	int oc;
	while ((oc = getopt(argc, argv, "k:t:f:hf")) >= 0) {
		switch (oc) {
			case 'k':
				K_value = atoi(optarg);
				//iskmer = true;
				break;
			case 't':
				t_FN = atoi(optarg);
				break;
			case 'f':
				f_FN = optarg;
				//is3 = true;
				break;
			case 'h':
				displayHelp(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'smallRNA_mix_fq -h'\n";
				exit(1);
				break;
		}
	}

	// std::ifstream f;

	// f.open(f_FN);
	// if (f.fail()) {
	// 	fprintf(stderr, "raw read info file '%s' does not exist.\n", f_FN.c_str());
	// 	exit(1);
	// }
	// f.close();
}

//string hash: transfer mer/read to integer
int str_hash_index(string s_seq)
{
    //upperCase
    transform(s_seq.begin(), s_seq.end(), s_seq.begin(), ::toupper);
    
    size_t b = str_hash(s_seq);

    int index = b % size_str_hash;

    return index;
}

//read fasta :ID_read_quality.txt(preprocess from row_fastq file)
void readrow_fastaFile(const char *refFile, vector<ReadNode> *readfreq) { // processing reference file

	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}

	char *id, *read;
	id = new char[MAX_CHAR_NUM];
	read = new char[MAX_CHAR_NUM];
	ReadNode node;
	int hash_index;
	int unique_num=0;
	int exist_flag=0;
	int number =0;

	char *quality;
	quality = new char[MAX_CHAR_NUM];
	
	std::ifstream inputfile(refFile); 
	std::string line;
	std::string line_strip;
	int count = 1;
	while(std::getline(inputfile, line)){
			//std::cout << line;
			//line.erase(std::remove(line.begin(),line.end(),'\n'),line.end());
			//line_strip = line;
			std::cout << line_strip;
			if (count % 40000000==0){
				cout << "read " << count/4 << "reads " << endl;
			}
			switch (count % 4){

				case 1:
					F_id.push_back(line);
					F_flag.push_back(0);
					count++;
					continue;
					
				case 2:
					F_read.push_back(line);
					strcpy(read,line.c_str());
					count++;
					continue;
				case 3:
					F_id2.push_back(line);
					count++;
					continue;
				case 0:
					F_quality.push_back(line);
					strcpy(quality,line.c_str());
					count++;
					break;

			}

	//while (fscanf(fp, "%s %*s %*s %s %*s %*s %*s %s", id, read, quality) != EOF) {
	
			// cout << "now, the "<< number << " entry is preprocess"<<endl;
			/*
			printf("%s %s %s",id,read, quality);

			F_id.push_back(id);
			F_read.push_back(read);
			F_flag.push_back(0);
			F_quality.push_back(quality);
			//*/
			node.read = read;
        	node.freq = 1;
 			node.quality = quality;
 			
        	exist_flag=0;

        	hash_index = str_hash_index(read);

        	// cout << "current read: " <<read <<" its'hash: "<< hash_index <<endl;

        	// cout << "hash_size: "<< readfreq[hash_index].size()<<endl;

		    for (unsigned int i =0; i < readfreq[hash_index].size(); i++)
		    {
		    	if( readfreq[hash_index][i].read == read )
		    	{
		    		// cout<<"find read "<<readfreq[hash_index][i].read<< "  "<< readfreq[hash_index][i].freq <<" flag: "<<exist_flag<<endl;
		    		//already exist
		    		readfreq[hash_index][i].freq ++;

		    		// when larger than k_value not be erronous reads anymore, no need to store
		    		if ( readfreq[hash_index][i].freq <= K_value )
		    		{
		    			readfreq[hash_index][i].group.push_back(number);
		    			// cout<< number << " entry is preprocess "<<readfreq[hash_index][i].read<< "  "<< readfreq[hash_index][i].freq <<" number "<<number<<endl;
		    			// cout << "existing total ID :" <<readfreq[hash_index][i].group.size()<<endl;

		    			// for (unsigned int k =0; k < readfreq[hash_index][i].group.size(); k++)
		    			// {
		    			// 	cout <<readfreq[hash_index][i].group[k]<< endl;
		    			// }

		    		}

		    		exist_flag=1;
		    		// cout<<"freq++ "<<readfreq[hash_index][i].read<< "  "<< readfreq[hash_index][i].freq <<" flag: "<<exist_flag<<endl;
		    		break;
		    	}	
		    }

		    //first detected
		    if (exist_flag == 0)
		    {
		    	node.group.push_back(number);

		    	readfreq[hash_index].push_back(node);
		    	unique_num++;
		    	// cout<< number << " entry is preprocess"<< "First time"<<node.read<< "  "<< node.freq <<" number "<<number<<endl;
		    	node.group.clear();

		    }

		    number ++;

	}
	fclose(fp);

	cout << "total unique: "<<unique_num<<endl;
}

//setting parameters, refering to read copy at least 6 in simu_datasets
inline void initial() {
	//nothing here
	readfreq = new vector<ReadNode>[MAX_CHAR_NUM];
}

//decensing order
bool GreaterSort (CandNode a,CandNode b) 
{ 
	return (a.freq > b.freq); 
}

//find candidates seq
int find_cand(string err_read, int frequency, vector<CandNode>& candidates) {

	set<string> set_exist;
	
	candidates.clear();
	candfreq.clear();

	// // cout <<"in function: "<<err_read << " - "<< frequency<<endl;
	// for (int i=0; i <frequency; i++)
	// {
	// 	candidates.push_back(err_read);
	// }
	// // cout <<"end function "<<endl<<endl;

	CandNode node;
	string read_cor;
	int total_push = 0;

	//search subs and ins candidates
	for(unsigned int i=0; i<err_read.size(); i++)
	{	
		//search subs candidate
		for(int j=0; j<4; j++)
		{	
			read_cor = err_read;
			if(err_read[i]  != invert_code_rule[j])
			{
				read_cor[i] = invert_code_rule[j];

				int hash_index = str_hash_index(read_cor);
			    for (unsigned int k =0; k < readfreq[hash_index].size(); k++)
			    {
			    	if( readfreq[hash_index][k].read == read_cor )
			    	{	
			    		// cout <<"sub_position: " << i << "potential: "<< j << "candidate: "<<readfreq[hash_index][k].read<<"freq: "<< readfreq[hash_index][k].freq<<endl;
				 		if( readfreq[hash_index][k].freq >=  K_value  )
				 		{
				 			if ( set_exist.count(read_cor) == 0 )
				 			{
					 			node.freq = readfreq[hash_index][k].freq;
					 			node.read = read_cor;
					 			node.quality = readfreq[hash_index][k].quality;
					 			candfreq.push_back(node);
					 			set_exist.insert(read_cor);
					 			total_push++;
				 			}

				 		}
				 	}	    			
			    }
			}
		}
		//search ins candidate
		read_cor = err_read;
		read_cor.erase(i,1);
		int hash_index = str_hash_index(read_cor);
		for (unsigned int k =0; k < readfreq[hash_index].size(); k++)
		{
	    	if( readfreq[hash_index][k].read == read_cor )
	    	{
	    		// cout <<"ins_position: " << i << "candidate: "<<readfreq[hash_index][k].read<<"freq: "<< readfreq[hash_index][k].freq<<endl;
		 		if( readfreq[hash_index][k].freq >= K_value )
		 		{
		 			if ( set_exist.count(read_cor) == 0 )
		 			{
			 			node.freq = readfreq[hash_index][k].freq;
			 			node.read = read_cor;
			 			node.quality = readfreq[hash_index][k].quality;
			 			candfreq.push_back(node);
			 			set_exist.insert(read_cor);
			 			total_push++;
		 			}
		 		}
		 	}  			
		}
	}
	//search delet candidate
	for(unsigned int i=0; i<=err_read.size(); i++)
	{	
		for(int j=0; j<4; j++)
		{
			read_cor = err_read;
		    read_cor.insert(i,&inbase[j],1);
	
			int hash_index = str_hash_index(read_cor);
			for (unsigned int k =0; k < readfreq[hash_index].size(); k++)
			{
			   	if( readfreq[hash_index][k].read == read_cor )
		    	{
		    		// cout <<"del_position: " << i << "potential: "<< j << "candidate: "<<readfreq[hash_index][k].read<<"freq: "<< readfreq[hash_index][k].freq<<endl;
			 		if( readfreq[hash_index][k].freq >=  K_value )
			 		{
			 			if ( set_exist.count(read_cor) == 0 )
			 			{
				 			node.freq = readfreq[hash_index][k].freq;
				 			node.read = read_cor;
				 			node.quality = readfreq[hash_index][k].quality;
				 			candfreq.push_back(node);
				 			set_exist.insert(read_cor);
				 			total_push++;
			 			}
			 		}
			 	}   			
			}
		}
	}

	if ( candfreq.size() > 0 )
	{
		// cout<<"total candidates: "<<candfreq.size() << "  total err_reads: "<<frequency<<endl;
		sort(candfreq.begin(),candfreq.end(),GreaterSort);

		if ( (frequency == 1) || (candfreq.size() == 1) )
		{	
			candidates.push_back(candfreq[0]);
			return 1;
		}
		else
		{			
			// cout<<"total candidates: "<<candfreq.size() << "  total err_reads: "<<frequency<<endl;
			//calculate proportion candidates and do proportion correction
			int sum_freq = 0;
			// int tmp_count = frequency;
		    for (unsigned int k =0; k < candfreq.size(); k++)
			{
				sum_freq = sum_freq + candfreq[k].freq;
				// cout<<candfreq[k].read<<" "<<candfreq[k].freq<<endl;
			}
       		// cout<<"sum : "<<sum_freq<<endl;

			float propor; 

			for (unsigned int k =0; k < candfreq.size(); k++)
			{
				propor = float(candfreq[k].freq) / sum_freq * frequency;
				// cout << propor << "propor float to int "<<round(double(propor))<<" "<<candfreq[k].read<<endl;

				for (int j = 0; j < round(double(propor)); j++)
				{
					candidates.push_back(candfreq[k]);
				}

				// tmp_count = tmp_count - round(double(propor));

				if ( int(candidates.size()) == frequency )
				{
					break;
				}
				else
				{
					if (round(float(candfreq[k+1].freq) / sum_freq * frequency) == 0)
					{
						candidates.push_back(candfreq[k]);
					}
				}
			}
			// //for print
			// for (unsigned int k =0; k < candidates.size(); k++)
			// {
			// 	cout <<"total candidates: "<<k<<" : "<<candidates[k]<<endl;
			// }
			return 2;
		}

		// cout << "candidates found;  function finish" <<endl<<endl;

	}
	else
	{
		//no candidates found
		// cout << "no candidates found;  function finish; read:" <<err_read<<endl;
		return 0;
	}
}




int main(int argc, char *argv[]) {
	
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	
	getPars(argc, argv);

	displayParams();
	initial();

	readrow_fastaFile(f_FN.c_str(), readfreq);
	//write correct reads
	if(!outfile){
    cout << "Unable to open output file";
        exit(1); // terminate with error
    }

    //write correct reads list
	if(!outlistfile){
    cout << "Unable to open correct list file";
        exit(1); // terminate with error

	}

	//write only changed reads
	if(!correctedfile){
    cout << "Unable to open only changed list file";
        exit(1); // terminate with error

	}

	auto start = std::chrono::system_clock::now();

	cout << "K_value:"<< K_value<<endl;

	cout << "Total reads: "<<F_read.size()<<endl;

	int cand_found;
	int no_found = 0;
	int cor_total = 0;
	int por_total = 0;
	//start correction
	for (unsigned int i = 0; i < F_read.size(); i++)
	{
		if (i % 600000 == 0){
			cout << i << "/" << F_read.size() << "reads finished, " << (float)(i-1)/F_read.size() <<endl;
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> runtime = end-start;
			cout << "running time: " << runtime.count() << endl;
		}
		int hash_index = str_hash_index(F_read.at(i));

		for (unsigned int j =0; j < readfreq[hash_index].size(); j++)
		{
			// cout<<"current read: "<<F_read.at(i) << " it's hash: "<< hash_index<<endl;
		   	if( readfreq[hash_index][j].read == F_read.at(i) )
	    	{
	    		// tmp_found++;
	    		// cout<< "found: " << readfreq[hash_index][j].read <<"it's freq: "<< readfreq[hash_index][j].freq<<endl;
		 		if( ( readfreq[hash_index][j].freq < K_value ) && ( readfreq[hash_index][j].freq > 0 ) )
		 		{
		 			//erronous read found
		 			// tmp_low ++;
		 			// cout << "erronous read found in subs corretion！！！ read—ID： " << F_id.at(i)<<" seq: "<<F_read.at(i) << " its'freq: "<< readfreq[hash_index][j].freq << endl;

		 			//check subs candidate
		 			cand_found = find_cand(F_read.at(i),readfreq[hash_index][j].freq, candidates);

		 			if (cand_found == 0)
		 			{
		 				no_found ++;
		 			}
		 			else if (cand_found == 1) 
		 			{
		 				for (int k =0; k < readfreq[hash_index][j].freq; k++ )
		 				{
		 					// cout <<k << " error reads index "<<readfreq[hash_index][j].group[k] << endl;
		 					// cout << F_id.at(readfreq[hash_index][j].group[k]) <<"  "<< F_read.at(readfreq[hash_index][j].group[k]) << endl;
		 					F_read.at(readfreq[hash_index][j].group[k]) = candidates[0].read; 
		 					F_flag.at(readfreq[hash_index][j].group[k]) = 1;
		 					F_quality.at(readfreq[hash_index][j].group[k]) = candidates[0].quality;
		 					// cout << F_id.at(readfreq[hash_index][j].group[k]) <<"  "<< F_read.at(readfreq[hash_index][j].group[k]) << endl;

		 				}
		 				// cout <<endl;
		 				cor_total = cor_total + readfreq[hash_index][j].freq;
		 			}
		 			else if (cand_found == 2)
		 			{	 				
		 				if ( int(candidates.size()) == readfreq[hash_index][j].freq )
		 				{
		 					for (int k =0; k < readfreq[hash_index][j].freq; k++ )
		 					{
		 						// cout << F_id.at(readfreq[hash_index][j].group[k]) <<"  "<< F_read.at(readfreq[hash_index][j].group[k]) << endl;
		 						F_read.at(readfreq[hash_index][j].group[k]) = candidates[k].read; 
		 						F_flag.at(readfreq[hash_index][j].group[k]) = 1;
		 						F_quality.at(readfreq[hash_index][j].group[k]) = candidates[k].quality;	
		 						// cout << F_id.at(readfreq[hash_index][j].group[k]) <<"  "<< F_read.at(readfreq[hash_index][j].group[k]) << " after correction" << endl;	 						
		 					}
		 					por_total = por_total + readfreq[hash_index][j].freq;
		 				}
		 				else
		 				{
		 					// cout << "proportion correction size doesn't match";
		 				}

		 			}


		 		}
		 	}   			
		}		
	}
	cout << "no_found: "<< no_found <<endl;
	cout << "total correction: "<< cor_total <<endl;
	cout << "total porpor correction: "<< por_total <<endl;



    //print output into files
	for (unsigned int i = 0; i < F_read.size(); ++i)
	{
		// //for fastq file
		// outfile<<F_id.at(i)<<endl<<F_read.at(i)<<endl<<"+"<<endl<<F_quality.at(i)<<endl; 
		// outlistfile<<F_id.at(i)<<" "<<F_read.at(i)<<" "<<F_quality.at(i)<<endl;

		//for fasta file
		outfile<<F_id.at(i)<<endl<<F_read.at(i)<<endl<<F_id2.at(i)<<endl<<F_quality.at(i)<<endl; 
		outlistfile<<F_id.at(i)<<" "<<F_read.at(i)<<endl;

		if(F_flag.at(i) == 1)
		{
			correctedfile<<F_id.at(i)<<endl<<F_read.at(i)<<endl<<F_quality.at(i)<<endl;
		}
	}

	cout << "code finished!" << endl;


	outfile.close();
	outlistfile.close();
	correctedfile.close();
}





