import sys
id1 = set({})
with open(sys.argv[1],"r") as if1:
        for line in if1:
                sline = line.strip().split(" ")
                for id in sline:
                        id1.add(id)
id2 = set({})
with open(sys.argv[2],"r") as if2:
        for line in if2:
                sline = line.strip().split(" ")
                for id in sline:
                        #if id not in id1:
                        id2.add(id)

print(len(id1),len(id2))
if len(id1) < len(id2):
	print(len(id2-id1))
	print(id2-id1)
else:
	print(len(id1-id2))
	print(id1-id2)
