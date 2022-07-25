import sys
IDs = set({})
with open(sys.argv[1],"r") as f:
    for line in f:
        IDs.add(line.split(" ")[0])
print(len(IDs))
with open(sys.argv[2],"w+") as wf:
    for id in list(IDs):
        wf.write(id+" ")
