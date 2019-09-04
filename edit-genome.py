import sys
import csv
import os

#takes csv file with
#ID, pos, posFromEnd, new, kindOfSNP

sequence_file = open(sys.argv[1])
sequence = ""
start = 0
for line in sequence_file:
    if line[0] == ">":
        ll = line.split(" ")
        start = int(ll[1])
    else:
        line = line.rstrip("\n")
        sequence = sequence + line

sequence = list(sequence)


edit_file = open(sys.argv[2])
edit = csv.reader(edit_file)



for line in edit:
    name = line[0]
    print(line)
    edit_type = line[4].lower()
    pos = line[1].split(":")[1]
    newfilename = name
    count = 0
    while os.path.isfile(newfilename+".fa"):
        count += 1
        newfilename = name + "-" + str(count)     
    newfile = open(newfilename+".fa", "w")
    newfile.write(">"+name+"\n")

    ll = sequence[:]

    if edit_type == "snp":
        pos = int(pos)-start
        new = line[3]
        if len(new) > 1:
            print(line)
            exit(0)

        #replace
        ll[pos] = new

    elif edit_type == "deletion":
        pos = line[2].split(":")[1]
        if "-" in pos:
            pos = pos.split("-")
            startpos = int(pos[0])-start
            endpos = int(pos[1])-start
        else:
            startpos = int(pos)-start
            endpos = startpos

        #delete
        del ll[startpos:endpos+1]
        
    elif edit_type == "insertion":
        pos = int(pos)-start
        new = line[3]

        #insert
        ll.insert(pos+1, new)

    newfile.write("".join(ll))
    newfile.close()


sequence_file.close()
edit_file.close()