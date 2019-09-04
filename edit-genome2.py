import sys
import csv

sequence_file = open(sys.argv[1])
start = 0
end = 0
linelength = 0
sequence = ""
#check file info
for line in sequence_file:
    if line[0] == ">":
        ll = line.split(" ")
        start = int(ll[1])
        end = int(ll[2])
    else:
        line = line.rstrip("\n")
        sequence = sequence + line

sequence = list(sequence)


edit_file = open(sys.argv[2])

newfile_name = sys.argv[3]

edit = csv.reader(edit_file)
newfile = open(newfile_name, "a")

number_per_seq = int(sys.argv[4])

if number_per_seq > 3:
    print("Error, only 3 possible")
    exit(-1)


##read in snps
snps = []

##snp
for line in edit:
    name = line[0]
    pos = line[2].split(":")
    chrom = int(pos[0])
    pos = int(pos[1])
    new = line[4]

    posInSeq = pos - start

    snp = [name,posInSeq,new]
    snps.append(snp)

if number_per_seq == 1:
    for item in snps:
        name = item[0]
        pos = item[1]
        new = item[2]

        newfile.write(">"+name+"\n")
        ll = sequence[:]
        ll[pos] = new

        newfile.write("".join(ll))
        newfile.write("\n")
        
elif number_per_seq == 2:
    for i in range(len(snps)):
        item = snps[i]
        name1 = item[0]
        pos1 = item[1]
        new1 = item[2]

        for j in range(i+1, len(snps)):
            item2 = snps[j]
            pos2 = item2[1]
            
            if pos2 == pos1:
                continue
            name2 = item[0]
            new2 = item[2]

            newfile.write(">"+name1+name2+"\n")
            ll = sequence[:]
            ll[pos1] = new1
            ll[pos2] = new2


            newfile.write("".join(ll))
            newfile.write("\n")

elif number_per_seq == 3:
    for i in range(len(snps)):
        item = snps[i]
        name1 = item[0]
        pos1 = item[1]
        new1 = item[2]

        for j in range(i+1, len(snps)):
            item2 = snps[j]
            pos2 = item2[1]
            name2 = item2[0]
            new2 = item2[2]

            for k in range(j+1, len(snps)):
                item3 = snps[k]
                pos3 = item3[1]
                name3 = item3[0]
                new3 = item3[2]
              #  print(pos1, pos2, pos3)


                if pos2 == pos3:
                    new3 = new2
                if pos1 == pos3:
                    new3 = new1
                if pos2 == pos1:
                    new2 = new1

                newfile.write(">"+name1+name2+name3+"\n")
                ll = sequence[:]
                ll[pos1] = new1
                ll[pos2] = new2
                ll[pos3] = new3
    

                newfile.write("".join(ll))
                newfile.write("\n")
newfile.close()
sequence_file.close()
edit_file.close()