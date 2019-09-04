import sys

down = int(sys.argv[3])
up = int(sys.argv[4])
if(down > up):
    down, up = up, down

f = open(sys.argv[1])
seqs = 0
#print ("Processing", f)
for line in f:
    if line[0] != "@":
        sl = line.split()
        ch = sl[3]
        p = int(sl[4])
        if ch == sys.argv[2]:
            if (p >= down) & (p <= up):
                seqs += 1
                if "bow" in sys.argv[1]: 
                    print(sl[3], sl[4], sl[2], sl[5], sl[-1])
                else:
                    print(sl[2], sl[3], sl[9])

#print ("Finished - found:",seqs,"sequences")
