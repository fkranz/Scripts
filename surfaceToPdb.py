import sys

if (len(sys.argv) != 6) & (len(sys.argv) != 3):
    print ("Usage: python surfaceToPdb.py inputFile outputPdbFile")
    print ("or:    python surfaceToPdb.py inputFile offsetX offsetY offsetZ outputPdbFile")
    print ("Sometimes pymol has troubles with the coordinates, because they are out of")
    print ("some range, then use offsets to correct")
    exit(-1)

f = open(sys.argv[1])
if len(sys.argv) == 3:
    ox = 0.0
    oy = 0.0
    oz = 0.0
    g = open(sys.argv[2], "w")
else:
    ox = float(sys.argv[2])
    oy = float(sys.argv[3])
    oz = float(sys.argv[4])
    g = open(sys.argv[5], "w")
c = 0
for line in f:
    pos = line[:-2]
    pos = pos.replace("(","")
    pos = pos.replace(")","")
    sl = pos.split(",")
    s = '{:<6}{:>5} {:<4}{:<1}{:>3} {:<1}{:>4}{:<1}   {:>8,.2f}{:>8,.3f}{:>8,.3f}{:>6,.3f}{:>6,.2f}{:>10}{:>2}{:>2}'.format("ATOM",c+1, "PS", "",
            "PSD", "P", 1, "",float(sl[0])+ox, float(sl[1])+oy,float(sl[2])+oz, 1.00, 0.00, "", "P", 0)
    c = (c+1) % 99999
    g.write(s)
    g.write("\n")

f.close()
g.close()
