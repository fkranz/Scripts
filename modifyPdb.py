import sys

def createThreeLetterTable():
    return {
        'trp':   -2.0285,
        'leu':   -1.9873,
        'met':   -1.6858,
        'ile':   -1.5691,
        'val':   -1.0727,
        'cys':   -0.9818,
        'phe':   -0.5965,
        'tyr':   -0.4582,
        'ala':   -0.3651,
        'pro':   -0.3085,
        'lys':   -0.0800,
        'arg':    0.3048,
        'thr':    0.6953,
        'gly':    0.9377,
        'ser':    1.0977,
        'his':    1.4126,
        'asn':    1.6949,
        'gln':    1.7191,
        'glu':    5.0083,
        'asp':    5.3240
    }

def modifypdb(infile, outfile, windowsize):
    table = createThreeLetterTable()
    lines = []
    bFactor = 0
    if windowsize == 1:
        for line in infile:
            if (line[0:4] != "ATOM"):
                lines.append(line)
                continue

            residuum = line[17:20]
            lowerResiduum = residuum.lower()
            if lowerResiduum not in table:
                outfile.write("Error: unknown residuum " + residuum)
                return

            bFactor = float(table[lowerResiduum])
            bFactorS = "{:>5}".format(str(round(bFactor, 2)))
            lines.append(line[0:61] + bFactorS + line[66:81])
    else:
        numberOfRes = 0
        linesBuffer = []
        lastResSeq = -1

        for line in infile:
            if (line[0:3] == "TER"):
                if numberOfRes > 1:
                    bFactor /= (numberOfRes - 1)
                    for l in linesBuffer:
                        bFactorS = "{:>5}".format(str(round(bFactor,2)))
                        lines.append(l[0:61] + bFactorS + l[66:81])
                    bFactor = 0
                    linesBuffer = []
                    numberOfRes = 0
            elif (line[0:4] != "ATOM"):
                lines.append(line)
                continue

            linesBuffer.append(line)

            residuum = line[17:20]
            resSeq = line[22:26]

            if (resSeq != lastResSeq):
                lowerResiduum = residuum.lower()
                if lowerResiduum not in table:
                    print("Error: unknown residuum " + residuum)
                    return

                bFactor += float(table[lowerResiduum])
                lastResSeq = resSeq
                numberOfRes += 1

            if (numberOfRes == windowsize+1):
                bFactor -= float(table[lowerResiduum])
                bFactor /= (numberOfRes - 1)
                numberOfRes = 1
                for l in linesBuffer:
                    bFactorS = "{:>5}".format(str(round(bFactor,2)))
                    lines.append(l[0:61] + bFactorS + l[66:81])

                tmp = linesBuffer[-1]
                linesBuffer = [tmp]
                bFactor = float(table[lowerResiduum])


    for line in lines:
        g.write(line)



if len(sys.argv) < 2:
    print("Usage: python3 program.py {in.pdb|in.txt} [windowsize]")
    print("\tInput can be a single pdb-file or a txt-file containing")
    print("\tmultiple pdb-files (one file path per line)")
    print("\tWindowsize is optional and by default 1.")



infile = sys.argv[1]
windowsize = 1
if len(sys.argv) > 2:
    windowsize = int(sys.argv[2])


if ".txt" in infile:
    inf = open(infile)
    for line in inf:
        f = open(line)
        of = line.replace(".pdb", "_modified_ws"+str(windowsize)+".pdb")
        g = open(of, 'w')
        modifypdb(f, g, windowsize)
else:
    f = open(infile)
    of = infile.replace(".pdb", "_modified_ws"+str(windowsize)+".pdb")
    g = open(of, 'w')

    modifypdb(f, g, windowsize)