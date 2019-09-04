from math import log, exp, log2
import sys, getopt, os

def print_help():
    print("\nUsage\n==========")
    print(arguments[0], "-d  directory -s sample -o output_path -f transcriptToGene-file [-t]")
    print(arguments[0], "-t  output transcripts (default: False)")
    print(arguments[0], "-h  prints this help")
    sys.exit()

def countToTpm(counts, effLen):
    tpms = {}
    rates = {}
    denom  = 0
    for item in counts.keys():
        t = counts[item] / effLen[item]
        rates[item] = t
        denom += t

    for item in rates.keys():
        tpms[item] = (rates[item]/denom * 1e6)

    return tpms


def countToFpkm(counts, effLen):
    N = sum(counts.values())
    fpkms = {}
    for item in counts.keys():
        if counts[item] == 0.0:
            fpkms[item] = 0.0
        else:
            t = exp(log(counts[item]) + log(1e9) - log(effLen[item]) - log(N))
            fpkms[item] = t

    return fpkms

def mylog(fpkms):
    lf = {}
    for item in fpkms.keys:
        if fpkms[item] == 0:
            lf[item] = float("Nan")
        else:
            lf[item] = log2(fpkms[item])
    return lf

def mylog2(f):
    if f == 0: 
        return float("Nan")
    else:
        return log2(f)

arguments = sys.argv
directory = ""
output = ""
sample = ""
control = ""
transcripts = False
transcriptToGene = "/home/snfrbert/Data/genomes/kallisto/hg38/transcriptToGene.txt"

try:
    opts, args = getopt.getopt(arguments[1:],"hd:s:o:t")
except getopt.GetoptError as errorInfo:
    print(errorInfo)
    print_help()
else:
    for opt, arg in opts:
        if (opt == "-h"):
            print_help()
        elif (opt == "-s"):
            sample = arg
        elif (opt == "-o"):
            output = arg
        elif (opt == "-d"):
            directory = arg
        elif(opt == "-t"):
            transcripts = True
        elif(opt == "-f"):
            transcriptToGene = arg



if (directory == "") or (output == "") or (sample == "") :
    print("Missing arguments!")
    print_help()


gene_trans = {}
mart = open(transcriptToGene)
for line in mart:
    if line[0] == "#": 
        continue
    ls = line.split("\t")
    g = ls[2].replace("\n", "")
    if g not in gene_trans:
        gene_trans[g] = [ls[1]]
    else:
        gene_trans[g].append(ls[1])


counts = {}
lens = {}

ilist = [f for f in os.listdir(directory) if sample in f]
for filename in ilist:
    path = directory + "/" + filename + "/abundance.tsv"
    f = open(path)
    for line in f:
        if line.startswith("target"): continue
        ls = line.split("\t")
        trans = ls[0]
        transk = trans.split(".")[0]
        if transk not in counts:
            counts[transk] = float(ls[3])
            lens[transk] = float(ls[2])
        else:
            counts[transk] += float(ls[3])
    f.close()


fpkm = countToFpkm(counts, lens)
tpm = countToTpm(counts, lens)

if transcripts:
    g = open(output+"-transcript", "w")
    print("transcript_id\tcount\tfpkm\ttpm", file=g)
    for item in counts.keys():
        print(item, counts[item], fpkm[item], tpm[item], sep="\t", file = g)
    g.close()
  # sys.exit(0)

res = {}
for item in gene_trans.keys():
    cs = 0
    ts = 0
    fs = 0
    for t in gene_trans[item]:
        if t in counts.keys():
            cs += counts[t]
            ts += tpm[t]
            fs += fpkm[t]
    res[item] = [cs, ts, fs]


g = open(output+"-gene", "w")
print("transcript_id\tcount\tfpkm\ttpm", file=g)
for item in res.keys():
    c, t, f = res[item]  
    print(item, c, f, t, sep="\t", file = g)
g.close()
sys.exit(0)
