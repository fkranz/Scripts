#! /usr/bin/env python3
import sys, getopt

chrom = {
    "GL000191.1":"1",
    "GL000192.1":"1",
    "GL000193.1":"4",
    "GL000194.1":"4",
    "GL000195.1":"7",
    "GL000196.1":"8",
    "GL000197.1":"8",
    "GL000198.1":"9",
    "GL000199.1":"9",
    "GL000200.1":"9",
    "GL000201.1":"9",
    "GL000202.1":"11",
    "GL000203.1":"17",
    "GL000204.1":"17",
    "GL000205.1":"17",
    "GL000206.1":"17",
    "GL000207.1":"18",
    "GL000208.1":"19",
    "GL000209.1":"19",
    "GL000210.1":"21"
}

def print_help():
    print("\nUsage\n==========")
    print(arguments[0], "-i input -a annotationfile [-o output -d directory]")
    print("Optional arguments:")
    print("\t-d directory:  output directory (default same as input)")
    print("\t-o ouptut:  output name (default inputname_anno)")
    print(arguments[0], "-h  prints this help")
    sys.exit()

arguments = sys.argv
directory = ""
output = ""
input = ""
anno = ""
try:
    opts, args = getopt.getopt(arguments[1:],"hi:o:d:a:")
except getopt.GetoptError as errorInfo:
    print(errorInfo)
    print_help()
else:
    for opt, arg in opts:
        if (opt == "-h"):
            print_help()
        elif (opt == "-i"):
            input = arg
        elif (opt == "-o"):
            output = arg
        elif (opt == "-d"):
            directory = arg
        elif (opt == "-a"):
            anno = arg

if (anno == "") or (input == ""):
    print("Error: input or annotationfile missing!")
    print_help()

if (output == ""):
    if(directory == ""):
        output = input.replace(".cod", "") + "_anno.cod"
    elif (directory == "."):
        output = input.split("/")[-1].replace(".cod", "") + "_anno.cod"
    else:
        output = directory + "/" + input.split("/")[-1].replace(".cod", "") + "_anno.cod"

annolist = {"1":[], "2":[], "3":[], "4":[], "5":[], "6":[], "7":[], "8":[], "9":[], "10":[], "11":[], "12":[], "13":[], "14":[], "15":[], "16":[], "17":[] ,"18":[], "19":[], "20":[], "21":[], "22":[], "X":[], "Y":[], "M":[]}
try:
    f = open(anno, "r")
except:
    print("Annoation file not found!")
    sys.exit()
else:
    for line in f:
        if line.startswith("#"): continue
        ls = line.split(",")
        chr = ls[0]
        chr = chr.upper()
        if chr in annolist:
            info = [int(ls[1]), ls[2], ls[3]]
            annolist[chr].append(info)
        elif chr == "MT":
            info = [int(ls[1]), ls[2], ls[3]]
            annolist["M"].append(info)
        elif chr in chrom:
            info = [int(ls[1]), ls[2], ls[3]]
            annolist[chrom[chr]].append(info)
        elif chr.startswith("HS"):
            info = [int(ls[1]), ls[2], ls[3]]
            chr = chr[5:7]
            chr = chr.replace("_", "")
            chr = chr.replace("L", "")
            annolist[chr].append(info)
    f.close()

try:
    f = open(input, "r")
except:
    print("Inputfile not found!")
    sys.exit()
else:
    g = open(output, "w")
    g.write("#seq_id\tchrom\tstart\t\tend\t\tstrand\tannotation\n")
    for line in f:
        if line.startswith("#"): continue
        ls = line.split()
        sid = ls[0]
        chr = ls[1]
        chr = chr.replace("chr", "")
        chr = chr.upper()
        start = int(ls[2])
        end = int(ls[3])
        strand = ls[4]
        mid = int(start + (end-start)/2)
        dist = abs(mid - annolist[chr][0][0])
        gene = ""
        for item in annolist[chr]:
            if True: #(strand == "1" and item[1] == "+") or (strand == "-1" and item[1] == "-"):
                ddist = abs(mid - item[0])
                if ddist > dist: 
                    break
                gene = item[2]
                dist = ddist
        g.write("{}\t{}\t{}\t{}\t{}\t{}".format(sid, ls[1], start, end, strand, gene))

    f.close()
    g.close()