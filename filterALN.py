import sys

if len(sys.argv) == 1: 
    filename = input("Please enter input filename: ")
    filename2 = input("Please enter output filename: ")
else:
    filename = sys.argv[1]
    filename2 = sys.argv[2]
if len(sys.argv) == 4:
    n = int(sys.argv[3])
else:
    n = int(input("Please enter maximum duplicates: "))

old = 0
new = 0
with open(filename2, "w") as g:
    with open(filename, "r") as f:
        for line in f:
            lst = line.split()
            newlist = lst[1:]
            m = int(lst[0])
            old += m
            if m > n:
                for i in range(n):
                    new += 1
                    g.write("\t".join(newlist))
                    g.write("\n")
            else:
                for i in range(m):
                    new += 1
                    g.write("\t".join(newlist))
                    g.write("\n")

print("Old number: ", old)
print("New number: ", new)
print("Percentage: ", (new/old) * 100, "%")
