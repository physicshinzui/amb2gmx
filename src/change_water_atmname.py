import sys 

def replace_water_atomname(top):
    fout = open("modified.top", "w")
    with open(top, "r") as fin:
        for line in fin:
            if "moleculetype" in line.split():
                fout.write(line)
                for line2 in fin:
                    if "SOL" in line2.split():
                        line2 = line2.replace(" O ", "OW ").replace(" H1", "HW1").replace(" H2", "HW2")
                        print("Replaced line: ", line2.strip())
                        fout.write(line2)
                    else:
                        fout.write(line2)
            else:  
                fout.write(line)

replace_water_atomname(sys.argv[1])
