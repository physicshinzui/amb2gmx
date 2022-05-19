import sys

def top2itp(top, out):
    with open(top, 'r') as fin:
        with open(out, 'w') as out: 
            for line in fin:

                if "defaults" in line:
                    for line in fin: 
                        if "[" in line:
                            out.write(line)
                            break

                elif "system" in line:
                    break

                else:
                    out.write(line)

top2itp(sys.argv[1], sys.argv[2])
