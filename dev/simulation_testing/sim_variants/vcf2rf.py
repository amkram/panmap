import sys

with open(sys.argv[1]) as vcf:
    vcf = vcf.readlines()
    vcf = vcf[1:]

    rf = open(sys.argv[2],"w")

    for line in vcf:
        line = line.split()

        new_line = line[1] + "\t" + line[3] + "\t" + line[4] + "\n"
        rf.write(new_line)