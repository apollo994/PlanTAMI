import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Tool for up and down extraction from edgeR results')
    parser.add_argument('--input', metavar='input', help='The input file')
    args = parser.parse_args()

################################################################################

    file_name=args.input
    file_name=file_name.split("\t")
    file_name=file_name[len(file_name)-1]

    down_genes=[]
    down_genes_ID=[]
    up_genes=[]
    up_genes_ID=[]

    with open (args.input) as my_file:
        salta=my_file.readline()
        for line in my_file:
            line=line.split("\t")
            if float(line[1])<0:
                down_genes.append(line)
                down_genes_ID.append(line[0]+"\n")
            if float(line[1])>0:
                up_genes.append(line)
                up_genes_ID.append(line[0]+"\n")


    sys.stdout.write("\n")
    sys.stdout.write("My sample = " + str(file_name) +"\n")
    sys.stdout.write("Significant DE genes = " + str(len(down_genes)+len(up_genes))+"\n")
    sys.stdout.write("Down regulated genes = " + str(len(down_genes))+"\n")
    sys.stdout.write("Up regulated genes = " + str(len(up_genes))+"\n")
    sys.stdout.write("\n")

    with open ("DW_"+file_name , "w") as out_1:
        out_1.write(salta)
        for line in down_genes:
            line="\t".join(line)
            out_1.write(line)

    with open ("UP_"+file_name , "w") as out_2:
        out_2.write(salta)
        for line in up_genes:
            line="\t".join(line)
            out_2.write(line)

    with open ("DW_ID_"+file_name , "w") as out_3:
        for line in down_genes_ID:
            out_3.write(line)

    with open ("UP_ID_"+file_name , "w") as out_4:
        for line in up_genes_ID:
            out_4.write(line)

    with open ("stat_"+file_name, "w") as out_5:
        out_5.write("My sample = " + str(file_name) +"\n")
        out_5.write("Significant DE genes = " + str(len(down_genes)+len(up_genes))+"\n")
        out_5.write("Down regulated genes = " + str(len(down_genes))+"\n")
        out_5.write("Up regulated genes = " + str(len(up_genes))+"\n")


if __name__ == "__main__":
    main()
