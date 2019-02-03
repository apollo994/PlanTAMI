import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Tool for conversion from geneID into otrho_family_ID (ath,sly,osa,hvu)')
    parser.add_argument('--input', metavar='input', help='The input file')
    parser.add_argument('--conv_tab', metavar='conv_tab', default="/Users/Apollo199/Tirocinio/HOMO_GENES/formatted_genefamily/gene_family_unique_gene.txt", help='Table for conversion')
    parser.add_argument('--sp', metavar='sp', help='species to convert (ath,sly,osa,hvu)')
    parser.add_argument('--sep', metavar='sep', default="tab",help='separator of the columns (tab, comma), def=tab')
    parser.add_argument('--out', metavar='out', default="coveted_ID_tab.txt",help='output name')
    args = parser.parse_args()

    unique_gene=args.conv_tab
    tab_to_convert=args.input
    separator=args.sep

    ################################################################################
    diz_ID={}
    line_counter=0
    ################################################################################

    if args.sp=="ath":
        index_sp=1
    if args.sp=="sly":
        index_sp=2
    if args.sp=="osa":
        index_sp=3
    if args.sp=="hvu":
        index_sp=4


    with open (unique_gene) as f:
        for line in f:
            row=line.strip().split("\t")
            diz_ID[row[index_sp]]=row[0]

    ###############################################################################

    with open(tab_to_convert) as l:
        if separator=="comma":
            head=l.readline().strip("\n").split(",")
            head="\t".join(head)
            with open (args.out,"w") as out:
                out.write(head+"\n")
        if separator=="tab":
            head=l.readline().strip("\n")
            with open (args.out,"w") as out:
                out.write(head+"\n")

        for line in l:
            if separator=="tab":
                row=line.strip().split("\t")
                if row[0] in diz_ID:
                    row[0]=diz_ID[row[0]]
                    line_counter=line_counter+1
                    with open (args.out,"a") as out:
                        out.write(("\t".join(str(i) for i in row))+"\n")
            if separator=="comma":
                row=line.strip().split(",")
                if row[0] in diz_ID1:
                    row[0]=diz_ID1[row[0]]
                    line_counter=line_counter+1
                    with open (args.out,"a") as out:
                        out.write(("\t".join(str(i) for i in row))+"\n")

        to_print="line converted = "+ str(line_counter)+"\n"
        sys.stderr.write(to_print)


if __name__ == "__main__":
    main()
