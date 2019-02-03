import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Get a conversion table with geneID and the relative plazaID')
    parser.add_argument('--input', metavar='input', default="/Users/Apollo199/Tirocinio/HOMO_GENES/genefamily_data.orth.txt", help='table from plaza')
    parser.add_argument('--sp', metavar='sp', help='species to get the conversion table (ath,sly,hvu,osa)')
    #parser.add_argument('--output', metavar='OUTPUTFILE', default="/dev/stdout", help='The output file.')
    #parser.add_argument('--fast', action='store_true', help="Fast parsing.")
    #parser.add_argument('--max', metavar='MAX', default="1024", type=int, help='Maximum.')
    #parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
    args = parser.parse_args()

    ############################################################################

    diz_ID={}

    output_file=args.sp+"_geneID_to_plazaID.txt"
    ############################################################################

    with open (args.input) as plaza_tab:
        for line in plaza_tab:
            line=line.split(";")
            if line[1]==args.sp:
                geneID=line[2]
                geneID=geneID.strip("\r\n")
                geneID=geneID.split(".")
                geneID=geneID[0]
                plazaID=line[0]
                diz_ID[geneID]=plazaID

    with open (output_file, "w") as output:
        for key in diz_ID:
            geneID=key
            palzaID=diz_ID[key]
            row_to_print=geneID+"\t"+palzaID+"\n"
            output.write(row_to_print)




if __name__ == "__main__":
    main()
