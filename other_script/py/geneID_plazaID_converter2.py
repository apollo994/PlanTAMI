#questo convertitore rispetto alla versione 1 restituisce sia il plazaID che il
#relativo speciesID che viene letto per la conversione 


import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='this tool convert from geneID to plazaID using a conversion table')
    parser.add_argument('--input', metavar='input', help='gene list to convert')
    parser.add_argument('--conv_tab', metavar='conv_tab', help='table with geneID and relative plaza ID')
    #parser.add_argument('--output', metavar='OUTPUTFILE', default="/dev/stdout", help='The output file.')
    #parser.add_argument('--fast', action='store_true', help="Fast parsing.")
    #parser.add_argument('--max', metavar='MAX', default="1024", type=int, help='Maximum.')
    #parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
    args = parser.parse_args()

    ###########################################################################

    diz_ID={}
    geneID_conv=[]
    geneID_not_conv=[]

    ###########################################################################

    input_name=args.input.split("/")
    input_name=input_name[len(input_name)-1]
    output_name="plazaID_tairID_"+input_name

    with open (args.conv_tab) as conv_tab:
        for line in conv_tab:
            line=line.split("\t")
            geneID=line[0]
            plazaID=line[1].strip("\n")
            diz_ID[geneID]=plazaID

    with open (args.input) as gene_to_convert:
        for line in gene_to_convert:
            geneID=line.strip("\n")
            if geneID in diz_ID:
                geneID_conv.append([diz_ID[geneID],geneID])
            else:
                geneID_not_conv.append(geneID)


    with open (output_name, "w") as out:
        for plazaID in geneID_conv:
            out.write(plazaID[0]+"\t"+plazaID[1]+"\n")


    sys.stderr.write("GeneID converted = "+str(len(geneID_conv))+"\n")
    sys.stderr.write("GeneID not converted = "+str(len(geneID_not_conv))+"\n")
    for item in geneID_not_conv:
        sys.stderr.write(item+"\n")


if __name__ == "__main__":
    main()
