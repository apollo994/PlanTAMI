#questo script prende come input un numero variabile di tabelle di risultati di deseq2 (sample.res.txt) e
#restituisce una tablla con solo i logFC e gli FDR delle tabelle in input
#e una tabella con solo i geni significativi in almeno un confronto

import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Tool for significant DE genes extraction from edgeR results')
    parser.add_argument('--input', metavar='input', nargs='+', help='The input file')
    parser.add_argument('--th', metavar='th', default="0.05", type=float, help='FDR threshold')
    parser.add_argument('--print_all', metavar='print_all', default="n", help='y to print table with all genes, def=n')
    args = parser.parse_args()

################################################################################

    logFC_FDR_allgenes=defaultdict(list)
    logFC_FDR_significant={}
    header=[]
    sample_number=len(args.input)

################################################################################

#funzione per estrarre solo logFC e FDR

    def extract_logFC_FDR(sample):
        with open (sample) as l:
            first_line=l.readline()
            for line in l:
                splitted_line=line.split("\t")
                geneID=splitted_line[0]
                logFC=splitted_line[1]
                FDR=splitted_line[5].strip("\n")
                logFC_FDR_allgenes[geneID].append(logFC)
                logFC_FDR_allgenes[geneID].append(FDR)

#funzione per estrarre solo i geni che passano la solgli in almeno un confronto

    def extract_significant():
        FDR_index=1
        counter=1
        while counter<=sample_number:
            for key in logFC_FDR_allgenes:
                if float(logFC_FDR_allgenes[key][FDR_index])<args.th:
                    logFC_FDR_significant[key]=logFC_FDR_allgenes[key]
            FDR_index=FDR_index+2
            counter=counter+1


#genero un header basandomi sui nomi dei campioni

    for pattern in args.input:
        sample=pattern.split("/")
        sample=sample[len(sample)-1]
        sample=sample.strip(".res.txt")
        header.append(sample+"_logFC")
        header.append(sample+"_FDR")

#genero il dizionario con tutti i geni e i relativi logFC e FDR

    for sample in args.input:
        extract_logFC_FDR(sample)

#genero il dizionario con solo i geni significativi

    extract_significant()

#scrivo la tabella con tutti i geni come primo file di output

    if args.print_all=="y": 
        with open ("logFC_FDR_allgenes.txt","w") as out:
            out.write("\t".join(header)+"\n")
            for key in logFC_FDR_allgenes:
                gene_row=key+"\t"+"\t".join(logFC_FDR_allgenes[key])+"\n"
                out.write(gene_row)

#scrivo la tabella con solo i geni significativi come secondo file di output

    with open ("logFC_FDR_%s_significant.txt" %str(args.th), "w") as out2:
        out2.write("\t".join(header)+"\n")
        for key in logFC_FDR_significant:
            gene_row=key+"\t"+"\t".join(logFC_FDR_significant[key])+"\n"
            out2.write(gene_row)

#stampo come stderr le statistiche

    stat_logFC_FDR_allgenes="Genes in total = "+str(len(logFC_FDR_allgenes))+"\n"
    sys.stderr.write(stat_logFC_FDR_allgenes)

    threshold_used="threshold to call significant = "+str(args.th)+"\n"
    sys.stderr.write(threshold_used)

    stat_logFC_FDR_significant="Genes significant in at least one condition = "+ str(len(logFC_FDR_significant))+"\n"
    sys.stderr.write(stat_logFC_FDR_significant)




if __name__ == "__main__":
    main()
