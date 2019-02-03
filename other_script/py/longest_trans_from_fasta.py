#questo tool legge un multi fasta (transcript ID) e ne restituisce uno nuovo con solo le
#sequenze del transcritto piu lungo per ogni gene, inoltre genera anche un file con
#gene e relativi trascritti con lunghezza

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='Tool for get the longest transcript from fasta')
    parser.add_argument('--input', metavar='input', help='The input multi-fasta file')
    parser.add_argument('--output', metavar='output', default="longest_transcript.fa", help='The output file.')

    args = parser.parse_args()

    my_new_fasta={}
    for_fasta_stat={}

    with open (args.input) as fasta:

        trans_ID=""
        gene_ID=""
        fasta_seq=""
        fasta_seq_len=0

        for line in fasta:

            if line[0]==">":
                trans_ID=line[1:len(line)].strip("\n")
                gene_ID=line[1:17]
                fasta_seq=""
                fasta_seq_len=0

            else:
                fasta_seq=line.strip("\n")
                fasta_seq_len=len(line.strip("\n"))

            if gene_ID in my_new_fasta:
                if fasta_seq_len>len(my_new_fasta[gene_ID][1]):
                    my_new_fasta[gene_ID]=[trans_ID,fasta_seq]

            else:
                my_new_fasta[gene_ID]=[trans_ID,fasta_seq]

            if gene_ID in for_fasta_stat:
                if fasta_seq_len!=0:
                    for_fasta_stat[gene_ID].append([trans_ID,fasta_seq_len])

            else:
                if fasta_seq_len!=0:
                    for_fasta_stat[gene_ID]=[[trans_ID,fasta_seq_len]]


    with open ("stat_"+args.output+".txt", "w") as out_stat:
        for gene in for_fasta_stat:
            out_stat.write(gene+"\n")
            for tr in for_fasta_stat[gene]:
                out_stat.write("\t"+tr[0]+"\t"+str(tr[1])+"\n")
            out_stat.write("\n")


    with open (args.output , "w") as out:
        for gene in my_new_fasta:
            out.write(">"+gene+"\n")
            out.write(my_new_fasta[gene][1]+"\n")


if __name__ == "__main__":
    main()
