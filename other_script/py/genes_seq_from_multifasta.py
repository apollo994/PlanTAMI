#questo tool legge un multi fasta e restituisce solo le sequenze di quelli elencati in un altro file input

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='Tool for get the longest transcript from fasta')
    parser.add_argument('--in_fa', metavar='input_fasta', help='The input multi-fasta file')
    parser.add_argument('--in_gene', metavar='input_gene_list', help='genes list to be extract')
    parser.add_argument('--output', metavar='output', default="multi_fa_gene_list.fa", help='The output file.')

    args = parser.parse_args()

    genes_list=[]

    with open (args.in_gene) as genes_file:
        for line in genes_file:
            genes_list.append(line.strip("\n"))

    all_my_fasta={}

    with open (args.in_fa) as fasta:

        gene_name=""
        gene_seq=""

        for line in fasta:
            if line[0]==">":
                #per pomodoro devo togliere "mrna:" e il .ntrans
                #gene_name=">"+line[6:20]

                gene_name=line.strip("\n")

            else:
                gene_seq=line.strip("\n")

            all_my_fasta[gene_name]=gene_seq

            if len(gene_seq)!=0:
                gene_name=""
                gene_seq=""


    for gene in genes_list:
        fasta_style_gene=">"+gene
        if fasta_style_gene in all_my_fasta:
            print (fasta_style_gene)
            print (all_my_fasta[fasta_style_gene])


if __name__ == "__main__":
    main()
