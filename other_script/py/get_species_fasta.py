#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Extract from multifasta just the sequence of ath,sly,osa and hvu')
    parser.add_argument('--input', metavar='INPUTFILE', default="/dev/stdin", help='The multifasta file.')

    args = parser.parse_args()


    multi_fasta=[]

    file_name=args.input.split("/")
    file_name=file_name[len(file_name)-1]

    print (file_name + "......done")

    with open (args.input) as my_input:
        for line in my_input:
            multi_fasta.append(line)
        with open ("my_sp_"+file_name , "w") as my_output:
            for i in range(0,len(multi_fasta)):
                #trova arabidopsis
                if multi_fasta[i][0:3]==">AT" and multi_fasta[i][4]=="G":
                    my_output.write(multi_fasta[i])
                    my_output.write(multi_fasta[i+1])
                #trova tomato
                if multi_fasta[i][0:4]==">Sol":
                    my_output.write(multi_fasta[i])
                    my_output.write(multi_fasta[i+1])
                #trova rice
                if multi_fasta[i][0:4]==">LOC":
                    my_output.write(multi_fasta[i])
                    my_output.write(multi_fasta[i+1])
                #trova barley
                if multi_fasta[i][0:4]==">HVU":
                    my_output.write(multi_fasta[i])
                    my_output.write(multi_fasta[i+1])



if __name__ == "__main__":
  main()
