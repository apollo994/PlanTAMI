#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Toll for conversion of rice IDs')
    parser.add_argument('--input', metavar='INPUTFILE' ,default="/dev/stdin", help='List of IDs to be converted')
    parser.add_argument('--all', metavar='ALL' ,default= "/Users/Apollo199/MEGA/MEGAsync Uploads/tirocinio/Riso/RAP-MSU_2018-03-29.txt", help='Complete IDs tuple')
    parser.add_argument('--output', metavar='OUTPUTFILE',default="con", help='"con" for converted, "tup" for IDs tuple (def="con")')
    parser.add_argument('--sample', metavar='sample',default="/dev/stdout" ,help='name for the output')


    args = parser.parse_args()

    IDs_to_be_converted=[]

    with open (args.input) as input_file:
        for line in input_file:
            IDs_to_be_converted.append(line.strip("\n"))


    ID_type="RAP"

    if IDs_to_be_converted[0]=="L":
        ID_type="MUS"


    all_IDs={}

    with open (args.all) as input_all:

        if ID_type=="RAP":
            for line in input_all:

                line=line.split("\t")
                RAP_ID=line[0]
                MUS_ID=line[1].strip("\n")

                if RAP_ID!="None" and MUS_ID!="None":
                    MUS_ID=MUS_ID[0:14]
                    all_IDs[RAP_ID]=MUS_ID

        if ID_type=="MUS":
            for line in input_file:

                line=line.split("\t")
                RAP_ID=line[0]
                MUS_ID=line[1].strip("\n")

                if RAP_ID!="None" and MUS_ID!="None":
                    MUS_ID=MUS_ID[0:14]
                    all_IDs[MUS_ID]=RAP_ID


    with open (args.sample, "w") as out:
        if args.output=="con":
            for ID in IDs_to_be_converted:
                if ID in all_IDs:
                    out.write(all_IDs[ID]+"\n")
                else:
                    sys.stderr.write(ID+" NOT FOUND"+"\n")

        if args.output=="tup":
            for ID in IDs_to_be_converted:
                if ID in all_IDs:
                    out.write(ID+"\t"+all_IDs[ID]+"\n")
                else:
                    sys.stderr.write(ID+" NOT FOUND"+"\n")


if __name__ == "__main__":
  main()
