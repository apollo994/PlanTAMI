import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='My nice tool.')
    parser.add_argument('--input', metavar="input", help='The input file.')
    parser.add_argument('--head', metavar="head", default="t",help='specify if the header is present, (t or f, def=t)')
    parser.add_argument('--col', metavar="col", default=0, type=int, help='specify the index of the column (def=0)')
    parser.add_argument('--sep', metavar="sep", default="tab",help='separator of the column(tab, comma or semi, def=tab)')
    parser.add_argument('--unique', metavar="unique", default="f",help='option to get also un unique element of the desired comun (def=f)')
    parser.add_argument('--out', metavar='out', default="def", help='The output file')

    args = parser.parse_args()

################################################################################

    if args.out=="def":
        input_file=args.input.split("/")
        input_file=input_file[len(input_file)-1]
        output_file="selected_column_"+input_file
    else:
        output_file=args.out

    element_to_print=[]

    with open (args.input) as my_file:
        if args.head=="t":
            head=my_file.readline()
        for line in my_file:
            if args.sep=="tab":
                line=line.split("\t")
                line=line[args.col].strip("\n")
                element_to_print.append(line+"\n")
            if args.sep=="comma":
                line=line.split(",")
                line=line[args.col].strip("\n")
                element_to_print.append(line+"\n")
            if args.sep=="semi":
                line=line.split(";")
                line=line[args.col].strip("\n")
                element_to_print.append(line+"\n")


    with open (output_file, "w") as out:
        for element in element_to_print:
            out.write(element)

    if args.unique=="t":
        with open ("unique_"+output_file, "w") as out:
            for element in set(element_to_print):
                out.write(element)

if __name__ == "__main__":
    main()
