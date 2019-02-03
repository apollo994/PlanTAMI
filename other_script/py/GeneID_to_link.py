import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='tool that conver a .txt gene list in .html with links to reference website')
    parser.add_argument('--input', metavar='input', help='gene list with one ID per line')
    parser.add_argument('--sp', metavar='sp', help='the species of the gene list (ath, sly, hvu, osa)')

    args = parser.parse_args()

    gene_list=[]

    with open (args.input) as my_list:
        for line in my_list:
            gene_list.append(line.strip("\n"))

    #salvo in una stringa il nome del file di input
    file_name=args.input.split("/")
    file_name=file_name[len(file_name)-1]
    file_name=file_name.strip(".txt")

    final_file_name=file_name+".html"

    if args.sp=="ath":
        with open (final_file_name, "w") as output:
            output.write(file_name+"<br>")
            for gene in gene_list:
                output.write("<a href="+"https://www.arabidopsis.org/servlets/TairObject?type=locus&name="+gene+">"+gene+"</a>"+"<br>")
                output.write("\n")

    if args.sp=="osa":
        with open (final_file_name, "w") as output:
            output.write(file_name+"<br>")

    if args.sp=="sly":
        with open (final_file_name, "w") as output:
            output.write(file_name+"<br>")

    if args.sp=="hvu":
        with open (final_file_name, "w") as output:
            output.write(file_name+"<br>")

    else:
        sys.stderr.write("Species NOT detected, specify with --sp parameter ")

if __name__ == "__main__":
    main()
