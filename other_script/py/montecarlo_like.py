import sys
import argparse
from collections import defaultdict
from collections import Counter
import re
import random


def main():
    parser = argparse.ArgumentParser(description='tool for plazaID counting and statistical analysis with montecarlo like test')
    parser.add_argument('--input', metavar='input', help='plazaID list to validate')
    parser.add_argument('--all_plaza', metavar='all_plaza', help='all plazaID for the species in analysis')
    parser.add_argument('--random', metavar='random',type=int ,default="10000" ,help='number of list of random plazaID to generate (def=10000)')
    parser.add_argument('--output', metavar='output',default="montecarlo_like_test_results.txt" ,help='output name')
    args = parser.parse_args()

################################################################################

    #genero una lista contente i plazaID del mio set
    my_plaza_ID=[]

    with open (args.input) as my_plaza:
        for line in my_plaza:
            my_plaza_ID.append(line.strip("\n"))

    #genero una lista contente tutti i plazaID
    all_plaza_ID=[]

    with open (args.all_plaza) as all_plaza:
        for line in all_plaza:
            all_plaza_ID.append(line.strip("\n"))

################################################################################

    #funzione per calcolare le volte che vedo una famiglia nel set analizzato

    def compute_count_list(plazaID_list):
        plazaID_freq=dict(Counter(plazaID_list))

        return plazaID_freq

################################################################################

    #genero la tabella di conti del mio dataset

    my_count=compute_count_list(my_plaza_ID)

################################################################################

    #genero un dizionario con chiavi uguali a i miei plazaID in cui vado ad
    #inserire i vari conti ogni volta che genero un dataset random

    random_count_diz={}

    for key in my_count:
        random_count_diz[key]=[]

################################################################################

    #ciclo per generare n (n=random) dataset random e inserire i valori nel
    #dizionario random count se tra i plazaID pescati a caso trovo uno dei
    #plazaID del mio dataset

    counter=0

    while counter<args.random:
        random_list=random.sample(all_plaza_ID,len(my_plaza_ID))
        random_count=compute_count_list(random_list)
        for key in random_count_diz:
            if key in random_count:
                random_count_diz[key].append(float(random_count[key]))
            else:
                random_count_diz[key].append(float(0))

        counter=counter+1

################################################################################

    #conto le volte i valori di una chiave sono piu alti del valore della chiave
    #nel dizionario my count

    result_diz={}

    for key in my_count:
        observed=float(my_count[key])
        more_or_equal=float(0)
        for element in random_count_diz[key]:
            if element>=observed:
                more_or_equal=more_or_equal+1
        result_diz[key]=[observed,more_or_equal,(more_or_equal/args.random)]

################################################################################

    #stampo il file di output

    with open (args.output,"w") as out:
        header="plazaID"+"\t"+"my_count"+"\t"+"count_in_"+str(args.random)+"_random_sample"+"\t"+"pvalue"+"\n"
        out.write(header)
        for key in result_diz:
            plazaID=str(key)
            my_count=str(result_diz[key][0])
            count_in_random=str(result_diz[key][1])
            pvalue=str(result_diz[key][2])
            the_line=plazaID+"\t"+my_count+"\t"+count_in_random+"\t"+pvalue+"\n"
            out.write(the_line)


if __name__ == "__main__":
    main()
