import sys
import argparse
from collections import Counter
from scipy import stats


def main():
    parser = argparse.ArgumentParser(description='tool for plazaID count and statistical test with fisher method')
    parser.add_argument('--input', metavar='input', help='plazaID list to validate')
    parser.add_argument('--all_plaza', metavar='all_plaza', help='all plazaID for the species in analysis')
    parser.add_argument('--output', metavar='output',default="fisher_test_results.txt" ,help='output name')
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

    #calcolo la frequenza per il mio set e per tutti i geni

    my_count=compute_count_list(my_plaza_ID)
    all_count=compute_count_list(all_plaza_ID)

################################################################################

    #genero un dizionario con le famiglie del mio set, i relativi conti
    #e i conti nel set totale (all_count)

    my_vs_all_count_enr_pv={}

    for fam in my_count:
        if fam in all_count:
            my_array=[]
            my_array.append([float(my_count[fam]),len(my_plaza_ID)])
            my_array.append([float(all_count[fam]),len(all_plaza_ID)])

            oddsratio, pvalue = stats.fisher_exact(my_array, alternative="greater")
            my_vs_all_count_enr_pv[fam]=[my_count[fam],all_count[fam],round(oddsratio,4),round(pvalue,4)]


################################################################################

    #genero l'output

    with open (args.output, "w") as out:
        header="plazaID"+"\t"+"my_count"+"\t"+"all_count"+"\t"+"enrichment"+"\t"+"pvalue"+"\n"
        out.write(header)
        for key in my_vs_all_count_enr_pv:
            plazaID=str(key)
            my_count=str(my_vs_all_count_enr_pv[key][0])
            all_count=str(my_vs_all_count_enr_pv[key][1])
            enrichment=str(my_vs_all_count_enr_pv[key][2])
            pvalue=str(my_vs_all_count_enr_pv[key][3])
            the_line=plazaID+"\t"+my_count+"\t"+all_count+"\t"+enrichment+"\t"+pvalue+"\n"
            out.write(the_line)


if __name__ == "__main__":
    main()
