#versione due, prende come input, non piu solo la lista di plazaID ma la lista di plazaID
#con reativo speciesID, in questo modo nell'output salvo anche la family con
#i relativi geni appartenenti alla famiglia nelle due specie analizzate

import sys
import argparse
from collections import defaultdict
import re
import random
from collections import Counter
import math
from decimal import Decimal
import os
from scipy.stats import binom

def main():
    parser = argparse.ArgumentParser(description='Tool for normalised pointwise mutal information computation across two DE gene list')
    parser.add_argument('--in_sp1', metavar='in_sp1', help='The significant plazaID of species 1')
    parser.add_argument('--in_sp2', metavar='in_sp2', help='The significant plazaID of species 2')
    parser.add_argument('--all_sp1', metavar='all_sp1', help='All plazaID of species 1')
    parser.add_argument('--all_sp2', metavar='all_sp2', help='All plazaID of species 2')
    parser.add_argument('--random', metavar='random',type=int ,default="10000" ,help='number of list of random plazaID to generate (def=10000)')
    parser.add_argument('--sample', metavar='sample', default="my_sample_result", help='sample name')

    args = parser.parse_args()

################################################################################

    #funzione per creare cartelle se non ancora esistenti
    def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)

    #creo la cartella dove salvare i risultati

    createFolder("./"+args.sample+"_results/")

################################################################################
    #liste degli input

    #genero lista contente tutti gli id significativi delle prima specie
    sig_ID_sp1=[]

    with open (args.in_sp1) as sig_plaza_sp1:
        for line in sig_plaza_sp1:
            line=line.split("\t")
            sig_ID_sp1.append(line[0].strip("\n"))

    #genero lista contente tutti gli id significativi della seconda specie
    sig_ID_sp2=[]

    with open (args.in_sp2) as sig_plaza_sp2:
        for line in sig_plaza_sp2:
            line=line.split("\t")
            sig_ID_sp2.append(line[0].strip("\n"))

    #genero una lista contente tutti i plazaID
    all_ID_sp1=[]

    with open (args.all_sp1) as all_plaza_sp1:
        for line in all_plaza_sp1:
            line=line.split("\t")
            all_ID_sp1.append(line[0].strip("\n"))

    #genero una lista contente tutti i plazaID
    all_ID_sp2=[]

    with open (args.all_sp2) as all_plaza_sp2:
        for line in all_plaza_sp2:
            line=line.split("\t")
            all_ID_sp2.append(line[0].strip("\n"))

    sys.stdout.write("significant ID in species 1 = "+str(len(sig_ID_sp1))+"\n")
    sys.stdout.write("significant ID in species 2 = "+str(len(sig_ID_sp2))+"\n")
    sys.stdout.write("total ID in species 1 = "+str(len(all_ID_sp1))+"\n")
    sys.stdout.write("total ID in species 2 = "+str(len(all_ID_sp2))+"\n")

################################################################################

    def common_elements(list1, list2):
        return list(set(list1) & set(list2))

    intersection_sig=common_elements(sig_ID_sp1,sig_ID_sp2)

    sys.stdout.write("significant common ID in sp1 and sp2 (unique) = "+str(len(intersection_sig))+"\n")

################################################################################

    #genero delle nuove liste di significativi eliminando quelli non presenti
    #nella lista degli ID comuni

    filt_sig_ID_sp1=[]
    filt_sig_ID_sp2=[]

    for element in sig_ID_sp1:
        if element in intersection_sig:
            filt_sig_ID_sp1.append(element)

    sys.stdout.write("ID after filter (also repeted) in sp1 = "+str(len(filt_sig_ID_sp1))+"\n")

    for element in sig_ID_sp2:
        if element in intersection_sig:
            filt_sig_ID_sp2.append(element)

    sys.stdout.write("ID after filter (also repeted) in sp2 = "+str(len(filt_sig_ID_sp2))+"\n")

################################################################################

    #funzione per calcolare le volte che vedo una famiglia nel set analizzato

    def compute_count_list(plazaID_list):
        plazaID_freq=dict(Counter(plazaID_list))

        return plazaID_freq

################################################################################

    #funzione per calcolare la npmi

    def compute_npmi(p_my,p_all):
        npmi= (math.log(p_my/p_all, 2.0))/(-(math.log(p_all,2.0)))

        return npmi

################################################################################

    #funzione per calcolare 1 - distribuzione binomiale cumulativa

    def bin_dist(n_obs,n_my_gene,p_all):
        k=n_obs-1
        n=n_my_gene
        p=p_all

        return round(1-binom.cdf(k,n,p),6)

################################################################################

    #dizionari in cui metto come chiavi gli ID risultanti dall'intersezione dei
    #significativi e come chiavi una lista dei valori di npmi generati dalla lista random

    random_npmi_res_sp1={}
    random_npmi_res_sp2={}

    for key in intersection_sig:
        random_npmi_res_sp1[key]=[]
        random_npmi_res_sp2[key]=[]

    #dizionari con conti di tutti gli delle due specie
    count_all_sp1=compute_count_list(all_ID_sp1)
    count_all_sp2=compute_count_list(all_ID_sp2)


    #ciclo per generare n liste random con elementi pari alle liste di partenza
    #di sp1 e sp2

    counter=1

    while counter<=args.random:

        #genero liste random con la stessa lunghezza della lista dei significativa
        random_ID_sp1=random.sample(all_ID_sp1,len(sig_ID_sp1))

        random_ID_sp2=random.sample(all_ID_sp2,len(sig_ID_sp2))

        #intersezione delle due liste rando
        intersection_random=common_elements(random_ID_sp1,random_ID_sp2)

        #cicli per eliminare dalle liste random gli ID non presenti nell'altra lista
        #e generare un dizionari con i conti degli ID rimasti
        filt_random_ID_sp1=[]
        filt_random_ID_sp2=[]

        for element in random_ID_sp1:
            if element in intersection_random:
                filt_random_ID_sp1.append(element)

        count_random_sp1=compute_count_list(filt_random_ID_sp1)

        for element in random_ID_sp2:
            if element in intersection_random:
                filt_random_ID_sp2.append(element)

        count_random_sp2=compute_count_list(filt_random_ID_sp2)

        #cicli che calcolano la npmi degli ID del mio set di ID partenza ma su conti random
        #separatamente di sp1 e sp2

        for key in random_npmi_res_sp1:
            if key in count_random_sp1:
                p_random_ID_sp1=float(count_random_sp1[key])/len(sig_ID_sp1)
                p_all_ID_sp1=float(count_all_sp1[key])/len(all_ID_sp1)

                npmi_random_ID_sp1=compute_npmi(p_random_ID_sp1,p_all_ID_sp1)

                random_npmi_res_sp1[key].append(npmi_random_ID_sp1)

            else:
                random_npmi_res_sp1[key].append(float(-1))

        for key in random_npmi_res_sp2:
            if key in count_random_sp2:
                p_random_ID_sp2=float(count_random_sp2[key])/len(sig_ID_sp2)
                p_all_ID_sp2=float(count_all_sp2[key])/len(all_ID_sp2)

                npmi_random_ID_sp2=compute_npmi(p_random_ID_sp2,p_all_ID_sp2)

                random_npmi_res_sp2[key].append(npmi_random_ID_sp2)

            else:
                random_npmi_res_sp2[key].append(float(-1))

        sys.stderr.write("generated "+str(counter)+" random list out of "+str(args.random)+"\n")

        counter=counter+1

################################################################################

    count_sig_sp1=compute_count_list(filt_sig_ID_sp1)
    count_sig_sp2=compute_count_list(filt_sig_ID_sp2)
    count_all_sp1=compute_count_list(all_ID_sp1)
    count_all_sp2=compute_count_list(all_ID_sp2)

    #calcolo la probabilita di vedere un ID nell insieme analizzato, nel caso
    #degli ID filtrati divido per il numero di ID filtrati, nel caso della lista
    #di tutti gli ID divido per la lunghezza della lista all_ID_sp

    result_diz={}

    for element in intersection_sig:
        result_diz[element]=[]

    for element in intersection_sig:
        result_diz[element]=[]
        result_diz[element].append((float(count_sig_sp1[element]))/(float(len(sig_ID_sp1))))
        result_diz[element].append((float(count_sig_sp2[element]))/(float(len(sig_ID_sp2))))
        result_diz[element].append((float(count_all_sp1[element]))/(float(len(all_ID_sp1))))
        result_diz[element].append((float(count_all_sp2[element]))/(float(len(all_ID_sp2))))

    #aggiungo un elemnto alla lista di valoi che corriponde al coefficenti di
    #mutual information

    for element in result_diz:
        p_my_sp1=result_diz[element][0]
        p_my_sp2=result_diz[element][1]
        p_all_sp1=result_diz[element][2]
        p_all_sp2=result_diz[element][3]

        ###############################

        #calcolo npmi sp1
        npmi_sp1=compute_npmi(p_my_sp1,p_all_sp1)
        result_diz[element].append(npmi_sp1)

        #calcolo pvalue sp1
        higher=0

        for value in random_npmi_res_sp1[element]:
            if value>=npmi_sp1:
                higher=higher+1


        p_value_sp1=float(higher)/len(random_npmi_res_sp1[element])
        result_diz[element].append(p_value_sp1)

        #calcolo FDR sp1
        FDR_sp1=p_value_sp1*len(result_diz)
        result_diz[element].append(FDR_sp1)

        #calcolo numero geni nella familgia per sp1
        prop_sp1=str(count_sig_sp1[element])+"/"+str(count_all_sp1[element])
        result_diz[element].append(prop_sp1)

        #calcolo valore binomial distribution della familgia per sp1
        p_bin_dist_sp1=bin_dist(float(count_sig_sp1[element]),float(len(sig_ID_sp1)),float(p_all_sp1))
        result_diz[element].append(p_bin_dist_sp1)

        ##################################

        #calcolo npmi sp2
        npmi_sp2=compute_npmi(p_my_sp2,p_all_sp2)
        result_diz[element].append(npmi_sp2)

        #calcolo pvalue sp2
        higher=0

        for value in random_npmi_res_sp2[element]:
            if value>=npmi_sp2:
                higher=higher+1

        p_value_sp2=float(higher)/len(random_npmi_res_sp2[element])
        result_diz[element].append(p_value_sp2)

        #calcolo FDR sp2
        FDR_sp2=p_value_sp2*len(result_diz)
        result_diz[element].append(FDR_sp2)

        #calcolo numero geni nella familgia per sp2
        prop_sp2=str(count_sig_sp2[element])+"/"+str(count_all_sp2[element])
        result_diz[element].append(prop_sp2)

        #calcolo valore binomial distribution della familgia per sp2
        p_bin_dist_sp2=bin_dist(float(count_sig_sp2[element]),float(len(sig_ID_sp2)),float(p_all_sp2))
        result_diz[element].append(p_bin_dist_sp2)

################################################################################

    #leggo i file di input ma in questo caso mi salvo anche i geni relativi a ogni
    #famiglia in modo da poterli contare e annotare sul file family degli output

    #liste degli input

    #genero lista contente tutti gli plazaID e spID significativi della prima specie
    sig_plaza_spID_sp1={}

    with open (args.in_sp1) as file_sig_plaza_spID_sp1:
        for line in file_sig_plaza_spID_sp1:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in sig_plaza_spID_sp1:
                sig_plaza_spID_sp1[plazaID].append(spID)
            else:
                sig_plaza_spID_sp1[plazaID]=[spID]

    #genero lista contente tutti gli plazaID e spID significativi della seconda specie
    sig_plaza_spID_sp2={}

    with open (args.in_sp2) as file_sig_plaza_spID_sp2:
        for line in file_sig_plaza_spID_sp2:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in sig_plaza_spID_sp2:
                sig_plaza_spID_sp2[plazaID].append(spID)
            else:
                sig_plaza_spID_sp2[plazaID]=[spID]


    #genero una lista contente tutti i plazaID, spID di sp1
    all_plaza_spID_sp1={}

    with open (args.all_sp1) as file_all_plaza_spID_sp1:
        for line in file_all_plaza_spID_sp1:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in all_plaza_spID_sp1:
                all_plaza_spID_sp1[plazaID].append(spID)
            else:
                all_plaza_spID_sp1[plazaID]=[spID]

    #genero una lista contente tutti i plazaID, spID di sp2
    all_plaza_spID_sp2={}

    with open (args.all_sp2) as file_all_plaza_spID_sp2:
        for line in file_all_plaza_spID_sp2:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in all_plaza_spID_sp2:
                all_plaza_spID_sp2[plazaID].append(spID)
            else:
                all_plaza_spID_sp2[plazaID]=[spID]

    #genero un file contente le le famiglie in comune, i geni e la proporzione
    #rispetto al set completo, per entrambe le specie

    family_result={}

    for family in intersection_sig:
        gene_list_sp1=sig_plaza_spID_sp1[family]
        prop_gene_sp1=str(len(sig_plaza_spID_sp1[family]))+"/"+str(len(all_plaza_spID_sp1[family]))
        gene_list_sp2=sig_plaza_spID_sp2[family]
        prop_gene_sp2=str(len(sig_plaza_spID_sp2[family]))+"/"+str(len(all_plaza_spID_sp2[family]))

        family_result[family]=[gene_list_sp1, prop_gene_sp1, gene_list_sp2, prop_gene_sp2]

################################################################################

    #genero due nuovi dizionari per le due specie in modo da poter scrivere
    #file di out con le specifiche gene per gene

    result_diz_sp1={}
    result_diz_sp2={}

    with open (args.in_sp1) as file_sig_plaza_spID_sp1:
        for line in file_sig_plaza_spID_sp1:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in intersection_sig:
                result_diz_sp1[spID]=[plazaID]
                #npmi_value
                result_diz_sp1[spID].append(str(round(result_diz[plazaID][4],7)))
                #pv_npmi
                result_diz_sp1[spID].append(str(round(result_diz[plazaID][5],7)))
                #number of gene
                result_diz_sp1[spID].append(str(result_diz[plazaID][7]))
                #p_bin_dist
                result_diz_sp1[spID].append(str(round(result_diz[plazaID][8],7)))


    with open (args.in_sp2) as file_sig_plaza_spID_sp2:
        for line in file_sig_plaza_spID_sp2:
            line=line.split("\t")
            plazaID=line[0]
            spID=line[1].strip("\n")
            if plazaID in intersection_sig:
                result_diz_sp2[spID]=[plazaID]
                result_diz_sp2[spID].append(str(round(result_diz[plazaID][9],7)))
                result_diz_sp2[spID].append(str(round(result_diz[plazaID][10],7)))
                result_diz_sp2[spID].append(str(result_diz[plazaID][12]))
                result_diz_sp2[spID].append(str(round(result_diz[plazaID][13],7)))

################################################################################

    #scrivo la tabella dei risultati, i nomi dell'intersezione e le statistiche

    with open ("./"+args.sample+"_results/"+args.sample+"_family.txt", "w") as out:
        out.write("plazaID"+"\t"+"p_my_sp1"+"\t"+"p_my_sp2"+"\t"+"p_all_sp1"+"\t"+"p_all_sp2"+"\t"+"npmi_sp1"+"\t"+"pv_npmi_sp1"+"\t"+"FDR_sp1"+"\t"+"n_gene_sp1"+"\t"+"p_bin_dist_sp1"+"\t"+"npmi_sp2"+"\t"+"pv_npmi_sp2"+"\t"+"FDR_sp2"+"\t"+"n_gene_sp2"+"\t"+"p_bin_dist_sp2"+"\n")
        for family in result_diz:
            p_my_sp1=str(round(result_diz[family][0],7))
            p_my_sp2=str(round(result_diz[family][1],7))
            p_all_sp1=str(round(result_diz[family][2],7))
            p_all_sp2=str(round(result_diz[family][3],7))
            npmi_sp1=str(round(result_diz[family][4],7))
            pv_npmi_sp1=str(round(result_diz[family][5],7))
            FDR_sp1=str(round(result_diz[family][6],7))
            n_gene_sp1=str(result_diz[family][7])
            p_bin_dist_sp1=str(round(result_diz[family][8],7))
            npmi_sp2=str(round(result_diz[family][9],7))
            pv_npmi_sp2=str(round(result_diz[family][10],7))
            FDR_sp2=str(round(result_diz[family][11],7))
            n_gene_sp2=str(result_diz[family][12])
            p_bin_dist_sp2=str(round(result_diz[family][13],7))


            line_to_print=family+"\t"+p_my_sp1+"\t"+p_my_sp2+"\t"+p_all_sp1+"\t"+p_all_sp2+"\t"+npmi_sp1+"\t"+pv_npmi_sp1+"\t"+FDR_sp1+"\t"+n_gene_sp1+"\t"+p_bin_dist_sp1+"\t"+npmi_sp2+"\t"+pv_npmi_sp2+"\t"+FDR_sp2+"\t"+n_gene_sp2+"\t"+p_bin_dist_sp2+"\n"

            out.write(line_to_print)

    with open ("./"+args.sample+"_results/"+args.sample+"_family_genes.txt", "w") as out2:
        for family in family_result:
            out2.write(family)
            out2.write("\t")
            out2.write(str(family_result[family][0]))
            out2.write("\t")
            out2.write(family_result[family][1])
            out2.write("\t")
            out2.write(str(family_result[family][2]))
            out2.write("\t")
            out2.write(family_result[family][3])
            out2.write("\n")


    with open ("./"+args.sample+"_results/"+args.sample+"_stat.txt", "w") as out3:
        out3.write("significant ID in species 1 = "+str(len(sig_ID_sp1))+"\n")
        out3.write("significant ID in species 2 = "+str(len(sig_ID_sp2))+"\n")
        out3.write("total ID in species 1 = "+str(len(all_ID_sp1))+"\n")
        out3.write("total ID in species 2 = "+str(len(all_ID_sp2))+"\n")
        out3.write("significant common ID in sp1 and sp2 (unique) = "+str(len(intersection_sig))+"\n")
        out3.write("ID after filter (also repeted) in sp1 = "+str(len(filt_sig_ID_sp1))+"\n")
        out3.write("ID after filter (also repeted) in sp2 = "+str(len(filt_sig_ID_sp2))+"\n")

    with open ("./"+args.sample+"_results/"+args.sample+"_genes_sp1.txt", "w") as out4:
        out4.write("geneID"+"\t"+"plazaID"+"\t"+"npmi_sp1"+"\t"+"pv_npmi_sp1"+"\t"+"n_gene_sp1"+"\t"+"p_bin_dist_sp1"+"\n")
        for gene in result_diz_sp1:
            plazaID=result_diz_sp1[gene][0]
            npmi_sp1=result_diz_sp1[gene][1]
            pv_npmi_sp1=result_diz_sp1[gene][2]
            n_gene_sp1=result_diz_sp1[gene][3]
            p_bin_dist_sp1=result_diz_sp1[gene][4]

            line_to_print=gene+"\t"+plazaID+"\t"+npmi_sp1+"\t"+pv_npmi_sp1+"\t"+n_gene_sp1+"\t"+p_bin_dist_sp1+"\n"

            out4.write(line_to_print)

    with open ("./"+args.sample+"_results/"+args.sample+"_genes_sp2.txt", "w") as out5:
        out5.write("geneID"+"\t"+"plazaID"+"\t"+"npmi_sp2"+"\t"+"pv_npmi_sp2"+"\t"+"n_gene_sp2"+"\t"+"p_bin_dist_sp2"+"\n")
        for gene in result_diz_sp2:
            plazaID=result_diz_sp2[gene][0]
            npmi_sp2=result_diz_sp2[gene][1]
            pv_npmi_sp2=result_diz_sp2[gene][2]
            n_gene_sp2=result_diz_sp2[gene][3]
            p_bin_dist_sp2=result_diz_sp2[gene][4]

            line_to_print=gene+"\t"+plazaID+"\t"+npmi_sp2+"\t"+pv_npmi_sp2+"\t"+n_gene_sp2+"\t"+p_bin_dist_sp2+"\n"

            out5.write(line_to_print)
    
################################################################################


if __name__ == "__main__":
    main()
