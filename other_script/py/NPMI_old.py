import sys
import argparse
from collections import defaultdict
import re
import random
from collections import Counter
import math
from decimal import Decimal
#pachhetto per parlare con il sistema operativo
import os
#pacchetto per la statistica in python
import scipy
from scipy.stats import binom
#pacchetto per dare comandi R
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')



def main():
    parser = argparse.ArgumentParser(description='Tool for normalised pointwise mutal information computation across two DE gene list')
    parser.add_argument('--plaza', metavar='plaza',default="/Users/Apollo199/Tirocinio/HOMO_GENES/genefamily_data.orth.txt" ,help='the table from plaza with all family, species and genes')
    parser.add_argument('--sp1', metavar='sp1',help='species 1, es. ath,sly,osa,hvu')
    parser.add_argument('--sp2', metavar='sp2',help='species 2, es. ath,sly,osa,hvu')
    parser.add_argument('--in_sp1', metavar='in_sp1', help='The significant plazaID of species 1')
    parser.add_argument('--in_sp2', metavar='in_sp2', help='The significant plazaID of species 2')
    parser.add_argument('--random', metavar='random',type=int ,default="10000" ,help='number of list of random plazaID to generate (def=10000)')
    parser.add_argument('--th', metavar='th',type=float ,default="0.05" ,help='BY correction threshold (def=0.05)')
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
            sys.exit(1)

    #creo la cartella dove salvare i risultati

    createFolder("./"+args.sample+"_results/")

    #salvo in una stringa il nome dei due campioni partendo da file di input
    sample_sp1=args.in_sp1.split("/")
    sample_sp1=sample_sp1[len(sample_sp1)-1]

    sample_sp2=args.in_sp2.split("/")
    sample_sp2=sample_sp2[len(sample_sp2)-1]

################################################################################

    #funzione per trovare intersezione tra due liste

    def common_elements(list1, list2):
        return list(set(list1) & set(list2))

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

    #funzione per ottenere gli elementi unici di una lista

    def unique(list1):
        # unique_list = []
        # for x in list1:
        #     if x not in unique_list:
        #         unique_list.append(x)
        # return unique_list
        return(list(set(list1)))

################################################################################

    #funzione per calcolare FDR con metodo BH

    def BH_FDR(serie):
        #indice per ciclare sui valori al contrario
        index=sorted(list(range(0,len(serie))), reverse=True)
        #lista in cui annotare i risultari
        BH_res=[]
        #numero di test
        n=len(serie)
        #variabili per controllare se il pv è uguale al precedente o FDR diventa piu piccolo del precedente
        ex_pv=1
        ex_adj_pv=1
        for idx in index:
            #aggiungo 1 all'indice che parte da 0 per ottenere il rank
            rank=idx+1
            #annoto il pv
            pv=serie[idx]
            #se il pv analizzato è uguale al pv precedente annoto il pv corretto come il precedente
            if pv==ex_pv:
                adj_pv=ex_adj_pv
                BH_res.append(adj_pv)
            #se invece il pv è piu basso del precedente aggiiorno ex_pv e calcolo il nuovo pv corretto
            if pv<ex_pv:
                ex_pv=pv
                #calcolo il il pv corretto
                adj_pv=pv*(n/rank)
                #verifico se adj_pv è piu piccolo del precedente
                if adj_pv<ex_adj_pv:
                    #annoto il pv corretto nella lista e nell ex_BH
                    ex_adj_pv=adj_pv
                    BH_res.append(adj_pv)
                #se non lo è uso il precedente
                else:
                    adj_pv=ex_adj_pv
                    BH_res.append(ex_adj_pv)


        #riordino la lista dei risultati in ordine crescente
        BH_res=sorted(BH_res)

        return (BH_res)

################################################################################

    #funzione per calcolare FDR con metodo BY

    def BY_FDR(serie):
        #indice per ciclare sui valori al contrario
        index=sorted(list(range(0,len(serie))), reverse=True)
        #lista in cui annotare i risultari
        BY_res=[]
        #numero di test
        n=len(serie)
        #variabili per controllare se il pv è uguale al precedente o FDR diventa piu piccolo del precedente
        ex_pv=1
        ex_adj_pv=1
        #calcolo il valore q (sommatoria di 1/rank) dell'ultima posizione da cui po sottraggo ogni volta
        q=0
        #print(q)
        for i in index:
            rank=(i+1)
            q=q+1/rank


        for idx in index:
            #aggiungo 1 all'indice che parte da 0 per ottenere il rank
            rank=idx+1
            #annoto il pv
            pv=serie[idx]
            #se il pv analizzato è uguale al pv precedente annoto il pv corretto come il precedente
            if pv==ex_pv:
                adj_pv=ex_adj_pv
                BY_res.append(adj_pv)
            #se invece il pv è piu basso del precedente aggiiorno ex_pv e calcolo il nuovo pv corretto
            if pv<ex_pv:
                ex_pv=pv
                #calcolo il il pv corretto
                adj_pv=pv*(n*q/rank)
                #verifico se adj_pv è piu piccolo del precedente
                if adj_pv<ex_adj_pv:
                    #annoto il pv corretto nella lista e nell ex_BY
                    ex_adj_pv=adj_pv
                    BY_res.append(adj_pv)
                #se non lo è uso il precedente
                else:
                    adj_pv=ex_adj_pv
                    BY_res.append(ex_adj_pv)


        #riordino la lista dei risultati in ordine crescente
        BY_res=sorted(BY_res)

        return (BY_res)

################################################################################

    #genero due dizionari per le due specie contenti le famiglie come chiave e il
    #relativo gene

    plaza_diz_sp1={}
    all_spID_sp1=[]
    all_plazaID_sp1=[]

    plaza_diz_sp2={}
    all_spID_sp2=[]
    all_plazaID_sp2=[]

    with open (args.plaza) as plaza_tab:
        for line in plaza_tab:
            line=line.split(";")
            fam=line[0]
            species=line[1]
            spID=line[2].strip("\n")
            spID=spID.strip("\r")
            #elimino il .1 o .2 in id di pomodoro
            if species=="sly":
                spID=spID[0:14]

            if species==args.sp1:
                plaza_diz_sp1[spID]=fam
                all_spID_sp1.append(spID)
                all_plazaID_sp1.append(fam)

            if species==args.sp2:
                plaza_diz_sp2[spID]=fam
                all_spID_sp2.append(spID)
                all_plazaID_sp2.append(fam)


################################################################################

    #genero un dizionario con spID e plazaID, la lista degli spID e la lista plazaID
    #il tutto per entrambe le specie

    sig_spID_plazaID_sp1={}
    sig_spID_sp1=[]
    sig_plazaID_sp1=[]

    with open (args.in_sp1) as file_sig_ID_sp1:
        for line in file_sig_ID_sp1:
            line=line.strip("\n")
            if line in plaza_diz_sp1:
                sig_spID_plazaID_sp1[line]=plaza_diz_sp1[line]
                sig_spID_sp1.append(line)
                sig_plazaID_sp1.append(plaza_diz_sp1[line])

    sig_spID_plazaID_sp2={}
    sig_spID_sp2=[]
    sig_plazaID_sp2=[]

    with open (args.in_sp2) as file_sig_ID_sp2:
        for line in file_sig_ID_sp2:
            line=line.strip("\n")
            if line in plaza_diz_sp2:
                sig_spID_plazaID_sp2[line]=plaza_diz_sp2[line]
                sig_spID_sp2.append(line)
                sig_plazaID_sp2.append(plaza_diz_sp2[line])


################################################################################

    #trovo le famiglie in comune tra sp1 e sp2

    intersection_sig=common_elements(sig_plazaID_sp1,sig_plazaID_sp2)

    #genero delle nuove liste di significativi eliminando quelli non presenti
    #nella lista dei plazaID comuni

    filt_sig_plazaID_sp1=[]
    filt_sig_plazaID_sp2=[]

    for element in sig_plazaID_sp1:
        if element in intersection_sig:
            filt_sig_plazaID_sp1.append(element)


    for element in sig_plazaID_sp2:
        if element in intersection_sig:
            filt_sig_plazaID_sp2.append(element)


################################################################################

    #dizionari in cui metto come chiavi gli ID risultanti dall'intersezione dei
    #significativi e come chiavi una lista dei valori di npmi generati dalla lista random

    random_npmi_res_sp1={}
    random_npmi_res_sp2={}

    for key in intersection_sig:
        random_npmi_res_sp1[key]=[]
        random_npmi_res_sp2[key]=[]

    #dizionari con conti di tutti gli delle due specie
    count_all_sp1=compute_count_list(all_plazaID_sp1)
    count_all_sp2=compute_count_list(all_plazaID_sp2)


    #ciclo per generare n liste random con elementi pari alle liste di partenza
    #di sp1 e sp2

    counter=1

    while counter<=args.random:

        #genero liste random con la stessa lunghezza della lista dei significativa
        random_ID_sp1=random.sample(all_plazaID_sp1,len(sig_plazaID_sp1))

        random_ID_sp2=random.sample(all_plazaID_sp2,len(sig_plazaID_sp2))

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
                p_random_ID_sp1=float(count_random_sp1[key])/len(sig_plazaID_sp1)
                p_all_ID_sp1=float(count_all_sp1[key])/len(all_plazaID_sp1)

                npmi_random_ID_sp1=compute_npmi(p_random_ID_sp1,p_all_ID_sp1)

                random_npmi_res_sp1[key].append(npmi_random_ID_sp1)

            else:
                random_npmi_res_sp1[key].append(float(-1))

        for key in random_npmi_res_sp2:
            if key in count_random_sp2:
                p_random_ID_sp2=float(count_random_sp2[key])/len(sig_plazaID_sp2)
                p_all_ID_sp2=float(count_all_sp2[key])/len(all_plazaID_sp2)

                npmi_random_ID_sp2=compute_npmi(p_random_ID_sp2,p_all_ID_sp2)

                random_npmi_res_sp2[key].append(npmi_random_ID_sp2)

            else:
                random_npmi_res_sp2[key].append(float(-1))

        sys.stdout.write("generated "+str(counter)+" random list out of "+str(args.random)+"\n")

        counter=counter+1

################################################################################

    count_sig_sp1=compute_count_list(filt_sig_plazaID_sp1)
    count_sig_sp2=compute_count_list(filt_sig_plazaID_sp2)
    count_all_sp1=compute_count_list(all_plazaID_sp1)
    count_all_sp2=compute_count_list(all_plazaID_sp2)

    result_diz_sp1={}
    result_diz_sp2={}

    #aggiungo i valori di p_my e p_all ai due dizionari

    for element in intersection_sig:
        result_diz_sp1[element]=[]
        result_diz_sp1[element].append(round((float(count_sig_sp1[element]))/(float(len(sig_plazaID_sp1))),7))
        result_diz_sp1[element].append(round((float(count_all_sp1[element]))/(float(len(all_plazaID_sp1))),7))

        result_diz_sp2[element]=[]
        result_diz_sp2[element].append(round((float(count_sig_sp2[element]))/(float(len(sig_plazaID_sp2))),7))
        result_diz_sp2[element].append(round((float(count_all_sp2[element]))/(float(len(all_plazaID_sp2))),7))

################################################################################

    #aggiungo i valori di npmi, prop_gene, p_value_bin, p_value_ran , bonf_cor alla lista di sp1

    for element in result_diz_sp1:
        p_my_sp1=result_diz_sp1[element][0]
        p_all_sp1=result_diz_sp1[element][1]

        #calcolo npmi sp1
        npmi_sp1=compute_npmi(p_my_sp1,p_all_sp1)
        result_diz_sp1[element].append(round(npmi_sp1,7))

        #calcolo prop_gene
        prop_sp1=str(count_sig_sp1[element])+"/"+str(count_all_sp1[element])
        result_diz_sp1[element].append(prop_sp1)

        #calcolo p_value_bin
        p_bin_dist_sp1=bin_dist(float(count_sig_sp1[element]),float(len(sig_plazaID_sp1)),float(p_all_sp1))
        result_diz_sp1[element].append(round(p_bin_dist_sp1,7))

        #calcolo p_value_ran
        higher=0

        for value in random_npmi_res_sp1[element]:
            if value>=npmi_sp1:
                higher=higher+1


        p_value_sp1=float(higher)/len(random_npmi_res_sp1[element])
        result_diz_sp1[element].append(p_value_sp1)

        #calcolo bonf_cor
        bonf_cor_sp1=p_value_sp1*len(result_diz_sp1)
        if bonf_cor_sp1>1:
            bonf_cor_sp1=1
        result_diz_sp1[element].append(round(bonf_cor_sp1,7))

    # for family in result_diz_sp1:
    #     print (family,result_diz_sp1[family])

    #ordino result_diz_sp1 basandomi sul p_value_ran (sesto elemento della lista)
    sorted_result_diz_sp1=sorted(result_diz_sp1.items(), key=lambda e: e[1][5])

    #salvo in una lista i p_value_ran ordinati
    ord_p_value_ran_sp1=[]
    for element in sorted_result_diz_sp1:
        ord_p_value_ran_sp1.append(element[1][5])

    #calcolo i p_value_adj con BY method con R
    BY_corr = list(stats.p_adjust(FloatVector(ord_p_value_ran_sp1), method = 'BY'))

    #calcolo i p_value_adj con BY method con mia funzione
    #BY_corr = BY_FDR(ord_p_value_ran_sp1)

    #ciclo per aggiungere il BY_corr alla lista di valori di sorted_result_diz_sp1
    #per comodita poi riscrivo il result diz
    index=range(0,len(BY_corr))
    for idx in index:
        sorted_result_diz_sp1[idx][1].append(round(BY_corr[idx],7))
        result_diz_sp1[sorted_result_diz_sp1[idx][0]]=sorted_result_diz_sp1[idx][1]

    # print ()
    # for family in result_diz_sp1:
    #     print (family,result_diz_sp1[family])
    # print ()
    # for family in sorted_result_diz_sp1:
    #     print (family[0],family[1])

################################################################################

    #aggiungo i valori di npmi, prop_gene, p_value_bin, p_value_ran , bonf_cor alla lista di sp2

    for element in result_diz_sp2:
        p_my_sp2=result_diz_sp2[element][0]
        p_all_sp2=result_diz_sp2[element][1]

        #calcolo npmi sp2
        npmi_sp2=compute_npmi(p_my_sp2,p_all_sp2)
        result_diz_sp2[element].append(round(npmi_sp2,7))

        #calcolo prop_gene
        prop_sp2=str(count_sig_sp2[element])+"/"+str(count_all_sp2[element])
        result_diz_sp2[element].append(prop_sp2)

        #calcolo p_value_bin
        p_bin_dist_sp2=bin_dist(float(count_sig_sp2[element]),float(len(sig_plazaID_sp2)),float(p_all_sp2))
        result_diz_sp2[element].append(round(p_bin_dist_sp2,7))

        #calcolo p_value_ran
        higher=0

        for value in random_npmi_res_sp2[element]:
            if value>=npmi_sp2:
                higher=higher+1


        p_value_sp2=float(higher)/len(random_npmi_res_sp2[element])
        result_diz_sp2[element].append(p_value_sp2)

        #calcolo bonf_cor
        bonf_cor_sp2=p_value_sp2*len(result_diz_sp2)
        if bonf_cor_sp2>1:
            bonf_cor_sp2=1
        result_diz_sp2[element].append(round(bonf_cor_sp2,7))

    #ordino result_diz_sp2 basandomi sul p_value_ran (sesto elemento della lista)
    sorted_result_diz_sp2=sorted(result_diz_sp2.items(), key=lambda e: e[1][5])

    #salvo in una lista i p_value_ran ordinati
    ord_p_value_ran_sp2=[]
    for element in sorted_result_diz_sp2:
        ord_p_value_ran_sp2.append(element[1][5])

    #calcolo i p_value_adj con BY method
    BY_corr = list(stats.p_adjust(FloatVector(ord_p_value_ran_sp2), method = 'BY'))

    #calcolo i p_value_adj con BY method con mia funzione
    #BY_corr = BY_FDR(ord_p_value_ran_sp2)

    #ciclo per aggiungere il BY_corr alla lista di valori di sorted_result_diz_sp2
    index=range(0,len(BY_corr))
    for idx in index:
        sorted_result_diz_sp2[idx][1].append(round(BY_corr[idx],7))
        result_diz_sp2[sorted_result_diz_sp2[idx][0]]=sorted_result_diz_sp2[idx][1]

    # print ()
    # for family in sorted_result_diz_sp2:
    #     print (family[0],family[1])

################################################################################

    #genero dizionario con plazaID e tutti i relativi geni, sia per sig che all

    sig_plazaID_to_spID_sp1={}
    for gene in sig_spID_plazaID_sp1:
        spID=gene
        plazaID=sig_spID_plazaID_sp1[gene]
        if plazaID in sig_plazaID_to_spID_sp1:
            sig_plazaID_to_spID_sp1[plazaID].append(spID)
        else:
            sig_plazaID_to_spID_sp1[plazaID]=[spID]

    #print (len(sig_plazaID_to_spID_sp1))

    sig_plazaID_to_spID_sp2={}
    for gene in sig_spID_plazaID_sp2:
        spID=gene
        plazaID=sig_spID_plazaID_sp2[gene]
        if plazaID in sig_plazaID_to_spID_sp2:
            sig_plazaID_to_spID_sp2[plazaID].append(spID)
        else:
            sig_plazaID_to_spID_sp2[plazaID]=[spID]

    #print (len(sig_plazaID_to_spID_sp2))

    all_plazaID_to_spID_sp1={}
    for gene in plaza_diz_sp1:
        spID=gene
        plazaID=plaza_diz_sp1[gene]
        if plazaID in all_plazaID_to_spID_sp1:
            all_plazaID_to_spID_sp1[plazaID].append(spID)
        else:
            all_plazaID_to_spID_sp1[plazaID]=[spID]

    #print (len(all_plazaID_to_spID_sp1))

    all_plazaID_to_spID_sp2={}
    for gene in plaza_diz_sp2:
        spID=gene
        plazaID=plaza_diz_sp2[gene]
        if plazaID in all_plazaID_to_spID_sp2:
            all_plazaID_to_spID_sp2[plazaID].append(spID)
        else:
            all_plazaID_to_spID_sp2[plazaID]=[spID]

    #print(len(all_plazaID_to_spID_sp2))

################################################################################

    #genero due dizionari con le specifiche per gene di sp1 e sp2 inoltre salvo
    #in una lista tutti i geni che hanno BY_cor < args.th

    gene_result_sp1={}
    sig_com_gene_sp1=[]

    for key in sig_spID_plazaID_sp1:
        spID=key
        plazaID=sig_spID_plazaID_sp1[key]
        if plazaID in result_diz_sp1:
            gene_result_sp1[spID]=[plazaID, result_diz_sp1[plazaID][3],result_diz_sp1[plazaID][4],result_diz_sp1[plazaID][5],result_diz_sp1[plazaID][6],result_diz_sp1[plazaID][7]]
            if result_diz_sp1[plazaID][7]<=args.th:
                sig_com_gene_sp1.append(spID)

    sorted_gene_result_sp1=sorted(gene_result_sp1.items(), key=lambda e: e[1][5])


    gene_result_sp2={}
    sig_com_gene_sp2=[]

    for key in sig_spID_plazaID_sp2:
        spID=key
        plazaID=sig_spID_plazaID_sp2[key]
        if plazaID in result_diz_sp2:
            gene_result_sp2[spID]=[plazaID, result_diz_sp2[plazaID][3],result_diz_sp2[plazaID][4],result_diz_sp2[plazaID][5],result_diz_sp2[plazaID][6],result_diz_sp2[plazaID][7]]
            if result_diz_sp2[plazaID][7]<=args.th:
                sig_com_gene_sp2.append(spID)

    sorted_gene_result_sp2=sorted(gene_result_sp2.items(), key=lambda e: e[1][5])

################################################################################

    #genero un dizionario contente le le famiglie in comune, i geni e la proporzione
    #rispetto al set completo, per entrambe le specie

    family_result={}

    for family in intersection_sig:
        gene_list_sp1=sig_plazaID_to_spID_sp1[family]
        prop_gene_sp1=str(len(sig_plazaID_to_spID_sp1[family]))+"/"+str(len(all_plazaID_to_spID_sp1[family]))
        gene_list_sp2=sig_plazaID_to_spID_sp2[family]
        prop_gene_sp2=str(len(sig_plazaID_to_spID_sp2[family]))+"/"+str(len(all_plazaID_to_spID_sp2[family]))

        family_result[family]=[gene_list_sp1, prop_gene_sp1, gene_list_sp2, prop_gene_sp2]

################################################################################

    #scrivo output

    with open ("./"+args.sample+"_results/"+args.sample+"_stat.txt", "w") as stat_out:
        stat_out.write(args.sample+"\n")
        stat_out.write("sample_sp1 = "+sample_sp1+"\n")
        stat_out.write("sample_sp2 = "+sample_sp2+"\n")
        stat_out.write("random lists generated = "+str(args.random)+"\n")
        stat_out.write("\n")
        stat_out.write("species 1 = "+args.sp1+"\n")
        stat_out.write("total genes in species 1 = "+str(len(all_spID_sp1))+"\n")
        stat_out.write("total families in species 1 = "+str(len(unique(all_plazaID_sp1)))+"\n")
        stat_out.write("significant genes in species 1 = "+str(len(sig_spID_sp1))+"\n")
        stat_out.write("significant families in species 1 = "+str(len(unique(sig_plazaID_sp1)))+"\n")
        stat_out.write("\n")
        stat_out.write("species 2 = "+args.sp2+"\n")
        stat_out.write("total genes in species 2 = "+str(len(all_spID_sp2))+"\n")
        stat_out.write("total families in species 2 = "+str(len(unique(all_plazaID_sp2)))+"\n")
        stat_out.write("significant genes in species 2 = "+str(len(sig_spID_sp2))+"\n")
        stat_out.write("significant families in species 2 = "+str(len(unique(sig_plazaID_sp2)))+"\n")
        stat_out.write("\n")
        stat_out.write("significant families in sp1 and sp2 = "+str(len(intersection_sig))+"\n")
        stat_out.write("Genes in common in sp1 = "+str(len(filt_sig_plazaID_sp1))+"\n")
        stat_out.write("Significant (th= "+str(args.th)+") genes in common in sp1 = "+str(len(sig_com_gene_sp1))+"\n")
        stat_out.write("Genes in common in sp2 = "+str(len(filt_sig_plazaID_sp2))+"\n")
        stat_out.write("Significant (th= "+str(args.th)+") genes in common in sp2 = "+str(len(sig_com_gene_sp2))+"\n")

    #risutati famiglie e relativi geni e proporzioni di sp1 e sp2

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

    #risultati famiglia sp1

    with open ("./"+args.sample+"_results/"+args.sample+"_family_sp1.txt", "w") as out:
        out.write("plazaID"+"\t"+"p_my_sp1"+"\t"+"p_all_sp1"+"\t"+"npmi"+"\t"+"n_genes"+"\t"+"pv_bin"+"\t"+"pv_ran"+"\t"+"Bonf_cor"+"\t"+"BY_cor"+"\n")
        for family in sorted_result_diz_sp1:
            plazaID=family[0]
            p_my_sp1=str(family[1][0])
            p_all_sp1=str(family[1][1])
            npmi=str(family[1][2])
            n_genes=str(family[1][3])
            pv_bin=str(family[1][4])
            pv_ran=str(family[1][5])
            Bonf_cor=str(family[1][6])
            BY_cor=str(family[1][7])

            line_to_print=plazaID+"\t"+p_my_sp1+"\t"+p_all_sp1+"\t"+npmi+"\t"+n_genes+"\t"+pv_bin+"\t"+pv_ran+"\t"+Bonf_cor+"\t"+BY_cor+"\n"

            out.write(line_to_print)

    #risultati famiglia sp2

    with open ("./"+args.sample+"_results/"+args.sample+"_family_sp2.txt", "w") as out3:
        out3.write("plazaID"+"\t"+"p_my_sp2"+"\t"+"p_all_sp2"+"\t"+"npmi"+"\t"+"n_genes"+"\t"+"pv_bin"+"\t"+"pv_ran"+"\t"+"Bonf_cor"+"\t"+"BY_cor"+"\n")
        for family in sorted_result_diz_sp2:
            plazaID=family[0]
            p_my_sp2=str(family[1][0])
            p_all_sp2=str(family[1][1])
            npmi=str(family[1][2])
            n_genes=str(family[1][3])
            pv_bin=str(family[1][4])
            pv_ran=str(family[1][5])
            Bonf_cor=str(family[1][6])
            BY_cor=str(family[1][7])

            line_to_print=plazaID+"\t"+p_my_sp2+"\t"+p_all_sp2+"\t"+npmi+"\t"+n_genes+"\t"+pv_bin+"\t"+pv_ran+"\t"+Bonf_cor+"\t"+BY_cor+"\n"

            out3.write(line_to_print)

    #risultati geni sp1

    with open ("./"+args.sample+"_results/"+args.sample+"_genes_sp1.txt", "w") as out4:
        out4.write("geneID"+"\t"+"plazaID"+"\t"+"n_genes"+"\t"+"pv_bin"+"\t"+"pv_ran"+"\t"+"Bonf_cor"+"\t"+"BY_cor"+"\n")
        for family in sorted_gene_result_sp1:
            geneID=family[0]
            plazaID=str(family[1][0])
            n_genes=str(family[1][1])
            pv_bin=str(family[1][2])
            pv_ran=str(family[1][3])
            Bonf_cor=str(family[1][4])
            BY_cor=str(family[1][5])

            line_to_print=geneID+"\t"+plazaID+"\t"+n_genes+"\t"+pv_bin+"\t"+pv_ran+"\t"+Bonf_cor+"\t"+BY_cor+"\n"

            out4.write(line_to_print)

    #risultati geni sp2

    with open ("./"+args.sample+"_results/"+args.sample+"_genes_sp2.txt", "w") as out5:
        out5.write("geneID"+"\t"+"plazaID"+"\t"+"n_genes"+"\t"+"pv_bin"+"\t"+"pv_ran"+"\t"+"Bonf_cor"+"\t"+"BY_cor"+"\n")
        for family in sorted_gene_result_sp2:
            geneID=family[0]
            plazaID=str(family[1][0])
            n_genes=str(family[1][1])
            pv_bin=str(family[1][2])
            pv_ran=str(family[1][3])
            Bonf_cor=str(family[1][4])
            BY_cor=str(family[1][5])

            line_to_print=geneID+"\t"+plazaID+"\t"+n_genes+"\t"+pv_bin+"\t"+pv_ran+"\t"+Bonf_cor+"\t"+BY_cor+"\n"

            out5.write(line_to_print)

    #lista geni significativi sp1
    with open ("./"+args.sample+"_results/"+args.sample+"_significant_genes_sp1.txt", "w") as out6:
        for gene in sig_com_gene_sp1:
            out6.write(gene+"\n")

    #lista geni significativi sp2
    with open ("./"+args.sample+"_results/"+args.sample+"_significant_genes_sp2.txt", "w") as out7:
        for gene in sig_com_gene_sp2:
            out7.write(gene+"\n")

################################################################################

if __name__ == "__main__":
    main()
