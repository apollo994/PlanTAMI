import sys
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='tool for mutiple test correction of a series of p_value')
    parser.add_argument('--input', metavar='INPUTFILE', help='Pv value formatted line by line')
    parser.add_argument('--output', metavar='OUTPUTFILE', default="adj_pv_table", help='Name of the output file')
    #parser.add_argument('--fast', action='store_true', help="Fast parsing.")
    #parser.add_argument('--max', metavar='MAX', default="1024", type=int, help='Maximum.')
    #parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
    args = parser.parse_args()

    ###########################################################################

    #importo la serie

    pv_test=[]

    with open (args.input) as my_file:
        for line in my_file:
            pv_test.append(float(line.strip("\n")))


    pv_test=sorted(pv_test)

    ###########################################################################

    #funzione per calcolare FDR con bonferroni
    def BON_corr(serie):
        #lista in cui annotare i risultari
        BON_res=[]
        #numero di test
        n=len(serie)
        #ciclo in cui moltiplico ogni pv per in numero di osservazioni
        for pv in serie:
            adj_pv=pv*n
            if adj_pv>1:
                adj_pv=1
                BON_res.append(adj_pv)
            else:
                BON_res.append(adj_pv)

        return(BON_res)


    ###########################################################################

    #funzione per calcolare FDR con metodo BH

    def BH_corr(serie):
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


    ###########################################################################

    #funzione per calcolare FDR con metodo BY

    def BY_corr(serie):
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

    ###########################################################################

    #caslcolo i pv adj con i diversi metodi
    BON_results=BON_corr(pv_test)
    BH_results=BH_corr(pv_test)
    BY_results=BY_corr(pv_test)

    ###########################################################################

    cicli=list(range(0,len(pv_test)))

    with open (args.output + ".txt", "w") as out:
        #scrivo l'header
        out.write("pv"+"\t"+"BON"+"\t"+"BH"+"\t"+"BY"+"\n")
        #ciclo per stampare sull'output i risultai
        for idx in cicli:
            out.write(str(pv_test[idx])+"\t")
            out.write(str(BON_results[idx])+"\t")
            out.write(str(BH_results[idx])+"\t")
            out.write(str(BY_results[idx])+"\n")


if __name__ == "__main__":
    main()
