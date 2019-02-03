import sys
import argparse
import re
import random
import os

def main():
    parser = argparse.ArgumentParser(description='Tool for random ID list generator')
    parser.add_argument('--all_sp1', metavar='all_sp1', help='All plazaID of species 1')
    parser.add_argument('--sp1', metavar='sp1', default="sp1", help='species 1 name (es. ath)')
    parser.add_argument('--all_sp2', metavar='all_sp2', help='All plazaID of species 2')
    parser.add_argument('--sp2', metavar='sp2', default="sp2", help='species 2 name (es. ath)')
    parser.add_argument('--random', metavar='random',type=int ,default="10000" ,help='number of list of random plazaID to generate (def=10000)')
    parser.add_argument('--small', metavar='small',type=int ,default="100" ,help='size of the smallest randomID list size (def=100)')
    parser.add_argument('--big', metavar='big',type=int ,default="1000" ,help='size of the biggest randomID list size (def=1000)')
    parser.add_argument('--inc', metavar='inc',type=int ,default="50" ,help='increment of size to go from small to big randomID list (def=50)')



    args = parser.parse_args()

################################################################################
    #liste degli input

    #genero una lista contente tutti i plazaID
    all_ID_sp1=[]

    with open (args.all_sp1) as all_plaza_sp1:
        for line in all_plaza_sp1:
            all_ID_sp1.append(line.strip("\n"))

    #genero una lista contente tutti i plazaID
    all_ID_sp2=[]

    with open (args.all_sp2) as all_plaza_sp2:
        for line in all_plaza_sp2:
            all_ID_sp2.append(line.strip("\n"))


    sys.stderr.write("total ID in species 1 = "+str(len(all_ID_sp1))+"\n")
    sys.stderr.write("total ID in species 2 = "+str(len(all_ID_sp2))+"\n")

################################################################################

    #funzione per creare cartelle se non ancora esistenti
    def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)

################################################################################


    #genero una lista contenent le grandezze delle liste di randomID che poi vado
    #a generare

    lists_sizes=[]

    smallest_list_size=args.small
    biggest_list_size=args.big

    while smallest_list_size<=biggest_list_size:
        lists_sizes.append(smallest_list_size)
        smallest_list_size=smallest_list_size+args.inc


    #ciclo che crea una nuova cartella per ogni grandezza di lista richiesta

    for size in lists_sizes:

        createFolder("./"+args.sp1+"_"+args.sp2+"/"+str(size)+"_"+args.sp1+"/")
        createFolder("./"+args.sp1+"_"+args.sp2+"/"+str(size)+"_"+args.sp2+"/")

        counter=1

        while counter<= args.random:

            #ciclo che genera un numero pari a random di liste della grandezza "size"

            random_ID_sp1=random.sample(all_ID_sp1,size)
            with open ("./"+args.sp1+"_"+args.sp2+"/"+str(size)+"_"+args.sp1+"/"+str(size)+"_random_"+args.sp1+"_ID_set"+str(counter)+".txt", "w") as out:
                for element in random_ID_sp1:
                    out.write(element+"\n")

            random_ID_sp2=random.sample(all_ID_sp2,size)
            with open ("./"+args.sp1+"_"+args.sp2+"/"+str(size)+"_"+args.sp2+"/"+str(size)+"_random_"+args.sp2+"_ID_set"+str(counter)+".txt", "w") as out:
                for element in random_ID_sp2:
                    out.write(element+"\n")

            counter=counter+1

        sys.stderr.write("generated "+str(args.random)+" list of "+str(size)+" random element"+"\n")

################################################################################


    #ciclo sulle liste generate per calcolare le npmi dei random ID
    #per annotarli uso una tabbela per ogn confronto posiibile (es 100sp1 vs 100sp2)

    






if __name__ == "__main__":
    main()
