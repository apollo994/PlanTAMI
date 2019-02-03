import sys
import argparse
from collections import defaultdict
import random


def main():
    parser = argparse.ArgumentParser(description='tool for extract random element from list')
    parser.add_argument('--list', metavar='list', help='list')
    parser.add_argument('--random', metavar='random', default="100", type=int, help='number of random element to extract(def=100)')

    args = parser.parse_args()

################################################################################

    starting_list=[]

    with open (args.list) as my_list:
        for line in my_list:
            starting_list.append(line.strip("\n"))

    random_list=random.sample(starting_list, args.random)

    sample_name=args.list.split("/")
    sample_name=sample_name[len(sample_name)-1]

    
    with open (str(args.random)+"_random_from_"+sample_name, "w") as out:
        for line in random_list:
            out.write(line+"\n")

if __name__ == "__main__":
    main()
