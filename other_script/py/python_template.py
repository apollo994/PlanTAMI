import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='My nice tool.')
    #parser.add_argument('--input', metavar='INPUTFILE', default="/dev/stdin", help='The input file.')
    #parser.add_argument('--output', metavar='OUTPUTFILE', default="/dev/stdout", help='The output file.')
    #parser.add_argument('--fast', action='store_true', help="Fast parsing.")
    #parser.add_argument('--max', metavar='MAX', default="1024", type=int, help='Maximum.')
    #parser.add_argument('integers', metavar='N', type=int, nargs='+', help='an integer for the accumulator')
    args = parser.parse_args()

if __name__ == "__main__":
    main()
