#!/usr/bin/python3
# ExtractingLines.py by Xin Chen
# 07/16/2015 
import sys, getopt
def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print ("test.py -i <inputfile> -o <outputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ("test.py -i <inputfile> -o <outputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg 
    data=[]
    flag=False
    outfile = open(outputfile, 'w')
    with open(inputfile, 'r') as f:
        for line in f:
            if line.startswith('Substitution'):
                flag=True
            if flag:
                data.append(line)
            if line.startswith('Nodes and Times'):
                flag=False
    for line in data:  
        print(line, file = outfile, end = '')
    print('Done')                
if __name__ == "__main__": main(sys.argv[1:])
