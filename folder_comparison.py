#!/usr/bin/python3
# listcomparison.py by Xin Chen
# 07/29/2015
# this is to compare and find out which files are missing in two folders with
# sequential files

import re
import os
def main():
    filelist1 = []
    filelist2 = []
    missingfiles = list()
    string = []
    path1 = os.getcwd() + "/group5"
    #print(path1 + "\n")
    filelist1 = os.listdir(path1) 
    print(len(filelist1))
    path2 = os.getcwd() + "/raxml_output/output_g5/results/besttree"
    #print(path2 + "\n")
    filelist2 = os.listdir(path2) 
    print(len(filelist2))
    for x in range(0, 404):
        flag = False
        string = "gene_" + str(x) + "."
        if any(string in s for s in filelist1) and not any(string in t for t in filelist2):
            flag = True
        if flag:
            missingfiles.append(string)
    print("missing files are:" + "\n")
    print(missingfiles)
if __name__ == "__main__": main() 
