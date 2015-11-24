#!/usr/bin/python3
# ExtractInfo.py by Xin Chen
# 07/24/2015


import re
import os
def main():
    outfile = open('geneinfo3.txt', 'w')
    print('locus\tac_num\tgene_name\tscore\te\tid\tgap\tnum_seq_hits\tqseq_coverage\tsseq_coverage\n')
    print('locus\tac_num\tgene_name\tscore\te\tid\tgap\tnum_seq_hits\n', file = outfile, end = '')
    for x in range(1, 404):
        data = []
        string = []
        tmp = []
        flag = True
        config_found = False
        endQ = 0
        endS = 0
        filename = "./blastoutput/result_L" + str(x) + ".out"
        print(filename)     
        try:
            f = open(filename, 'r')  
            print("successful open!")
            for line in f:
                if line.startswith('Length='):
                    tmp = re.search('([\d]+)',line)
                    #print(qlength)
                    qlength = tmp.group(1)
                    qlength = int(qlength)
                    staQ = qlength
                    break
            for line in f:
                if line.startswith('>'):
                    string = re.search('\|(.+?)\|(.+)\n', line)
                    acNum = string.group(1)
                    #print(acNum)
                    geneName = string.group(2)
                    #print(geneName)
                    flag = False
                    break
            if flag:
                print ("no hit has been found!")
                continue  
            for line in f:
                if line.find('Length=') == -1:
                    string = re.search('(.+)\n', line)
                    geneName = geneName + string.group(1)
                    #print(geneName)
                elif line.find('Length=') != -1:
                    tmp = re.search('([\d]+)',line)           
                    slength = tmp.group(1)
                    slength = int(slength)
                    staS = slength
                    #print(slength)
                    break
            for line in f:   
                if line.find('Score =') != -1:
                    #print("found it")
                    config_found = True    
                if config_found:
                    data.append(line)  
                    #print(data)  
                if line.find('Strand=') != -1:                   
                    config_found=False                      
                    break
            string = ' '.join(line.rstrip('\n') for line in data)   
            #print(string)                        
            tmp = re.search('Score.+?(\d+) bits', string)
            #print(tmp)
            score = tmp.group(1)
            #print(score)
            tmp = re.search('Expect = (.+?) Identities', string)  
            e = tmp.group(1)
            #print(e)
            tmp = re.search('Identities.+?\((\d+\%)\)', string)      
            identities = tmp.group(1)
            #print(identities)
            tmp = re.search('Gaps.+?\((\d+\%)\)', string)
            gaps = tmp.group(1)   
            #print(gaps)
            for line in f:                
                if line.startswith('Query '):
                    tmp = re.search('Query[\s]+([\d]+)[\s]+[\-\w]+[\s]+([\d]+)$', line)
                    tmp1 = tmp.group(1)
                    #print("tmp1=",tmp1)
                    tmp2 = tmp.group(2)
                    #print("tmp2=",tmp2)
                    if endQ < int(tmp2):
                        endQ = int(tmp2) 
                    if endQ > int(tmp2):
                        #print(staQ)
                        #print(endQ)
                        #print(staS)
                        #print(endS)
                        break
                    if staQ > int(tmp1):
                        staQ = int(tmp1)
                        #print("staNum",staNum)   
                if line.startswith('Sbjct '):
                    tmp = re.search('Sbjct[\s]+([\d]+)[\s]+[\-\w]+[\s]+([\d]+)$', line)
                    tmp1 = tmp.group(1)
                    #print("tmp1=",tmp1)
                    tmp2 = tmp.group(2)
                    #print("tmp2=",tmp2)
                    if endS < int(tmp2):
                        endS = int(tmp2) 
                    if staS > int(tmp1):
                        staS = int(tmp1) 
            f.close()              
            qseq_percent=((endQ-staQ+1)/qlength)*100   
            sseq_percent=((endS-staS+1)/slength)*100        
            f1 = open(filename)
            count = 0
            for line1 in f1:
                if acNum in line1:
                    count += 1
            f1.close()                                                                      
            print(x,'\t',acNum,'\t',geneName,'\t', score,'\t',e,'\t',identities,'\t',gaps,'\t',count/2,'\t',qseq_percent,"%",'\t',sseq_percent,"%",'\n',  file = outfile, end = '')

        except IOError:
            print ("no such file!")
if __name__ == "__main__": main() 