#!/usr/bin/python3
# ExtractingLines.py by Xin Chen
# 07/16/2015 

def main():
    data=[]
    flag=False
    outfile = open('subrate8.txt', 'w')
    with open("mlb_8", 'r') as f:
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
if __name__ == "__main__": main()
