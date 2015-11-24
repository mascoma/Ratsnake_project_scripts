for i in $(seq 1 1 403); do
cat ./g6/g6_gene_$i.nex mbsets.txt >> ./g6/g6_inputgene_$i.nex
 
 
done 