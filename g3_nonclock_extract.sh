for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g3_nonclock_gene_$i.txt >> g3_nonclock.txt
sed -n "/lnL(ntime/p" g3_nonclock_gene_$i.txt >> g3_nonclock.txt

done