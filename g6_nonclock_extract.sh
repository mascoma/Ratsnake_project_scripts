for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g6_nonclock_gene_$i.txt >> g6_nonclock.txt
sed -n "/lnL(ntime/p" g6_nonclock_gene_$i.txt >> g6_nonclock.txt

done