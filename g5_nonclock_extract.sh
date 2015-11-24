for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g5_nonclock_gene_$i.txt >> g5_nonclock.txt
sed -n "/lnL(ntime/p" g5_nonclock_gene_$i.txt >> g5_nonclock.txt

done