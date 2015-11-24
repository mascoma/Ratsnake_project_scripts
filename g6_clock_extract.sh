for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g6_clock_gene_$i.txt >> g6_clock.txt
sed -n "/lnL(ntime/p" g6_clock_gene_$i.txt >> g6_clock.txt

done