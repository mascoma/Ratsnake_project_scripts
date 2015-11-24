for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g3_clock_gene_$i.txt >> g3_clock.txt
sed -n "/lnL(ntime/p" g3_clock_gene_$i.txt >> g3_clock.txt

done