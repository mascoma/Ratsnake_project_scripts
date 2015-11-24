for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g2_clock_gene_$i.txt >> g2_clock.txt
sed -n "/lnL(ntime/p" g2_clock_gene_$i.txt >> g2_clock.txt

done