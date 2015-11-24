for i in $(seq 1 1 403)
do
sed -n "/BASEML/p" g5_clock_gene_$i.txt >> g5_clock.txt
sed -n "/lnL(ntime/p" g5_clock_gene_$i.txt >> g5_clock.txt

done