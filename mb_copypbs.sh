for a in $(seq 2 1 403); do
cp mb_gene1.pbs mb_gene$a.pbs
sed -i -e "s/mb_gene1"/"mb_gene$a/" mb_gene$a.pbs
sed -i -e "s/inputgene_1"/"inputgene_$a/" mb_gene$a.pbs
sed -i -e "s/gene_1"/"gene_$a/" mb_gene$a.pbs
done