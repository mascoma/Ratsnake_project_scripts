for i in $(seq 1 1 1000)
do 

cp control_g2 control_g2_gene_$i
sed -i -e "s/g2_BStreeline_0"/"g2_BStreeline_$i/" control_g2_gene_$i

done
