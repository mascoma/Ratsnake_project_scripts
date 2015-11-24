for i in $(seq 1 1 403); do

java -jar jModelTest.jar -d group2/group2.txt.gene_$i.phy -s 3 -i -g 4 -f -BIC -tr 7 -o group2out/g2outgene_$i.txt

done
