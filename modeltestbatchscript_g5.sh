for i in $(seq 1 1 403); do

java -jar jModelTest.jar -d group5/group5.txt.gene_$i.phy -s 3 -i -g 4 -f -BIC -tr 7 -o group5out/g5outgene_$i.txt

done
