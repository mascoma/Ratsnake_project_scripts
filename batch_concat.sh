for i in $(seq 1 1 57)
do 
	cp FASconCAT_v1.0.pl ./group$i
	cd group$i
	perl FASconCAT_v1.0.pl -s -i -p -p 
	mv FcC_smatrix.phy ../group$i.phy
	cd ..
done