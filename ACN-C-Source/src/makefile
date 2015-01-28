# simple version 1.1
default: clean

gmb:
	gcc -o gmb -O3 gmb.c -I../cfitsio -L../cfitsio -lcfitsio -lm
gmf:
	gcc -o gmf gmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm
bmf:
	gcc :q
:U-boat-o bmf bmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm
nmf:
	gcc -o nmf nmf.c -I../cfitsio -L../cfitsio -lcfitsio -lm
compare:
	gcc -o compare compare.c -I../cfitsio -L../cfitsio -lcfitsio -lm
showdata:
	gcc -o showdata showdata.c -I:../cfitsio -L../cfitsio -lcfitsio -lm
cleanobjectfile:
	gcc -o cleanobjectfile cleanobjectfile.c -I../cfitsio -L../cfitsio -lcfitsio -lm

centroid:
	gcc -o centroid centroid.c -I../cfitsio -L../cfitsio -lcfitsio -lm
acn-aphot:
	gcc -o acn-aphot acn-aphot.c -I../cfitsio -L../cfitsio -lcfitsio -lm

listdir:
	gcc -o listdir listdir.c -lm -lnsl


	
gmbdouble:
	gcc -o gmb-double gmb-double.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
	

	
rrffast:
	gcc -o rrf -O3 rrf.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl


gmf-test:
	gcc -o gmf-test gmf-test.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl

subtractbias:
	gcc -o subtractbias subtractbias.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
reduction:
	gcc -o reduction reduction.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
pdra:
	gcc -o pdra pdra.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
masterbias:
	gcc -o masterbias masterbias.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl

imarith2:
	gcc -o imarith2 imarith2.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
imarith:
	gcc -o imarith imarith.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl
	
imlist:
	gcc -o imlist imlist.c -I../cfitsio -L../cfitsio -lcfitsio -lm -lnsl


listhead: 
	gcc -o listhead listhead.c -I../cfitsio -L../cfitsio -lcfitsio -lm
	
clean:
	rm masterbias
	rm mb.fits
	
