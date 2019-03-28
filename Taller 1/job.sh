#!/bin/bash
#gfortran -g -fcheck=all -Wall tuprograma.f90


gfortran -g -fcheck=all -Wall punto3_RK2.f90
#for j in {31,41,51,61,71,81,91,101,111,121,131,141,151,161,171,181,191,201,211,221,231,241,251,261,271,281,291,301,311,321}
for j in {21..23}
do
 	echo "$j" | ./a.out >> nodo"$j".dat
 	#echo "$i" | ./FEesfera
	#echo "$j $i" | ./a.out 
done

 #python error2pto.py >> error.dat
#done 





##mkdir punto_3
#gfortran -g -fcheck=all -Wall Punto3_euler.f90
#./a.out
#echo "$1" | python Graficas_punto2.py
#mkdir punto_3/Euler_Forward
#mv punto3_euler.eps punto_3/Euler_Forward
#mv DCT_CDS.dat DCT_CDS_t.dat DCT_CDS_NU.dat malla_CDS.dat punto_3/Euler_Forward

#gfortran -g -fcheck=all -Wall punto3_RK2.f90
#./a.out
#mkdir punto_3/RK2
#echo "$2" | python Graficas_punto2.py
#mv punto3_RK2.eps punto_3/RK2
#mv DCT_RK2_t.dat DCT_RK2.dat DCT_RK2_NU.dat malla_RK2.dat punto_3/RK2



#for i in {1,2,3}
 #do 
 #	echo "$i" | ./a.out
 #done
 




# esta bien para difusion conveccion 
#gfortran difusion_conveccion.f90

#for j in {31,41,51,61,71,81,91,101,111,121,131,141,151,161,171,181,191,201,211,221,231,241,251,261,271,281,291,301,311,321}
#do
 #for i in {1,3}
 #do 
	#echo "$i" | ./FEesfera
#	echo "$j $i" | ./a.out 
 #done
 #python error2pto.py >> error.dat
#done 

