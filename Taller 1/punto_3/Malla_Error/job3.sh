#!/bin/bash

gfortran punto3_RK2.f90
for j in {21..61}
 do 
 	echo "$j" "nodo$j" | ./a.out
 done
 