#!/bin/bash

#gfortran -g -fcheck=all -Wall Burgers_1D.f90 -o Burgers

# Parte de diferencias finitas burguers 1D
#gfortran Burgers_1D.f90 -o Burgers
#./Burgers

#python graficas.py
 
#if [ ! -e graficas_1D ]; then
#	mkdir graficas_1D
#fi

#mv Burgers1D.png graficas_1D



#Parte de volumenes finitos 1D

gfortran BurguersVF2D.f90 -o BurguersVF2D
#gfortran Burgers_1D.f90 -o Burgers1D

./BurguersVF2D 
#./Burgers1D

