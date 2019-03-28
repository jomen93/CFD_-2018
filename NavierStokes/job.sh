clear 
k=1000
name1=Files
name2=Animaciones

if [ ! -e "${name1}" ]; then
	mkdir "${name1}"
	echo "No existe la carpeta ${name1}, se crea "
	sleep 1
	clear
else 
	rm "${name1}"/*.txt
	echo "Se borran los elementos antiguos en ${name1}"
	sleep 1
	clear
fi

gfortran -g -fcheck=all -Wall vort_stre_NE.f90 -o vort_stre_NE
echo "$k" | ./vort_stre_NE 
sleep 3


mv *.txt "${name1}"
echo "Datos listos para ser graficados"
sleep 1
clear
echo "Comienzo de las graficas"
echo "$k" | python graficas.py 

mv *.png "${name1}"
clear
sleep 1

echo "Se empiezan las animaciones "
clear 
sleep 1

echo "$k" | python animacion.py

if [ ! -e "${name2}" ]; then
	mkdir "${name2}"
	echo "No existe la carpeta de ${name2}, se crea "
	sleep 1
	clear
else 
	rm "${name2}"/*.png
	echo "Se borran los elementos antiguos en ${name2}"
	sleep 1
	clear
fi

mv *.gif "${name2}"
mv *.mp4 "${name2}"

echo "Listo!"









