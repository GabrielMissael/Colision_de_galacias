#Este archivo corresponde al Makefile del programa

#Escribiendo "make" en la terminal, se compila automaticamente, incluyendo los archivos que ocupa el proyecto
main: funciones.c main.c libs
	
	gcc funciones.c main.c -o main.out -lm -I ./libs
	mkdir Graficas
	mkdir Graficas/Graficas3D
	mkdir Graficas/GraficasXY
	mkdir Graficas/GraficasXZ
	mkdir puntos

#Al escribir el comando "make clean" se borra el archivo ejectuable y la carpeta con los archivos generados previamente.
#Se crea otra carpeta con el mismo nombre, para guardar futuros archivos
clean:
	
	rm main.out
	rm -r puntos
	rm -r Graficas
	

