#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "galaxia.h"

//En este archivo se encuentra definidas todas las funciones que utiliza el programa

//Esta función imprime la matriz en la terminal
void imprime(double **aux, int n){
	for(int i=0 ; i<n ; i++){
		for(int j=0 ; j<9 ; j++){
			printf("%0.3lf \t", aux[i][j]);
		}	
		printf("\n");
	}
}

//Esta función crea números aleatorios postivos o negativos, con 5 decimales, no mas grandes o mas pequeños que un valor M
double crea(double M){	
	//Si pos es par, num será negativo	
	int pos = rand();
	pos%=2;

	//Generamos num y lo limitamos a +-M, y le genereamos 5 decimales
	double num = rand();
	int rest = num/(M*10000);
	num -= rest*(M*10000);
	num /=10000;
	if(pos==0){
		return num*-1;
	}
	if(pos==1){
		return num;	
	}
}

//Esta función crea esferas de radio r, n puntos y con un desplazamiento "t" en el eje "x", con velocidad maxima "v" para cada punto
double **esfera(double **aux, int n, double r, double t, double v, double rel_x, double rel_y){
	double count =0;
	
	//Realizamos n iteraciones	
	for (int i=0 ; i<n ; i++){
		for(int j=0 ; j<6 ; j++){
			//para las posiciones, se crean numeros aleatorios x tal que -r<x<r			
			if(j<3){		
				aux[i][j] = crea(r);
			}

			//para las velocidades, se crean numeros aleatorios x tal que -v<x<v	
			if(j>=3){
				aux[i][j] = crea(v);	
			}

			//Se revisa si el punto esta dentro de la esfera, de caso contrario, se regresa a j=0 para generar otro
			if(j==3){				
				count = sqrt(aux[i][0]*aux[i][0]+aux[i][1]*aux[i][1]+aux[i][2]*aux[i][2]);
				if(count >r){
					j=0;				
				}
			}
		}
		
		//Desplazamos el punto generado una distancia "t" en el eje x
		aux[i][0]+=t;
		//Ponemos la velocidad relativa en x y en y
		aux[i][3]+=rel_x;
		aux[i][4]+=rel_y;		
	}
	return aux;
}

//Esta función guarda los puntos de ambas esferas en un archivo
void archivo(double **es1, double **es2, int n, int N, int pos, char *supp, int len){	
	//A lo mas, se generarán 99,999 archivos, y se guarda unicamente cada 10 iteraciones
	if(pos>99999 || pos%10!=0){
		return;
	}	
	//Generamos el nombre del archivo, que tendra la forma "numero.txt", donde numero es el numero de iteración que representa
	char name[9] = {'0','0', '0', '0', '0', '.', 't', 'x', 't'};
	name[0]=(pos/100000)+'0';
	pos%=100000;	
	name[1]=(pos/10000)+'0';
	pos%=10000;	
	name[2]=(pos/1000)+'0';
	pos%=1000;
	name[3]=(pos/100)+'0';
	pos%=100;
	name[4	]=(pos/10)+'0';

	//Copiamos el nombre del archivo a la cadena supp, que contiene previamente la dirección donde se guardará el archivo
	for(int i=0 ; i<9 ; i++){
		supp[i+len-1]=name[i];
	}

	//Creamos el archivo de escritura, en la dirección indicada con el nombre que se ha generado
	FILE *crea;
	crea = fopen(supp, "w");

	fprintf(crea, "## Las primeras %i lineas corresponden a la esfera 1 \n##las siguienets %i lineas a la esfera 2\n", n, N);
	fprintf(crea, "## Las primeras tres columnas corresponden a la posición en 'x', 'y' y 'z', dadas en kpc, las siguientes tres a las velocidades, dadas en km/s\n");
	
	//Imprimimos n lineas correspondientes a la esfera 1 con 6 datos (posiciones y velocidades en x, y y z)
	for(int i=0 ; i<n ; i++){
		for(int j=0 ; j<6 ; j++){
			fprintf(crea, "%0.3lf \t", es1[i][j]);
		}	
		fprintf(crea, "\n");
	}

	//Imprimimos N lineas correspondientes a la esfera 2
	for(int i=0 ; i<N ; i++){
		for(int j=0 ; j<6 ; j++){
			fprintf(crea, "%0.3lf \t", es2[i][j]);
		}	
		fprintf(crea, "\n");
	}
	
	//Cerramos el archivo de escritura
	fclose(crea);
	return;
}

//Esta función obtiene la fuerza neta actuante en "x", "y" y "z" sobre una particula
double *fuerza(double **es1, double **es2, int mat, int pos, double *pun, int n, int N){
	
	double x, y, z, G, r;
	//G es la constante gravitacional
	G = -4300;

	//Para que la aceleración quede en kPc/gAño²
	double factor = 1.0443;
	
	//Inicializamos las fuerzas en 0
	pun[6]=0;
	pun[7]=0;
	pun[8]=0;

	//Guardamos la posicion de la particula
	x = pun[0];
	y = pun[1];
	z = pun[2];

	//Recorremos las n particulas de la esfera 1
	for(int i=0 ; i<n ; i++){
		//Excluimos a la particula del conteo		
		if(i != pos || mat==1){
			
			//Obtenemos la distancia entre el punto analizado y el punto i de la esfera 1
			r = sqrt((x-es1[i][0])*(x-es1[i][0])+(y-es1[i][1])*(y-es1[i][1])+(z-es1[i][2])*(z-es1[i][2]));	

			//Obtenemos la fuerza que genera el punto i en el punto analizado		
			pun[6] += factor*(G*(x-es1[i][0]))/(r*r*r);
			pun[7] += factor*(G*(y-es1[i][1]))/(r*r*r);
			pun[8] += factor*(G*(z-es1[i][2]))/(r*r*r);
		}	
	}

	//Realizamos lo mismo, pero con los puntos de la esfera 2
	for(int i=0 ; i<N ; i++){
		if(i != pos || mat==0){			
			r = sqrt((x-es2[i][0])*(x-es2[i][0])+(y-es2[i][1])*(y-es2[i][1])+(z-es2[i][2])*(z-es2[i][2]));			
			pun[6] += factor*(G*(x-es2[i][0]))/(r*r*r);
			pun[7] += factor*(G*(y-es2[i][1]))/(r*r*r);
			pun[8] += factor*(G*(z-es2[i][2]))/(r*r*r);	
		}
	}

	//Regresamos el resultado obtenido
	return pun;
}

//Esta función obtiene la energia potencial total del sistema 
double potencial(double **es1, double **es2, int n, int N){
	double total=0, r, G=4300;
	double factor = 1.0443;
	for(int i=0 ; i<n ; i++){
		for(int j=i+1 ; j<n ; j++){
			if(i!=j){			
				r = sqrt((es1[i][0]-es1[j][0])*(es1[i][0]-es1[j][0])+(es1[i][1]-es1[j][1])*(es1[i][1]-es1[j][1])+(es1[i][2]-es1[j][2])*(es1[i][2]-es1[j][2]));
				if(r>0){			
					total += (factor*G/r);
				}	
			}	
		}	
		for(int j=0 ; j<N ; j++){			
			r = sqrt((es1[i][0]-es2[j][0])*(es1[i][0]-es2[j][0])+(es1[i][1]-es2[j][1])*(es1[i][1]-es2[j][1])+(es1[i][2]-es2[j][2])*(es1[i][2]-es2[j][2]));
			if(r>0){			
				total += (factor*G/r);
			}			
		}
			
	}
	for(int i=0 ; i<N ; i++){
		for(int j=i+1 ; j<N ; j++){			
			if(j!=i){				
				r = sqrt((es2[i][0]-es2[j][0])*(es2[i][0]-es2[j][0])+(es2[i][1]-es2[j][1])*(es2[i][1]-es2[j][1])+(es2[i][2]-es2[j][2])*(es2[i][2]-es2[j][2]));
				if(r>0){			
					total += (factor*G/r);
				}	
			}	
		}		
	}
	return total;
}

