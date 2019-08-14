/*
El presente proyecto consiste en simular la colisión de dos galaxias, usando para ello el metodo de Leap-frog. 
Cada galaxia se idealiza como un espacio esferico ocupado por n puntos distribuidos de forma aleatoria, con velocidades iniciales también aleatorias y masas iguales
*/

//Incluimos las librerias a ocupar, incluyendo la creada para este proyecto, "galaxia.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "galaxia.h"

//En este archivo se encuentra el cuerpo principal del programa

int main(){		
	
	//Declaramos las variables que ocuparemos de condición inicial
	//n y N son el numero de puntos de cada esfera, pos es el paso de tiempo en el que vamos
	int n, N, pos=0, gen, len;

	//r_1 y r_2 es el radio de cada esfera, R la distancia que las separa, vo la velocidad maxima para un punto dado
	//h el paso de tiempo y T el tiempo total a analizar
	double r_1, r_2, R, vo, h, T, rel_x, rel_y, pot, cin;

	//La cadena de caracteres "dir" es para guardar la dirección donde se guardaran los archivos de los puntos a cada paso de tiempo
	//arch es una cadena para, en caso de ser necesario, guardar el nombre del archivo con posiciones y velocidades inciciales
	char dir[300], arch[40];	

	//Leeremos las condiciones iniciales de este archivo
	FILE *lee;
	lee = fopen("inicio.txt", "r");

	//Leemos los datos
	fscanf(lee, "%i", &n);
	fscanf(lee, "%lf", &r_1);
	fscanf(lee, "%i", &N);
	fscanf(lee, "%lf", &r_2);
	fscanf(lee, "%lf", &vo);
	fscanf(lee, "%lf", &R);
	fscanf(lee, "%lf", &T);
	fscanf(lee, "%lf ", &h);
	fgets(dir, 300, lee);
	fscanf(lee, "%i \n", &gen);
	fgets(arch, 40, lee);
	fscanf(lee, "%lf", &rel_x);
	fscanf(lee, "%lf", &rel_y);

	//Cerramos el archivo de lectura
	fclose(lee);

	len = strlen(dir);
	int itera;
	itera = T/h;
	itera+=2;
		
	//Creamos y reservamos el espacio para 4 matrices, 2 para cada esfera (para condiciones actuales y anteriores), ademas de una matriz para la energia mecanica total
	double **es1, **es2, **nes1, **nes2, **energia;
	energia = (double **) malloc(itera*sizeof(double));
	es1 = (double **) malloc(n*sizeof(double));	
	es2 = (double **) malloc(N*sizeof(double));
	nes1 = (double **) malloc(n*sizeof(double));	
	nes2 = (double **) malloc(N*sizeof(double));

	//Guardaran n registros para cada punto, con 9 datos cada uno; posicion, velocidad y fuerza neta, en "x", "y" y "z"
	//Matrices para esfera 1
	for (int i=0 ; i<n ; i++){
		es1[i] = (double*) malloc (9*sizeof(double));
		nes1[i] = (double*) malloc (9*sizeof(double));
	}

	//Matrices para esfera2
	for (int i=0 ; i<N ; i++){
		es2[i] = (double*) malloc (9*sizeof(double));
		nes2[i] = (double*) malloc (9*sizeof(double));
	}

	//Matriz para energia mecanica
	for(int i=0 ; i<itera ; i++){
		energia[i] = (double*) malloc(4*sizeof(double));
		energia[i][0]=i+1;
	}

	//Se generan las esferas si así se indica
	if(gen==1){	
		es1 = esfera(es1, n, r_1, 0, vo, 0.0, 0.0);
		es2 = esfera(es2, N, r_2, R, vo, rel_x, rel_y);
	}

	//Se leen los datos iniciales para las esferas, de un archivo, si asi se indica
	if(gen==0){
		//El archivo debe tener un numero n, y luego n lineas correspondientes a la esfera 1, y luego un numero N, seguido de N lineas correspondientes a la esfera 2		
		FILE *galax;
		char *ex;
		int term;
		term = strlen(arch);
		term-=1;
		ex = (char *) malloc(term*sizeof(char));	
		for(int i=0 ; i<term ; i++){		
			ex[i] = arch[i];
		}
		galax = fopen(ex, "r");
		fscanf(galax, "%i", &n);
		for(int i=0 ; i<n ; i++){
			for(int j=0 ; j<6 ; j++){
				fscanf(galax, "%lf ", &es1[i][j]);
			}	
		}
		fscanf(galax, "%i", &N);
		for(int i=0 ; i<N ; i++){
			for(int j=0 ; j<6 ; j++){
				fscanf(galax, "%lf", &es2[i][j]);
			}	
		}
	}

	//Inicializamos las fuerzas para la esfera 1 y 2, usando la función fuerza()
	for(int i=0 ; i<n ; i++){	
		es1[i] = fuerza(es1, es2, 0, i, es1[i], n, N);
		
	}
	for(int i=0 ; i<N ; i++){
		es2[i] = fuerza(es1, es2, 1, i, es2[i], n, N);
	}

	//Creamos el primer archivo, con todos los datos iniciales, con la función archivo()
	
	archivo(es1, es2, n, N, pos, dir, len);
	pos++;	
	
	//En este ciclo se realiza la simulación, desde i hasta T, dando pasos de tamaño h	
	for(double i=0 ; i<=T; i+=h){	
		cin=0.0;
		pot=0.0;
		//El factor es la conversion para pasar de km/s a kPc/gAño, que son las unidades con las que medimos distancia y tiempo
		double factor = 1.0219053;
		//Primero, generamos las nuevas posiciones para cada punto de las dos esferas		
		for(int j=0 ; j<n ; j++){			
			cin+=(0.5)*pow((sqrt(es1[j][3]*es1[j][3]+es1[j][4]*es1[j][4]+es1[j][5]*es1[j][5])*factor), 2);		
			nes1[j][0] = es1[j][0] + h*es1[j][3]*factor + (0.5)*es1[j][6]*(h*h);
			nes1[j][1] = es1[j][1] + h*es1[j][4]*factor + (0.5)*es1[j][7]*(h*h);
			nes1[j][2] = es1[j][2] + h*es1[j][5]*factor + (0.5)*es1[j][8]*(h*h);
		}
		for(int j=0 ; j<N ; j++){
			cin+=(0.5)*pow((sqrt(es2[j][3]*es2[j][3]+es2[j][4]*es2[j][4]+es2[j][5]*es2[j][5])*factor), 2);		
			nes2[j][0] = es2[j][0] + h*es2[j][3]*factor + (0.5)*es2[j][6]*(h*h);
			nes2[j][1] = es2[j][1] + h*es2[j][4]*factor + (0.5)*es2[j][7]*(h*h);
			nes2[j][2] = es2[j][2] + h*es2[j][5]*factor + (0.5)*es2[j][8]*(h*h);	
		}	

		//Calculamos las nuevas fuerzas y velocidades para cada punto en las dos esferas
		for(int j=0 ; j<n ; j++){
			nes1[j] = fuerza(nes1, nes2, 0, j, nes1[j], n, N);
	
			nes1[j][3] = es1[j][3] + (0.5)*(es1[j][6]+nes1[j][6])*h;
			nes1[j][4] = es1[j][4] + (0.5)*(es1[j][7]+nes1[j][7])*h;
			nes1[j][5] = es1[j][5] + (0.5)*(es1[j][8]+nes1[j][8])*h;		
		}
		for(int j=0 ; j<N ; j++){			
			nes2[j] = fuerza(nes1, nes2, 1, j, nes2[j], n, N);

			nes2[j][3] = es2[j][3] + (0.5)*(es2[j][6]+nes2[j][6])*h;
			nes2[j][4] = es2[j][4] + (0.5)*(es2[j][7]+nes2[j][7])*h;
			nes2[j][5] = es2[j][5] + (0.5)*(es2[j][8]+nes2[j][8])*h;	
		}
		pot = potencial(es1, es2, n, N);
		//printf("ok %i\n", pos-1);
		energia[pos-1][1] = pot;
		energia[pos-1][2] = cin;
		energia[pos-1][3] = cin+pot;

		//Igualamos las esferas anteriores con las nuevas, para despues generar un archivo con esos datos
		es1 = nes1;
		es2 = nes2;
		archivo(es1, es2, n, N, pos, dir, len);

		//Avanzamos el contador
		pos++;
	}

	//creamos el archivp "size.txt" que contiene el numero de puntos por esfera para poder graficar los datos

	FILE *ini;
	char *extra;
	extra = dir;
	char name[8] = {'s','i','z','e','.','t','x', 't'};
	for(int i=0 ; i<8 ; i++){
		extra[len+i-1]=name[i];
	}
	ini = fopen(extra, "w");
	fprintf(ini, "%i\t%i\t%i\n", n, N, pos);
	fclose(ini);

	//Creamos el archivo "energia.txt" que contiene la energia potencial, cinetica y mecanica a cada paso de tiempo
	
	extra = dir;
	char name1[11] = {'e','n','e','r','g','i','a', '.', 't', 'x', 't'};
	for(int i=0 ; i<11 ; i++){
		extra[len+i-1]=name1[i];
	}
	ini = fopen(extra, "w");
	fprintf(ini, "Las columnas son: \nIteración \tEnergía Potencial\tEnergia Cínetica\tEnergia Mecanica total\n");
	for(int i=0 ; i<itera ; i++){
		fprintf(ini, "%0.0lf \t", energia[i][0]);		
		for(int j=1 ; j<4 ; j++){			
			fprintf(ini, "%0.1lf \t", energia[i][j]);
		}
		fprintf(ini, "\n");
	}
	fclose(ini);
	return 0;
}
