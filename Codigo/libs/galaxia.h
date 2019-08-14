#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Este archivo corresponde a la libreria creada para este proyecto


//Se incluyen las funciones usadas, mismas que estan definidas en el archivo "funciones.c"
void imprime(double **aux, int n);
double **esfera(double **aux, int n, double r, double t, double v, double rel_x, double rel_y);
double limit(double aux, double r);
double crea(double M);
void archivo(double **es1, double **es2, int n, int N, int pos, char *dir, int len);
double *fuerza(double **es1, double **es2, int mat, int pos, double *pun, int n, int N);
double potencial(double **es1, double **es2, int n, int N);

