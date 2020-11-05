/****************************************************************************
 * Tarea #1  - [ILN250] Gestión de información y operación                  *
 *																																					*
 * Este archivo es parte de la entrega del código de la seción de alogorit- *
 * mos de minimización programados para el problema SSWB.										*
  *																																					*
 ****************************************************************************/
 /**
 * @file tarea_1.c
 * @author Mauricio Aravena - Diego Badillo
 * @date 4 Nov 2020
 * @brief Main de la tarea #1, llama algoritmos y calcula el tiempo de ejecución
 *
 * Este archivo es el encargado tanto de generar los sets de valores para
 * entregar a los algoritmos llamados, como también contiene la lógica para
 * calcular el tiempo de ejecución de cada uno, asi como entregar los prints,
 * con esta información.
 *
 */

 /**
 * [Bibliotecas estandar importadas]
 */
#include <stdio.h>		/**< Estandar de I/O en C */
#include <stdlib.h>		/**< Estandar de C */
#include <time.h>			/**< Funciones relacionadas a la medición de tiempo */
#include <math.h>			/**< Biblioteca de matematicas */
#include <sys/time.h> /**< System Time */

/**
* [Biblioteca definida por usuario]
*/
#include "utils.h"		/**< Archivo con las declaraciones de las funciones algoritmo */


// #define EPSILON 0.001
// #define KMAX 30
// #define I_SIZE 1000
#define INSTANCES 2 // TEMPORAL
// #define ALPHA 0.3 // 0<ALPHA<0.5
// #define BETA 0.6 // 0<BETA<1


//DECLARACIONES
// Set* getNewSet(int size);
// void generateDarray(Set* ptrSet);
// void printValues(Set* ptrSet);
// void deleteArray(Set* ptrSet);
// void initPoint(Set* ptrSet);
// void weiszfeld(Set* ptrSet, double epsilon,  int kmax, double Zk1);
// void distance(Set* ptrSet, double* ptrDistances, double x0, double y0);
// double sum(int i, Set* ptrSet, double* distances);
// void newPoint(Set* ptrSet, double* ptrDistances);
// void Dfunction(Set* ptrSet, double x, double y, double* returnArray);
// double backtracking(Set* ptrSet, double* difArray, double* nabArray, double alpha, double beta);
// double function(Set* ptrSet, double x, double y);
// double multArray12x21(double* Array1, double* Array2);
// void multArray22x21(double arrayA[][2], double* arrayB, double* arrayReturn);
// void gradient(Set* ptrSet, double epsilon);
// void invHfunction(Set* ptrSet, double x, double y, double returnArray[][2]);
// void newton(Set* ptrSet, double epsilon);

// Variables globales

// instances, max, avg (se almacena max y avg en los ultimos elementos

/*
[0] - [INSTANCES - 1] Tiempos
[INSTANCES] - [INSTANCES +1] Datos temporales {max, avg}
[INSTANCES + 2] - [INSTANCES + 3] Datos de iteracion {max, avg}
[INSTANCES + 4] - [INSTANCES + 5] Datos de F.O. {max, avg}
*/
unsigned long avgArrayW[INSTANCES + 2]; //Weiszfeld
unsigned long avgArrayG[INSTANCES + 2]; //Gradient
unsigned long avgArrayN[INSTANCES + 2]; //Newton

double avgIteration[3]  = {0., 0., 0.};
int maxIteration[3]  = {0, 0, 0};
double avgFc[3] = {0., 0., 0.};
double maxFc[3] = {0., 0., 0.};
/*--------------------------------------------------------------*/
int main(int argc, char const *argv[])
{
	srand48(time(0));
	//INICIALIZACION EN CERO
	avgArrayW[INSTANCES + 1] = avgArrayG[INSTANCES + 1] = avgArrayN[INSTANCES + 1] = 0.;
	struct timeval stop, start;
	unsigned long timeTaken;
	unsigned long prev;
	int nIteration[3] = {0, 0, 0}; // N° iteraciones
	double Fc[3] = {0., 0., 0.};

	for(int i=0; i<INSTANCES; i++) {
		Set* test  = getNewSet(I_SIZE);
		//printValues(test);
		printf("\nInstancia nº %i \n\n", i);
		printf("Soluciones iniciales X0 = %f, Y0 = %f\n", test->x0, test->y0);
		for(int j=0; j<3; j++) { // para cada algoritmo
			initPoint(test);
			switch (j)
			{
				case 0:
					printf("Weiszfeld\n");
					gettimeofday(&start, NULL);
					nIteration[0] = weiszfeld(test, EPSILON, KMAX, 0);
					gettimeofday(&stop, NULL);
					printf("Soluciones finales X* = %f, Y* = %f\n", test->x0, test->y0);
					printf("Nº iteraciones: %i\n", nIteration[0]);
					break;
				case 1:
					printf("Gradiente\n");
					gettimeofday(&start, NULL);
					nIteration[1] = gradient(test, EPSILON);
					printf("Soluciones finales X* = %f, Y* = %f\n", test->x0, test->y0);
					gettimeofday(&stop, NULL);
					break;
				case 2:
					printf("Newton\n");
					gettimeofday(&start, NULL);
					nIteration[2] = newton(test, EPSILON, 1);
					printf("Soluciones finales X* = %f, Y* = %f\n", test->x0, test->y0);
					gettimeofday(&stop, NULL);
					break;
			}
			timeTaken = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;

			switch(j)
			{
				case 0:
					//Analisis temporal
					avgArrayW[i] = timeTaken;
					printf("tiempo: %lu [us]\n", timeTaken);
					printf("nIteration: %d\n", nIteration[0]);
					prev = (i == 0)?0:avgArrayW[i-1];
					if (timeTaken > prev)
						avgArrayW[INSTANCES] = timeTaken; //se guarda max
					avgArrayW[INSTANCES + 1] += timeTaken; //se guarda sum
					//Analisis iteracion
					if (nIteration[0] > maxIteration[0])
						maxIteration[0] = nIteration[0];
					avgIteration[0] += nIteration[0];
					//Funcion objetivo
					Fc[0] = function(test, test->x0, test->y0);
					printf("Funcion objetivo: %lf\n", Fc[0]);
					if (Fc[0] > maxFc[0])
						maxFc[0] = Fc[0];
					avgFc[0] += Fc[0];
					break;
				case 1:
					//Analisis temporal
					avgArrayG[i] = timeTaken;
					printf("tiempo: %lu [us]\n", timeTaken);
					printf("nIteration: %d\n", nIteration[1]);
					prev = (i == 0)?0:avgArrayG[i-1];
					if (timeTaken > prev)
						avgArrayG[INSTANCES] = timeTaken; //se guarda max
					avgArrayG[INSTANCES + 1] += timeTaken; //se guarda sum
					//Analisis iteracion
					if (nIteration[1] > maxIteration[1])
						maxIteration[1] = nIteration[1];
					avgIteration[1] += nIteration[1];
					//Funcion objetivo
					Fc[1] = function(test, test->x0, test->y0);
					printf("Funcion objetivo: %lf\n", Fc[1]);
					if (Fc[1] > maxFc[1])
						maxFc[1] = Fc[1];
					avgFc[1] += Fc[1];
					break;

				case 2:
					//Analisis temporal
					avgArrayN[i] = timeTaken;
					printf("tiempo: %lu [us]\n", timeTaken);
					printf("nIteration: %d\n", nIteration[2]);
					prev = (i == 0)?0:avgArrayN[i-1];

					if (timeTaken > prev)
						avgArrayN[INSTANCES] = timeTaken; //se guarda max
					avgArrayN[INSTANCES + 1] += timeTaken; //se guarda sum
					//Analisis iteracion
					if (nIteration[2] > maxIteration[2])
						maxIteration[2] = nIteration[2];
					avgIteration[2] += nIteration[2];
					//Funcion objetivo
					Fc[2] = function(test, test->x0, test->y0);
					printf("Funcion objetivo: %lf\n", Fc[0]);
					if (Fc[2] > maxFc[2])
						maxFc[2] = Fc[2];
					avgFc[2] += Fc[2];
					break;
			}
		}
		deleteArray(test);
	}
	//calculo promedios
	avgArrayW[INSTANCES + 1] /= INSTANCES;
	avgArrayG[INSTANCES + 1] /= INSTANCES;
	avgArrayN[INSTANCES + 1] /= INSTANCES;

	avgIteration[0] /= INSTANCES;
	avgIteration[1] /= INSTANCES;
	avgIteration[2] /= INSTANCES;

	avgFc[0] /= INSTANCES;
	avgFc[1] /= INSTANCES;
	avgFc[2] /= INSTANCES;

	printf("\n\n -------- MAXIMOS Y PROMEDIOS --------\n\n");
	printf("Weiszfeld\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", avgArrayW[INSTANCES], avgArrayW[INSTANCES + 1] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[0], avgIteration[0]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[0], avgFc[0]);


	printf("Gradiente\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", avgArrayG[INSTANCES], avgArrayG[INSTANCES + 1] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[1], avgIteration[1]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[1], avgFc[1]);


	printf("Newton\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", avgArrayN[INSTANCES], avgArrayN[INSTANCES + 1] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[2], avgIteration[2]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[2], avgFc[2]);

	return 0;
}
/*----------------------------------------------------------------------------*/
/*funciones*/
// Set* getNewSet(int size)
// {
//     Set* newSet = NULL;
//     if(size > 0) {
//         newSet = (Set*)malloc(sizeof(Set));
//     }
//     newSet->size = size;
//     generateDarray(newSet);
// 	initPoint(newSet);
//     return newSet;
//
// }
//
// void generateDarray(Set* ptrSet)
// {
// 	int size = ptrSet->size;
// 	double* arrTempX = malloc(sizeof(double)*size);
// 	double* arrTempY = malloc(sizeof(double)*size);
// 	double* arrTempW = malloc(sizeof(double)*size);
//
// 	for( int i = 0; i < size; i++){
// 		arrTempX[i] = lrand48() % (size+1);
// 		arrTempY[i] = lrand48() % (size+1);
// 		//The lrand48() function returns non-negative, long integers, uniformly distributed over the interval [0,RAND_MAX].
//         int lower = size/10;
//         int upper = lower*5;
// 		arrTempW[i] = ((double)lrand48()/(double)(RAND_MAX))*(upper-lower) + lower;
// 	}
// 	ptrSet->ArrayX = arrTempX;
// 	ptrSet->ArrayY = arrTempY;
// 	ptrSet->ArrayW = arrTempW;
// 	initPoint(ptrSet);
// }
//
// void printValues(Set* ptrSet) {
// 	int size = ptrSet->size;
// 	for(int i = 0; i < size; i++){
// 		printf("[%i] Valores Xi = %f, Yi = %f, Wi = %f\n", (i+1),ptrSet->ArrayX[i],ptrSet->ArrayY[i], ptrSet->ArrayW[i]);
// 	}
// }
//
// void deleteArray(Set* ptrSet)
// {
// 	free(ptrSet->ArrayX);
// 	free(ptrSet->ArrayY);
// 	free(ptrSet->ArrayW);
// 	free(ptrSet);
// }
// void initPoint(Set* ptrSet)
// {
// 	int size = ptrSet->size;
// 	double tempSumWX = 0;
// 	double tempSumWY = 0;
// 	double tempSumW = 0;
// 	double tempX  = 0;
// 	double tempY  = 0;
// 	double tempW = 0;
// 	for(int i = 0; i < size; i++){
// 		tempX  = ptrSet->ArrayX[i];
// 		tempY  = ptrSet->ArrayY[i];
// 		tempW = ptrSet->ArrayW[i];
// 		tempSumWX += tempW*tempX;
// 		tempSumWY += tempW*tempY;
// 		tempSumW  += tempW;
// 	}
// 	ptrSet->x0 = tempSumWX/tempSumW;
//     ptrSet->y0 = tempSumWY/tempSumW;
// }
//
//
// void weiszfeld(Set* ptrSet, double epsilon, int kmax, double Zk1)
// {
// 	double* distances = malloc(sizeof(double)*(ptrSet->size));
// 	distance(ptrSet, distances, ptrSet->x0, ptrSet->y0);
// 	double Zk = sum(ptrSet->size, ptrSet, distances);
// 	//Esto se ejecuta la primera vez
// 	newPoint(ptrSet, distances); // Seteo de nuevo punto inicial
// 	free(distances);
// 	if ( (kmax <= 0) || fabs(Zk - Zk1) < epsilon ){
// 		//Terminar y guardar
// 		//printf("Weiszfeld listo! \n");
// 	}
// 	else{
// 		//printf("yip yip! \n");
// 		weiszfeld(ptrSet, epsilon, kmax-1, Zk);
// 	}
// }
//
// void distance(Set* ptrSet, double* ptrDistances, double x0, double y0)
// {
// 	 double* Xs = ptrSet->ArrayX;
// 	 double* Ys = ptrSet->ArrayY;
// 	for(int i = 0; i < ptrSet->size; i++) {
// 		ptrDistances[i] = sqrt(pow((double)(Xs[i]-x0), (double)(2)) + pow((double)(Ys[i]-y0), (double)(2)));
// 	}
// }
//
// double sum(int i, Set* ptrSet, double* distances)
// {
// 	if(i == 1) {
// 		return ( (ptrSet->ArrayW[0]) * distances[0] );
// 	}
// 	return sum(i-1, ptrSet, distances) + (ptrSet->ArrayW[i-1]) * distances[i-1];
// }
//
// void newPoint(Set* ptrSet, double* ptrDistances)
// {
// 	int size = ptrSet->size;
// 	double tempSumWX = 0;
// 	double tempSumWY = 0;
// 	double tempSumW = 0;
// 	double tempX  = 0;
// 	double tempY  = 0;
// 	double tempD = 0;
// 	double tempW = 0;
// 	for(int i = 0; i < size; i++){
// 		tempX  = ptrSet->ArrayX[i];
// 		tempY  = ptrSet->ArrayY[i];
// 		tempD  = ptrDistances[i];
// 		tempW  = ptrSet->ArrayW[i];
// 		tempSumWX += tempW*tempX/tempD;
// 		tempSumWY += tempW*tempY/tempD;
// 		tempSumW  += tempW/tempD;
// 	}
// 	ptrSet->x0 = tempSumWX/tempSumW;
//     ptrSet->y0 = tempSumWY/tempSumW;
// }
//
// /*----------------------------------------------------------*/
//
// void gradient(Set* ptrSet, double epsilon)
// {
// 	//calculo deltax
// 	double difArray[2];
// 	Dfunction(ptrSet, ptrSet->x0, ptrSet->y0, difArray);
// 	double nabArray[2];
// 	for (int i=0; i<2; i++) {
// 		nabArray[i] = difArray[i];
// 	}
// 	for(int i=0; i<2; i++) {
// 		difArray[i] = -difArray[i];
// 	}
//
// 	double t = backtracking(ptrSet, difArray, nabArray, ALPHA, BETA);
// 	ptrSet->x0 = ptrSet-> x0 + (t*difArray[0]);
// 	ptrSet->y0 = ptrSet-> y0 + (t*difArray[1]);
// 	double norma = pow(nabArray[0], 2) + pow(nabArray[1],2);
// 	if(norma >= pow(epsilon,2)  && t >= 0.0000000005){
// 		//printf("%lf \n", norma);
// 		gradient(ptrSet, epsilon);
// 	}
//
// }
//
//
// double backtracking(Set* ptrSet, double* difArray, double* nabArray, double alpha, double beta)
// {
// 	double xk = ptrSet->x0;
// 	double yk = ptrSet->y0;
// 	double t = 1;
// 	double a = function(ptrSet, xk + t*difArray[0], yk + t*difArray[1]);
// 	double b = function(ptrSet, xk, yk) + alpha*t*multArray12x21(nabArray, difArray);
// 	while (a > b) {
// 		t = beta*t;
// 		a = function(ptrSet, xk + t*difArray[0], yk + t*difArray[1]);
// 		b = function(ptrSet, xk, yk) + alpha*t*multArray12x21(nabArray, difArray);
// 	}
// 	//printf("t = %.10lf\n", t);
// 	return t;
// }
//
// double function(Set* ptrSet, double x, double y)
// {
// 	double* ArrayX = ptrSet->ArrayX;
// 	double* ArrayY = ptrSet->ArrayY;
// 	double* ArrayW = ptrSet->ArrayW;
// 	double retorno = 0;
// 	for(int i = 0; i < ptrSet->size; i++){
// 		retorno += ArrayW[i]*sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));
// 	}
// 	return retorno;
// }
//
// void Dfunction(Set* ptrSet, double x, double y, double* returnArray)
// {
// 	double* ArrayX = ptrSet->ArrayX;
// 	double* ArrayY = ptrSet->ArrayY;
// 	double* ArrayW = ptrSet->ArrayW;
//
// 	returnArray[0] = 0;
// 	returnArray[1] = 0;
//
// 	for(int i = 0; i < ptrSet->size; i++){
// 		returnArray[0] += ArrayW[i]*(1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2)))*(x - ArrayX[i]);
// 		returnArray[1] += ArrayW[i]*(1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2)))*(y - ArrayY[i]);
// 	}
// }
//
// void invHfunction(Set* ptrSet, double x, double y, double returnArray[][2])
// {
// 	double* ArrayX = ptrSet->ArrayX;
// 	double* ArrayY = ptrSet->ArrayY;
// 	double* ArrayW = ptrSet->ArrayW;
// 	returnArray[0][0] = 0;
// 	returnArray[0][1] = 0;
// 	returnArray[1][0] = 0;
// 	returnArray[1][1] = 0;
// 	// Calculo H
// 	for(int i = 0; i < ptrSet->size; i++) {
// 		double d = 1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));
// 		double d2 = 1/(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));
//
// 		returnArray[0][0] += ArrayW[i]*d*(1-pow(ArrayX[i] - x, (double) 2)*d2);
// 		returnArray[1][1] += ArrayW[i]*d*(1-pow(ArrayY[i] - y, (double) 2)*d2);
// 		double temp = -ArrayW[i]*(ArrayX[i] - x)*(ArrayY[i] - y)*pow(d, (double) 3);
// 		returnArray[0][1] += temp;
// 		returnArray[1][0] += temp;
// 	}
// 	// Calculo invH
// 	double	a = returnArray[0][0],
// 			b = returnArray[0][1],
// 			c = returnArray[1][0],
// 			d = returnArray[1][1];
// 	double det = 1/(a*d - b*c);
// 	returnArray[0][0] = d*det;
// 	returnArray[1][1] = a*det;
// 	returnArray[0][1] = -b*det;
// 	returnArray[1][0] = -d*det;
// }
//
// double multArray12x21(double* Array1, double* Array2)
// {
// 	return (Array1[0]*Array2[0] + Array1[1]*Array2[1]);
// }
//
// void multArray22x21(double arrayA[][2], double* arrayB, double* arrayReturn)
// {
// 	arrayReturn[0] = arrayA[0][0]*arrayB[0] + arrayA[0][1]*arrayB[1];
// 	arrayReturn[1] = arrayA[1][0]*arrayB[0] + arrayA[1][1]*arrayB[1];
// }
//
// void newton(Set* ptrSet, double epsilon)
// {
// 	//calculo deltax
// 	double difArray[2];
// 	double nabArray[2];
// 	double HArray[2][2];
// 	Dfunction(ptrSet, ptrSet->x0, ptrSet->y0, nabArray);
// 	invHfunction(ptrSet, ptrSet->x0, ptrSet->y0, HArray);
// 	multArray22x21(HArray, nabArray, difArray);
//
// 	double tempArray[2];
// 	for (int i=0; i<2; i++) {
// 		tempArray[i] = difArray[i];
// 	}
//
// 	for(int i=0; i<2; i++) {
// 		difArray[i] = -difArray[i];
// 	}
// 	//calculo lambda^2
// 	double lambda2 = multArray12x21(nabArray, tempArray);
// 	double t = 1;
// 	//criterio de parada
// 	if(lambda2/2 >= epsilon  && t >= 0.0000000005){
// 		//printf("%lf \n", lambda2);
// 		t = backtracking(ptrSet, difArray, nabArray, ALPHA, BETA);
// 		//actualizar x e y
// 		ptrSet->x0 = ptrSet-> x0 + (t*difArray[0]);
// 		ptrSet->y0 = ptrSet-> y0 + (t*difArray[1]);
// 		gradient(ptrSet, epsilon);
// 	}
// }
//
// // https://www.wolframalpha.com/input/?i=derivative+w*%28%28a+-+x%29%5E2++%2B+%28b+-+y%29%5E2%29%5E%28-0.5%29+*+%28x-a%29
// // falta: evaluacion de fc objetivo en cada instancia
// // imprimir numero de iteraciones para cada instancia y guardarlos para calcular promedio y maximo
// Como iteracion es distinta, no tiene mucho sentido comparar el promedio entre las iteraciones, por lo que haremos será comparar en una iteración la función objetivo converge a lo mismo, ya que se está comparando los tres algoritmos para el mismo set. De todas manera, por cumplir los requirimientos será de todas maneras comparar el máximo y promedio ALMACENAR
// Comentar diferencias entre eficiencia de iteraciones dy de tiempo
