#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"
/*----------------------------------------------------------------------------*/
/*funciones*/
Set* getNewSet(int size)
{
    Set* newSet = NULL;
    if(size > 0) {
        newSet = (Set*)malloc(sizeof(Set));
    }
    newSet->size = size;
    generateDarray(newSet);
	initPoint(newSet);
    return newSet;

}

void generateDarray(Set* ptrSet)
{
	int size = ptrSet->size;
	double* arrTempX = malloc(sizeof(double)*size);
	double* arrTempY = malloc(sizeof(double)*size);
	double* arrTempW = malloc(sizeof(double)*size);

	for( int i = 0; i < size; i++){
		arrTempX[i] = lrand48() % (size+1);
		arrTempY[i] = lrand48() % (size+1);
		//The lrand48() function returns non-negative, long integers, uniformly distributed over the interval [0,RAND_MAX].
        int lower = size/10;
        int upper = lower*5;
		arrTempW[i] = ((double)lrand48()/(double)(RAND_MAX))*(upper-lower) + lower;
	}
	ptrSet->ArrayX = arrTempX;
	ptrSet->ArrayY = arrTempY;
	ptrSet->ArrayW = arrTempW;
	initPoint(ptrSet);
}

void printValues(Set* ptrSet) {
	int size = ptrSet->size;
	for(int i = 0; i < size; i++){
		printf("[%i] Valores Xi = %f, Yi = %f, Wi = %f\n", (i+1),ptrSet->ArrayX[i],ptrSet->ArrayY[i], ptrSet->ArrayW[i]);
	}
}

void deleteArray(Set* ptrSet)
{
	free(ptrSet->ArrayX);
	free(ptrSet->ArrayY);
	free(ptrSet->ArrayW);
	free(ptrSet);
}
void initPoint(Set* ptrSet)
{
	int size = ptrSet->size;
	double tempSumWX = 0;
	double tempSumWY = 0;
	double tempSumW = 0;
	double tempX  = 0;
	double tempY  = 0;
	double tempW = 0;
	for(int i = 0; i < size; i++){
		tempX  = ptrSet->ArrayX[i];
		tempY  = ptrSet->ArrayY[i];
		tempW = ptrSet->ArrayW[i];
		tempSumWX += tempW*tempX;
		tempSumWY += tempW*tempY;
		tempSumW  += tempW;
	}
	ptrSet->x0 = tempSumWX/tempSumW;
    ptrSet->y0 = tempSumWY/tempSumW;
}


void weiszfeld(Set* ptrSet, double epsilon, int kmax, double Zk1)
{
	double* distances = malloc(sizeof(double)*(ptrSet->size));
	distance(ptrSet, distances, ptrSet->x0, ptrSet->y0);
	double Zk = sum(ptrSet->size, ptrSet, distances);
	//Esto se ejecuta la primera vez
	newPoint(ptrSet, distances); // Seteo de nuevo punto inicial
	free(distances);
	if ( (kmax <= 0) || fabs(Zk - Zk1) < epsilon ){
		//Terminar y guardar
		//printf("Weiszfeld listo! \n");
	}
	else{
		//printf("yip yip! \n");
		weiszfeld(ptrSet, epsilon, kmax-1, Zk);
	}
}

void distance(Set* ptrSet, double* ptrDistances, double x0, double y0)
{
	 double* Xs = ptrSet->ArrayX;
	 double* Ys = ptrSet->ArrayY;
	for(int i = 0; i < ptrSet->size; i++) {
		ptrDistances[i] = sqrt(pow((double)(Xs[i]-x0), (double)(2)) + pow((double)(Ys[i]-y0), (double)(2)));
	}
}

double sum(int i, Set* ptrSet, double* distances)
{
	if(i == 1) {
		return ( (ptrSet->ArrayW[0]) * distances[0] );
	}
	return sum(i-1, ptrSet, distances) + (ptrSet->ArrayW[i-1]) * distances[i-1];
}

void newPoint(Set* ptrSet, double* ptrDistances)
{
	int size = ptrSet->size;
	double tempSumWX = 0;
	double tempSumWY = 0;
	double tempSumW = 0;
	double tempX  = 0;
	double tempY  = 0;
	double tempD = 0;
	double tempW = 0;
	for(int i = 0; i < size; i++){
		tempX  = ptrSet->ArrayX[i];
		tempY  = ptrSet->ArrayY[i];
		tempD  = ptrDistances[i];
		tempW  = ptrSet->ArrayW[i];
		tempSumWX += tempW*tempX/tempD;
		tempSumWY += tempW*tempY/tempD;
		tempSumW  += tempW/tempD;
	}
	ptrSet->x0 = tempSumWX/tempSumW;
    ptrSet->y0 = tempSumWY/tempSumW;
}

/*----------------------------------------------------------*/

void gradient(Set* ptrSet, double epsilon)
{
	//calculo deltax
	double difArray[2];
	Dfunction(ptrSet, ptrSet->x0, ptrSet->y0, difArray);
	double nabArray[2];
	for (int i=0; i<2; i++) {
		nabArray[i] = difArray[i];
	}
	for(int i=0; i<2; i++) {
		difArray[i] = -difArray[i];
	}

	double t = backtracking(ptrSet, difArray, nabArray, ALPHA, BETA);
	ptrSet->x0 = ptrSet-> x0 + (t*difArray[0]);
	ptrSet->y0 = ptrSet-> y0 + (t*difArray[1]);
	double norma = pow(nabArray[0], 2) + pow(nabArray[1],2);
	if(norma >= pow(epsilon,2)  && t >= 0.0000000005){
		//printf("%lf \n", norma);
		gradient(ptrSet, epsilon);
	}

}


double backtracking(Set* ptrSet, double* difArray, double* nabArray, double alpha, double beta)
{
	double xk = ptrSet->x0;
	double yk = ptrSet->y0;
	double t = 1;
	double a = function(ptrSet, xk + t*difArray[0], yk + t*difArray[1]);
	double b = function(ptrSet, xk, yk) + alpha*t*multArray12x21(nabArray, difArray);
	while (a > b) {
		t = beta*t;
		a = function(ptrSet, xk + t*difArray[0], yk + t*difArray[1]);
		b = function(ptrSet, xk, yk) + alpha*t*multArray12x21(nabArray, difArray);
	}
	//printf("t = %.10lf\n", t);
	return t;
}

double function(Set* ptrSet, double x, double y)
{
	double* ArrayX = ptrSet->ArrayX;
	double* ArrayY = ptrSet->ArrayY;
	double* ArrayW = ptrSet->ArrayW;
	double retorno = 0;
	for(int i = 0; i < ptrSet->size; i++){
		retorno += ArrayW[i]*sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));
	}
	return retorno;
}

void Dfunction(Set* ptrSet, double x, double y, double* returnArray)
{
	double* ArrayX = ptrSet->ArrayX;
	double* ArrayY = ptrSet->ArrayY;
	double* ArrayW = ptrSet->ArrayW;

	returnArray[0] = 0;
	returnArray[1] = 0;

	for(int i = 0; i < ptrSet->size; i++){
		returnArray[0] += ArrayW[i]*(1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2)))*(x - ArrayX[i]);
		returnArray[1] += ArrayW[i]*(1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2)))*(y - ArrayY[i]);
	}
}

void invHfunction(Set* ptrSet, double x, double y, double returnArray[][2])
{
	double* ArrayX = ptrSet->ArrayX;
	double* ArrayY = ptrSet->ArrayY;
	double* ArrayW = ptrSet->ArrayW;
	returnArray[0][0] = 0;
	returnArray[0][1] = 0;
	returnArray[1][0] = 0;
	returnArray[1][1] = 0;
	// Calculo H
	for(int i = 0; i < ptrSet->size; i++) {
		double d = 1/sqrt(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));
		double d2 = 1/(pow(ArrayX[i] - x, (double) 2) + pow(ArrayY[i] - y, (double) 2));

		returnArray[0][0] += ArrayW[i]*d*(1-pow(ArrayX[i] - x, (double) 2)*d2);
		returnArray[1][1] += ArrayW[i]*d*(1-pow(ArrayY[i] - y, (double) 2)*d2);
		double temp = -ArrayW[i]*(ArrayX[i] - x)*(ArrayY[i] - y)*pow(d, (double) 3);
		returnArray[0][1] += temp;
		returnArray[1][0] += temp;
	}
	// Calculo invH
	double	a = returnArray[0][0],
			b = returnArray[0][1],
			c = returnArray[1][0],
			d = returnArray[1][1];
	double det = 1/(a*d - b*c);
	returnArray[0][0] = d*det;
	returnArray[1][1] = a*det;
	returnArray[0][1] = -b*det;
	returnArray[1][0] = -d*det;
}

double multArray12x21(double* Array1, double* Array2)
{
	return (Array1[0]*Array2[0] + Array1[1]*Array2[1]);
}

void multArray22x21(double arrayA[][2], double* arrayB, double* arrayReturn)
{
	arrayReturn[0] = arrayA[0][0]*arrayB[0] + arrayA[0][1]*arrayB[1];
	arrayReturn[1] = arrayA[1][0]*arrayB[0] + arrayA[1][1]*arrayB[1];
}

void newton(Set* ptrSet, double epsilon)
{
	//calculo deltax
	double difArray[2];
	double nabArray[2];
	double HArray[2][2];
	Dfunction(ptrSet, ptrSet->x0, ptrSet->y0, nabArray);
	invHfunction(ptrSet, ptrSet->x0, ptrSet->y0, HArray);
	multArray22x21(HArray, nabArray, difArray);

	double tempArray[2];
	for (int i=0; i<2; i++) {
		tempArray[i] = difArray[i];
	}

	for(int i=0; i<2; i++) {
		difArray[i] = -difArray[i];
	}
	//calculo lambda^2
	double lambda2 = multArray12x21(nabArray, tempArray);
	double t = 1;
	//criterio de parada
	if(lambda2/2 >= epsilon  && t >= 0.0000000005){
		//printf("%lf \n", lambda2);
		t = backtracking(ptrSet, difArray, nabArray, ALPHA, BETA);
		//actualizar x e y
		ptrSet->x0 = ptrSet-> x0 + (t*difArray[0]);
		ptrSet->y0 = ptrSet-> y0 + (t*difArray[1]);
		gradient(ptrSet, epsilon);
	}
}

// https://www.wolframalpha.com/input/?i=derivative+w*%28%28a+-+x%29%5E2++%2B+%28b+-+y%29%5E2%29%5E%28-0.5%29+*+%28x-a%29
// falta: evaluacion de fc objetivo en cada instancia
// imprimir numero de iteraciones para cada instancia y guardarlos para calcular promedio y maximo
