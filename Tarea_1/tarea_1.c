#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// TEMP
typedef struct set {
	double* ArrayX;
	double* ArrayY;
	double* ArrayW;
	int size;
    double x0;
    double y0;
}Set;

#define EPSILON 0.0001
#define KMAX 100
#define I_SIZE 100

//DECLARACIONES
Set* getNewSet(int size);
void generateDarray(Set* ptrSet);
void printValues(Set* ptrSet);
void deleteArray(Set* ptrSet);
void initPoint(Set* ptrSet);
void printInit(Set* ptrSet);
void weiszfeld(Set* ptrSet, double epsilon,  int kmax, double Zk1);
void distance(Set* ptrSet, double* ptrDistances, double x0, double y0);
double sum(int i, Set* ptrSet, double* distances);
void newPoint(Set* ptrSet, double* ptrDistances);

/*----------------------------------------------------------------------------*/
int main(int argc, char const *argv[]) {
	printf("Hola mundo, soy la tarea de GIO! \n");
	Set* test  = getNewSet(I_SIZE);
	/*
	int xs[] = {0, 0, 2, 2};
	int ys[] = {0, 2, 0, 2};
	double ws[] = {1.0, 1.0, 1.0, 1.0};
	Set* test = malloc(sizeof(Set));
	test->size = 4;
	test->ArrayX = xs;
	test->ArrayY = ys;
	test->ArrayW = ws;
	*/
	initPoint(test);
	printValues(test);
	printInit(test);
	weiszfeld(test, EPSILON, KMAX, 0);
	printInit(test);
	//deleteArray(test);
	return 0;
}
/*------------------------------------------------------------*/
// for (1 -- 100){}
//     genrarset (10)
//     for (1 3){
//
//         marca tiemo incio
//         algoritmo i
//         guardar en arreglo el tiempo
//         marca tiempo final
//     }
//     freeset();
// }
// calcular promedio


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
    srand48(time(0));
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

void printInit(Set* ptrSet)
{
	printf("Soluciones iniciales Xo = %f, Yo = %f\n", ptrSet->x0, ptrSet->y0);
}

void weiszfeld(Set* ptrSet, double epsilon,  int kmax, double Zk1)
{
	double* distances = malloc(sizeof(double)*(ptrSet->size));
	distance(ptrSet, distances, ptrSet->x0, ptrSet->y0);

	double Zk = sum(ptrSet->size, ptrSet, distances);
	//Esto se ejecuta la primera vez
	newPoint(ptrSet, distances); // Seteo de nuevo punto inicial
	free(distances);
	if ( (kmax <= 0) || abs(Zk - Zk1) < epsilon ){
		//Terminar y guardar
		printf("Weiszfeld listo! \n");
	}
	else{
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
