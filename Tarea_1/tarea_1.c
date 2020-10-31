#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// TEMP
typedef struct set {
	int* ArrayX;
	int* ArrayY;
	float* ArrayW;
	int size;
    float x0;
    float y0;
}Set;

#define EPSILON 0.01
#define KMAX 2
//DECLARACIONES
Set* getNewSet( int );
void generateDarray(Set* );
void printValues(Set* );
int randGenerator();
void deleteArray(Set* );
void initPoint(Set* );
void printInit(Set* );
// ESTAS SON NUEVAS
void distance(Set* ptrSet, float* ptrDistances, int x0, int y0);
void weiszfeld(Set* ptrSet, float epsilon,  int kmax, float Zk1);
float sum(int i, Set* ptrSet, float* distances);
void newPoint(Set* ptrSet, float* ptrDistances);

/*----------------------------------------------------------------------------*/
int main(int argc, char const *argv[]) {
	printf("Hola mundo, soy la tarea de GIO! \n");
	//Set* test  = getNewSet(10);
	int xs[] = {0, 0, 2, 2};
	int ys[] = {0, 2, 0, 2};
	float ws[] = {1.0, 1.0, 1.0, 1.0};
	Set* test = malloc(sizeof(Set));
	test->size = 4;
	test->ArrayX = xs;
	test->ArrayY = ys;
	test->ArrayW = ws;
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
Set* getNewSet( int size){
    Set* newSet = NULL;
    if(size > 0){
        newSet = (Set*)malloc(sizeof(Set));
    }
    newSet->size = size;
    generateDarray(newSet);
	initPoint(newSet);
    return newSet;

}

void generateDarray(Set* ptrSet) {
	int size = ptrSet->size;
	int* arrTempX = malloc(sizeof( int)*size);
	int* arrTempY = malloc(sizeof( int)*size);
	float* arrTempW = malloc(sizeof(float)*size);
    srand48(time(0));
	for( int i = 0; i < size; i++){
		arrTempX[i] = lrand48() % (size+1);
		arrTempY[i] = lrand48() % (size+1);
		//The lrand48() function returns non-negative, long integers, uniformly distributed over the interval [0,RAND_MAX].
        int lower = size/10;
        int upper = lower*5;
		arrTempW[i] = ((float)lrand48()/(float)(RAND_MAX))*(upper-lower) + lower;
	}
	ptrSet->ArrayX = arrTempX;
	ptrSet->ArrayY = arrTempY;
	ptrSet->ArrayW = arrTempW;
	initPoint(ptrSet);
}

void printValues(Set* ptrSet) {
	int size = ptrSet->size;
	for(int i = 0; i < size; i++){
		printf("[%u] Valores Xi = %u, Yi = %u, Wi = %f\n", (i+1),ptrSet->ArrayX[i],ptrSet->ArrayY[i], ptrSet->ArrayW[i]);
	}
}

void deleteArray(Set* ptrSet){
	free(ptrSet->ArrayX);
	free(ptrSet->ArrayY);
	free(ptrSet->ArrayW);
	free(ptrSet);
}
void initPoint(Set* ptrSet) {
	int size = ptrSet->size;
	float tempSumWX = 0;
	float tempSumWY = 0;
	float tempSumW = 0;
	int tempX  = 0;
	int tempY  = 0;
	float tempW = 0;
	for(int i = 0; i < size; i++){
		tempX  = ptrSet->ArrayX[i];
		tempY  = ptrSet->ArrayY[i];
		tempW = ptrSet->ArrayW[i];
		tempSumWX += tempW*tempX;
		tempSumWY += tempW*tempY;
		tempSumW  += tempW;
	}
	ptrSet->x0 = round(tempSumWX/tempSumW);
    ptrSet->y0 = round(tempSumWY/tempSumW);
}

void printInit(Set* ptrSet){
	printf("Soluciones iniciales Xo = %u, Yo = %u\n", ptrSet->x0, ptrSet->y0);
}

void weiszfeld(Set* ptrSet, float epsilon,  int kmax, float Zk1) {
	float* distances = malloc(sizeof(float)*(ptrSet->size));
	distance(ptrSet, distances, ptrSet->x0, ptrSet->y0);

	float Zk = sum(ptrSet->size, ptrSet, distances);
	//Esto se ejecuta la primera vez
	newPoint(ptrSet, distances); // Seteo de nuevo punto inicial
	free(distances);
	if ( (kmax <= 0) || abs(Zk - Zk1) < epsilon ){
		//Terminar y guardar
		printf("Termine ! \n");
	}
	else{
		weiszfeld(ptrSet, epsilon, kmax-1, Zk);
	}
}

void distance(Set* ptrSet, float* ptrDistances, int x0, int y0) {
	 int* Xs = ptrSet->ArrayX;
	 int* Ys = ptrSet->ArrayY;
	for(int i = 0; i < ptrSet->size; i++) {
		ptrDistances[i] = pow(Xs[i]-x0,2) + pow(Ys[i]-y0,2);
		ptrDistances[i] = sqrt(ptrDistances[i]);
	}
}

float sum(int i, Set* ptrSet, float* distances) {
	if(i == 1) {
		return ( (ptrSet->ArrayW[0]) * distances[0] );
	}
	return sum(i-1, ptrSet, distances) + (ptrSet->ArrayW[i-1]) * distances[i-1];
}

void newPoint(Set* ptrSet, float* ptrDistances) {
	int size = ptrSet->size;
	float tempSumWX = 0;
	float tempSumWY = 0;
	float tempSumW = 0;
	int tempX  = 0;
	int tempY  = 0;
	float tempD = 0;
	float tempW = 0;
	for(int i = 0; i < size; i++){
		tempX  = ptrSet->ArrayX[i];
		tempY  = ptrSet->ArrayY[i];
		tempD  = ptrDistances[i];
		tempW  = ptrSet->ArrayW[i];
		tempSumWX += tempW*tempX/tempD;
		tempSumWY += tempW*tempY/tempD;
		tempSumW  += tempW/tempD;
	}
	ptrSet->x0 = round(tempSumWX/tempSumW);
    ptrSet->y0 = round(tempSumWY/tempSumW);
}
