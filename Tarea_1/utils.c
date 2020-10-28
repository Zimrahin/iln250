#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"
/*----------------------------------------------------------------------------*/
/*funciones*/
Set* getNewSet(unsigned int size){
    Set* newSet = NULL;
    if(size > 0){
        newSet = (Set*)malloc(sizeof(Set));
    }
    newSet->size = size;
    generateDarray(newSet);
	initPoint(newSet);
    return newSet;

}

void generateDarray(Set* ptrSet){
	unsigned int size = ptrSet->size;
	unsigned int* arrTempX = malloc(sizeof(unsigned int)*size);
	unsigned int* arrTempY = malloc(sizeof(unsigned int)*size);
	float* arrTempW = malloc(sizeof(float)*size);
    srand(time(0));
	for(unsigned int i = 0; i < size; i++){
		arrTempX[i] = rand() % (size+1);
		arrTempY[i] = rand() % (size+1);
        int lower = size/10;
        int upper = lower*5;
		arrTempW[i] = ((float)rand()/(float)(RAND_MAX))*(upper-lower) + lower;
	}
	ptrSet->ArrayX = arrTempX;
	ptrSet->ArrayY = arrTempY;
	ptrSet->ArrayW = arrTempW;
	initPoint(ptrSet);
}

void printValues(Set* ptrSet){
	unsigned int size = ptrSet->size;
	for(unsigned int i = 0; i < size; i++){
		printf("[%u] Valores Xi = %u, Yi = %u, Wi = %f\n", (i+1),ptrSet->ArrayX[i],ptrSet->ArrayY[i], ptrSet->ArrayW[i]);
	}
}

void deleteArray(Set* ptrSet){
	free(ptrSet->ArrayX);
	free(ptrSet->ArrayY);
	free(ptrSet->ArrayW);
	free(ptrSet);
}

void initPoint(Set* ptrSet){
	unsigned int size = ptrSet->size;
	float tempSumWX = 0;
	float tempSumWY = 0;
	float tempSumW = 0;
	unsigned int tempX  = 0;
	unsigned int tempY  = 0;
	float tempW = 0;
	for(unsigned int i = 0; i < size; i++){
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

/*------------*/
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
