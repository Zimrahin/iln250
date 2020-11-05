/****************************************************************************
 * Tarea #1  - [ILN250] Gestión de información y operación                  *
 *                                                                          *
 * Este archivo es parte de la entrega del código de la seción de alogorit- *
 * mos de minimización programados para el problema SSWB.                   *
 *                                                                          *
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
#include <time.h>		/**< Funciones relacionadas a la medición de tiempo */
#include <math.h>		/**< Biblioteca de matematicas */
#include <sys/time.h>	/**< System Time */
/**
* [Biblioteca definida por usuario]
*/
#include "utils.h"		/**< Archivo con las declaraciones de las funciones */

/**
 * [ Constantes de funcionamientos ]
*/
#define INSTANCES 100   /**< Numero de instancias por defecto a ejecutar */

/**
 * [ Variables globales]
*/
unsigned long avgTime[3] = {0., 0., 0.};  /**< arreglo tiempos promedios*/
unsigned long maxTime[3] = {0., 0., 0.};  /**< arreglo tiempos maximos */
double avgIteration[3]  = {0., 0., 0.};   /**< arreglo promedio iteraciones */
int maxIteration[3]  = {0, 0, 0};         /**< arreglo maximo iteraciones */
double avgFc[3] = {0., 0., 0.};           /**< arreglo promedio FO */
double maxFc[3] = {0., 0., 0.};			  /**< Arreglo maximo FO */

/****************************************************************************
									MAIN
 ****************************************************************************/
int main()
{
	srand48(time(0));
	struct timeval stop, start;
	unsigned long timeTaken;
	int nIteration[3] = {0, 0, 0}; // N° iteraciones
	double Fc[3] = {0., 0., 0.}; // Fc objetivo

	for(int i=0; i<INSTANCES; i++) {
		Set* test  = getNewSet(I_SIZE);
		printf("-------------------------------------------------------");
		printf("\nInstancia nº %i \n\n", i);
		printf("Soluciones iniciales: X0 = %f, Y0 = %f\n\n", test->x0, test->y0);
		for(int j=0; j<3; j++) { // para cada algoritmo
			initPoint(test);
			switch (j)
			{
				case 0:
					printf("Weiszfeld\n");
					gettimeofday(&start, NULL);
					nIteration[0] = weiszfeld(test, EPSILON, KMAX, 0);
					gettimeofday(&stop, NULL);
					printf("Soluciones finales: X* = %f, Y* = %f\n", test->x0, test->y0);
					break;
				case 1:
					printf("Gradiente\n");
					gettimeofday(&start, NULL);
					nIteration[1] = gradient(test, EPSILON);
					gettimeofday(&stop, NULL);
					printf("Soluciones finales: X* = %f, Y* = %f\n", test->x0, test->y0);
					break;
				case 2:
					printf("Newton\n");
					gettimeofday(&start, NULL);
					nIteration[2] = newton(test, EPSILON, 1);
					gettimeofday(&stop, NULL);
					printf("Soluciones finales: X* = %f, Y* = %f\n", test->x0, test->y0);
					break;
			}

			timeTaken = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;

			//prints
			printf("tiempo: %lu [us]\n", timeTaken);
			printf("Nº iteraciones: %i\n", nIteration[j]);
			//Analisis temporal
			if (timeTaken > maxTime[j])
				maxTime[j] = timeTaken; //se guarda max
			avgTime[j] += timeTaken; //se guarda sum
			//Analisis iteracion
			if (nIteration[j] > maxIteration[j])
				maxIteration[j] = nIteration[j];
			avgIteration[j] += nIteration[j];
			//Funcion objetivo
			Fc[j] = function(test, test->x0, test->y0);
			printf("Funcion objetivo: %lf\n\n", Fc[j]);
			if (Fc[j] > maxFc[j])
				maxFc[j] = Fc[j];
			avgFc[j] += Fc[j];
		}
		deleteArray(test);
	}

	//calculo promedios
	for (int i = 0; i < 3; i++) {
		avgTime[i] /= INSTANCES;
		avgIteration[i] /= INSTANCES;
		avgFc[i] /= INSTANCES;
	}

	printf("\n\n ---------------- MAXIMOS Y PROMEDIOS ----------------\n\n");
	printf("Weiszfeld\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", maxTime[0], avgTime[0] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[0], avgIteration[0]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[0], avgFc[0]);

	printf("Gradiente\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", maxTime[1], avgTime[1] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[1], avgIteration[1]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[1], avgFc[1]);

	printf("Newton\n");
	printf("tmax = %lu [us], tavg = %lu [us]\n", maxTime[2], avgTime[2] );
	printf("nIterationMax = %i, nIterationAvg = %.2lf \n", maxIteration[2], avgIteration[2]);
	printf("maxFc = %.4lf, avgFc = %.4lf\n\n", maxFc[2], avgFc[2]);

	return 0;
}
