/****************************************************************************
 * Tarea #1  - [ILN250] Gestión de información y operación                  *
 *                                                                          *
 * Este archivo es parte de la entrega del código de la seción de alogorit- *
 * mos de minimización programados para el problema SSWB.                   *
 *                                                                          *
 ****************************************************************************/
 /**
 * @file utils.h
 * @author Mauricio Aravena - Diego Badillo
 * @date 4 Nov 2020
 * @brief Archivo header, de las funciones implementadas en esta tarea.
 *
 * Este archivo contiene las declaraciones de las funciones implementadas para
 * cumplir con los requisitos pedidos en la tarea, contiene los dos algoritmo
 * bases, Gradiente y Newton, asi como todas las funciones relacionadas, como
 * la función para estimar el paso t, por el metodo de backtracking. Tambien, se
 * encuentra la declaracion de una estructura, la cual se convierte permite
 * agrupar los parametros del problema, como sus soluciones iniciales.
 *
 */
#ifndef UTILS_H
#define UTILS_H

// TEMPORAL POR COMENTAR

#define EPSILON 0.001
#define KMAX 30
#define I_SIZE 1000
#define ALPHA 0.3 // 0<ALPHA<0.5
#define BETA 0.6 // 0<BETA<1

/**
 * @brief Set: Estructura que agrupa los distintos conjuntos del problema
 *
 * La estructura Set, permite agrupar los conjuntos que están definidos para el
 * SSWB, permite agrupar los arreglos que contienen las coordenadas X,Y de los
 * distintos nodos, como tambien el conjunto de pesos asociado a cada uno.
 * También contiene ciertos parámetros de interés que facilitan la resolución
 * mediante los algoritmos implementados.
 */
typedef struct set {
	double* ArrayX;		/**< Puntero a un arreglo que contiene los puntos Xi   */
	double* ArrayY;		/**< Puntero a un arreglo que contiene los puntos Yi   */
	double* ArrayW;		/**< Puntero a un arreglo que contiene los pesos Wi    */
	int size;					/**< Tamaño de los arreglos, asociado a cardinalidad I */
  double x0;				/**< Solucion inicial en X */
  double y0;				/**< Solucion inicial en Y */
}Set;

/**
 * @brief getNewSet inicializa una nueva estructura, dado un tamaño
 *
 * La función se encarga, de inicializar la estructura, generando los arreglos
 * que contienen los los puntos Xi/Yi y los pesos asociados Wi y asignando los
 * punteros a las variables correspondientes de la estructura. Para esto, se le
 * debe entregar el tamaño, que esta asociado a la cardinalidad del conjunto I
 * deseado. Además se calcula y asigna los asociados a las soluciones iniciales
 * para X/Y a utilizar en los algoritmos
 *
 * @param[in] size La cardinalida del conjunto I deseado.
 * @param[out] newSet El puntero que contiene el set generado
*/
Set* getNewSet(int size);
/**
 * @brief generateDarray genera los arrays que contienen Xi, Yi y Wi.
 *
 * Esta funcion es la encargada de generar los arreglos que contienen los puntos
 * Xi, Yi de forma uniformemente distribuida discreta, asi como también el
 * arreglo de los pesos asociados Wi, de forma uniformemente distribuida.
 * Luego de generar los arreglos asociada los punteros a los punteros del Set
 * pasado. Para esto, al menos debe estar ya asignado el tamaño, en la
 * estructura Set pasada como parametro.
 *
 * @param[in] ptrSet Estructura a la cual se le generara Xi, Yi y Wi
*/
void generateDarray(Set* ptrSet);

void printValues(Set* ptrSet);
void deleteArray(Set* ptrSet);
void initPoint(Set* ptrSet);
void weiszfeld(Set* ptrSet, double epsilon,  int kmax, double Zk1);
void distance(Set* ptrSet, double* ptrDistances, double x0, double y0);
double sum(int i, Set* ptrSet, double* distances);
void newPoint(Set* ptrSet, double* ptrDistances);
void Dfunction(Set* ptrSet, double x, double y, double* returnArray);
double backtracking(Set* ptrSet, double* difArray, double* nabArray, double alpha, double beta);
double function(Set* ptrSet, double x, double y);
double multArray12x21(double* Array1, double* Array2);
void multArray22x21(double arrayA[][2], double* arrayB, double* arrayReturn);
void gradient(Set* ptrSet, double epsilon);
void invHfunction(Set* ptrSet, double x, double y, double returnArray[][2]);
void newton(Set* ptrSet, double epsilon);



#endif
