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

#define EPSILON 0.0001
#define KMAX 30
#define I_SIZE 10000
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
	int size;			/**< Tamaño de los arreglos, asociado a cardinalidad I */
	double x0;			/**< Solucion inicial en X */
	double y0;			/**< Solucion inicial en Y */
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
/**
 * @brief deleteArray procede a dealocar los arreglos y estructuras creadas
 *
 * Elimina los arreglos asociados al conjunto set pasado, para finalmente
 * eliminar este ultimo.
 *
 * @param[in] ptrSet Puntero, a la estructura Set a dealocar.
*/
void deleteArray(Set* ptrSet);
/**
 * @brief initPoint procede a calcular los puntos iniciales de un Set.
 *
 * La funcion debe ser llamada luego de generar los arreglos asociados a la
 * estructura Set, esta procede a utilizar la formula entregada para los puntos
 * inciales dada en el documento de trabajo. La funcion luego de calcular los
 * valores de estos puntos, los asocia a los elementos X0 e Y0 de la estructura.
 *
 * @param[in] ptrSet Puntero a la estructura a operar sobre.
*/
void initPoint(Set* ptrSet);
/**
 * @brief weiszfeld aplicar el algoritmo del mismo nombre sobre un Set.
 *
 * Esta funcion aplica el algoritmo de Weiszfeld de forma recursiva sobre un Set
 * pasado, dado una condición de termino definida por un epsilon, un número
 * máximo de iteraciones dado por kmax, y un valor Zk1, el cual debe ser
 * inicializado a cero al llamar la funcion, ya que es utilizado por el
 * el algoritmo para la recursividad.
 * Funcionamiento interno: En primera instancia, se calcula un vector de
 * distancias, el cual permite calcular el valor Zk(0) y nuevos puntos
 * iniciales para el set. Luego se verifica la condición de termino, si
 * fabs(Zk - Zk1) < epsilon. Si esto se cumple, el algoritmo termina, en caso
 * contrario se llama de forma recursiva, pasando como argumento Zk1 el
 * Zk(0) calculado en la presente iteracion. De forma analoga tambien, para cada
 * iteracion se verifica si se llego al numero maximo de iteraciones kmax, si es
 * asi, el algoritmo se detiene. Los resultados se alamacenan como los nuevos
 * puntos X0 e Y0 del set pasado.
 *
 * @param[in] ptrSet Puntero al set sobre el cual se aplicara el algoritmo
 * @param[in] epsilon Valor para la condición de termino, generalmente 10^-4
 * @param[in] kmax Numero máximo de iteraciones para el algoritmo
 * @param[in] Zk1 Parametro para llamadas recursiva, inicializar a 0
*/
int weiszfeld(Set* ptrSet, double epsilon,  int kmax, double Zk1);
/**
 * @brief distance inicializa un arreglo con las distancias di
 *
 * La función inicializa un arreglo, con el cálculo de las distancias
 * eucliadianas entre los puntos Xi e Yi, con los puntos iniciales X0, Y0, que
 * se obtienen del set pasado. La funcion retorna un entero que indica el
 * número total de iteraciones realizadas
 *
 * @param[in] ptrSet Puntero al set sobre el cual se calculara las distancias
 * @param[in] ptrDistances Puntero a un arreglo donde se guarda el resultado
 * @param[out] cInterno contador interno de las iteraciones del algoritmo.
*/
void distance(Set* ptrSet, double* ptrDistances);
/**
 * @brief sum  Calcula el valor de la sumatoria que define a Z(k)
 *
 * La función se llama de manera recursiva para calcular el valor de la
 * sumatoria que está definida por Wi*di(k). Retorna el valor del resultado.
 *
 * @param[in] i Corresponde a la cardinalidad del conjunto I
 * @param[in] ptrSet Puntero a la estructura sobre la cual se operara.
 * @param[in] distances Puntero a arreglo de distancias de la iteracion k.
 * @param[out] sum Resultado que contiene el valor de Zk.
*/
double sum(int i, Set* ptrSet, double* distances);
/**
 * @brief newPoint Recalcula los puntos iniciales considerando la distancia d(k)
 *
 * La función realiza el recálculo de los puntos iniciales, considerando la
 * distancia calculada en la iteración k-ésima - 1, para poder estimar los
 * puntos iniciales X0(k) e Y0(k), estos se reescriben en la estructura.
 *
 * @param[in] ptrSet Puntero al Set sobre el cual se operara.
 * @param[in] distances Puntero al arreglo de distancias de iteracion k-1
*/
void newPoint(Set* ptrSet, double* ptrDistances);
/**
 * @brief gradient Aplica el algoritmo de gradiente con backtracking sobre Set
 *
 * La función aplica el algoritmo resolutivo del método de gradiente con
 * backtracking sobre la estructura set pasado, de forma recursiva.
*/
int gradient(Set* ptrSet, double epsilon);
void Dfunction(Set* ptrSet, double x, double y, double* returnArray);
double backtracking(Set* ptrSet, double* difArray, double* nabArray, double alpha, double beta);
double function(Set* ptrSet, double x, double y);
double multArray12x21(double* Array1, double* Array2);
void multArray22x21(double arrayA[][2], double* arrayB, double* arrayReturn);
void invHfunction(Set* ptrSet, double x, double y, double returnArray[][2]);
int newton(Set* ptrSet, double epsilon);


//Revision para eliminacion
void printValues(Set* ptrSet);



#endif
