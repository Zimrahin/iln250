#ifndef UTILS_H
#define UTILS_H
// TEMP
typedef struct set {
	double* ArrayX;
	double* ArrayY;
	double* ArrayW;
	int size;
    double x0;
    double y0;
}Set;

#define EPSILON 0.001
#define KMAX 30
#define I_SIZE 1000
#define INSTANCES 1
#define ALPHA 0.3 // 0<ALPHA<0.5
#define BETA 0.6 // 0<BETA<1



// DECLARACIONES
Set* getNewSet(int size);
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
