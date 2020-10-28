#ifndef UTILS_H
#define UTILS_H
// TEMP
typedef struct set {
	unsigned int* ArrayX;
	unsigned int* ArrayY;
	float* ArrayW;
	unsigned int size;
    unsigned int x0;
    unsigned int y0;
}Set;

//DECLARACIONES
Set* getNewSet(unsigned int );
void generateDarray(Set* );
void printValues(Set* );
int randGenerator();
void deleteArray(Set* );
void initPoint(Set* );
void printInit(Set* );


/*----------------------------------------------------------------------------*/



#endif
