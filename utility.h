#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED
#include "complex.h"
#include <stdlib.h>
#include <stdio.h>
#include "lora_phy.h"
#include <stdint.h>
#include <inttypes.h>

/// Remember to free the memory of new vectors after the use with the command free();

//A FUNCTION THAT DYNAMICALLY ALLOCATES AN INITIALIZED VECTOR OF TYPE DOUBLE (ENDPOINTS INCLUDED)//
double* newVetRange(int64_t inizio, int64_t fine);

//FUNCTION THAT DYNAMICALLY ALLOCATES A DOUBLE-COMPLEX VECTOR [a+ I*b] OF DIMENSION = size//
double complex* NewVetCmplx(int64_t size);

//A FUNCTION THAT PRINTS ALL ELEMENTS OF A COMPLEX VECTOR SUCH AS A+iB//
void stampaVetCmplx(double complex v[],int64_t size);

//FUNCTION THAT PRINTS ALL ELEMENTS OF A DOUBLE VECTOR//
void stampaVetDouble(double v[], int64_t size);

//FUNCTION THAT PRINTS ALL ELEMENTS OF AN INTEGERS VECTOR//
void stampaVetInt(int64_t v[], int64_t size);

#endif // UTILITY_H_INCLUDED
