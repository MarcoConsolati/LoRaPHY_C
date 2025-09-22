#include "utility.h"
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <inttypes.h>

double* newVetRange(int64_t inizio, int64_t fine) {
	int64_t lunghezza = fine - inizio;
	double* vettore = (double*)malloc(lunghezza * sizeof(double));
	if ( vettore == NULL ) {
        printf("\nError: Memory not allocated correctly for the vector!\n");
        exit(1);
	}
	for (int64_t i = 0; i < lunghezza; i++) {			/// Array initialized with values
		vettore[i] = i + inizio;						/// v[0]=0.0,v[1]=1.0,...,v[10]=10.0;
	}
return vettore;
}

double complex* NewVetCmplx(int64_t size) {
	double complex *vettore = (double complex *) malloc (size * sizeof(double complex));
	if (vettore == NULL ) {
		printf("\nError: Memory not allocated correctly for the vector!\n");
        exit(1);
	}
	else{
        for (int64_t i = 0; i < size; i++) {
            vettore[i] = 0.0 + 0.0*I;					/// Array initialized with zero values
            }                                           /// composed of a real part and an imaginary part
    }
return vettore;
}

void stampaVetCmplx(double complex v[], int64_t size){
	if(size <= 0 ){
		printf("\nArray size can't be negative or null!\n");

	}
	else{
		for(int64_t i = 0 ; i < size ; i++){
			printf("\nv[%"PRId64"]= %f + (%f*I)\t",i, creal(v[i]), cimag(v[i]) );
		}
	}printf("\n");
}

void stampaVetDouble( double v[], int64_t size){
	if(size <= 0 ){
		printf("\nArray size can't be negative or null!\n");

	}
	else{
		for(int64_t i = 0 ; i < size ; i++){
			printf("\nv[%"PRId64"]= %f\t",i,v[i]);
		}
	}printf("\n");
}

void stampaVetInt(int64_t v[], int64_t size){
	if(size <= 0 ){
		printf("\nArray size can't be negative or null!\n");

	}
	else{
		for(int64_t i = 0 ; i < size ; i++){
			printf("\nv[%"PRId64"]= %"PRId64"\t",i,v[i]);
		}
	}printf("\n");
}

