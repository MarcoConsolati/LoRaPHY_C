#ifndef LORA_PHY_H_INCLUDED
#define LORA_PHY_H_INCLUDED
#include <complex.h>
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <math.h>
#include <fftw3.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

///Structures for storing the output values of some functions

typedef struct {
    int64_t   index;
    double value;
} Peak;

typedef struct {
	double complex* vet;
	int64_t size_chirp;
} ResVetCmplx;

typedef struct {
    double complex* vet;
    int64_t size_file;
} ResFileCmplx;

typedef struct {
    double* vet;
    int64_t size_file;
} ResFileDouble;

typedef struct {
    int64_t* vet;
    int64_t length_vet;
    double cfo;
} ResVetInt;


//Main struct LoRaPHY//

typedef struct{
	int64_t	    rf_freq;                   	/// carrier frequency
	int64_t	    sf;			            	/// spreading factor (7,8,9,10,11,12)
	int64_t	    bw;							/// bandwidth (125kHz 250kHz 500kHz)
	int64_t	    fs;							/// sampling frequency
	int64_t	    cr;							/// code rate: (1:4/5 2:4/6 3:4/7 4:4/8)
	int64_t     payload_len;			    /// payload length
	int64_t	    has_header;					/// explicit header: 1, implicit header: 0
	int64_t	    crc;        				/// crc = 1 if CRC Check is enabled else 0
	int64_t	    ldr;       					/// ldr = 1 if Low Data Rate Optimization is enabled else 0
	int64_t	    whitening_seq[255];         /// whitening sequence inizialized in the buildLora function
	int64_t	    header_checksum_matrix[5][12];///we use a 5 x 12 matrix to calculate header checksum
	int64_t		preamble_len;				/// preamble length
	double complex* sig;                    /// input baseband signal
	int64_t         sig_size;               /// number of floats elements in file, read with readCmplx before demodulate!
	double complex*	downchirp;				/// ideal chirp with decreasing frequency from B/2 to -B/2
	double complex*	upchirp;				/// ideal chirp with increasing frequency from -B/2 to B/2
	int64_t		size_downchirp;				/// size of downchirp complex vector
	int64_t		size_upchirp;				/// size of upchirp complex vector
	int64_t		sample_num;					/// number of sample points per symbol
	int64_t		bin_num;					/// number of bins after FFT (with zero padding)
	int64_t		zero_padding_ratio;			/// FFT zero padding ratio
	int64_t     fft_len;					/// FFT size
	int64_t		preamble_bin;				/// reference bin in current decoding window, used to eliminate CFO
    double	    cfo;						/// carrier frequency offset
	int64_t		is_debug;					/// set 'true'=1 for debug information set 'false'=0 for no info
	int64_t		hamming_decoding_en;		/// enable hamming decoding set 'true'=1 or 'false'=0
} LoRaPHY;


//FUNCTION BUILDLORA CONSTRUCTS AND INITIALIZES SOME PARAMETERS OF THE LoRaPHY STRUCT AND INVOKES INIT FUNCTION//
void buildLora(LoRaPHY* Lora , int64_t rf_freq, int64_t  sf, int64_t  bw, int64_t  fs);

//FUNCTION INIT THAT INITIALIZES THE SUBSEQUENT LoRaPHY PARAMETERS//
void init(LoRaPHY* Lora);

//FUNCTION CALC_SYM_NUM CALCULATE THE NUMBER OF SYMBOLS.
//ARGUMENTS: LORA, PLEN(Payload length)
//OUTPUT: SYM_NUM (numero di simboli)
double calc_sym_num(LoRaPHY* Lora, int64_t plen);

//FUNCTION MAX_CALC CALCULATE THE MAX NUMBER FROM TWO VAR//
double max_calc(double a , double b);

//FUNCTION CHIRP GENERATES UPCHIRP_DOWNCHIRP SIGNAL//
ResVetCmplx chirp(bool is_up, int64_t sf, int64_t bw, int64_t fs, int64_t h , double cfo , int64_t tdelta );

//FUNCTION WRITE_FILE WRITE DOUBLE COMPLEX DATA INTO A FILE CALLED FILENAME//
int16_t writeCmplxFile(double complex* vet, int64_t size, const char* filename);

void readCmplxFile(double complex* vet, int64_t size, const char* filename);

//FUNCTION READ_FILE READ THE DATA FROM A FILE OF DOUBLE COMPLEX DATA//
double complex* read_file(char* filename,int64_t* size);

//FUNCTION READCMPLX THAT READS FROM THE SIG. CFILE AND RETURNS ITS CONTENTS IN A COMPLEX FLOAT VECTOR AND ITS SIZE//
ResFileCmplx readCmplx(const char* filePath);

//FUNCTION READDOUBLE READS FROM THE FILE INCLUDED IN THE PATH AND RETURNS ITS CONTENTS IN A DOUBLE VECTOR AND ITS SIZE (USED FOR THE FIR FILTER AND ANTIALIASING COEFFICIENTS)
ResFileDouble readDouble(const char* filePath);

//FUNCTION INFO_LORA RETURNS VALUES INSIDE THE LoRaPHY STRUCT//
void info_Lora(LoRaPHY* Lora);

//FUNCTION THAT RETURNS THE CONVOLUTION OF TWO DOUBLE COMPLEX SIGNALS INTO A DOUBLE COMPLEX SIGNAL//
void ConvolveCmplx(double complex* input, int64_t inLength, double* h_real, double complex* h_imag, int64_t h_Length, double complex* output);

//FUNCTION THAT CONVERTS AN INTEGER TO 10-BIT BINARY //
int64_t* convertiInt10Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 9-BIT BINARY //
int64_t* convertiInt9Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 8-BIT BINARY //
int64_t* convertiInt8Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 7-BIT BINARY //
int64_t* convertiInt7Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 6-BIT BINARY //
int64_t* convertiInt6Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 5-BIT BINARY //
int64_t* convertiInt5Bit(int64_t value);

//FUNCTION THAT CONVERTS AN INTEGER TO 4-BIT BINARY //
int64_t* convertiInt4Bit(int64_t value);

//FUNCTION THAT CONVERTS A BINARY VECTOR TO A DECIMAL NUMBER//
int64_t binToDec(int64_t* vetBin, int64_t length);

//FUNCTION MOD_ SIMILAR TO MATLAB BECAUSE MOD(%) IN C BEHAVES DIFFERENTLY//
double mod_double(double a , double b);

//FUNCTION BIT_REDUCE IT EXTRACTS THE BITS OF A SERIES OF PREFIXED INDEXES FROM THE WORDS AND MAKES XOR BIT BY BIT//
int64_t* bit_reduce(int64_t* words, int64_t* index, int64_t lengthwords, int64_t sizeindex);

//FUNCTION WORD_REDUCE ENCODES THE BITS OF A SERIES OF PREFIXED INDEXES IN WORDS AND MAKES OR BIT BY BIT//
int64_t word_reduce(int64_t* words, int64_t lengthwords);

Peak dechirp(LoRaPHY* Lora, int64_t  x , bool is_up);

int64_t detectnew(LoRaPHY* Lora, int64_t startIndx);

int64_t syncnew(LoRaPHY* Lora, int64_t x);

double calc_sym_num(LoRaPHY* Lora, int64_t plen);

int64_t parse_header(LoRaPHY* Lora, int64_t* data , int64_t lengthdata);

int64_t calc_payload_len(LoRaPHY* Lora, int64_t slen, bool no_redundant_bytes);

ResVetInt demodulate(LoRaPHY* Lora, double complex* vet_from_file, int64_t vetsize_from_file );

ResVetCmplx modulate(LoRaPHY* Lora , int64_t* symbols , int64_t length_symbols);

int64_t* dynamic_compensation(LoRaPHY* Lora, int64_t* data , int64_t length_data);

int64_t* gray_coding(LoRaPHY* Lora, int64_t* din , int64_t length_din);

int64_t* gray_decoding(LoRaPHY* Lora, int64_t* symbols_i , int64_t lengthsymbols_i);

int64_t* diag_deinterleave(LoRaPHY* Lora, int64_t* symbols , int64_t lengthsymbol , int64_t ppm);

int64_t* diag_interleave(LoRaPHY* Lora, int64_t* codewords, int64_t lengthwords, int64_t rdd);

int64_t* hamming_decode(LoRaPHY* Lora,int64_t* words,int64_t lengthwords,int64_t rdd);

int64_t* hamming_encode(LoRaPHY* Lora , int64_t* nibbles , int64_t length_nibbles);

ResVetInt encode(LoRaPHY* Lora, int64_t* payload , int64_t length_payload);

ResVetInt decode(LoRaPHY* Lora, int64_t* symbols_m ,  int64_t number_cols , int64_t number_rows);

int64_t** whiten(LoRaPHY* Lora, int64_t* data , int64_t lengthdata);

int64_t* dewhiten(LoRaPHY* Lora, int64_t* bytes , int64_t lengthbytes);

int64_t* gen_header(LoRaPHY* Lora, int64_t payload_length);

ResVetInt symbols_to_bytes(LoRaPHY* Lora, int64_t* symbols , int64_t length_symbols);

double time_on_air(LoRaPHY* Lora, int64_t plen);


#endif // LORA_PHY_H_INCLUDED
