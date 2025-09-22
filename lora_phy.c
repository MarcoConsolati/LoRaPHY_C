#include "lora_phy.h"
#include <tgmath.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <complex.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

void buildLora(LoRaPHY* Lora , int64_t rf_freq, int64_t  sf, int64_t  bw, int64_t  fs){
	Lora->rf_freq = rf_freq;
	Lora->sf = sf;
	Lora->bw = bw;
	Lora->fs = fs;
	Lora->cr = -1;               ///initialized to -1 as a sentinel value, since it will be set by the user ///
	Lora->payload_len = -1;      ///initialized to -1 as a sentinel value, since it will be set by the calc_payload_len function///
	Lora->preamble_bin = -1;     ///initialized to -1 as a sentinel value, since it will be set by the syncnew function///
	Lora->has_header = 1;
	Lora->crc = 0;
	Lora->is_debug = 0;
	Lora->hamming_decoding_en = 1;
	Lora->zero_padding_ratio = 10;
	Lora->cfo = 0.0;
	unsigned char data[5][12]={{1,1,1,1,0,0,0,0,0,0,0,0},{1,0,0,0,1,1,1,0,0,0,0,1},{0,1,0,0,1,0,0,1,1,0,1,0},{0,0,1,0,0,1,0,1,0,1,1,1},{0,0,0,1,0,0,1,0,1,1,1,1}};
	for (int64_t i = 0; i < 5; i++) {
        for (int64_t j = 0; j < 12; j++) {
            Lora->header_checksum_matrix[i][j] = data[i][j];
        }
    }
    unsigned char reg = 0xFF;										///setting the whitening sequence
    unsigned char mask8 = 1 << 7;									///
    unsigned char mask6 = 1 << 5;									///
    unsigned char mask5 = 1 << 4;									///
    unsigned char mask4 = 1 << 3;									///
    unsigned char bit8 , bit6 , bit5 , bit4 , bitshift;				///
//    printf("\nSequenza di whitening:\n");                         ///
    for(int64_t i = 0 ; i < 255 ; i++){								///
    Lora->whitening_seq[i] = reg;									///
    	bit8 = (reg & mask8) >> 7;									///
    	bit6 = (reg & mask6) >> 5;									///
    	bit5 = (reg & mask5) >> 4;									///
    	bit4 = (reg & mask4) >> 3;									///
    	bitshift = reg << 1;										///
    	reg = bitshift ^ (bit8 ^ (bit6 ^ (bit5 ^ bit4)));
//        printf("%3d ",Lora->whitening_seq[i]);
	}
	Lora->preamble_len = 8;
	init(Lora);
}

void init(LoRaPHY* Lora){
	ResVetCmplx result_chirp1;
	ResVetCmplx result_chirp2;
	Lora->bin_num = (pow(2, Lora->sf) * Lora->zero_padding_ratio);
	Lora->sample_num = 2 * pow(2, Lora->sf);
	Lora->fft_len = Lora->sample_num * Lora->zero_padding_ratio;
	result_chirp1 = chirp(false , Lora->sf , Lora->bw , 2*Lora->bw , 0 , Lora->cfo , 0 );
	Lora->downchirp = result_chirp1.vet;
	Lora->size_downchirp = result_chirp1.size_chirp;
	result_chirp2 = chirp(true  , Lora->sf , Lora->bw , 2*Lora->bw , 0  , Lora->cfo , 0 );
	Lora->upchirp = result_chirp2.vet;
	Lora->size_upchirp = result_chirp2.size_chirp;
	Lora->sig = NULL;
	if (pow(2,Lora->sf) / (double)Lora->bw > 16e-3)
		Lora->ldr = 1;
	else
		Lora->ldr = 0;
}

double max_calc(double a , double b){
double result = 0.0;
result = (a > b) ? a : b;
return result;
}

#define M_PI 3.14159265358979323846
ResVetCmplx chirp(bool is_up, int64_t sf, int64_t bw, int64_t fs, int64_t h , double cfo , int64_t tdelta ) {
	int64_t tscale = 1;
	int64_t k = 0;
	long long b_w = (long long) bw * bw;                ///If I don't cast the product operation bw*bw I lose significant amounts of resolving digits
	float fzero = 0.0, phi = 0.0;						///is_up:`true` if constructing an up-chirp signal `false` if constructing a down-chirp	signal										//sf: Spreading Factor	//bw: Bandwidth	//h: Start frequency offset (0 to 2^SF-1)
	int64_t snum1 = 0,snum2 = 0,snum = 0, size_res =0;	///fs: Sampling Frequency //tdelta: Time offset (0 to 1/fs)
	int64_t N = pow(2,sf);                              ///tscale: Scaling the sampling frequency
//	printf("\nN=%"PRId64" , sf=%"PRId64" , bw=%"PRId64" , fs=%"PRId64" , h=%"PRId64"",N , sf , bw , fs , h);
	double T = (double) N / bw;
	int64_t samp_per_sym =  round(((double)fs / bw)* N);
	int64_t h_orig = h;
	h = round(h);
	cfo = (cfo + (h_orig - h)) / (double) N * bw;
	if (is_up == true){
		k = (b_w) / (double) N;
		fzero = - bw / (double)(2.0 + cfo);
	}
		else {
			k = - ( b_w ) / (double) N;
			fzero = bw /(double)(2.0 + cfo);
		}
	snum = (double) samp_per_sym * (N - h) / (double)N;
//	printf("\n\tsnum=%"PRId64"",snum);
	snum1 = snum;
	double* t0 = newVetRange(0, snum + 1);
	for (int64_t i = 0; i <= snum ; i++) {

		t0[i] = (t0[i] /(double)fs) * tscale + tdelta;
//        printf("\nt[%"PRId64"]=%f\t",i,t0[i]);
	}
		double complex* c1 = NewVetCmplx(snum + 1);
		for (int64_t i = 0; i <= snum ; i++) {

			c1[i] = cexp(I * 2 * M_PI *(t0[i]*(fzero + k * T * h /  N + 0.5 * k * t0[i])));

		}
//		stampaVetCmplx(c1 , snum);
		if (snum == 0) {
			phi=0.0;
		}else{
            phi = carg(c1[snum]);						///phase in radiant
            }
//        printf("\nc1(snum)= %f+%f*I",creal(c1[snum]),cimag(c1[snum]));
        free(t0);
        t0 = NULL;
        snum = (double) (samp_per_sym * h )/ N; //modified da N - 1 a N//
//        printf("\n\tsnum2=%"PRId64"",snum);
        snum2 = snum;
        double* t_ = newVetRange(0, snum + 1);
        for (int64_t i = 0; i <= snum ; i++) {

            t_[i] = (t_[i] / (double)fs) + tdelta;
//            printf("\nt[%"PRId64"]=%f\t",i,t_[i]);
        }
        double complex* c2 = NewVetCmplx(snum);
        for (int64_t i = 0; i < snum; i++) {

            c2[i] = cexp(I*(phi + 2 * M_PI *( t_[i] * (fzero + 0.5 * k * t_[i]))));
        }
//        stampaVetCmplx(c2 , snum);
        size_res = snum1 + snum2;
        double complex* y = NewVetCmplx(size_res);              ///Copy the first snum1 elements from c1
        for (int64_t i = 0; i < snum1 ; i++) {
            y[i] = c1[i];										///the first snum1 elements from c1 are copied, then from index 0 to snum1
        }
                                                                ///Copy the remaining elements of c2
        for (int64_t i = 0 ; i < snum2 ; i++) {
            y[snum1 + i] = c2[i];								///the next snum1+i elements from c2 are copied, from index snum1+1 to snum1+snum2
        }
//	stampaVetCmplx(y,size_res);
	ResVetCmplx resultchirp;
	resultchirp.size_chirp = size_res;
    resultchirp.vet = (double complex*)calloc(size_res, sizeof(double complex));
    for(int64_t j = 0 ; j < size_res ; j++){
        resultchirp.vet[j] = y[j];
    }

    free(t_);
    t_ = NULL;
	free(c1);
	c1 = NULL;
	free(c2);
	c2 = NULL;
	free(y);
	y = NULL;
	return resultchirp;
}

int16_t writeCmplxFile(double complex* vet, int64_t size, const char* filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening the file %s\n", filename);
        return -1;
    }
    int16_t v = 0;
    int16_t* IQ = malloc(2 * size * sizeof(int16_t));
    if (IQ == NULL) {
        printf("Error: memory allocation failed\n");
        fclose(file);
        return -1;
    }
    for (int64_t i = 0; i < size; i++) {
        // Scrivi la parte reale e immaginaria come byte nel file
        int16_t real = (int16_t) (creal(vet[i]) * 32767);
        int16_t imag = (int16_t) (cimag(vet[i]) * 32767);
        // Assegnare i valori reali e immaginari al vettore interleaved
        IQ[2 * i] = real;
        IQ[2 * i + 1] = imag;
//        printf("\n%"PRId64"] --> IQ =%"PRId16" %"PRId16"", i ,real,imag);
    }
        // Scrivere il vettore interleaved nel file
    v = fwrite(IQ, sizeof(int16_t), 2 * size, file);
    if (ferror(file)) {
        printf("\n%"PRId16"",v);
        printf("\nError writing the file\n");
    }

    free(IQ);
    fclose(file);
return v;
}

void readCmplxFile(double complex* vet, int64_t size, const char* filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        printf("Error opening the file %s\n", filename);
        return;
    }
    for (int64_t i = 0; i < size; i++) {
        /// Read the REAL PART and COMPLEX PART as a byte from the file
        float real, imag;
        fread(&real, sizeof(float), 1, file);
        fread(&imag, sizeof(float), 1, file);
        vet[i] = real + imag * I;
    }
    fclose(file);
}

double complex* read_file(char* filename,int64_t* size) {
    FILE* f = fopen(filename, "rb");
    if(f == NULL){
        printf("\nERROR: Failed to open the file: %s\n", filename);
        *size = 0;
        f = NULL;
        return NULL;
    }
    // Determine the size of the file
    fseek(f, 0, SEEK_END);
    int64_t file_size = ftell(f);
    rewind(f);

    // Check if the file size is compatible with complex double data
    if(file_size % (2 * sizeof(double)) != 0){
        printf("\nERROR: File to open is not compatible with double complex data\n");
        *size = 0;
        f = NULL;
        return NULL;  // Invalid file size
    }

    *size = file_size / (2.0 * sizeof(double));
    double* file_data = malloc(2 * (*size) * sizeof(double));

    if(file_data == NULL){
        printf("\nERROR: Failed to allocate memory in file reading!");
        f = NULL;
        *size = 0;
        return NULL;  // Memory allocation failed
    }
    // Read the file content
    fread(file_data, sizeof(double), 2 * (*size), f);
    // Convert the data back to double complex
    double complex* result = malloc(*size * sizeof(double complex));
    if (result == NULL) {
        printf("\nERROR: Failed to allocate memory in file reading!");
        f = NULL;
        free(file_data);
        file_data = NULL;
        *size = 0;
        return NULL;  // Memory allocation failed
    }

    for (int64_t i = 0; i < *size; i++) {
        result[i] = file_data[i] + I * file_data[*size + i];
    }
    fclose(f);
    f = NULL;
    free(file_data);
    file_data = NULL;
    return result;
}

ResFileCmplx readCmplx(const char* filePath) {
    FILE* file = fopen(filePath, "rb");
    if (file == NULL) {
        printf("\nError: Failed to open the path of file: %s\n", filePath);
        file = NULL;
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    int64_t dimTotFile = ftell(file);
    fseek(file, 0, SEEK_SET);
    int64_t numElementi = dimTotFile / sizeof(float complex);      //number of elements in file = total dimension / type of variable of the data(float complex)
//    printf("\ndimTot=%"PRId64"\tnumEle=%"PRId64"\tsizeof(floatcomplex)=%"PRId64"",dimTotFile , numElementi , sizeof(float complex));
    float complex* vettore = (float complex *)malloc(numElementi * sizeof(float complex));
    if (vettore == NULL) {
        printf("\nError: memory allocation failed\n");
        file = NULL;
        exit(1);
    }

    fread(vettore, sizeof(float complex), numElementi, file);                     /// Read the content of the file in array
/*    printf("Il vettore complesso contiene %"PRId64" elementi.\n", numElementi); ///Print number of elements and content of the array
    printf("Contenuto del vettore complesso:\n");
    for (int64_t i = 0; i < numElementi; i++) {
        printf("%.3f+ %.3fi\t", creal(vettore[i]), cimag(vettore[i]));
    }*/
    double complex* temp = (double complex *)malloc(numElementi * sizeof(double complex));
    if (temp == NULL) {
        printf("\nError: memory allocation failed\n");
        file = NULL;
        exit(1);
    }
    for(int64_t i = 0 ; i < numElementi ; i++ ){
           temp[i] = vettore[i];
//        printf("%.4f+%.4f*I  ",creal(temp[i]), cimag(temp[i]));
    }
    ResFileCmplx result;
    result.size_file = numElementi;
    result.vet = temp;
    fclose(file);
    file = NULL;
    free(vettore);
    vettore = NULL;
    return result;
}

ResFileDouble readDouble(const char* filePath){
FILE* file = fopen(filePath, "r");
    if(file == NULL){
        printf("\nERROR: Failed to open the path of file: %s\n", filePath);
        file = NULL;
        exit(1);
    }
    int64_t numElementi = 0;
    float valore;
    while(fscanf(file, "%f", &valore) == 1) {
         numElementi++;
    }
    float* vettore = (float*)malloc(numElementi * sizeof(float));
    if (vettore == NULL) {
        printf("\nError: memory allocation failed\n");
        file = NULL;
        exit(1);
    }
    double* temp = (double*)malloc(numElementi * sizeof(double));
    if (temp == NULL) {
        printf("\nError: memory allocation failed\n");
        file = NULL;
        exit(1);
    }
    fseek(file, 0, SEEK_SET);
    for(int64_t i = 0; i < numElementi ; i++){
        fscanf(file, "%f", &vettore[i]);
        temp[i] = vettore[i];
    }
    ResFileDouble result;
    result.size_file = numElementi;
    result.vet = temp;
    fclose(file);
    file = NULL;
    free(vettore);
    vettore = NULL;
    return result;
}

void ConvolveCmplx(double complex* input, int64_t inLength, double* h_real, double complex* h_imag, int64_t h_Length, double complex* output) {
    for (int64_t i = 0 ; i < inLength ; i++) {
        output[i] = 0.0 + I * (0.0);
        for (int64_t j = 0; j < h_Length ; j++) {
            if (i - j >= 0 && i - j < inLength) {
                output[i] += input[i - j] * (h_real[j] + I * h_imag[j]);
//                printf("\ninput=%f+%f*I\threal=%f\t+\thimag=%f*I", creal(input[i-j]) , cimag(input[i-j]), h_real[j] , h_imag[j]);
            }
        }
    }
}

int64_t* convertiInt10Bit(int64_t value){
int64_t* bitArray =(int64_t*)malloc(10 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 1023 || value < 0){
    printf("\nValue can't be converted with only 10bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 9; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[9 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt9Bit(int64_t value){
int64_t* bitArray =(int64_t*)malloc(9 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 511 || value < 0){
    printf("\nValue can't be converted with only 9bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 8; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[8 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt8Bit(int64_t value){
int64_t* bitArray =(int64_t*)malloc(8 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 255 || value < 0){
    printf("\nValue can't be converted with only 8bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 7; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[7 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt7Bit(int64_t value){
int64_t* bitArray = (int64_t*)malloc(7 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 127 || value < 0){
    printf("\nValue can't be converted with only 7bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 6; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[6 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt6Bit(int64_t value){
int64_t* bitArray = (int64_t*)malloc(6 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 63 || value < 0){
    printf("\nValue can't be converted with only 6bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 5; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[5 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt5Bit(int64_t value){
int64_t* bitArray = (int64_t*)malloc(5 * sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 31 || value < 0){
    printf("\nValue can't be converted with only 5bits!");
    bitArray = NULL;
    exit(1);
}
else{
    for (int64_t i = 4; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[4 - i] = bit;
    }
}
 return bitArray;
}

int64_t* convertiInt4Bit(int64_t value){
int64_t* bitArray =(int64_t*)calloc(4,sizeof(int64_t));
if (bitArray == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
}
if(value > 15 || value < 0){
    printf("\nValue can't be converted with only 4bits!");
    bitArray = NULL;
    exit(1);
}
else{
   for (int64_t i = 3; i >= 0 ; i--) {
        int64_t bit = ( value >> i) & 1;
        bitArray[3 - i] = bit;
    }
}
 return bitArray;
}

int64_t binToDec(int64_t* vetBin, int64_t length) {
    int64_t decimale = 0;
    for(int64_t i = 0 ; i < length ; i++) {
        decimale = decimale * 2 + vetBin[i];
    }
 return decimale;
}

double mod_double(double a , double b){
  if(b == 0)
    return a;
  double result = fmod(a, b);
    return ((result >= 0 && b > 0) || (a <= 0 && b < 0)) ? result : (result + b);
}

int64_t* bit_reduce(int64_t* words, int64_t* index, int64_t lengthwords, int64_t sizeindex){
 int64_t* b = calloc(lengthwords,sizeof(int64_t));
 if (b == NULL) {
    fprintf(stderr,"\nError: Memory allocation failed!");
    exit(1);
 }
 int64_t* bit_pos = (int64_t*)calloc(lengthwords, sizeof(int64_t));
 int64_t maskbitxor;
 for(int64_t i = 0 ; i < sizeindex ; i++ ){
//    printf("\nbitxor=");
    maskbitxor = 1 << index[i] ;
    for(int64_t j = 0; j < lengthwords ; j++){
        bit_pos[j] = (words[j] & maskbitxor) >> index[i];
//        printf("%"PRId64" ",bit_pos[j]);
    }
    if( i == 0){
        for(int64_t j = 0; j < lengthwords ; j++){
            b[j] = bit_pos[j];
        }
    }else{
        for(int64_t k = 0; k < lengthwords ; k++){
            b[k] = b[k] ^ bit_pos[k];
        }
    }
 }
 free(bit_pos);
 bit_pos = NULL;
 return b;
}

int64_t word_reduce(int64_t* words, int64_t lengthwords){
int64_t w;
w = words[0];
//printf("\nw=%"PRId64"",w);
for(int64_t i = 1 ; i < lengthwords ; i++){
    w = w | words[i];
//    printf("\nw=%"PRId64"",w);
}

 return w;
}

void info_Lora(LoRaPHY* Lora) {
	printf("\nInitialization of the struct LoRaPHY:");
	printf("\nrf_freq=%15"PRId64"    [Hz]\nsf=\t%15"PRId64"    --->Spreading factor(7,8,9,10,11,12)\nbw=\t%15"PRId64"    [Hz]    (125kHz 250kHz 500kHz)\nfs=\t%15"PRId64"    [Hz]", Lora->rf_freq, Lora->sf, Lora->bw, Lora->fs);
	printf("\ncr=\t%15"PRId64"    (1->4/5 2->4/6 3->4/7 4->4/8)    {If -1 (sentinel value), means the value is to be setted}\npayload_len=%11"PRId64"    --->Payload length        \t{If -1 (sentinel value), means the value is to be setted}\nhas_header=%12"PRId64"    {1:Enabled/0:Disabled}\ncrc=\t%15"PRId64"    {1:Enabled/0:Disabled}",Lora->cr,Lora->payload_len,Lora->has_header,Lora->crc);
	printf("\nldr=\t%15"PRId64"    {1:Enabled/0:Disabled}\npreamble_len=%10"PRId64"    --->Preamble length\nsample_num=%12"PRId64"    --->Sample points per symbol",Lora->ldr,Lora->preamble_len,Lora->sample_num);
	printf("\nbin_num=%15"PRId64"\nzero_padding_ratio=%4"PRId64"    --->Fft padding ratio\nfft_len=%15"PRId64"    --->Fft length",Lora->bin_num,Lora->zero_padding_ratio,Lora->fft_len);
	printf("\npreamble_bin=%10"PRId64"    --->Preamble bin    \t\t{If -1 (sentinel value), means the value is to be setted}\ncfo=%19.4f    [Hz]\nis_debug=%14"PRId64"    {1:Enabled/0:Disabled}\nhamming_decoding_en=%3"PRId64"    {1:Enabled/0:Disabled}\n",Lora->preamble_bin,Lora->cfo,Lora->is_debug,Lora->hamming_decoding_en);
    printf("\nHeader Checksum Matrix:\n");
	for (int64_t i = 0; i < 5; i++) {
        for (int64_t j = 0; j < 12; j++) {
            printf("%"PRId64" ",Lora->header_checksum_matrix[i][j]);
        }printf("\n");
    }
}

void freeLora(LoRaPHY* Lora){
	if(Lora->downchirp != NULL){
		free(Lora->downchirp);
		Lora->downchirp = NULL;
	}

	if(Lora->upchirp != NULL){
		free(Lora->upchirp);
		Lora->upchirp = NULL;
	}
}
