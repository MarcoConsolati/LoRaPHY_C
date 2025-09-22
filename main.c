#include "lora_phy.h"
#include "crc_ccitt.h"
#include "utility.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>
#include <complex.h>
#include <inttypes.h>
#include <string.h>


//FUNCTION THAT DECHIRP THE SIGNAL SIG. CFILE FROM THE INITIAL INDEX  X' OF A SYMBOL,
//MAKES THE FFT OF A PORTION OF THE SIGNAL AND THEN FINDS THE MAXIMUM VALUE (ABSOLUTE VALUE)
//ARGUMENTS: LORA, X (integer index of a symbol), IS_UP (true or false decide if downchirped or upchirped)
//OUTPUT: RESULTPK(is a struct described in lora_phy.h header file)
Peak dechirp(LoRaPHY* Lora, int64_t  x , bool is_up){
	int64_t N = Lora->fft_len;
    double complex* c = (double complex*)malloc(Lora->size_upchirp*sizeof(double complex));
    if (c == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
	double complex* input = (double complex*)calloc(N,sizeof(double complex));
	if (input == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
	double complex* output = (double complex*)calloc(N,sizeof(double complex));
	if (output == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
	double* ft_abs = (double*)calloc((N / 2.0),sizeof(double));
    if (ft_abs == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
	if(!is_up){
        for(int64_t i = 0; i < Lora->size_upchirp ; i++){
            c[i] = Lora->upchirp[i];
        }
//      stampaVetCmplx(c ,Lora->size_upchirp );
//      printf("\nsizeUp=%"PRId64"",Lora->size_upchirp);
	}
    else {
        for(int64_t i = 0; i < Lora->size_downchirp ; i++){
            c[i] = Lora->downchirp[i];
        }
//		stampaVetCmplx(c ,Lora->size_downchirp );
//		printf("\nsizeDown=%"PRId64"",Lora->size_downchirp);
        }
	for (int64_t i = 0 ; i <= Lora->sample_num - 1 ; i++) {         ///Array truncated from index X to X+sample_num+1
        if ((i + x) > Lora->sig_size - 1){
                input[i] = 0.0;
        } else {
            input[i] = Lora->sig[i + x] * c[i];
        }
//      printf("sig elem pos %"PRId64": (%f, %f) * (%f,%f) = (%f, %f)\n", i, creal(Lora->sig[i + x]), cimag(Lora->sig[i + x]), creal(c[i]), cimag(c[i]), creal(input[i]), cimag(input[i]));                                                                               //gli indici aggiunti sul finale si paddano con numeri complessi nulli
//      printf("\ninput[%"PRId64"]=%f +(%f*i)\t",i,creal(input[i]),cimag(input[i]));
    }
//    stampaVetCmplx(input, Lora->sample_num - 1);
    fftw_plan myplan;                                                                   //Il k-esimo elemento del vettore Output[k],corrisponde
 	myplan = fftw_plan_dft_1d( N , input , output  , FFTW_FORWARD , FFTW_ESTIMATE);     //alla fft discreta su N punti del vettore troncato trunc_input[i],
 	fftw_execute(myplan);                                                               //(dall'indice x all'indice x+sample_num-1),rappresentando la
    Peak resultPk;                                                                      //componente frequenziale di freq:k/T e T=periodo di campionamento
 	double pk_max = 0.0;
    for(int64_t i = 0; i < N / 2 ; i++){
        ft_abs[i] =  cabs(output[i]) + cabs(output[i + ( N / 2 )]);                     //Sommo il valore assoluto della prima metà
        if(ft_abs[i] > pk_max){
            pk_max = ft_abs[i];
            resultPk.index = i;
            resultPk.value = ft_abs[i];
        }
    }
//    stampaVetCmplx(input, N);
//    stampaVetCmplx(output, N);
//    stampaVetDouble(ft_abs,N/2);
    fftw_destroy_plan(myplan);
    free(c);
    c = NULL;
    free(input);
    input = NULL;
    free(output);
    output = NULL;
    free(ft_abs);
    ft_abs = NULL;
    fftw_cleanup();
    return resultPk;
}


//FUNCTION DETECTNEW LOOKING FOR THE PREAMBLE WITHIN THE SIGNAL SIG. CFILE STARTING FROM THE INITIAL INDEX startIndx.
//ARGUMENTS: LORA, STARTINDX.
//OUTPUT: X (INDEX X IF PREAMBLE FOUND, ELSE RETURN X=-1);
int64_t detectnew(LoRaPHY* Lora, int64_t startIndx){
    int64_t ii = startIndx;
    int64_t  max_bin_list_size = Lora->sig_size - (Lora->sample_num * Lora->preamble_len);
    int64_t* pk_bin_list = (int64_t *)malloc(max_bin_list_size * sizeof(int64_t));
    if (pk_bin_list == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    int64_t count_list_element = 0;
    Peak pkzero;
    int64_t x = 0;
    int64_t bindiff = -1;
    while( ii <= max_bin_list_size ){
        if( count_list_element == Lora->preamble_len - 1){
            x = ii - round((double)(pk_bin_list[count_list_element - 1])/ (Lora->zero_padding_ratio * 2));
            for(int64_t i = 0; i < count_list_element ; i++){
//                printf("\npkbinlist[%"PRId64"]=%"PRId64"",i,pk_bin_list[i]);
            }
            free(pk_bin_list);
            return x;
        }
        pkzero = dechirp(Lora, ii , true);
        printf("\ndetectnew\nii=%"PRId64"\npkzero=%.4f from the index=%5"PRId64"",ii, pkzero.value,pkzero.index);
        if( count_list_element != 0 ){
            bindiff = (pk_bin_list[ count_list_element - 1] - pkzero.index) - ((pk_bin_list[ count_list_element - 1] - pkzero.index) / Lora->bin_num) * Lora->bin_num;
            bindiff = mod_double((pk_bin_list[ count_list_element - 1] - pkzero.index) , Lora->bin_num);
            if( bindiff > (Lora->bin_num / 2)){
                bindiff = Lora->bin_num - bindiff;
            }
            if( bindiff <= Lora->zero_padding_ratio){
                pk_bin_list[count_list_element] = pkzero.index;
                count_list_element++;
                for(int64_t i=0; i < count_list_element ; i++){
//                    printf("\npkbinlist[%"PRId64"]=%"PRId64"",i,pk_bin_list[i]);
                }
            }else{
                pk_bin_list[0] = pkzero.index;
                count_list_element = 1;
            }
        }else{
            pk_bin_list[0] = pkzero.index;
            count_list_element = 1;
            for(int64_t i=0; i < count_list_element ; i++){
//                printf("\npkbinlist[%"PRId64"]=%"PRId64"",i,pk_bin_list[i]);
            }
        }
        ii += Lora->sample_num;
    }
    free(pk_bin_list);
    pk_bin_list = NULL;
    return -1;
}


//FUNCTION SYNCNEW TO SYNCHRONIZE SIG SIGNAL PACKETS FROM THE INITIAL INDEX X.
//ARGUMENTS: LORA , X(integer value of a index);
//OUTPUT: X_sync (integer value of the index after up-down alignment)
int64_t syncnew(LoRaPHY* Lora, int64_t x){
    // Packet synchronization find downchirp
    bool found = false;
//    printf("\nsigsize=%"PRId64" samplenum=%"PRId64" index_x=%"PRId64"", Lora->sig_size, Lora->sample_num , x);
    Peak up_peak;
    Peak down_peak;
    while(x < Lora->sig_size - Lora->sample_num){
        up_peak = dechirp(Lora, x , true);      // if arguments of dechirp < 2 ---> c=downchirp;
        down_peak = dechirp(Lora, x , false);   // if arguments of dechirp=2 and is_up=false ---> c=upchirp;
        if (fabs(down_peak.value) > fabs(up_peak.value)) {
//            printf("\ndownpeak=%f maggiore di uppeak=%f", down_peak.value , up_peak.value);
            // downchirp detected
            found = true;
        }
        x += Lora->sample_num;
//        printf("\nindex_x=%"PRId64"",x);
        if (found) {
            break;
        }
    }
    if (!found) {
        return -1;  // Return type is int64_t, so returning -1 to indicate no synchronization found
    }
    // Up-Down Alignment
    // NOTE preamble_len >= 6
    // NOTE there are two NETID chirps between preamble and SFD
    // NOTE `detect` has already shifted the up peak to position 0
    Peak pkd = dechirp(Lora, x, false);
    int64_t to;
    if (pkd.index > (Lora->bin_num / 2.0)) {
        to = round( (double)(pkd.index - 1 - Lora->bin_num) / Lora->zero_padding_ratio);
    } else {
        to = round( (double)(pkd.index - 1) / Lora->zero_padding_ratio);
    }
    x += to;

    // set preamble reference bin for CFO elimination
    Peak pku = dechirp(Lora, (x - 4 * Lora->sample_num) , true);  // if arguments of dechirp < 2 ---> c=downchirp;
    Lora->preamble_bin = pku.index;
//    printf("\npkuindex=%"PRId64"",pku.index);
    double temp_cfo;
    if (Lora->preamble_bin > (Lora->bin_num / 2.0)) {
        temp_cfo = fabs(Lora->preamble_bin - Lora->bin_num) * (Lora->bw /(double)Lora->bin_num);
                                                                            //In the original matlab version the index found was 1
    } else {                                                                        //in C version i add fabs because if index
        temp_cfo =  (Lora->preamble_bin + 1) * (Lora->bw / (double)Lora->bin_num);  //was 0 we've got negative cfo!
                                                                                    //In the original matlab version the index found was 1
                                                                                    //in C version is 0 so i had to add one!
    }
    // set x to the start of data symbols
    Peak pku2;
    Peak pkd2;
    pku2 = dechirp(Lora, x - Lora->sample_num , true);       // if arguments of dechirp < 2 ---> c=downchirp;
    pkd2 = dechirp(Lora, x - Lora->sample_num , false);      // if arguments of dechirp=2 and is_up=false ---> c=upchirp;
    int64_t x_sync;
    if ( fabs(pku2.value) > fabs(pkd2.value) ) {
        // current symbol is the first downchirp
        x_sync = x + round(2.25 * Lora->sample_num);
    } else {
        // current symbol is the second downchirp
        x_sync = x + round(1.25 * Lora->sample_num);
    }
    Lora->cfo = temp_cfo;

 return x_sync;
}


//FUNCTION PARSE_HEADER SET PARSE HEADER AND OTHER PARAMETERS.
//ARGUMENTS: LORA, DATA(Array of 8th elements containing the symbols of the header packet)
//OUTPUT: IS_VALID 'true'(1) if is a valid header 'false'(0) otherwise.
int64_t parse_header(LoRaPHY* Lora , int64_t* data , int64_t lengthdata){
int64_t is_valid = -1;
int64_t* symbols = NULL;
int64_t* symbols_g = NULL;
int64_t* codewords = NULL;
int64_t* nibbles = (int64_t*)malloc((Lora->sf - 2) * sizeof(int64_t));
if (nibbles == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
symbols = dynamic_compensation(Lora , data , lengthdata);
symbols_g = gray_coding(Lora , symbols , lengthdata);
codewords = diag_deinterleave(Lora , symbols_g , 8 , Lora->sf - 2);
printf("\nSymbols after parse header dyn_comp:");
for(int64_t i = 0; i < 8 ; i++){
    printf(" %"PRId64"", symbols[i]);
}
printf("\nSymbols_g after gray coding:");
for(int64_t i = 0; i < 8 ; i++){
    printf(" %"PRId64"", symbols_g[i]);
}
printf("\nCodewords:");
for(int64_t i = 0; i < Lora->sf - 2 ; i++){
    printf(" %"PRId64"", codewords[i]);
}
nibbles = hamming_decode(Lora , codewords , Lora->sf - 2 , 8);
Lora->payload_len = (nibbles[0]*16 + nibbles[1]);
Lora->crc = (nibbles[2] & 1);
Lora->cr = (nibbles[2] >> 1);
int64_t* header_checksum = (int64_t*)malloc(5*sizeof(int64_t));
if (header_checksum == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t* header_nibbles = (int64_t*)malloc(12 * sizeof(int64_t));
int64_t* nibbles_bits = (int64_t*)malloc(5 * sizeof(int64_t));;
int64_t* bitarray = convertiInt4Bit(nibbles[4]);
printf("\nHeader_checksum bits:");
for(int64_t i = 0 , j = 0 ; i < 5  ; i++ ){
    if( i == 0){
        header_checksum[i] = nibbles[3] & 1 ;
    }else{
        header_checksum[i] = bitarray[j];
        j++;
    }
    printf(" %"PRId64"", header_checksum[i]);
}
int64_t* nibbles_bin0 = convertiInt4Bit(nibbles[0]);
int64_t* nibbles_bin1 = convertiInt4Bit(nibbles[1]);
int64_t* nibbles_bin2 = convertiInt4Bit(nibbles[2]);
int64_t vet_bin_nibbles[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
printf("\nnibbles[0]: ");
for(int64_t i = 0; i < 4 ; i++){
    vet_bin_nibbles[i] = nibbles_bin0[i];
    printf("%"PRId64" ",nibbles_bin0[i]);
}
printf("\nnibbles[1]: ");
for(int64_t i = 0; i < 4 ; i++){
    vet_bin_nibbles[i+4] = nibbles_bin1[i];
    printf("%"PRId64" ",nibbles_bin1[i]);
}
printf("\nnibbles[2]: ");
for(int64_t i = 0; i < 4 ; i++){
    vet_bin_nibbles[i+8] = nibbles_bin2[i];
    printf("%"PRId64" ",nibbles_bin2[i]);
}
printf("\n\nHeader Checksum Matrix: \n");
for(int64_t i = 0; i < 5; i++) {
    for(int64_t j = 0; j < 12; j++) {
        printf("%"PRId64" ", Lora->header_checksum_matrix[i][j]);
    }
    printf("\n");
}
printf("\n\n(*)");
printf("\n\nVet_Nibbles: ");
for(int64_t i = 0; i < 12 ; i++){
    printf("%"PRId64" ",vet_bin_nibbles[i]);
}
printf("\n\n = Result\t\t\t\t\t\tnibbles_bits\n");
for(int64_t i = 0; i < 5; i++){
    int64_t total = 0;
    for(int64_t j = 0; j < 12; j++) {
        header_nibbles[j] = Lora->header_checksum_matrix[i][j] * vet_bin_nibbles[j];
        printf("%"PRId64" ", header_nibbles[j]);
        if(header_nibbles[j] % 2 != 0)
            total++;
    }
    if(total % 2 == 0)
        nibbles_bits[i] = 0;
    else
        nibbles_bits[i] = 1;
    printf("    (Binary Sum of the %"PRId64" row)------->%"PRId64"", i+1 ,nibbles_bits[i]);
    printf("\n");
}
for(int64_t i = 0; i < 5 ; i++){
    if(header_checksum[i] != nibbles_bits[i]){
        printf("\nInvalid Header Checksum!");
        is_valid = 0;
    }else{
        is_valid = 1;
    }
}

free(symbols);
symbols = NULL;
free(symbols_g);
symbols_g = NULL;
free(codewords);
codewords = NULL;
free(nibbles);
nibbles = NULL;
free(bitarray);
bitarray = NULL;
free(header_checksum);
header_checksum = NULL;
free(nibbles_bin0);
nibbles_bin0 = NULL;
free(nibbles_bin1);
nibbles_bin1 = NULL;
free(nibbles_bin2);
nibbles_bin2 = NULL;
free(nibbles_bits);
nibbles_bits = NULL;
free(header_nibbles);
header_nibbles = NULL;
return is_valid;
}


double calc_sym_num(LoRaPHY* Lora, int64_t plen){
double sym_num = 0.0;
if(Lora->cr == -1){
    printf("\n Code rate not setted!");
    exit(1);
}else{
//printf("\nLora->ldr = %"PRId64"", Lora->ldr);
//printf("\nplen = %"PRId64"" ,plen);
double num = (double)(2 * plen - Lora->sf + 7 + 4 * Lora->crc - 5 *(1 - Lora->has_header));
double den = (double)(Lora->sf - 2 * Lora->ldr);
//printf("\nnum=%f",num);
//printf("\nden=%f",den);
sym_num = (double) 8.0 + max_calc((4+Lora->cr)*ceil(num / den) , 0.0);
}
 return sym_num;
}

//FUNCTION CALC_PAYLOAD_LEN CALCULATE THE LENGTH OF THE PAYLOAD PACKET.
//ARGUMENTS: LORA, SLEN(number of symbols), NO_REDUNDANT_BYTES(boolean 'TRUE' or 'FALSE')
//OUTPUT: PLEN(payload packet length)
// plen_d possibly has fractional part 0.5, which means
// there would be 0.5 uncontrollable redundant byte in a packet.
// The 0.5 byte results in unexpected symbols when called by
// function `symbols_to_bytes`. To make all specified symbols
// controllable, we use `ceil` instead of `floor` when
// no_redundant_bytes is true.
int64_t calc_payload_len(LoRaPHY* Lora, int64_t slen, bool no_redundant_bytes){
int64_t plen;
double temp;
//Lora->cr = 4; //Lora->cr dovrebbe essere 1 - 2 - 3 - 4 in base al code rate!
//printf("\ncr=%"PRId64"",Lora->cr);
if(Lora->cr == -1){
    printf("\nCode Rate not setted!");
    exit(1);
}else if(Lora->cr != 1 && Lora->cr != 2 && Lora->cr != 3 && Lora->cr != 4){
    printf("\nWARNING! Code rate setted not correctly!");
    exit(1);
}
else{
double plen_d = (((double)Lora->sf - 2) / 2.0) - 2.5 * (double)Lora->has_header + (((double)Lora->sf - (double)Lora->ldr * 2) / 2.0) * ceil(((double)slen - 8.0) / (double)(Lora->cr + 4));
//printf("\nLora->has_header=%"PRId64" Lora->ldr=%"PRId64" Lora->cr=%"PRId64"",Lora->has_header,Lora->ldr,Lora->cr);
printf("\nplen_d = %f",plen_d);
if( no_redundant_bytes == true )
    temp = ceil(plen_d);
else
    temp = floor(plen_d);
plen = (int64_t)temp;
printf("\t--->plen = %"PRId64"",plen);
}
 return plen;
}


//FUNCTION DYNAMIC_COMPESATION TO COMPENSATE BIN DRIFT CAUSED BY FREQUENCY SAMPLING OFFSET.
//ARGUMENTS: LORA, DATA(array containing the symbols with BIN DRIFT)
//OUTPUT: SYMBOLS(array containing ricalibrated symbols)
int64_t* dynamic_compensation(LoRaPHY* Lora, int64_t* data , int64_t length_data){
double* sfo_drift = newVetRange(1 , length_data + 1);
int64_t* symbols = (int64_t*)malloc(length_data * sizeof(int64_t));
if (symbols == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
printf("\n");
for(int64_t i = 0 ; i < length_data ; i++){
    sfo_drift[i] = (sfo_drift[i] + 1) * (pow(2, Lora->sf) * Lora->cfo / (double)Lora->rf_freq);
//    printf("\nsfo_drift[%"PRId64"]=%f", i , sfo_drift[i]);
    symbols[i] = mod_double(round(data[i] - sfo_drift[i]) , pow(2, Lora->sf));     //Caso in cui LDR=0 l'output sarà questo di symbols
//    printf("\tmod(%4"PRId64" - %7f , %"PRId64")\t-->symbols=%"PRId64"",data[i] , sfo_drift[i] , result , symbols[i]);
//    printf("\nsymbols=%"PRId64"",symbols[i]);
}
double* v = (double*) malloc(length_data * sizeof(double));
if (v == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
if(Lora->ldr){
    double bin_delta = 0.0;
    double binoffset = 0.0;
    double v_last = 1.0;
    for(int64_t i = 0; i < length_data ; i++){
        v[i] = symbols[i];
        bin_delta = mod_double(round(v[i] - v_last) , 4);
//        printf("\nmod(%4.4f - %4.4f , %d)\t-->bin_delta=%.4f", v[i] , v_last , 4 , bin_delta);
        if(bin_delta < 2.0){
            binoffset -= bin_delta;
        }
        else{
            binoffset -= bin_delta + 4.0;
        }
        v_last = v[i];
        symbols[i] = mod_double(round(v[i] + binoffset) , pow(2, Lora->sf));
//        printf("\n");                                 //Caso in cui LDR=1 l'output sarà questo di symbols
//        printf("\nsymbols=%"PRId64"",symbols[i]);     //dalla FUNCTION di Dynamic_compensation
    }
}

free(sfo_drift);
sfo_drift = NULL;
free(v);
v = NULL;
return symbols;
}


//FUNCTION GRAY_CODING DECODES THE SYMBOLS PASSED TO IT VIA A VECTOR AND RECALIBRATES THEM.
//ARGUMENTS: LORA, DIN(pointer to a vector containing the symbols to be decoded)
//OUTPUT: *SYMBOLS(vector containing the symbols after calibration)
int64_t* gray_coding(LoRaPHY* Lora, int64_t* din , int64_t length_din){
int64_t* symbols_g = (int64_t*) malloc (length_din * sizeof(int64_t));
if (symbols_g == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t* s = (int64_t*) malloc (length_din * sizeof(int64_t));
if (s == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t* symbols = (int64_t*)malloc(length_din * sizeof(int64_t));
if (symbols == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
//printf("\nlength of din=%"PRId64"", length_din);
printf("\n");
for(int64_t i = 0; i < length_din ; i++){
    if( i < 8){
        symbols_g[i] = floor(din[i] / 4.0);
    }
    else{
        if( Lora->ldr){
            symbols_g[i] = floor(din[i] / 4.0);
        }
        else{
            symbols_g[i] = mod_double((din[i] - 1) , pow(2,Lora->sf));
        }
    }
    printf("\ns[%2"PRIu64"]=%5"PRIu64"", i , symbols_g[i]);
    s[i] =  symbols_g[i];
    symbols[i] = s[i] ^ (s[i] >> 1);
    printf("\tafter Gray coding--->\tsymbols[%2"PRIu64"]=%5"PRIu64"", i , symbols[i]);
}
free(symbols_g);
symbols_g = NULL;
free(s);
s = NULL;
return symbols;
}


//FUNCTION GRAY_DECODING DECODE SYMBOLS PASSED TO IT USING A VECTOR.
//ARGUMENTS: LORA, SYMBOLS_I(pointer to a vector containing the symbols to be decoded)
//OUTPUT: *SYMBOLS(vector containing the new symbols )
int64_t* gray_decoding(LoRaPHY* Lora, int64_t* symbols_i , int64_t lengthsymbols_i){
int64_t* symbols = (int64_t*)malloc(lengthsymbols_i * sizeof(int64_t));
if (symbols == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
printf("\n");
for(int64_t i = 0; i < lengthsymbols_i ; i++){
    int64_t num = symbols_i[i];
    int64_t mask = symbols_i[i] >> 1;
    while( mask!= 0){
        num ^= mask;
        mask >>= 1;
    }
//    printf("\nnum after while: %"PRId64"\nLora->ldr=%"PRId64"",num,Lora->ldr);
    if( i < 8 || Lora->ldr ){
        symbols[i] = mod_double((num * 4 + 1) , pow(2,Lora->sf));
    }else{
        symbols[i] = mod_double((num + 1) , pow(2,Lora->sf));
    }
    printf("\nsymbols_i[%2"PRIu64"]= %3"PRIu64"\tafter Gray decoding--->\tsymbols[%2"PRIu64"]= %3"PRIu64"", i , symbols_i[i], i, symbols[i]);
}

 return symbols;
}


//FUNCTION DIAG_DEINTERLEAVING REGROUP SYMBOLS IN THEIR ORIGINAL ORDER BEFORE INTERLEAVING(after code gray).
//ARGUMENTS: LORA, SYMBOLS(pointer to a vector containing the symbols after gray coding), PPM(parity bit size)
//OUTPUT: CODEWORDS(symbols updated after deinterterleaving)
int64_t* diag_deinterleave(LoRaPHY* Lora, int64_t* symbols, int64_t lengthsymbol, int64_t ppm){
int64_t* codewords = (int64_t*) malloc(ppm * sizeof(int64_t));
if(codewords == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t** matrice = (int64_t**)malloc(lengthsymbol * sizeof(int64_t*)); // Allocazione dinamica delle righe, ogni riga rappresenta un simbolo
if(matrice == NULL) {
    printf("\nRows memory allocation error\n!");
    exit(1);
    }
for(int64_t i = 0 ; i < lengthsymbol ; i++) {                          // Allocazione dinamica delle colonne per ogni riga, ogni colonna rappresenta bits di ridondanza 'ppm'(es:code rate=4/5 -->rdd=5bits)
    matrice[i] = (int64_t*)malloc(ppm * sizeof(int64_t));
    if (matrice[i] == NULL) {
        printf("\nColumns memory allocation error\n!");
        exit(1);
    }
}
int64_t temp = 0;
int64_t tempMat = 0;
switch( ppm ){
    case 5:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){          //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt5Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 5 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){              //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];                     //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){          //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    case 6:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){         //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt6Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 6 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){             //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];                    //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){         //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    case 7:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){     //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt7Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 7 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];                //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){         //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    case 8:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){     //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt8Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 8 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];    //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    case 9:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){     //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt9Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 9 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];    //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    case 10:
        printf("\n");
        for(int64_t i = 0; i < lengthsymbol ; i++){     //seleziona le colonne
            temp = symbols[i];
            int64_t* bits = convertiInt10Bit(temp);
            printf("\nSymbols[%"PRId64"]=%"PRId64"\t---->converted in 10 bit: ",i,symbols[i]);
            for( int64_t j = 0; j < ppm ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[i][j] = bits[j];    //bit ordinati in modo L-MSB, scrivo le symbol in bit per riga, dalla riga in cima al fondo
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary symbols: ");
        printf("\n(symbols bit are inserted from the first row to the last,\nfrom left to right(MSB-->LSB)");
        for(int64_t i = 0; i < lengthsymbol ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < ppm ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < lengthsymbol; i++){         //shifto a sinistra ogni elemento di ogni riga della matrice di posizione
            int64_t posizioni_shift = i;                    //pari a n°riga - 1
            for (int64_t k = 0; k < posizioni_shift; k++) {
                tempMat = matrice[i][0];
                for (int64_t j = 0; j < ppm - 1; j++) {
                    matrice[i][j] = matrice[i][j + 1];
                }
                matrice[i][ppm - 1] = tempMat;
            }
        }
        break;
    }
printf("\nLeft shifted Matrix in binary:\n(every bit of the row is shifted of a position = #row - 1)");
for(int64_t i = 0; i < lengthsymbol ; i++){             //stampa la matrice
    printf("\n");
    for( int64_t j = 0; j < ppm ; j++){
        printf("%"PRId64" ", matrice[i][j]);
    }
}
                                                        //Le codewords non sono altro che le colonne binarie convertite dall'ultima alla prima
for(int64_t j = ppm - 1 , k = 0; j >= 0 ; j-- , k++){
    int64_t dec = 0;
    for(int64_t i = lengthsymbol - 1 ; i >= 0 ; i--){
        dec = dec * 2 + matrice[i][j];
    }
    codewords[k] = dec;
    printf("\n-->codewords[%"PRId64"]= %"PRId64"", k , codewords[k]);
}
for(int64_t i = 0; i < lengthsymbol; i++){
    free(matrice[i]);
}
free(matrice);
matrice = NULL;
 return codewords;
}


//FUNCTION DIAG_INTERLEAVE
//ARGUMENTS: LORA, CODEWORDS(pointer to a vector containing the words), LENGTHWORDS, RDD(number of redundancy bits)
//OUTPUT: SYMBOLS_I(pointer to a vector containing symbols after interleaving the words)
int64_t* diag_interleave(LoRaPHY* Lora, int64_t* codewords, int64_t lengthwords, int64_t rdd){
int64_t* symbols_i =(int64_t*)malloc(rdd * sizeof(int64_t));
if (symbols_i == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t** matrice = (int64_t**)malloc(rdd * sizeof(int64_t*));  //Allocazione dinamica delle righe in base ai bits di ridondanza(es:code rate=4/5 -->rdd=5bits)
if (matrice == NULL) {
    printf("\nRows memory allocation error!\n");
    exit(1);
    }
for (int64_t i = 0; i < rdd; i++) {                     // Allocazione dinamica delle colonne per ogni riga, ogni colonna rappresenta una codewords
    matrice[i] = (int64_t*)malloc(lengthwords * sizeof(int64_t));
    if (matrice[i] == NULL) {
        printf("\nColumns memory allocation error!");
        exit(1);
    }
}
int64_t temp = 0;
switch( rdd ){
    case 5:
        printf("\n");
        for(int64_t i = 0; i < lengthwords ; i++){      //seleziona le colonne
            temp = codewords[i];
            printf("\nCodeword[%"PRId64"]=%"PRId64"\t---->converted in 5bit: ",i,temp);
            int64_t* bits = convertiInt5Bit(temp);
            for( int64_t j = 0; j < rdd ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[j][lengthwords - i - 1] = bits[rdd - j - 1];//bit ordinati in modo R-MSB, scrivo le codewords dall'ultima colonna alla prima
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary codewords: ");
        printf("\n(Codewords inserted from the last column to the first,\nfrom top to bottom (LSB--->MSB)");
        for(int64_t i = 0; i < rdd ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < lengthwords ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < rdd; i++){              //shifto circolarmente a destra ogni bit di ogni riga pari ad una posizione (n°riga-1)
            int64_t posizioni_shift = i;
            for (int64_t k = 0; k < posizioni_shift; k++) {
                temp = matrice[i][lengthwords - 1];
            for (int64_t j = lengthwords - 1; j > 0; j--) {
                matrice[i][j] = matrice[i][j - 1];
            }
                matrice[i][0] = temp;
            }
        }
        break;
    case 6:
        printf("\n");
        for(int64_t i = 0; i < lengthwords ; i++){      //seleziona le colonne
            temp = codewords[i];
            printf("\nCodeword[%"PRId64"]=%"PRId64"\t---->converted in 6bit: ",i,temp);
            int64_t* bits = convertiInt6Bit(temp);
            for( int64_t j = 0; j < rdd ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[j][lengthwords - i - 1] = bits[rdd - j - 1];//bit ordinati in modo R-MSB, scrivo le codewords dall'ultima colonna alla prima
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary codewords: ");
        printf("\n(Codewords inserted from the last column to the first,\nfrom top to bottom (LSB--->MSB)");
        for(int64_t i = 0; i < rdd ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < lengthwords ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < rdd; i++){              //shifto circolarmente a destra ogni bit di ogni riga pari ad una posizione (n°riga-1)
            int64_t posizioni_shift = i;
            for (int64_t k = 0; k < posizioni_shift; k++) {
                temp = matrice[i][lengthwords - 1];
            for (int64_t j = lengthwords - 1; j > 0; j--) {
                matrice[i][j] = matrice[i][j - 1];
            }
                matrice[i][0] = temp;
            }
        }
        break;
    case 7:
        printf("\n");
        for(int64_t i = 0; i < lengthwords ; i++){      //seleziona le colonne
            temp = codewords[i];
            printf("\nCodeword[%"PRId64"]=%"PRId64"\t---->converted in 7bit: ",i,temp);
            int64_t* bits = convertiInt7Bit(temp);
            for( int64_t j = 0; j < rdd ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[j][lengthwords - i - 1] = bits[rdd - j - 1];//bit ordinati in modo R-MSB, scrivo le codewords dall'ultima colonna alla prima
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary codewords: ");
        printf("\n(Codewords inserted from the last column to the first,\nfrom top to bottom (LSB--->MSB)");
        for(int64_t i = 0; i < rdd ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < lengthwords ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < rdd; i++){              //shifto circolarmente a destra ogni bit di ogni riga pari ad una posizione (n°riga-1)
            int64_t posizioni_shift = i;
            for (int64_t k = 0; k < posizioni_shift; k++) {
                temp = matrice[i][lengthwords - 1];
            for (int64_t j = lengthwords - 1; j > 0; j--) {
                matrice[i][j] = matrice[i][j - 1];
            }
                matrice[i][0] = temp;
            }
        }
        break;
    case 8:
        printf("\n");
        for(int64_t i = 0; i < lengthwords ; i++){     //seleziona le colonne
            temp = codewords[i];
            printf("\nCodeword[%"PRId64"]=%"PRId64"\t---->converted in 8bit: ",i,temp);
            int64_t* bits = convertiInt8Bit(temp);
            for( int64_t j = 0; j < rdd ; j++){         //seleziona le righe
                printf("%"PRId64" ",bits[j]);
                matrice[j][lengthwords - i - 1] = bits[rdd - j - 1];    //bit ordinati in modo R-MSB, scrivo le codewords dall'ultima colonna alla prima
            }
            free(bits);
            bits = NULL;
        }
        printf("\nMatrix with binary codewords: ");
        printf("\n(Codewords inserted from the last column to the first,\nfrom top to bottom (LSB--->MSB)");
        for(int64_t i = 0; i < rdd ; i++){              //stampa la matrice
            printf("\n");
            for( int64_t j = 0; j < lengthwords ; j++){
                printf("%"PRId64" ", matrice[i][j]);
            }
        }
        for (int64_t i = 1; i < rdd; i++){              //shifto circolarmente a destra ogni bit di ogni riga pari ad una posizione (n°riga-1)
            int64_t posizioni_shift = i;
            for (int64_t k = 0; k < posizioni_shift; k++) {
                temp = matrice[i][lengthwords - 1];
            for (int64_t j = lengthwords - 1; j > 0; j--) {
                matrice[i][j] = matrice[i][j - 1];
            }
                matrice[i][0] = temp;
            }
        }
        break;
    }
printf("\nRight shifted Matrix in binary:\n(every bit of the row is shifted of a position = #row - 1)");
for(int64_t i = 0; i < rdd ; i++){              //stampa la matrice
    printf("\n");
    for( int64_t j = 0; j < lengthwords ; j++){
        printf("%"PRId64" ", matrice[i][j]);
    }
    symbols_i[i] = binToDec( matrice[i] , lengthwords);
    printf("\t-------->\tsymbol_i[%"PRId64"]=%"PRId64"", i , symbols_i[i]);
}

for (int64_t i = 0; i < rdd; i++) {
        free(matrice[i]);
    }
    free(matrice);
matrice = NULL;
 return symbols_i;
}


//FUNCTION HAMMING_DECODE implements Hamming decoding to correct any errors in a set of codewords.
//The function takes the codewords as input after the interleaving (codewords) and the number of redundancy bits (rdd).
//Returns nibbles after Hamming decoding.
//ARGUMENTS: LORA, CODEWORDS(passed by the DIAG_DEINTERLEAVE function) , RDD(number of redundancy bits)
//OUTPUT: NIBBLES (values ​​assumed by the codewords after decoding)
int64_t* hamming_decode(LoRaPHY* Lora , int64_t* words , int64_t lengthwords , int64_t rdd){
    int64_t* nibbles =(int64_t*)calloc(lengthwords,sizeof(int64_t));
    if (nibbles == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    int64_t* parity = (int64_t*)calloc(lengthwords,sizeof(int64_t));
    if (parity == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    int64_t* pf = (int64_t*)calloc(lengthwords,sizeof(int64_t));
    if (pf == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    int64_t* p1 = NULL;
    int64_t* p2 = NULL;
    int64_t* p3 = NULL;
    int64_t* p4 = NULL;
    int64_t* p5 = NULL;
    int64_t vet1[] = {7,3,2,0};
    int64_t vet2[] = {6,3,1,0};
    int64_t vet3[] = {4,2,1,0};
    int64_t vet4[] = {4,3,2,1,0};
    int64_t vet5[] = {5,3,2,1};
    p1 = bit_reduce(words,vet1,lengthwords,4);
    p2 = bit_reduce(words,vet2,lengthwords,4);
    p3 = bit_reduce(words,vet3,lengthwords,4);
    p4 = bit_reduce(words,vet4,lengthwords,5);
    p5 = bit_reduce(words,vet5,lengthwords,4);
    if(Lora->hamming_decoding_en == 1){
        switch( rdd ){
    case 5:

    case 6:
        for(int64_t i = 0; i < lengthwords ; i++){
            nibbles[i] = mod_double( words[i] , 16);
        }
        break;
    case 7:

    case 8:
        for(int64_t i = 0; i < lengthwords ; i++){
            parity[i] = p2[i] * 4 + p3[i] * 2 + p5[i];
        }
        for(int64_t i = 0; i < lengthwords ; i++){
        switch(parity[i]){
            case 3:
                pf[i] = 4;
                break;
            case 5:
                pf[i] = 8;
                break;
            case 6:
                pf[i] = 1;
                break;
            case 7:
                pf[i] = 2;
                break;
            default:
                pf[i] = 0;
            }
            words[i] = words[i] ^ (int64_t)pf[i];
            nibbles[i] = mod_double( words[i] , 16);
        }

        break;
    default:
        printf("\nInvalid code rate!");
        return;
        }
    }
    printf("\nNibbles after Hamming decode:");
    for(int64_t i = 0; i < lengthwords ; i++){
        printf("%"PRId64" ", nibbles[i]);

    }
 return nibbles;
}


//FUNCTION HAMMING_ENCODE implements Hamming encoding to correct any errors in a set of codewords.
//The function takes the Nibbles as input and returns the codewords after encoding.
//ARGUMENTS: LORA, NIBBLES(pointer to a data vector before encoding)
//OUTPUT: CODEWORDS(pointer to a vector containing the words after encoding)
int64_t* hamming_encode(LoRaPHY* Lora , int64_t* nibbles , int64_t length_nibbles){
int64_t* codewords = (int64_t*)malloc(length_nibbles * sizeof(int64_t));
if (codewords == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
int64_t cr_now = 0;
int64_t* p1 = NULL;
int64_t* p2 = NULL;
int64_t* p3 = NULL;
int64_t* p4 = NULL;
int64_t* p5 = NULL;
int64_t vet1[] = {0,2,3};
int64_t vet2[] = {0,1,3};
int64_t vet3[] = {0,1,2};
int64_t vet4[] = {0,1,2,3};
int64_t vet5[] = {1,2,3};
int64_t* ws = NULL;
for( int64_t i = 0; i < length_nibbles ; i++){
    p1 = bit_reduce(nibbles, vet1, length_nibbles, 3);
    p2 = bit_reduce(nibbles, vet2, length_nibbles, 3);
    p3 = bit_reduce(nibbles, vet3, length_nibbles, 3);
    p4 = bit_reduce(nibbles, vet4, length_nibbles, 4);
    p5 = bit_reduce(nibbles, vet5, length_nibbles, 3);
    if( i < Lora->sf - 2)                                    /// the first SF-2 nibbles use CR=4/8
        cr_now = 4;
    else
        cr_now = Lora->cr;
    switch( cr_now ){
case 1:
    ws = (int64_t*)malloc(1 * sizeof(int64_t));
    ws[0] = p4[i] * 16;
    codewords[i] = ws[0] | nibbles[i];
    break;
case 2:
    ws = (int64_t*)malloc(3 * sizeof(int64_t));
    ws[0] = p5[i] * 32;
    ws[1] = p3[i] * 16;
    ws[2] = nibbles[i];
    codewords[i] = word_reduce(ws , 3);
    break;
case 3:
    ws = (int64_t*)malloc(4 * sizeof(int64_t));
    ws[0] = p2[i] * 64;
    ws[1] = p5[i] * 32;
    ws[2] = p3[i] * 16;
    ws[3] = nibbles[i];
    codewords[i] = word_reduce(ws , 4);
    break;
case 4:
    ws = (int64_t*)malloc(5 * sizeof(int64_t));
    ws[0] = p1[i] * 128;
    ws[1] = p2[i] * 64;
    ws[2] = p5[i] * 32;
    ws[3] = p3[i] * 16;
    ws[4] = nibbles[i];
    codewords[i] = word_reduce(ws , 5);
    break;
default:
    printf("\nInvalid code rate!");
    free(p1);
    p1 = NULL;
    free(p2);
    p2 = NULL;
    free(p3);
    p3 = NULL;
    free(p4);
    p4 = NULL;
    free(p5);
    p5 = NULL;
    ws = NULL;
    return;
    }
}

 return codewords;
}


//FUNCTION DEMODULATE WHICH DEMODULATES THE SIGNAL IN SIG.CFILE FROM THE PACKAGE INDEX AFTER THE PREAMBLE,
//AFTER 5 CONSECUTIVE DOWNCHIRPS (DETECTNEW FINDS THE INDEX).
//ARGUMENTS: LORA
//OUTPUT: Struct SYM_BOL(which contains a vector of the demodulated symbols, the cfo, the length of the symbol vector)
ResVetInt demodulate(LoRaPHY* Lora, double complex* vet_from_file, int64_t vetsize_from_file ){
    init(Lora);
    Lora->cfo = 0.0;
    int64_t sig_sizeFile = vetsize_from_file;
    double complex* sigdem = vet_from_file;
/*    for(int64_t i = 0; i < sig_sizeFile ; i++){
        printf("\nsigdem[%"PRId64"]=%f+(%f*I)",i,creal(sigdem[i]),cimag(sigdem[i]));
    }
    printf("\nsiglength=%"PRId64"",sig_sizeFile);*/
    ResFileDouble Fir = readDouble("./FirlpCoeff.txt");
    int64_t sig1_size = Fir.size_file + sig_sizeFile;
/*    printf("\nFIRlength=%"PRId64"",Fir.size_file);
    for(int64_t i = 0; i < Fir.size_file ; i++){
        printf("\nCoeff.[%"PRId64"]=%f",i,Fir.vet[i]);
    }*/
    double complex* h_imag = (double complex*)malloc(Fir.size_file * sizeof(double complex));
    if (h_imag == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for(int64_t i = 0; i < Fir.size_file ; i++){
        h_imag[i] = 0.0*I;
    }
    double complex* sig1 = (double complex*)malloc(sig1_size * sizeof(double complex));
    if (sig1 == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for(int64_t i = 0; i < sig1_size ; i++ ){
        if( i < sig_sizeFile){
           sig1[i] = creal(sigdem[i]) + cimag(sigdem[i])*I;
        }
        else{
            sig1[i] = 0;
        }
//    printf("\nsig1[%"PRId64"]=%f+%f*I",i,creal(sig1[i]),cimag(sig1[i]));
    }
    double complex* sigFir = (double complex*)malloc((sig1_size) * sizeof(double complex));
    if (sigFir == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    ConvolveCmplx(sig1 , sig1_size , Fir.vet , h_imag , Fir.size_file , sigFir);
/*    for(int64_t i = 0; i < sig1_size ; i++ ){
    printf("\nsigFir[%"PRId64"]=%f+(%f*I)",i,creal(sigFir[i]),cimag(sigFir[i]));
    }*/
    int64_t w = floor(Fir.size_file / 2.0);
    double complex* sigTrunc = (double complex*)malloc(sig_sizeFile * sizeof(double complex));
    if (sigTrunc == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for (int64_t i = 0 ; i < sig_sizeFile ; i++) {
		sigTrunc[i] = sigFir[i + w];
//		printf("\nsigTrun[%"PRId64"]=%f+%f*I",i,creal(sigTrunc[i]),cimag(sigTrunc[i]));
    }
    double p = (double)Lora->bw * 2.0;
    double pNorm = p / p;
    double q = (double)Lora->fs;
    int64_t qNorm = q / p;
    double sum_h = 0.0;
    ResFileDouble Lpfilt = readDouble("./LpfCoeff.txt");
//    printf("\nLplength=%"PRId64"",Lpfilt.size_file);
    for(int64_t i = 0; i < Lpfilt.size_file ; i++){
        sum_h += Lpfilt.vet[i];
    }
//    printf("\nsum_h=%f",sum_h);
    for(int64_t i = 0; i < Lpfilt.size_file ; i++){
        Lpfilt.vet[i] *= ( pNorm / sum_h);
//        printf("\nLpfCoeff.[%"PRId64"]=%f",i,Lpfilt.vet[i]);
    }
    int64_t sigLpf_size = (sig_sizeFile + Lpfilt.size_file);
    double complex* sigLpf = (double complex*)malloc(sigLpf_size * sizeof(double complex));
    if (sigLpf == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    h_imag = (double complex*)realloc(h_imag, Lpfilt.size_file * sizeof(double complex));
    if (h_imag == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for(int64_t i = 0; i < Lpfilt.size_file ; i++){
        h_imag[i] = 0.0*I;
//        printf("\nhimag[%"PRId64"]=%f", i , cimag(h_imag[i]));
    }
    double complex* sigResamp = (double complex*)calloc(sigLpf_size,sizeof(double complex));
    if (sigResamp == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for(int64_t i = 0; i < sigLpf_size ; i++ ){
        if(i < sig_sizeFile){
            sigLpf[i] = sigTrunc[i];
        }else{
            sigLpf[i] = 0.0 + 0.0 * I;
        }
//    printf("\nsigTruLPF[%"PRId64"]=%f+(%f)*I",i,creal(sigLpf[i]),cimag(sigLpf[i]));
    }
    ConvolveCmplx(sigLpf , sigLpf_size , Lpfilt.vet , h_imag , Lpfilt.size_file , sigResamp);
    for(int64_t i = 0; i < sigLpf_size - 1 ; i++ ){
//        printf("\nsigResam[%"PRId64"]=%f+(%f)*I",i,creal(sigResamp[i]),cimag(sigResamp[i]));         ///errore sul sigResamp!
        }
    int64_t len_y = round(sig_sizeFile * (p / q));
//    printf("\nleny=%"PRId64" qnorm=%"PRId64"",len_y,qNorm);
    int64_t z = ceil((Lpfilt.size_file + 2) / 2.0);
    double complex* y = (double complex*)malloc(len_y * sizeof(double complex));
    if (y == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    for(int64_t i = z , j = 0 ; i < sig_sizeFile + z ; i += qNorm , j++ ){  //decimo i campioni di un fattore 4
        y[j] = sigResamp[i];                                                //ovvero prendo un campione ogni 4
                                                                            //ma prima scarto i primi 44 campioni dovuti
//        printf("\ny[%"PRId64"]=%f+(%f)*I",j,creal(y[j]),cimag(y[j]));     //al ritardo del filtro antialiasing (N/2)
    }
    Lora->sig = y;
    Lora->sig_size = len_y;
    int64_t x = 0;
    int64_t temp = 0;
    ResVetInt sym_bols;
    double** pk_list_total = NULL;
    while(x < Lora->sig_size){
        x = detectnew(Lora, x);
        if(x < 0){
            break;
        }
        // align symbols with SFD
        x = syncnew(Lora, x);
        if( x == -1) {
                printf("\nNo syncronization found!");
        } else {
                printf("\nafter syncronization x=%"PRId64"",x);
        }                      //the goal is to extract payload_len from PHY header
        int64_t ii;            //header is in the first 8 symbols
        double pk_list[8][2];  //pk_list dovrebbe avere [8+sym_num]righe x [2]colonne ma al primo utilizzo estrae solo i primi 8 simboli
        int64_t symbols[8];    //i quali poi ricopierò nella matrice totale pk_list_total
        for(int64_t i = 0; i < 8 ; i++){
            pk_list[i][0] = 0.0;
            pk_list[i][1] = 0.0;
//            printf("\npk_list[%2d][0]=%.3f\tpk_list[%2d][1]=%.3f",i,pk_list[i][0],i,pk_list[i][1]);
        }
        if(x > Lora->sig_size - (8 * Lora->sample_num + 1)){
            return;
        }
        Peak pk1;
        for(ii = 0; ii < 8; ii++){
            pk1 = dechirp(Lora, x + (ii * Lora->sample_num) , true); // if arguments of dechirp < 2 ---> c=downchirp;
            pk_list[ii][0] = pk1.value;
            pk_list[ii][1] = pk1.index;
            printf("\npklist\nvalue=%8.4f from the index=%8"PRId64" --->",pk1.value,pk1.index);
//            printf("\nLora->binnum=%"PRId64" Lora->preamblebin=%"PRId64"", Lora->bin_num , Lora->preamble_bin);
            symbols[ii] = mod_double(((round(pk1.index + Lora->bin_num - Lora->preamble_bin) / (double)Lora->zero_padding_ratio)) , pow(2,Lora->sf));
            printf("\tsymbols=%4"PRId64"",symbols[ii]);
        }
        int64_t lengthdata = sizeof(symbols) / sizeof(symbols[0]);
//        printf("\nlengthdata=%"PRId64"",lengthdata);
        if(Lora->has_header){
            int64_t is_valid = parse_header(Lora, symbols , lengthdata);
            printf("\nis_valid=%"PRId64" after parse header!", is_valid);
            if(!is_valid){
                x = x + 7 * Lora->sample_num;
                printf("\nx=%"PRId64"", x);
                continue;
            }
        }
        // number of symbols in the packet
        int64_t sym_num = calc_sym_num(Lora , Lora->payload_len);
        printf("\nsym_num=%"PRId64"",sym_num);
        temp = sym_num;
        // demodulate the rest LoRa data symbols
        if(x > Lora->sig_size - sym_num * Lora->sample_num + 1 ){
            return;
        }

        //pk_list_total viene allocata come matrice dinamica a seconda della quantità di symbols rimasti da decodificare//
        // Allocazione dinamica della matrice//
        pk_list_total = (double**)malloc(sym_num * sizeof(double*));
        if (pk_list_total == NULL) {
            printf("\nRows memory allocation error\n");
            exit(1);
        }
        // Allocazione dinamica delle colonne per ogni riga
        for (int64_t i = 0; i < sym_num; i++) {
            pk_list_total[i] = (double*)malloc(2 * sizeof(double));
            if (pk_list_total[i] == NULL) {
                printf("\nColumns memory allocation error\n");
                exit(1);
            }
        }
//        printf("\nMatrix pk_list_total initialized for cointaining peak and index:\nPEAK\tINDEX");
        for(int64_t i = 0; i < sym_num ; i++){
//            printf("\n");
            for( int64_t j = 0; j < 2 ; j++){
                pk_list_total[i][j] = 0.0;
//                printf("%.3f\t", pk_list_total[i][j]);
            }
        }
        int64_t* symbols_tot = (int64_t*)malloc(sym_num * sizeof(int64_t));   //ricopio i primi 8 simboli calcolati prima nella matrice pk_list_total
        if (symbols_tot == NULL) {
            printf("\nError: Memory allocation failed!\n");
            exit(1);
        }
        for(int64_t j = 0 ; j < 8 ; j++){
            pk_list_total[j][0] = pk_list[j][0];
            pk_list_total[j][1] = pk_list[j][1];
            symbols_tot[j] = symbols[j];
//            printf("\npreamble symbols=%"PRId64"",symbols_tot[j]);
        }
        Peak pk2;
        for(ii = 8; ii < sym_num ; ii++){
            pk2 = dechirp(Lora, x + ii * Lora->sample_num, true);
            pk_list_total[ii][0] = pk2.value;
            pk_list_total[ii][1] = pk2.index;
            printf("\npklist\nvalue=%8.4f from the index=%8"PRId64" --->",pk2.value,pk2.index);
            symbols_tot[ii] = mod_double(round(((pk2.index + Lora->bin_num - Lora->preamble_bin) / (double)Lora->zero_padding_ratio)), pow(2,Lora->sf));
            //printf("\n");
            printf("\tsymbols=%4"PRId64"",symbols_tot[ii]);
        }
        x = x + sym_num * Lora->sample_num;
        // compensate CFO drift
        int64_t* symbDynComp ;
        symbDynComp = dynamic_compensation(Lora, symbols_tot , sym_num);
//        printf("\n");
/*        for( int64_t i = 0; i < sym_num ; i++){
            printf("\nsymbDynC[%"PRId64"]=%.3f", i , (double)symbDynComp[i]);
        }*/
        sym_bols.vet = (int64_t*)malloc(sym_num * sizeof(int64_t));
        if (sym_bols.vet == NULL) {
            printf("\nError: Memory allocation failed!\n");
            exit(1);
        }
        for(int64_t i = 0; i < sym_num; i++){
            sym_bols.length_vet = sym_num;
            sym_bols.vet[i] = mod_double((double)symbDynComp[i] , pow(2,Lora->sf));
            sym_bols.cfo = Lora->cfo;
        }
        free(symbDynComp);
        symbDynComp = NULL;
        free(symbols_tot);
        symbols_tot = NULL;
    }

    if( temp == 0){
        printf("No preamble detected!\n");
    }
free(Lpfilt.vet);
Lpfilt.vet = NULL;
free(Fir.vet);
Fir.vet = NULL;
free(h_imag);
h_imag = NULL;
free(sig1);
sig1 = NULL;
free(sigFir);
sigFir = NULL;
free(sigTrunc);
sigTrunc = NULL;
free(sigLpf);
sigLpf = NULL;
free(sigResamp);
sigResamp = NULL;
free(y);
y = NULL;
return sym_bols;
}


//FUNCTION MODULATE THAT MODULATES A BASEBAND SIGNAL.
//ARGUMENTS: LORA,SYMBOLS(pointer to a vector containing chirped symbols which must be modulated, valid symbol range [0...2^sf-1])
//OUTPUT: Struct RESULTCHIRP(Struct which contains vector of double complex baseband signals for Lora, and the length of the vector)
ResVetCmplx modulate(LoRaPHY* Lora, int64_t* symbols , int64_t length_symbols){
ResVetCmplx uc = chirp(true , Lora->sf , Lora->bw , Lora->fs , 0 , Lora->cfo , 0);
ResVetCmplx dc = chirp(false , Lora->sf , Lora->bw , Lora->fs , 0 , Lora->cfo , 0);
int64_t chirp_len = uc.size_chirp;
//printf("\nsizechirp= %"PRId64"",uc.size_chirp);
int64_t new_length_preamble = Lora->preamble_len * uc.size_chirp;
//printf("\nlength_preamble=%"PRId64"", new_length_preamble);
//stampaVetCmplx( uc.vet , chirp_len);
//stampaVetCmplx( dc.vet , chirp_len);
double complex* preamble = (double complex*) malloc((new_length_preamble) * sizeof(double complex));
if (preamble == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
for(int64_t i = 0; i < Lora->preamble_len ; i++){
    for(int64_t j = 0; j < chirp_len; j++){
        preamble[i * chirp_len + j] = uc.vet[j];
    }
}
//stampaVetCmplx(preamble , Lora->preamble_len * uc.size_chirp);
ResVetCmplx netid1 = chirp(true , Lora->sf , Lora->bw , Lora->fs , 24 , Lora->cfo , 0);
ResVetCmplx netid2 = chirp(true , Lora->sf , Lora->bw , Lora->fs , 32 , Lora->cfo , 0);
double complex* netid = (double complex*)malloc((chirp_len * 2) * sizeof(double complex));
if (netid == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
for(int64_t i = 0; i < 2 ; i++){
    for(int64_t j = 0; j < chirp_len  ; j++){
        if( i == 0)
            netid[j] = netid1.vet[j];
        else
            netid[chirp_len + j] = netid2.vet[j];
    }
}
//printf("\nnetidlength=%"PRId64"",chirp_len * 2);
//stampaVetCmplx(netid , chirp_len*2);
int64_t sfd_length = chirp_len * 2 + round((double)chirp_len / 4.0);
//printf("\nsfdlength=%"PRId64"",sfd_length);
double complex* sfd = (double complex*)malloc(sfd_length * sizeof(double complex));
if (sfd == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
for(int64_t i = 0; i < 3  ; i++ ){
    if(i < 2){
        for(int64_t j = 0; j < chirp_len ; j++){
            sfd[i * chirp_len + j] = dc.vet[j];
        }
    }
    else{
        for(int64_t j = 0; j < round((double)chirp_len / 4.0) ; j++){
            sfd[i * chirp_len + j] = dc.vet[j];
        }
    }
}
for(int64_t i = 0; i < sfd_length; i++ ){
    sfd[i] = dc.vet[i%chirp_len];
}
//stampaVetCmplx(sfd,sfd_length);
int64_t length_data = length_symbols * chirp_len;
//printf("\ndatalength=%"PRId64"",length_data);
int64_t length_s = length_data + sfd_length + chirp_len * 2 + new_length_preamble;
double complex* s = (double complex*)malloc(length_s * sizeof(double complex));
if (s == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
//printf("\nlength_s=%"PRId64"",length_s);
double complex* datachirp = (double complex*)malloc(length_data * sizeof(double complex));
if (datachirp == NULL) {
    printf("\nError: Memory allocation failed!\n");
    exit(1);
}
ResVetCmplx data;
//data.vet = (double complex*)malloc(chirp_len * sizeof(double complex));
for(int64_t i = 0; i < length_symbols ; i++){
    data = chirp(true , Lora->sf , Lora->bw , Lora->fs , symbols[i] , Lora->cfo , 0);
    for(int64_t j = 0 ; j < chirp_len ; j++){
        datachirp[j + i * chirp_len] = data.vet[j];
    }
    free(data.vet);
    data.vet = NULL;
}

//stampaVetCmplx(datachirp, length_data);
int64_t j = 0;
for(int64_t i = 0 ; i < length_s ; i++) {
    if (i < new_length_preamble){
        s[i] = preamble[i];
    }
    else if(i < new_length_preamble + chirp_len * 2) {
        s[i] = netid[i - new_length_preamble];
    }
    else if(i < new_length_preamble + chirp_len * 2 + sfd_length) {
        s[i] = sfd[i - (new_length_preamble + chirp_len * 2)];
    }
    else {
        s[i] = datachirp[j];
        j++;
    }
}
//stampaVetCmplx(s , length_s);
ResVetCmplx res_modulate;
res_modulate.size_chirp = length_s;
res_modulate.vet = (double complex*)malloc(length_s * sizeof(double complex));
for(int64_t k = 0 ; k < length_s ; k++){
    res_modulate.vet[k] = s[k];
}

//stampaVetCmplx(res_modulate.vet, res_modulate.size_chirp);

///This part is for saving output from "modulate" to a binary file named output.bin.
/*char filename[] = "output.bin";
int64_t success = 0;
success = write_file(s, length_s , filename);
if(success){
    printf("\nDati scritti con successo nel file %s.\n", filename);
}else{
    printf("\nErrore durante la scrittura nel file %s.\n", filename);
}*/

free(uc.vet);
uc.vet = NULL;
free(dc.vet);
dc.vet = NULL;
free(netid1.vet);
netid1.vet = NULL;
free(netid2.vet);
netid2.vet = NULL;
free(preamble);
preamble = NULL;
free(netid);
netid = NULL;
free(sfd);
sfd = NULL;
free(s);
s = NULL;
free(datachirp);
datachirp = NULL;
return res_modulate;
}


//FUNCTION ENCODE ENCODE THE BYTES IN THE PAYLOAD INTO SYMBOLS.
//ARGUMENTS: LORA, PAYLOAD(pointer to a vector containing the payload)
//OUTPUT: SYMBOLS(pointer to a vector containing symbols encoded)
ResVetInt encode(LoRaPHY* Lora, int64_t* payload , int64_t length_payload){
ResVetInt symbols;
int64_t* checksum = NULL;
int64_t* data = NULL;
int64_t length_data;
int64_t* data_nibbles = NULL;
int64_t* header_nibbles = NULL;
int64_t* data_nibbles_hasheader = NULL;
//  printf("\npayload length = %"PRId64"", length_payload);
int64_t sym_num = calc_sym_num(Lora , length_payload);
//  printf("\nsym_num=%"PRId64"",sym_num);
int64_t nibble_num = (Lora->sf - 2 + (sym_num - 8) / (Lora->cr + 4)*(Lora->sf - 2 * Lora->ldr));
//  printf("\nnibblenum=%"PRId64"",nibble_num);
if(Lora->crc){
    checksum = calc_crc(Lora, payload , length_payload);
    data = (int64_t*)malloc((length_payload + 2)* sizeof(int64_t));
    if (data == NULL) {
        printf("\nError: Memory allocation failed!\n");
        exit(1);
    }
    printf("\nData: \n");
    for(int64_t i = 0; i < length_payload + 2 ; i++){
        if(i < length_payload)
            data[i] = payload[i];
        else if(i == length_payload)
            data[i] = checksum[0];
        else
            data[i] = checksum[1];
        printf("%"PRId64" ",data[i]);
    }
    length_data = length_payload + 2;
//    printf("\nlenghtdata=%"PRId64"",length_data);
    if(checksum != NULL){
        free(checksum);
        checksum = NULL;
    }
}else{
    data = (int64_t*)malloc( length_payload * sizeof(int64_t));
    if (data == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    printf("\nData: \n");
    for(int64_t i = 0; i < length_payload ; i++){
        data[i] = payload[i];
        printf("%"PRId64" ",data[i]);
    }
    length_data = length_payload;
//    printf("\nlenghtdata=%"PRId64"",length_data);
}
int64_t a = ceil((double)(nibble_num - 2 * length_data)/ 2.0);
//    printf("\na=%"PRId64"",a);
int64_t* data_W = (int64_t*)malloc((length_data + a) * sizeof(int64_t));
if (data_W == NULL) {
    printf("\nError: Memory allocation failed!");
    exit(1);
}
printf("\nData_W : \n");
for(int64_t i = 0; i < length_data + a ; i++){
    if( i < length_data ){
        data_W[i] = data[i];
   }else{
        data_W[i] = 255;
    }
printf("%"PRId64" ", data_W[i]);
}
int64_t** data_w_mat = whiten(Lora , data_W , length_payload);
int64_t* data_w_new = (int64_t*)malloc((length_data + a) * sizeof(int64_t));
if (data_w_new == NULL) {
    printf("\nError: Memory allocation failed!");
    exit(1);
}
printf("\nData after whitening: \n");
for(int64_t ii = 0; ii < length_data + a ; ii++){       ///salvo i valori della diagonale della matrice dopo il whitening
    if( ii < length_payload){
        data_w_new[ii] = data_w_mat[ii][ii];
   }else{
        data_w_new[ii] = data_W[ii];
    }
printf("%"PRId64" ",data_w_new[ii]);
}
data_nibbles = (int64_t*)calloc(nibble_num , sizeof(int64_t));
if (data_nibbles == NULL) {
    printf("\nError: Memory allocation failed!");
    exit(1);
}
//printf("\nData nibbles without header nibbles: \n");
for(int64_t i = 0 , j = 1 ; i < nibble_num ; i++ , j++){
    int64_t idx = ceil((double)j / 2.0);
//    printf("\n%"PRId64" ",idx);
    if( mod_double(j,2) == 1 ){
        data_nibbles[i] = data_w_new[idx - 1] & 15;
    }else{
        data_nibbles[i] = data_w_new[idx - 1] >> 4;
    }
//    printf("%"PRId64" ",data_nibbles[i]);
}
int64_t sizenibbles = 0;
if( Lora->has_header ){
    header_nibbles = gen_header(Lora , length_payload); ///sono sempre 5 i nibbles in uscita///
    printf("\nheader nibbles: ");
    for(int64_t i = 0; i < 5 ; i++){
        printf("%"PRId64" ", header_nibbles[i]);
    }
    sizenibbles = 5 + nibble_num;
//    printf("\nsizenibbles=%"PRId64"", sizenibbles);
    data_nibbles_hasheader = (int64_t*)calloc(sizenibbles,sizeof(int64_t));
    if (data_nibbles_hasheader == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    printf("\nHeader+Datanibbles: \n");
    int64_t j = 0;
    for(int64_t i = 0 ; i < sizenibbles ; i++ ){
        if(i < 5)
            data_nibbles_hasheader[i] = header_nibbles[i];
        else{
            data_nibbles_hasheader[i] = data_nibbles[j];
            j++;
        }
        printf("%"PRId64" ",data_nibbles_hasheader[i]);
    }
}else{
    sizenibbles = nibble_num;
//    printf("\nsizenibbles=%"PRId64"", sizenibbles);
    data_nibbles_hasheader = (int64_t*)malloc(sizenibbles * sizeof(int64_t));
    if (data_nibbles_hasheader == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    printf("\nDatanibbles without header: \n");
    for(int64_t i = 0; i < nibble_num ; i++){
        data_nibbles_hasheader[i] = data_nibbles[i];
    }
    for(int64_t i = 0; i < sizenibbles ; i++){
        printf("%"PRId64" ",data_nibbles_hasheader[i]);
    }
 }
    int64_t* codewords;
    codewords = hamming_encode(Lora, data_nibbles_hasheader ,sizenibbles);
    printf("\nCodewords after Hamming encode: \n");
    for(int64_t i = 0; i < sizenibbles ; i++){
        printf("%"PRId64" ", codewords[i]);
    }
///first 8 symbols use CR=4/8 ---> rdd = 8 setted in diag_interleave///
///from the first 5 codewords we have 8 symbols interleaved after diag_interleave///
    int64_t* symbols_i = diag_interleave(Lora, codewords , Lora->sf - 2 , 8); ///With the first (sf - 2) codewords i get always the first 8 symbols_i
    int64_t* temp = (int64_t*)malloc(8 * sizeof(int64_t));
    if (temp == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    for(int64_t k = 0; k < 8 ; k++){
        temp[k] = symbols_i[k];
    }
    int64_t ppm = Lora->sf - 2 * Lora->ldr;
    int64_t rdd = Lora->cr + 4;
    int64_t symbols_i_length = 0;
    for(int64_t i = Lora->sf - 2 ; i < sizenibbles - ppm + 1 ; i += ppm){
        int64_t* current_symbols = diag_interleave(Lora, &codewords[i] , ppm , rdd);
        for(int64_t j = 0 ; j < rdd ; j++) {
            symbols_i = (int64_t*)realloc(symbols_i, (symbols_i_length + rdd) * sizeof(int64_t));
            if (symbols_i == NULL) {
                printf("\nError: Memory allocation failed!");
                exit(1);
            }
            symbols_i[symbols_i_length + j] = current_symbols[j];
        }
        symbols_i_length += rdd;
        free(current_symbols);
        current_symbols = NULL;

    }
    printf("\nSymbols_i total after interleaving: \n");
    int64_t* symbols_i_tot = (int64_t*)malloc((8 + symbols_i_length) * sizeof(int64_t));
    if (symbols_i_tot == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    memcpy(symbols_i_tot, temp, 8 * sizeof(int64_t));
    memcpy(symbols_i_tot + 8, symbols_i, symbols_i_length * sizeof(int64_t));
    for(int64_t k = 0; k < 8 + symbols_i_length ; k++){
       printf("%"PRId64" ", symbols_i_tot[k]);
    }
    symbols.vet = (int64_t*)malloc((8 + symbols_i_length) * sizeof(int64_t));
    if (symbols.vet == NULL) {
        printf("\nError: Memory allocation failed!");
        exit(1);
    }
    symbols.vet = gray_decoding(Lora, symbols_i_tot, 8 + symbols_i_length);
    symbols.length_vet = 8 + symbols_i_length;

free(data);
data = NULL;
free(data_W);
data_W = NULL;
free(data_w_mat);
data_w_mat = NULL;
free(data_w_new);
data_w_new = NULL;
free(data_nibbles);
data_nibbles = NULL;
free(codewords);
codewords = NULL;
free(data_nibbles_hasheader);
data_nibbles_hasheader = NULL;
free(temp);
temp = NULL;
free(symbols_i_tot);
symbols_i_tot = NULL;
free(symbols_i);
symbols_i = NULL;
 return symbols;
}


//FUNCTION DECODE
//ARGUMENTS: LORA, SYMBOL_M(pointer to a vector containing symbols to be decoded),LENGTH_SYMBOLM(size of the vector)
//OUTPUT: Struct DATA(Struct containing decoded data, which in the last two elements there are the checksum value if enabled))
ResVetInt decode(LoRaPHY* Lora, int64_t* symbols_m , int64_t number_cols , int64_t number_rows){
int64_t* nibbles = NULL;
int64_t* codewords = NULL;
int64_t* symbols_g_column = NULL;
int64_t* bitarray = NULL;
int64_t* nibbles_bin0 = NULL;
int64_t* nibbles_bin1 = NULL;
int64_t* nibbles_bin2 = NULL;
int64_t* bytes = NULL;
int64_t* checksum = NULL;
ResVetInt data;
int64_t* nibbles_total= (int64_t*)calloc(number_rows*2,sizeof(int64_t));
int64_t* header_checksum = (int64_t*)calloc(5,sizeof(int64_t));
int64_t* header_nibbles = (int64_t*)calloc(12,sizeof(int64_t));
int64_t* nibbles_bits = (int64_t*)calloc(5,sizeof(int64_t));
if (nibbles_total == NULL || header_checksum == NULL || header_nibbles == NULL || nibbles_bits == NULL) {
    printf("\nError: Memory allocation failed in decode!");
    exit(1);
}
    symbols_g_column = gray_coding(Lora, symbols_m , number_rows);
    codewords = diag_deinterleave(Lora, symbols_g_column , 8 , Lora->sf - 2);
    if(!Lora->has_header){
        printf("\nThe data packet has not Header!\n");
        nibbles = hamming_decode(Lora, codewords , Lora->sf - 2 , 8);
        for(int64_t i = 0; i < Lora->sf - 2 ; i++){
//            printf("%"PRId64" ",nibbles[i]);
            nibbles_total[i] = nibbles[i];
        }
    }else{
        nibbles = hamming_decode(Lora, codewords , Lora->sf - 2 , 8);
        printf("\nNibbles with Header after hamming decode: ");
        for(int64_t i = 0 , j = 0 ; i < Lora->sf - 2 ; i++ ){       ///Save the last nibbles minus the header nibbles(5)///
            printf("%"PRId64" ",nibbles[i]);
            if( i > 4 ){
                nibbles_total[j] = nibbles[i];
                j++;
            }
        }
        Lora->payload_len = nibbles[0] * 16 + nibbles[1];
        Lora->crc = nibbles[2] & 1;
        Lora->cr = nibbles[2] >> 1;
        bitarray = convertiInt4Bit(nibbles[4]);
        printf("\nHeader_checksum bits:");
        for(int64_t i = 0 , j = 0 ; i < 5  ; i++ ){
            if( i == 0){
                header_checksum[i] = nibbles[3] & 1 ;
           }else{
                header_checksum[i] = bitarray[j];
                j++;
            }
            printf(" %"PRId64"", header_checksum[i]);
        }
        nibbles_bin0 = convertiInt4Bit(nibbles[0]);
        nibbles_bin1 = convertiInt4Bit(nibbles[1]);
        nibbles_bin2 = convertiInt4Bit(nibbles[2]);
        int64_t vet_bin_nibbles[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        printf("\nnibbles[0]: ");
        for(int64_t i = 0; i < 4 ; i++){
            vet_bin_nibbles[i] = nibbles_bin0[i];
            printf("%"PRId64" ",nibbles_bin0[i]);
        }
        printf("\nnibbles[1]: ");
        for(int64_t i = 0; i < 4 ; i++){
            vet_bin_nibbles[i+4] = nibbles_bin1[i];
            printf("%"PRId64" ",nibbles_bin1[i]);
        }
        printf("\nnibbles[2]: ");
        for(int64_t i = 0; i < 4 ; i++){
            vet_bin_nibbles[i+8] = nibbles_bin2[i];
        printf("%"PRId64" ",nibbles_bin2[i]);
        }

        printf("\n\nHeader Checksum Matrix: \n");
        for(int64_t i = 0; i < 5; i++) {
            for(int64_t j = 0; j < 12; j++) {
                printf("%"PRId64" ", Lora->header_checksum_matrix[i][j]);
            }
            printf("\n");
        }
        printf("\n\n(*)");
        printf("\n\nVet_Nibbles: ");
        for(int64_t i = 0; i < 12 ; i++){
            printf("%"PRId64" ",vet_bin_nibbles[i]);
        }
        printf("\n\n = Result\t\t\t\t\t\tnibbles_bits\n");
        int64_t total;
        for(int64_t i = 0; i < 5; i++) {
            total = 0;
            for(int64_t j = 0; j < 12; j++) {
                header_nibbles[j] = Lora->header_checksum_matrix[i][j] * vet_bin_nibbles[j];
                printf("%"PRId64" ", header_nibbles[j]);
                if(header_nibbles[j] % 2 != 0)
                    total++;
            }
            if(total % 2 == 0){
                nibbles_bits[i] = 0;
            }else{
                nibbles_bits[i] = 1;
            }
            printf("    (Binary Sum of the row)---------->%"PRId64"" ,nibbles_bits[i]);
            printf("\n");
            if(header_checksum[i] != nibbles_bits[i]){
                printf("\nHeader Checksum bit: %"PRId64" and Nibbles bit: %"PRId64" are not equal!", header_checksum[i] , nibbles_bits[i]);
                printf("\nInvalid Header Checksum!");
            }
        }
    }
    int64_t rdd = Lora->cr + 4;
    int64_t nibbles_length;
    if(!Lora->has_header){
        nibbles_length = Lora->sf - 2;
        Lora->payload_len = 0;
    }else{
        nibbles_length = Lora->sf - 2 - 5;
//    printf("\nnibbles_length=%d",nibbles_length);
    }
    int64_t* current_codewords = NULL;
    int64_t* current_nibbles = NULL;
    int64_t codewords_length = 0;
    for(int64_t ii = 8 ; ii < number_rows ; ii += rdd){
        current_codewords = diag_deinterleave(Lora, &symbols_g_column[ii] , rdd  , (Lora->sf - 2 * Lora->ldr) );
        current_nibbles = hamming_decode(Lora, current_codewords , (Lora->sf - 2 * Lora->ldr) , rdd);
        printf("\nCurrent nibbles in loop: \n");
        for (int64_t i = 0; i < (Lora->sf - 2 * Lora->ldr) ; i++) {
            printf("%"PRId64" ", current_nibbles[i]);
            nibbles_total[nibbles_length + i] = current_nibbles[i];
        }
        nibbles_length += (Lora->sf - 2 * Lora->ldr);
        codewords_length += (Lora->sf - 2 * Lora->ldr);
        free(current_codewords);
        current_codewords = NULL;
        free(current_nibbles);
        current_nibbles = NULL;
/*        printf("\nnibbles_length=%d",nibbles_length);
        printf("\ncodewords_length=%d",codewords_length);*/
    }

    int64_t bytesLength = fmin(255, floor(nibbles_length / 2));
    bytes = (int64_t*)calloc(bytesLength,sizeof(int64_t));
    if (bytes == NULL) {
        printf("\nMemory allocation error for Bytes");
        exit(1);
    }
//    printf("\nnibbleslength=%"PRId64" bytelength=%"PRId64"\n",nibbles_length , bytesLength);
    printf("\nBytes before dewhitening: \n");
    for (int64_t ii = 0; ii < bytesLength ; ii++) {
        bytes[ii] = (nibbles_total[2 * ii]) | (16 * nibbles_total[2 * ii + 1]);
        printf("%"PRId64" ",bytes[ii]);
    }
/*    for (int64_t ii = 0; ii < nibbles_length ; ii++) {
        printf("\nnibbles_total[%"PRId64"]=%"PRId64"", ii , nibbles_total[ii]);
    }*/
    printf("\n");
    int64_t len = Lora->payload_len;                ///WARNING: if the packet has not header (has_header=0) Lora->payload_len is setted to zero!!!
//    printf("\nLora->payload_len = %"PRId64"", len);
    if(Lora->crc && Lora->has_header){
        data.length_vet = len + 2;
        data.vet = dewhiten(Lora, bytes , len + 2);
        printf("\nData dewhiten: \n");
        for(int64_t k = 0; k < len ; k++){
            printf("%"PRId64" ", data.vet[k]);
        }
        data.vet[len] = bytes[len];
        data.vet[len + 1] = bytes[len + 1];
        printf("\nData decoded: \n");
        for(int64_t j = 0; j < len + 2 ; j++){
            printf("%"PRId64" ", data.vet[j]);
        }
        printf("\n");
        checksum = calc_crc(Lora, data.vet , len );
        if(checksum != NULL){
            free(checksum);
            checksum = NULL;
        }
    }else{
        if(Lora->has_header){
            data.length_vet = len;
            data.vet = dewhiten(Lora, bytes, len + 2);
            printf("\nData decoded without Checksum bits!\n");
            for(int64_t j = 0; j < len ; j++){
                printf("%"PRId64" ", data.vet[j]);
            }
        }else{
            printf("\nData can't be decoded without header!\n");
        }
    }
    symbols_g_column = NULL;
    codewords = NULL;
    bytes = NULL;
    nibbles = NULL;
    bitarray = NULL;
    nibbles_total = NULL;
    free(header_checksum);
    header_checksum = NULL;
    free(nibbles_bin0);
    nibbles_bin0 = NULL;
    free(nibbles_bin1);
    nibbles_bin1 = NULL;
    free(nibbles_bin2);
    nibbles_bin2 = NULL;
    nibbles_bits = NULL;
 return data;
}


//FUNCTION WHITEN implement a whitening process to make the data pseudo-casual using bitwise XOR operation.
//ARGUMENTS: LORA, DATA(pointer to a vector containing datas before the whitening)
//OUTPUT: DATA_W (pointer to a matrix containing datas after whitening)
int64_t** whiten(LoRaPHY* Lora , int64_t* data , int64_t lengthdata){
int64_t* data_w = (int64_t*)malloc(lengthdata * sizeof(int64_t));
if (data_w == NULL) {
    printf("\nRows memory allocation error\n");
    exit(1);
}
int64_t** matrice = (int64_t**)malloc(lengthdata * sizeof(int64_t*));
if (matrice == NULL) {
    printf("\nRows allocation error\n");
    exit(1);
    }
// Allocazione dinamica delle colonne per ogni riga
for (int64_t i = 0; i < lengthdata; i++) {
    matrice[i] = (int64_t*)malloc(lengthdata * sizeof(int64_t));
    if (matrice[i] == NULL) {
        printf("\nColumns memory allocation error\n");
        exit(1);
    }
}

for(int64_t i = 0; i < lengthdata ; i++){
//    printf("\n");
    for( int64_t j = 0; j < lengthdata ; j++){
        data_w[i] = Lora->whitening_seq[j] ^ data[i];
        matrice[i][j] = data_w[i];
//        printf("%3d ",data_w[i]);
    }
}
free(data_w);
data_w = NULL;
 return matrice;
}


//FUNCTION DEWHITEN Implement a reverse bleaching process to make the data correlated with the data sent.
//ARGUMENTS: LORA, BYTES(pointer to a vector containing bytes after the interleaving)
//OUTPUT: BYTES_W(pointer to a vector containing bytes after de-whitening)
int64_t* dewhiten(LoRaPHY* Lora , int64_t* bytes , int64_t lengthbytes){
int64_t* bytes_dew = (int64_t*)calloc(lengthbytes,sizeof(int64_t));
int64_t* bytes_w = (int64_t*)calloc(lengthbytes,sizeof(int64_t));
int64_t** matrice = (int64_t**)calloc(lengthbytes,sizeof(int64_t*));
if (bytes_dew == NULL || bytes_w == NULL || matrice == NULL) {
    printf("\nError: Memory allocation failed!");
    exit(1);
}
// Allocazione dinamica delle colonne per ogni riga
for (int64_t i = 0; i < lengthbytes; i++) {
    matrice[i] = (int64_t*)calloc(lengthbytes,sizeof(int64_t));
    if (matrice[i] == NULL) {
        printf("\nColumns memory allocation error!\n");
        exit(1);
    }
}
for(int64_t i = 0; i < lengthbytes ; i++){
//    printf("\n");
    for( int64_t j = 0; j < lengthbytes ; j++){
        bytes_w[i] = Lora->whitening_seq[j] ^ bytes[i];
        matrice[i][j] = bytes_w[i];
        }
}
for(int64_t i = 0; i < lengthbytes ; i++){
    bytes_dew[i] = matrice[i][i];
//    printf("%3d ",bytes_dew[i]);
}

free(bytes_w);
bytes_w = NULL;
for(int64_t i = 0; i < lengthbytes; i++) {
    free(matrice[i]);
}
free(matrice);
matrice = NULL;
 return bytes_dew;
}


//FUNCTION GEN_HEADER GENERATES THE HEADER INFORMATION BLOCK IN THE LORA PACKAGE.
//ARGUMENTS: LORA, PLEN(payload length)
//OUTPUT: HEADER_NIBBLES(pointer to a vector containing decimal nibble)
int64_t* gen_header(LoRaPHY* Lora, int64_t payload_length){         ///SONO SEMPRE 5 I NIBBLES IN USCITA DELL'HEADER///
int64_t* header_nibbles = (int64_t*)calloc(12,sizeof(int64_t));
int64_t* header_nibbles_bits = (int64_t*)calloc(5,sizeof(int64_t));
int64_t* nibbles = (int64_t*)calloc(5,sizeof(int64_t));
if(header_nibbles == NULL || header_nibbles_bits == NULL || nibbles == NULL ) {
    printf("\nMemory allocation error in gen_header\n");
    exit(1);
}
printf("\nLora->cr= %"PRId64"\tLora->crc=%"PRId64"",Lora->cr, Lora->crc);
nibbles[0] = payload_length >> 4;
nibbles[1] = payload_length & 15;
nibbles[2] = (2*Lora->cr) | Lora->crc;
//printf("\nnibble0=%"PRId64"\tnibble1=%"PRId64"\tnibble2=%"PRId64"",nibbles[0],nibbles[1],nibbles[2]);
int64_t* nibbles_bin0 = convertiInt4Bit(nibbles[0]);
int64_t* nibbles_bin1 = convertiInt4Bit(nibbles[1]);
int64_t* nibbles_bin2 = convertiInt4Bit(nibbles[2]);
int64_t vet_bin_nibbles[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
printf("\nnibbles[0]: ");
for(int64_t i = 0; i < 4 ; i++){
    vet_bin_nibbles[i] = nibbles_bin0[i];
    printf("%"PRId64" ",nibbles_bin0[i]);
}
printf("\nnibbles[1]: ");
for(int64_t i = 0 ; i < 4 ; i++){
    vet_bin_nibbles[i + 4] = nibbles_bin1[i];
    printf("%"PRId64" ",nibbles_bin1[i]);
}
printf("\nnibbles[2]: ");
for(int64_t i = 0 ; i < 4 ; i++){
    vet_bin_nibbles[i + 8] = nibbles_bin2[i];
    printf("%"PRId64" ",nibbles_bin2[i]);
}
printf("\nheader_checksum_matrix (*) header_nibbles[0,1,2]: \n");
for(int64_t i = 0; i < 5; i++) {
    int64_t total = 0;
    for(int64_t j = 0; j < 12; j++) {
        header_nibbles[j] = Lora->header_checksum_matrix[i][j] * vet_bin_nibbles[j];
        printf("%"PRId64" ", header_nibbles[j]);
        if(header_nibbles[j] % 2 != 0)
            total++;
    }
    if(total % 2 == 0)
        header_nibbles_bits[i] = 0;
    else
        header_nibbles_bits[i] = 1;
    printf("\t\t(Binary sum of the row)--->%"PRId64"\n",header_nibbles_bits[i]);
}
nibbles[3] = header_nibbles_bits[0];
for(int64_t i = 0; i < 4; i++) {
    nibbles[4] |= (int64_t) (header_nibbles_bits[ i + 1 ] * pow(2, 3 - i ));
}
printf("\nHeader Nibbles: ");
for(int64_t i = 0; i < 5; i++){
    printf("%"PRId64" ",nibbles[i]);
}

free(nibbles_bin0);
nibbles_bin0 = NULL;
free(nibbles_bin1);
nibbles_bin1 = NULL;
free(nibbles_bin2);
nibbles_bin2 = NULL;
free(header_nibbles);
header_nibbles = NULL;
free(header_nibbles_bits);
header_nibbles_bits = NULL;
 return nibbles;
}


//FUNCTION SYMBOLS_TO_BYTES DERIVE PAYLOAD BYTES FROM SYMBOLS.
//ARGUMENTS: LORA, SYMBOLS(pointer to a vector cointaining the symbols), LENGTH_SYMBOLS(size of the vector containing the symbols)
//OUTPUT: Struct BYTES_DEC(Struct of payload bytes derived from symbols)
ResVetInt symbols_to_bytes(LoRaPHY* Lora, int64_t* symbols , int64_t length_symbols){
init(Lora);
Lora->hamming_decoding_en = 0;
int64_t slen_tmp;
if(length_symbols <= 4){
    slen_tmp = 8 + Lora->has_header * (Lora->cr + 4);
}else{
    slen_tmp = 8 + ceil((length_symbols - 4 *(1 - Lora->has_header)) / 4.0 ) * (Lora->cr + 4);
}
int64_t plen = 0;
plen = calc_payload_len(Lora, slen_tmp , true);
Lora->payload_len = plen;
int64_t sym_num = 0;
sym_num = (int64_t)calc_sym_num(Lora, Lora->payload_len);
//printf("\nslen_tmp=%"PRId64"\tsym_num=%"PRId64"", slen_tmp , sym_num);
int64_t* symbols_ = (int64_t*)calloc(sym_num,sizeof(int64_t));
int64_t jj = 0;
if(Lora->has_header){
    jj = 8;
}else{
    jj = 0;
}
for(int64_t ii = 0; ii < length_symbols; ii += 4){
    if(ii + 3 <= length_symbols ){
        symbols_[jj] = symbols[ii];
        symbols_[jj+1] = symbols[ii+1];
        symbols_[jj+2] = symbols[ii+2];
        symbols_[jj+3] = symbols[ii+3];
    }else{
        symbols_[jj] = symbols[ii];
        if(ii + 1 < length_symbols){
            symbols_[jj+1] = symbols[ii+1];
        }else{
            symbols_[jj+1] = 0;
        }
        if(ii + 2 < length_symbols){
            symbols_[jj+2] = symbols[ii+2];
        }else{
            symbols_[jj+2] = 0;
        }
        if(ii + 3 < length_symbols){
            symbols_[jj+3] = symbols[ii+3];
        }else{
            symbols_[jj+3] = 0;
        }
    }
    if(jj == 1){
        jj = 9;
    }else{
        jj = jj + Lora->cr + 4;
    }
}
ResVetInt symbols_tmp;
symbols_tmp.vet = (int64_t*)calloc(sym_num, sizeof(int64_t));
int64_t* payload_zeros = (int64_t*)calloc(Lora->payload_len, sizeof(int64_t));
if(Lora->has_header){
    symbols_tmp = encode(Lora, payload_zeros , plen);
//    printf("\n Symbols tmp in symbols_to_byte: \n");
    for(int64_t z = 0 ; z < 8 ; z++){
        symbols_[z] = symbols_tmp.vet[z];
//        printf("%"PRId64" ", symbols_tmp[z]);
    }
}
printf("\n\nSymbols_ are: \n");
for(int64_t k = 0; k < sym_num ; k++){
    printf("%"PRId64" ", symbols_[k]);
}

ResVetInt bytes_dec;
bytes_dec.vet = (int64_t*)calloc(plen, sizeof(int64_t));
bytes_dec = decode(Lora, symbols_ , 1 , sym_num);
Lora->hamming_decoding_en = 1;
Lora->payload_len = plen;
free(symbols_);
symbols_ = NULL;
free(symbols_tmp.vet);
symbols_tmp.vet = NULL;
free(payload_zeros);
payload_zeros = NULL;
 return bytes_dec;
}


//FUNCTION TIME_ON_AIR CALCULATE THE FLYING TIME OF A LORA PACKET.
//ARGUMENTS: LORA, PLEN(payload length).
//OUTPUT: TIME_MS(flying time in millisec.)
double time_on_air(LoRaPHY* Lora, int64_t plen){
double time_ms = 0.0;
double sym_num = calc_sym_num(Lora, plen);
time_ms = ((sym_num + 4.25 + Lora->preamble_len) * ((double)pow(2,Lora->sf) / (double)(Lora->bw)) * 1000.0);
//printf("\n%f", time_ms);
 return time_ms;
}
