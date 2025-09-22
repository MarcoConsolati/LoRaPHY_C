#include "lora_phy.h"
#include "crc_ccitt.h"
#include "utility.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>
#include <complex.h>
#include <inttypes.h>
#include <string.h>


///TEST DEMODULATE, DECODE///

int main(int argc, char** argv) {
	LoRaPHY Lora;
    buildLora(&Lora, 868.1e6, 7 , 125e3, 2e6);
    ///These bytes should be: 1byte(MHDR)+7bytes(MAC Payload)+4bytes(MIC)
    int64_t phy_payload1[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x00 , 0x00 , 0x00 , 0x00 , 0x00 };
    int64_t data[] = {0x61 , 0x61 , 0x61 , 0x61, 0x61 , 0x61 , 0x61 , 0x61 , 0x61 , 0x61 };
    int64_t length_data = sizeof(data) / sizeof(data[0]);
    int64_t length_phy = sizeof(phy_payload1) / sizeof(phy_payload1[0]);
    int64_t* symbols = (int64_t*)malloc((length_data + length_phy)*sizeof(int64_t));
    Lora.cr = 4;
    Lora.crc = 1;
    Lora.has_header = 1;
    Lora.payload_len = length_data + length_phy;
    for(int64_t i = 0 ; i < length_data + length_phy ; i++){
        if( i < length_phy){
            symbols[i] = phy_payload1[i];
        }else{
            symbols[i] = data[i - length_phy];
             }
//        printf("symbols[%"PRIu64"] = 0x%2"PRIu64"\n" , i , symbols[i]);
    }
    info_Lora(&Lora);
    ///This is the Modulation part of the symbols data
    ResVetInt symbols_enc = encode(&Lora, symbols, Lora.payload_len);   ///TEST ENCODE, MODULATE///
    ResVetCmplx result;
    result = modulate(&Lora, symbols_enc.vet ,   symbols_enc.length_vet);

    ///This is the Demodulation part of the symbols
    ResVetInt symbols_dem;
    symbols_dem = demodulate(&Lora, result.vet , result.size_chirp);///Demodulate from vector "vet" in the struct "result"///
    ResVetInt symbols_dec;
    symbols_dec = decode(&Lora, symbols_dem.vet , 1 , symbols_dem.length_vet);

    ///free all the pointer allocated dynamically
    free(symbols_enc.vet);
    symbols_enc.vet = NULL;
    free(result.vet);
    result.vet = NULL;
    free(symbols);
    symbols = NULL;
    free(symbols_dem.vet);
    symbols_dem.vet = NULL;
    free(symbols_dec.vet);
    symbols_dec.vet = NULL;
    if(symbols_dem.vet == NULL && symbols_enc.vet == NULL && result.vet == NULL && symbols_dec.vet == NULL) {
        printf("\nMemory has been freed correctly.\n");
    } else {
        printf("\nMemory has not been freed correctly.\n");
    }
    freeLora(&Lora);
    return 0;
}



    ///TEST DEMODULATE, DECODE///
    ///Demodulate from a file of IQ samples///
/*    ResFileCmplx sig_complex;
    sig_complex = readCmplx("./sig.cfile");
    symbols_dem = demodulate(&Lora, sig_complex.vet , sig_complex.size_file);
    stampaVetCmplx(sig_complex.vet , sig_complex.size_file);
    printf("\n\nSymbols demodulated from sig.cfile: ");
    stampaVetInt(symbols_dem.vet,symbols_dem.length_vet);
    free(sig_complex.vet);
    free(symbols_dem.vet);*/


//    int64_t phy_payload2[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x20 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload3[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x40 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload4[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x60 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload5[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x80 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload6[] = {0x80, 0x11 , 0x22 , 0x33 , 0x44  , 0x00 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload7[] = {0x80, 0x11 , 0x22 , 0x33 , 0x44  , 0x20 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload8[] = {0x80, 0x11 , 0x22 , 0x33 , 0x44  , 0x40 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload9[] = {0x80, 0x11 , 0x22 , 0x33 , 0x44  , 0x60 , 0x00 , 0x00 , 0x00 , 0x00 };    /* + FRMPayload + MIC(4bytes/4ottects) */
//    int64_t phy_payload10[] = {0x80, 0x11 , 0x22 , 0x33 , 0x44  , 0x80 , 0x00 , 0x00 , 0x00 , 0x00 };   /* + FRMPayload + MIC(4bytes/4ottects) */

//    stampaVetCmplx( result.vet , result.size_chirp);

/*    int64_t v = 0;
    const char filename[] = "file_written_cmplx.bin";
    v = writeCmplxFile(result.vet, result.size_chirp , filename);
    if(v){
        printf("\nDati scritti con successo nel file %s.\n", filename);
    }else{
        printf("\nErrore durante la scrittura nel file %s.\n", filename);
    }
    double complex* vet = (double complex*)calloc( result.size_chirp, sizeof(double complex));
    readCmplxFile(vet, result.size_chirp, filename);*/
//    stampaVetCmplx(vet, result.size_chirp);




    ///TEST TIME_ON_AIR///
/*    int64_t plen = 10;
    double t = 0.0;
    Lora.cr = 4; //se non si imposta il code rate non ha senso!
    t = time_on_air(&Lora, plen);
    printf("\nt= %f (in milliseconds)",t);*/

    ///TEST DEWHITEN, WHITEN///
    /*    int64_t* bytes ;
    int64_t bytes_w[] = {254,252,255,252,245,231,197,141,2,29};
    bytes = bytes_w;
    int64_t** matrice1;
    matrice1 = whiten(&Lora, symbols_enc.vet,  symbols_enc.length_vet);
    int64_t* data_dew ;
    data_dew = dewhiten(&Lora, bytes,   symbols_enc.length_vet);
//    printf("\nMatrice di dati dopo il whitening: \n");
    for (int64_t i = 0; i < 10; i++) {
        for (int64_t j = 0; j < 10; j++) {
//            printf("%"PRId64"\t", matrice1[i][j]);
        }
//        printf("\n");
    }
//    printf("\nDati dopo il dewhitening: \n");
    for (int64_t i = 0; i < 10; i++) {
//        printf("%"PRId64"\t", data_dew[i]);
    }*/

    ///TEST CRC///
/*    int64_t data[]={1,2,3,4,5,6,7,8,9,10};
    int64_t* dataptr;
    dataptr = data;
	int64_t dataLength = (sizeof(data) / sizeof(data[0]));
	calc_crc(&Lora, dataptr , dataLength);*/


    ///TEST PARSE_HEADER///
/*    int64_t z;
    int64_t* data;
    int64_t symbols[] = {97,57,1,77,33,29,61,125,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,34,61,40,36,61,17,119,94,85,86,107,43,86,78,90,83};
    int64_t lengthdata = 40;
    data = symbols;
    z = parse_header(&Lora,data,lengthdata);
    if(z == 1)
        printf("\nValid z=%"PRId64"",z);
    else
        printf("\nInvalid z=%"PRId64"",z);*/


    ///TEST GRAY CODING, HAMMING_DECODE///
/*    int64_t dinint[] = {97,57,1,77,33,29,61,125,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,34,61,40,36,61,17,119,94,85,86,107,43,86,78,90,83};
    int64_t* din;
    din = dinint;
    int64_t lengthdin = sizeof(dinint) / sizeof(dinint[0]);
    int64_t* symbols = gray_coding(&Lora , din , lengthdin );
    printf("\n\nSymbols after Gray coding: ");
    for(int64_t i = 0; i < lengthdin; i++){
        printf("\nSymbols[%2"PRIu64"]=%3"PRIu64"", i , symbols[i]);
    }*/


    ///TEST DIAG_DEINTERLEAVING, HAMMING_DECODE, HAMMING_ENCODE///
    ///Before Hamming_encode you must set the code rate!///
/*    int64_t* codewords;
    int64_t symb[] = {58,93,127,127,119,109,125,58};
    int64_t* symbolsde;
    symbolsde = symb;
    codewords = diag_deinterleave(&Lora, symbolsde , 8 ,7);
    int64_t* nibbles = hamming_decode(&Lora, codewords , Lora.sf - 2 ,8);
    printf("\n\nNibbles vale:\n");
    for(int64_t i = 0; i < Lora.sf - 2 ; i++){
        printf(" %"PRId64"",nibbles[i]);
    }
    int64_t data_nibbles[] = {0,10,9,0,11,14,15,12,15,15,15,12,15,5,15,7,14,5,12,13,8,2,0,13,1,6,10,15,7,15,15,15,15,15,15,15,15,15};
    int64_t* head_nibbles;
    head_nibbles = data_nibbles;
    int64_t lengthnib = sizeof(data_nibbles) / sizeof(data_nibbles[0]);
    Lora.cr = 4;
    int64_t* wordscode = hamming_encode(&Lora, head_nibbles, lengthnib );
    printf("\n\nWords after Hamming_encode with Code rate = %"PRId64", are: \n", Lora.cr);
    for(int64_t i = 0 ; i < lengthnib ; i++){
        printf("%"PRId64" ", wordscode[i]);
    }*/


    ///TEST GRAY_DECODING, CALC_PAYLOAD_LEN///
/*    int64_t symbols_i[] = {20,9,0,26,12,4,8,16,58,93,127,127,119,109,125,58,47,78,127,90,108,93,79,5,49,34,52,50,34,24,77,115,126,127,95,63,127,107,117,123};
    int64_t* symboli;
    symboli = symbols_i;
    int64_t lengthsymbols = sizeof(symbols_i) / sizeof(symbols_i[0]);
    int64_t* packet = gray_decoding(&Lora, symboli, lengthsymbols);
    ///Before calc_payload_len you must set the code rate properly!///
    Lora.cr = 1;
    int64_t len = calc_payload_len(&Lora, lengthsymbols , true);*/


    ///TEST DIAG_INTERLEAVE///
/*    int64_t words[] = {6,9,5,25,5,10,254};
    int64_t* codewords;
    codewords = words;
    int64_t* symbols = diag_interleave(&Lora, codewords, 7,8);*/


    ///Per usare correttamente ENCODE bisogna///
    ///impostare i valori di Lora.cr,Lora.crc,Lora.has_header!
    ///TEST ENCODE///
/*    int64_t nibbles[] = {1,2,3,4,5,6,7,8,9,10};
    int64_t* symbols;
    symbols = nibbles;
    int64_t lengthsymb = sizeof(nibbles) / sizeof(nibbles[0]);
    Lora.sf = 8;
    Lora.cr = 4;
    Lora.crc = 1;
    Lora.has_header = 1;
    int64_t* symbols_enc = encode(&Lora, symbols, lengthsymb);
    for(int64_t k = 0; k < Lora.vet_size ; k++){
        printf("\nsymbols_enc[%"PRId64"] = %"PRId64"", k , symbols_enc[k]);
    }*/


    ///Per usare correttamente DECODE bisogna///
    ///impostare i valori di Lora.cr,Lora.crc,Lora.has_header!
    ///TEST DECODE///
/*    int64_t* symbols_m;
//14    int64_t symb_m[] = {97,57,1,77,33,29,61,125,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,34,61,40,36,61,17,119,94,85,86,107,43,86,78,90,83}; ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=4)///
//04    int64_t symb_m[] = {1,5,1,65,29,29,13,121,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,95,61,40,45,60,20,119,94};    ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=4)///
//13    int64_t symb_m[] = {29,53,5,49,37,1,1,97,45,106,86,86,91,74,87,54,117,86,109,73,106,118,34,61,40,36,61,17,119,85,86,107,43,86,78,90};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=3)///
//03    int64_t symb_m[] = {125,9,5,61,25,1,49,101,45,106,86,86,91,74,87,54,117,86,109,73,106,118,95,61,40,45,60,20,119};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=3)///
//12    int64_t symb_m[] = {97,5,25,49,29,29,13,97,45,106,86,86,91,74,54,117,86,109,73,106,34,61,40,36,61,17,85,86,107,43,86,78};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=2)///
//02    int64_t symb_m[] = {1,57,25,61,33,29,61,101,45,106,86,86,91,74,54,117,86,109,73,106,95,61,40,45,60,20}; ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=2)///
//11   int64_t symb_m[] = {97,9,1,49,25,97,1,121,45,106,86,86,16,54,117,86,109,127,34,61,40,36,46,85,86,107,43,64};    //DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=1)///
//01    int64_t symb_m[] = {1,53,1,61,37,97,49,125,45,106,86,86,16,54,117,86,109,127,95,61,40,45,46};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=1)///
//004    int64_t symb_m[] = {77,89,85,85,89,73,41,109,42,84,86,84,46,92,44,26,122,16,95,79,31,72,37,120,86,43,107,75,86,83,86,86};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,hasheader=0,cr=4)///
    int64_t lengthcols = 1;
    symbols_m = symb_m;
    int64_t length_symb_m = sizeof(symb_m) / sizeof(symb_m[0]);
    int64_t** matrice = (int64_t**)malloc(length_symb_m * sizeof(int64_t*));
    if (matrice == NULL) {
        perror("Errore nell'allocazione delle righe");
        exit(EXIT_FAILURE);
    }
    for (int64_t i = 0; i < length_symb_m; i++) {
        matrice[i] = (int64_t*)malloc(lengthcols * sizeof(int64_t));
        if (matrice[i] == NULL) {
            perror("Errore nell'allocazione delle colonne");
            exit(EXIT_FAILURE);
        }
    }
    for(int64_t i = 0; i < length_symb_m ; i++){
        for( int64_t j = 0; j < lengthcols ; j++){
            matrice[i][j] = symb_m[i];
//            printf("%3d ", matrice[i][j]);
        }
//         printf("\n");
    }
//    Lora.sf = 12;
    Lora.cr = 4;
    Lora.crc = 1;
    Lora.has_header = 1;
    decode(&Lora, matrice , lengthcols , length_symb_m);*/



    ///Per usare correttamente TEST SYMBOLS_TO_BYTES bisogna
    ///impostare i valori di Lora.cr,Lora.has_header,Lora.crc!
    ///ovviamente i simboli da decodif./codific. devono corrispondere
    ///ai parametri di Lora con cui sono stati spediti dal TX!
    ///TEST SYMBOLS_TO_BYTES///
//    int64_t* symbols;
/*14*///    int64_t symb[] = {97,57,1,77,33,29,61,125,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,34,61,40,36,61,17,119,94,85,86,107,43,86,78,90,83}; ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=4)///
/*04*///    int64_t symb[] = {1,5,1,65,29,29,13,121,45,106,86,86,91,74,87,45,54,117,86,109,73,106,118,7,95,61,40,45,60,20,119,94};    ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=4)///
/*13*///    int64_t symb[] = {29,53,5,49,37,1,1,97,45,106,86,86,91,74,87,54,117,86,109,73,106,118,34,61,40,36,61,17,119,85,86,107,43,86,78,90};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=3)///
/*03*///    int64_t symb[] = {125,9,5,61,25,1,49,101,45,106,86,86,91,74,87,54,117,86,109,73,106,118,95,61,40,45,60,20,119};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=3)///
/*12*///    int64_t symb[] = {97,5,25,49,29,29,13,97,45,106,86,86,91,74,54,117,86,109,73,106,34,61,40,36,61,17,85,86,107,43,86,78};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=2)///
/*02*///    int64_t symb[] = {1,57,25,61,33,29,61,101,45,106,86,86,91,74,54,117,86,109,73,106,95,61,40,45,60,20}; ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=2)///
/*11*///    int64_t symb[] = {97,9,1,49,25,97,1,121,45,106,86,86,16,54,117,86,109,127,34,61,40,36,46,85,86,107,43,64};    //DATA=1,2,3,4,5,6,7,8,9,10(crc=1,cr=1)///
/*01*///    int64_t symb[] = {1,53,1,61,37,97,49,125,45,106,86,86,16,54,117,86,109,127,95,61,40,45,46};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,cr=1)///
/*004*///   int64_t symb[] = {77,89,85,85,89,73,41,109,42,84,86,84,46,92,44,26,122,16,95,79,31,72,37,120,86,43,107,75,86,83,86,86};   ///DATA=1,2,3,4,5,6,7,8,9,10(crc=0,hasheader=0,cr=4)///
/*    symbols = symb;
    int64_t length_symb = sizeof(symb) / sizeof(symb[0]);
    int64_t plen = Lora.payload_len;
    int64_t* bytes_dec = (int64_t*)calloc(plen, sizeof(int64_t));
    ///prima di usare symbol_to_bytes bisogna impostare i valori per i parametri: Cr,Crc,Has_header!!!
    Lora.has_header = 1;
    Lora.crc = 1;
    Lora.cr = 4;
    bytes_dec = symbols_to_bytes(&Lora, symbols, length_symb);*/


    ///Per usare correttamente GEN_HEADER bisogna
    ///impostare i valori di Lora.cr e Lora.crc!
    ///TEST GEN_HEADER///
/*    int64_t plen = 10;
    Lora.crc = 1;
    Lora.cr = 1;
    int64_t* nibbles = (int64_t*)calloc(5,sizeof(int64_t));
    nibbles = gen_header(&Lora, plen);
    printf("\nNibbles_bits:");
    for(int64_t i = 0; i < 5 ; i++){
        printf("%"PRId64" ",nibbles[i]);
    }*/


    ///TEST WORD_REDUCE///
/*    int64_t* ws;
    int64_t words[] = {128,0,32,16,4};
    ws = words;
    int64_t result = word_reduce(ws, 5);
    printf("\nresult=%"PRId64"",result);*/



///TCP/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*    If we are creating a connection between client and server using TCP then it has a few functionalities like,
    TCP is suited for applications that require high reliability, and transmission time is relatively less critical.
    It is used by other protocols like HTTP, HTTPs, FTP, SMTP, Telnet. TCP rearranges data packets in the order specified.
    There is absolute guarantee that the data transferred remains intact and arrives in the same order in which it was sent.
    TCP does Flow Control and requires three packets to set up a socket connection before any user data can be sent.
    TCP handles reliability and congestion control. It also does error checking and error recovery.
    Erroneous packets are retransmitted from the source to the destination.

    TCP Server –

    1-using create(), Create TCP socket.
    2-using bind(), Bind the socket to server address.
    3-using listen(), put the server socket in a passive mode, where it waits for the client to approach the server to make a connection
    4-using accept(), At this point, connection is established between client and server, and they are ready to transfer data.
    5-Go back to Step 3.

    TCP Client –

    1-Create TCP socket.
    2-connect newly created client socket to server.


///TCP CLIENT CODE///////////////////////////////////////////////////////////////////////////////////////


#include <arpa/inet.h> // inet_addr()
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> // bzero()
#include <sys/socket.h>
#include <unistd.h> // read(), write(), close()
#define MAX 80
#define PORT 8080
#define SA struct sockaddr
void func(int sockfd)
{
	char buff[MAX];
	int n;
	for (;;) {
		bzero(buff, sizeof(buff));
		printf("Enter the string : ");
		n = 0;
		while ((buff[n++] = getchar()) != '\n')
			;
		write(sockfd, buff, sizeof(buff));
		bzero(buff, sizeof(buff));
		read(sockfd, buff, sizeof(buff));
		printf("From Server : %s", buff);
		if ((strncmp(buff, "exit", 4)) == 0) {
			printf("Client Exit...\n");
			break;
		}
	}
}

int main()
{
	int sockfd, connfd;
	struct sockaddr_in servaddr, cli;

	// socket create and verification
	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	if (sockfd == -1) {
		printf("socket creation failed...\n");
		exit(0);
	}
	else
		printf("Socket successfully created..\n");
	bzero(&servaddr, sizeof(servaddr));

	// assign IP, PORT
	servaddr.sin_family = AF_INET;
	servaddr.sin_addr.s_addr = inet_addr("127.0.0.1");
	servaddr.sin_port = htons(PORT);

	// connect the client socket to server socket
	if (connect(sockfd, (SA*)&servaddr, sizeof(servaddr))
		!= 0) {
		printf("connection with the server failed...\n");
		exit(0);
	}
	else
		printf("connected to the server..\n");

	// function for chat
	func(sockfd);

	// close the socket
	close(sockfd);
}

*/

