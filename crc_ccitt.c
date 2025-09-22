#include "crc_ccitt.h"
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

static bool             crc_tabccitt_init       = false;
static uint16_t         crc_tabccitt[256];

/*
 * uint16_t crc_xmodem( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_xmodem() performs a one-pass calculation of an X-Modem CRC
 * for a byte string that has been passed as a parameter.
 */

uint16_t crc_xmodem( const unsigned char *input_str, size_t num_bytes ) {

	return crc_ccitt_generic( input_str, num_bytes, CRC_START_XMODEM );

}  /* crc_xmodem */

/*
 * uint16_t crc_ccitt_1d0f( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_ccitt_1d0f() performs a one-pass calculation of the CCITT
 * CRC for a byte string that has been passed as a parameter. The initial value
 * 0x1d0f is used for the CRC.
 */

uint16_t crc_ccitt_1d0f( const unsigned char *input_str, size_t num_bytes ) {

	return crc_ccitt_generic( input_str, num_bytes, CRC_START_CCITT_1D0F );

}  /* crc_ccitt_1d0f */

/*
 * uint16_t crc_ccitt_ffff( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_ccitt_ffff() performs a one-pass calculation of the CCITT
 * CRC for a byte string that has been passed as a parameter. The initial value
 * 0xffff is used for the CRC.
 */

uint16_t crc_ccitt_ffff( const unsigned char *input_str, size_t num_bytes ) {

	return crc_ccitt_generic( input_str, num_bytes, CRC_START_CCITT_FFFF );

}  /* crc_ccitt_ffff */

/*
 * static uint16_t crc_ccitt_generic( const unsigned char *input_str, size_t num_bytes, uint16_t start_value );
 *
 * The function crc_ccitt_generic() is a generic implementation of the CCITT
 * algorithm for a one-pass calculation of the CRC for a byte string. The
 * function accepts an initial start value for the crc.
 */

uint16_t crc_ccitt_generic( const unsigned char *input_str, size_t num_bytes, uint16_t start_value ) {

	uint16_t crc;
	const unsigned char *ptr;
	size_t a;

	if ( ! crc_tabccitt_init ) init_crcccitt_tab();

	crc = start_value;
	ptr = input_str;

	if ( ptr != NULL ) for (a=0; a<num_bytes; a++) {

		crc = (crc << 8) ^ crc_tabccitt[ ((crc >> 8) ^ (uint16_t) *ptr++) & 0x00FF ];
	}

	return crc;

}  /* crc_ccitt_generic */

/*
 * uint16_t update_crc_ccitt( uint16_t crc, unsigned char c );
 *
 * The function update_crc_ccitt() calculates a new CRC-CCITT value based on
 * the previous value of the CRC and the next byte of the data to be checked.
 */

uint16_t update_crc_ccitt( uint16_t crc, unsigned char c ) {

	if ( ! crc_tabccitt_init ) init_crcccitt_tab();

	return (crc << 8) ^ crc_tabccitt[ ((crc >> 8) ^ (uint16_t) c) & 0x00FF ];

}  /* update_crc_ccitt */

/*
 * static void init_crcccitt_tab( void );
 *
 * For optimal performance, the routine to calculate the CRC-CCITT uses a
 * lookup table with pre-compiled values that can be directly applied in the
 * XOR action. This table is created at the first call of the function by the
 * init_crcccitt_tab() routine.
 */

void init_crcccitt_tab( void ) {
	uint16_t i;
	uint16_t j;
	uint16_t crc;
	uint16_t c;
	for (i=0; i<256; i++) {
		crc = 0;
		c   = i << 8;
		for (j=0; j<8; j++) {
			if ( (crc ^ c) & 0x8000 ) crc = ( crc << 1 ) ^ CRC_POLY_CCITT;
			else                      crc =   crc << 1;
		c = c << 1;
		}
	crc_tabccitt[i] = crc;
	}
	crc_tabccitt_init = true;

}  /* init_crcccitt_tab */


int64_t* calc_crc(LoRaPHY* Lora, int64_t* data , int64_t dataLength){
    int64_t* check = (int64_t*)malloc(2 * sizeof(int64_t));
	uint64_t Crc_int = 0;
	int64_t len = dataLength - 2 ;
	unsigned char* input = (unsigned char*)malloc(len* sizeof(unsigned char));
	for(int64_t i = 0 ; i < len ; i++){
		input[i] = data[i];
	}
	Crc_int = crc_xmodem(input, len);
	unsigned char CrcByte1 = (unsigned char)(Crc_int & 0xFF);			/// Takes the least significant byte of the CRC_int
    unsigned char CrcByte2 = (unsigned char)((Crc_int >> 8) & 0xFF);	/// Takes the most significant byte of the CRC_int
	check[0] = CrcByte1 ^ data[len + 1];
	check[1] = CrcByte2 ^ data[len];
	printf("\nCRC: 0x%"PRIX64"\nChecksum: %"PRId64"\t%"PRId64"\n", Crc_int , check[0] , check[1]);
    free(input);
    input = NULL;
 return check;
}
