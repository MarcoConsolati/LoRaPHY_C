#ifndef CRC_CCITT_H_INCLUDED
#define CRC_CCITT_H_INCLUDED

/* Library: libcrc
 * File:    src/crcccitt.c
 * Author:  Lammert Bies
 *
 * This file is licensed under the MIT License as stated below
 *
 * Copyright (c) 1999-2016 Lammert Bies
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Description
 * -----------
 * The module src/crcccitt.c contains routines which are used to calculate the
 * CCITT CRC values of a string of bytes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include "lora_phy.h"
#include <inttypes.h>
/*
 * #define CRC_POLY_xxxx
 *
 * The constants of the form CRC_POLY_xxxx define the polynomials for some well
 * known CRC calculations.
 */

#define		CRC_POLY_16			0xA001
#define		CRC_POLY_32			0xEDB88320ul
#define		CRC_POLY_64			0x42F0E1EBA9EA3693ull
#define		CRC_POLY_CCITT		0x1021
#define		CRC_POLY_DNP		0xA6BC
#define		CRC_POLY_KERMIT		0x8408
#define		CRC_POLY_SICK		0x8005

/*
 * #define CRC_START_xxxx
 *
 * The constants of the form CRC_START_xxxx define the values that are used for
 * initialization of a CRC value for common used calculation methods.
 */

#define		CRC_START_8			0x00
#define		CRC_START_16		0x0000
#define		CRC_START_MODBUS	0xFFFF
#define		CRC_START_XMODEM	0x0000
#define		CRC_START_CCITT_1D0F	0x1D0F
#define		CRC_START_CCITT_FFFF	0xFFFF
#define		CRC_START_KERMIT	0x0000
#define		CRC_START_SICK		0x0000
#define		CRC_START_DNP		0x0000
#define		CRC_START_32		0xFFFFFFFFul
#define		CRC_START_64_ECMA	0x0000000000000000ull
#define		CRC_START_64_WE		0xFFFFFFFFFFFFFFFFull


/*
 * uint16_t crc_xmodem( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_xmodem() performs a one-pass calculation of an X-Modem CRC
 * for a byte string that has been passed as a parameter.
 */

uint16_t crc_xmodem( const unsigned char *input_str, size_t num_bytes );

/* crc_xmodem */

/*
 * uint16_t crc_ccitt_1d0f( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_ccitt_1d0f() performs a one-pass calculation of the CCITT
 * CRC for a byte string that has been passed as a parameter. The initial value
 * 0x1d0f is used for the CRC.
 */

uint16_t crc_ccitt_1d0f( const unsigned char *input_str, size_t num_bytes );
/* crc_ccitt_1d0f */
/*
 * uint16_t crc_ccitt_ffff( const unsigned char *input_str, size_t num_bytes );
 *
 * The function crc_ccitt_ffff() performs a one-pass calculation of the CCITT
 * CRC for a byte string that has been passed as a parameter. The initial value
 * 0xffff is used for the CRC.
 */

uint16_t crc_ccitt_ffff( const unsigned char *input_str, size_t num_bytes );

uint16_t crc_ccitt_generic( const unsigned char *input_str, size_t num_bytes, uint16_t start_value );

uint16_t update_crc_ccitt( uint16_t crc, unsigned char c );

void init_crcccitt_tab( void );

//FUNCTION THAT CALCULATE CHECKSUM FROM DATA SEQUENCE//
int64_t* calc_crc(LoRaPHY* Lora, int64_t* data, int64_t dataLength);
#endif // CRC_CCITT_H_INCLUDED
