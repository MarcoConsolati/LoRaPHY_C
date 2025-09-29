## LoRaPHY_C
LoRaPHY_C is a complete C language implementation of LoRa physical layer, including baseband modulation, baseband demodulation, encoding and decoding.
LoRaPHY_C is organized as some file.c, file.h and a linked FFTW library for Fast Fourier Transform.

This repository is the implementation of the following paper:

Zhenqiang Xu, Pengjin Xie, Shuai Tong, Jiliang Wang. From Demodulation to Decoding: Towards Complete LoRa PHY Understanding and Implementation. ACM Transactions on Sensor Networks 2022.  
Avaible at: https://dl.acm.org/doi/10.1145/3546869

And this repository is a complete translation from the Matlab Project [LoRaPHY] created by [jkadbear] avaible at:  https://github.com/jkadbear/LoRaPHY

## Components
- LoRa Modulator
- LoRa Demodulator
- LoRa Encoder
- LoRa Decoder

## Supported features
- Extremely low SNR demodulation (-20 dB)
- Clock drift correction
- All spreading factors (SF = 7,8,9,10,11,12)
- All code rates (CR = 4/5,4/6,4/7,4/8)
- Explicit/Implicit PHY header mode
- PHY header/payload CRC check
- Low Data Rate Optimization (LDRO)

## How to compile and execute
Git clone this repo or download all the file.c and the file.h in the same directory, create a new project with your preferred IDE and link them all togheter.
For the Fast Fourier Transform I used the FFTW avaible for different O.S. at the link: https://fftw.org/  
You'll need to download the file.h of the library(fftw3.h) compatible with your O.S. in the same directory and link it to the project.  
For calculating the CRC checksum I used another library already included in this repo called "crc_ccitt" created by Lammert Bies and avaible at: https://github.com/lammertb/libcrc
 
## Prerequisites
- GCC compiler
- Code::Blocks (or others IDE)
- Git

## Example of usage
Below is a complete example showing how to build a LoRaPHY object, encode and modulate a payload, then demodulate and decode it. This demonstrates the full transmission and reception chain using the provided API.

```c
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

int main(int argc, char** argv) {
    LoRaPHY Lora;
    buildLora(&Lora, 868.1e6, 7 , 125e3, 2e6); /*This function build Lora object struct with four parameters: Carrier frequency, Spreading factor, Bandwidth, Sampling frequency*/
    int64_t phy_payload1[] = {0x40, 0x11 , 0x22 , 0x33 , 0x44  , 0x00 , 0x00 , 0x00 , 0x00 , 0x00 };
    int64_t data[] = {0x61 , 0x61 , 0x61 , 0x61, 0x61 , 0x61 , 0x61 , 0x61 , 0x61 , 0x61 };
    int64_t length_data = sizeof(data) / sizeof(data[0]);
    int64_t length_phy = sizeof(phy_payload1) / sizeof(phy_payload1[0]);
    int64_t* symbols = (int64_t*)malloc((length_data + length_phy)*sizeof(int64_t));
    Lora.cr = 4; /*To use the encode,decode,modulate,demodulate function you'll need to set before this parameters: Code rate, Crc enabled/disabled, Header enabled*/
    Lora.crc = 1;
    Lora.has_header = 1;
    Lora.payload_len = length_data + length_phy;
    for(int64_t i = 0 ; i < length_data + length_phy ; i++){
        if( i < length_phy){
            symbols[i] = phy_payload1[i];
        }else{
            symbols[i] = data[i - length_phy];
        }
    }
    info_Lora(&Lora);
    ResVetInt symbols_enc = encode(&Lora, symbols, Lora.payload_len);
    ResVetCmplx result = modulate(&Lora, symbols_enc.vet , symbols_enc.length_vet);
    ResVetInt symbols_dem = demodulate(&Lora, result.vet , result.size_chirp);
    ResVetInt symbols_dec = decode(&Lora, symbols_dem.vet , 1 , symbols_dem.length_vet);

    free(symbols_enc.vet); symbols_enc.vet = NULL;
    free(result.vet); result.vet = NULL;
    free(symbols); symbols = NULL;
    free(symbols_dem.vet); symbols_dem.vet = NULL;
    free(symbols_dec.vet); symbols_dec.vet = NULL;

    if(symbols_dem.vet == NULL && symbols_enc.vet == NULL && result.vet == NULL && symbols_dec.vet == NULL) {
        printf("\nMemory has been freed correctly.\n");
    } else {
        printf("\nMemory has not been freed correctly.\n");
    }
    freeLora(&Lora);
    return 0;
}
```

## Manual compiling from Linux terminal
```bash
gcc yourcode.c utility.c main.c lora_phy.c crc_ccitt.c -lm -lfftw3 -o yourcode
```

## Credits

This project is a translation of [LoRaPHY] by [jkadbear], originally written in MATLAB.  
Original repository: https://github.com/jkadbear/LoRaPHY

## License

This project is licensed under the MIT License. See the LICENSE file for details.
