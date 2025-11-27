# LoRaPHY C Implementation

![C Language](https://img.shields.io/badge/C-ANSI_C-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20Embedded-orange)

**LoRa / LoRaWAN Physical Layer Implementation in C** - Complete baseband modulation / demodulation with FFT for Software Defined Radio (SDR) and embedded systems.

## 📖 Overview

A portable C implementation of the LoRa physical layer, featuring baseband modulation, demodulation, encoding and decoding. Designed for **LoRaWAN networks**, **SDR platforms**, and **embedded IoT devices**.

## 🔬 Academic Reference

This repository implements the research from:
> Zhenqiang Xu, Pengjin Xie, Shuai Tong, Jiliang Wang. *"From Demodulation to Decoding: Towards Complete LoRa PHY Understanding and Implementation"*. ACM Transactions on Sensor Networks 2022.

Available at: https://dl.acm.org/doi/10.1145/3546869

## 🔄 MATLAB Translation

This is a complete C translation of the MATLAB project **[LoRaPHY](https://github.com/jkadbear/LoRaPHY)** by [jkadbear](https://github.com/jkadbear), organized into modular `.c`/`.h` files with FFTW library integration for Fast Fourier Transform.

**Original MATLAB Repository:** https://github.com/jkadbear/LoRaPHY

## ✨ Features & Components

### 🎛️ Core Components
- **LoRa Modulator** - Baseband signal generation
- **LoRa Demodulator** - Robust signal reception
- **LoRa Encoder** - Forward error correction
- **LoRa Decoder** - Data recovery

### 🚀 Advanced Capabilities
- **Extremely low SNR demodulation** (-20 dB)
- **Clock drift correction** for timing recovery
- **Full spreading factor support** (SF7-SF12)
- **Multiple code rates** (4/5 to 4/8)
- **PHY header modes** (Explicit/Implicit)
- **CRC verification** for header and payload
- **Low Data Rate Optimization** (LDRO)


## 🛠️ Installation & Compilation

### Prerequisites
- **FFTW3 library** for Fast Fourier Transform ([download here](https://fftw.org/))
- **C compiler** (GCC recommended)
- **CRC library** (included via `crc_ccitt` from [Lammert Bies' libcrc](https://github.com/lammertb/libcrc))
- Ubuntu 20.04+ or compatible Linux
- MATLAB 2023+ (for reference comparison)
- Code::Blocks for Windows10 (or others IDE)


### Using Command Line (Linux/Unix)
```bash
# Clone the repository
git clone https://github.com/MarcoConsolati/LoRaPHY_C.git
cd LoRaPHY_C

# Install FFTW3 (Ubuntu/Debian example)
sudo apt-get install libfftw3-dev

# Compile the project
gcc -o lora_phy *.c -lfftw3 -lm

# Run the executable
./lora_phy
```

## Manual compiling from Linux terminal
```bash
gcc main.c lora_phy.c utility.c crc_ccitt.c -lfftw3 -lm -o lora_phy_demo
```

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

## Credits

This project is a translation of [LoRaPHY](https://github.com/jkadbear/LoRaPHY) by [jkadbear](https://github.com/jkadbear), originally written in MATLAB.  
Original repository: https://github.com/jkadbear/LoRaPHY

## License

This project is licensed under the MIT License. See the LICENSE file for details.
