#LoRaPHY_C
LoRaPHY_C is a C language complete implementation of LoRa physical layer, including baseband modulation, baseband demodulation, encoding and decoding.
LoRaPHY_C is organized as some file.c and  file.h and a linked FFTW library for the Fast Fourier Transform.

This repo is the implementation of the following paper:

Zhenqiang Xu, Pengjin Xie, Shuai Tong, Jiliang Wang. From Demodulation to Decoding: Towards Complete LoRa PHY Understanding and Implementation. ACM Transactions on Sensor Networks 2022. avaible at: https://dl.acm.org/doi/10.1145/3546869

And this repo is a complete translation from the Matlab Project LoRaPHY created by jkadbear avaible at: https://github.com/jkadbear/LoRaPHY

## Components
*LoRa Modulator
*LoRa Demodulator
*LoRa Encoder
*LoRa Decoder

## Supported feature
*Extremely low SNR demodulation (-20 dB)
*Clock drift correction
*All spreading factors (SF = 7,8,9,10,11,12)
*All code rates (CR = 4/5,4/6,4/7,4/8)
*Explicit/Implicit PHY header mode
*PHY header/payload CRC check
*Low Data Rate Optimization (LDRO)

## How compile and execute
You'll need to download all the file.c and the file.h in the same directory, create a new project with your IDE and link them all togheter.
For the Fast Fourier Transform I used the FFTW avaible for different O.S. at the link: https://fftw.org/
You'll need to download the file.h of the library(fftw3.h) in the same directory and link it to the project or the FFT in the Modulte/Demodulate functions won't work.

### Prerequisites
- GCC compiler
- CodeBlocks (or others IDE)
- Git

### Manual compiling
```bash
gcc -o main.c lora_phy.c utility.c crcccitt.c