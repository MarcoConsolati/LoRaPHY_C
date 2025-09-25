## LoRaPHY_C
LoRaPHY_C is a complete C language implementation of LoRa physical layer, including baseband modulation, baseband demodulation, encoding and decoding.
LoRaPHY_C is organized as some file.c, file.h and a linked FFTW library for Fast Fourier Transform.

This repo is the implementation of the following paper:

Zhenqiang Xu, Pengjin Xie, Shuai Tong, Jiliang Wang. From Demodulation to Decoding: Towards Complete LoRa PHY Understanding and Implementation. ACM Transactions on Sensor Networks 2022.
avaible at: https://dl.acm.org/doi/10.1145/3546869

And this repo is a complete translation from the Matlab Project LoRaPHY created by jkadbear avaible at:  https://github.com/jkadbear/LoRaPHY

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
You'll need to download all the file.c and the file.h in the same directory, create a new project with your preferred IDE and link them all togheter.
For the Fast Fourier Transform I used the FFTW avaible for different O.S. at the link: https://fftw.org/
You'll need to download the file.h of the library(fftw3.h) in the same directory and link it to the project.
For calculating the CRC checksum I used another library ready to use called "crc_ccitt" created by Lammert Bies and avaible at: https://github.com/lammertb/libcrc
 
### Prerequisites
- GCC compiler
- Code::Blocks (or others IDE)
- Git

### Manual compiling
```bash
gcc yourcode.c utility.c main.c lora_phy.c crc_ccitt.c -lm -lfftw3 -o yourcode
