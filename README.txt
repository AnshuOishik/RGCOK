************************************ MtK-mer & RGCOK ***************************************************************************
					MtK-mer: Multi-Threaded Optimal k-mer length
		RGCOK: Reference Genome Compression Algorithm using Optimal k-mer length
							https://github.com/AnshuOishik/RGCOK
								Copyright (C) 2023 
*****************************************************************************************************************************
Information
1. Code is written in Notepad++ 
2. Implemented using java
3. Run on Linux operating system
4. The executable "bsc" should be put in the same directory of main class file
5. Set "bsc" mode "chmod 777"
*****************************************************************************************************************************
:::MtK-mer:::
Compilation Command:
javac Main.java

Execution Command:
java -Xmx20g Main chr.fa 8 

Testing kMer length from 9 to 21

Where, chr.fa is the input file contains chromosomes
8 is the number of threads
-Xmx20g to increase heap memory size
*****************************************************************************************************************************
RGCOK:
Compilation Command:
javac -d . *.java

:::Using BSC or Huffman compressor in the final stage:::
Execution Command:
Compression:
java -Xmx20g rgcok.RGCOK chr.fa comp 8

Decompression:
java -Xmx20g  rgcok.RGCOK chr.fa decomp 8

Where, "chr.fa" is the input file contains chromosomes file names
Argument "comp" for compression, "decomp" for decompression
"BscC.bsc" is the final compressed file produced using BSC compressor
"HCchr.id" and "HCchr.rgcok" are the final compressed files produced using Huffman compressor
2 or 8 is the number of threads
-Xmx20g to increase heap memory size

*************************************************************************************************************************
Commands for "bsc" executable file generation from available code at
https://github.com/IlyaGrebnov/libbsc
Compilation commands:
g++ -c libbsc/adler32/adler32.cpp
g++ -c libbsc/bwt/libsais/libsais.c
g++ -c libbsc/bwt/bwt.cpp
g++ -c libbsc/coder/coder.cpp
g++ -c libbsc/coder/qlfc/qlfc.cpp
g++ -c libbsc/coder/qlfc/qlfc_model.cpp
g++ -c libbsc/filters/detectors.cpp
g++ -c libbsc/filters/preprocessing.cpp
g++ -c libbsc/libbsc/libbsc.cpp
g++ -c libbsc/lzp/lzp.cpp
g++ -c libbsc/platform/platform.cpp
g++ -c libbsc/st/st.cpp
g++ -c bsc.cpp

Linking command:
g++ -o bsc bsc.o adler32.o bwt.o coder.o detectors.o libbsc.o libsais.o lzp.o platform.o preprocessing.o qlfc.o qlfc_model.o st.o
*****************************************************************************************************************************
