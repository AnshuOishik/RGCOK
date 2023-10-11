***********************************RGCOK********************************************
	Reference Genome Compression Algorithm using Optimal k-mer length
			https://github.com/AnshuOishik/RGCOK
				Copyright (C) 2023
************************************************************************************
Information
1. Code is written in Notepad++ 
2. Implemented using java
3. Run on Linux operating system
4. The executable "bsc" should be put in the same directory of main class file
5. Set "bsc" mode "chmod 777"
************************************************************************************
Compilation Command:
javac -d . *.java

For RGCOK_Ubuntu_32_BH/RGCOK_Ubuntu_Free_BH/RGCOK_Windows_BH (BSC or (BSC+Huffman)):
Execution Command:
Compression:
Ubuntu:
Macro:
java -Xmx20g rgcok.RGCOK chr.fa comp comp.fa 8
Micro/Windows:
java rgcok.RGCOK chr.fa comp comp.fa 2

Decompression:
Ubuntu:
Macro:
java -Xmx20g  rgcok.RGCOK chr.fa decomp comp.fa 8
Micro/Windows:
java rgcok.RGCOK chr.fa decomp comp.fa 2

Where, "chr.fa" is the input file contains chromosomes file name
Argument "comp" for compression, "decomp" for decompression
"comp.fa" is the final compressed file 
2 or 8 is the number of threads
-Xmx20g to increase heap memory size

****************************************************************************************
For RGCOK_Ubuntu_32_HB/RGCOK_Ubuntu_Free_HB/RGCOK_Windows_HB (Huffman or (Huffman+BSC)):
Execution Command:
Compression:
Ubuntu:
java -Xmx20g rgcok.RGCOK chr.fa comp chr.id chr.rgcok 8
Windows:
java rgcok.RGCOK chr.fa comp chr.id chr.rgcok 2

Decompression:
Ubuntu:
java -Xmx20g  rgcok.RGCOK chr.fa decomp HCchr.id HCchr.rgcok
Windows:
java  rgcok.RGCOK chr.fa decomp HCchr.id HCchr.rgcok

Where, "chr.fa" is the input file contains chromosomes file name
Argument "comp" for compression, "decomp" for decompression
"chr.id" and "chr.rgcok" are the input files to Huffman encoding 
"HCchr.id" and "HCchr.rgcok" are Huffman encoded files
Final compressed file is "FinalBsc.bsc"
2 or 8 is the number of threads
-Xmx20g to increase heap memory size

*************************************************************************************
Execution Command (By Creating jar File):
Creating a manifest file:
echo Main-Class: RGCOK >manifest.txt
(echo Main-Class: MainClassName >manifest.txt)

Creating our jar file:
jar cvfm RGCOKTest.jar manifest.txt *.class
(jar cvfm JarFileName.jar manifest.txt *.class)

Executing jar File:
java -jar RGCOKTest.jar chr.fa comp 8

**********************************************************************
For bsc executable file generation from available code:
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
