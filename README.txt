				       ****************************** MtK-mer & RGCOK ****************************
							MtK-mer: Multi-Threaded Optimal k-mer Length
					RGCOK: Reference Genome Compression Algorithm using Optimal k-mer Length
							https://github.com/AnshuOishik/RGCOK
								Copyright (C) 2023 
=============================================================================================================================
Introduction
To utilize the code, please use the Notepad++ editor.
Java has been utilized by us in the implementation.
Please use Linux as your operating system.
Please confirm that the physical memory on your computer is larger than 10 GB.
=============================================================================================================================
MtK-mer: Multi-Threaded Optimal k-mer Length
The MtK-mer algorithm used the randomization method to determine the ideal k-mer length. RGCOK uses the ideal k-mer length discovered using MtK-mer to compress the particular sequence.

# Compilation Command:
> javac -d *.java

# Execution Command:
> java -Xms10240m mtkmer.Main chr.fa 9 8 

Notice:
# -Xms10240m is the initial allocation of memory (MB)
# The list of target file directories and the reference file path (the first line) are both found in chr.fa
# Nine is k's value. This value should be changed in order to check for various values of k (9 through 21)
# The number of threads is eight (4, by default, is the optional value)
=============================================================================================================================
RGCOK: Reference Genome Compression Algorithm using Optimal k-mer Length
The MtK-mer technique yields the ideal k-mer length, which is used by RGCOK to compress sequences.

Compilation Command:
> javac -d . *.java

Execution Command:
Compression:
> java -Xms10240m rgcok.RGCOK chr.fa comp 8
Decompression:
> java -Xms10240m rgcok.RGCOK chr.fa decomp 8

Notice:
# The list of target file directories and the reference file path (the first line) are both found in chr.fa
# "decomp" is the argument for decompression, and "comp" for compression
# The final compressed file created by the BSC compressor is called "BscC.bsc".
# The number of threads is eight (4, by default, is the optional value)
# -Xms10240m is the initial allocation of memory (MB)
4. Please place the executable "bsc" in the main class file's directory.
5. Kindly set "chmod 0777" for "bsc" mode.
=============================================================================================================================
Commands for "bsc" executable file generation from available code at https://github.com/IlyaGrebnov/libbsc
Compilation commands:
> g++ -c libbsc/adler32/adler32.cpp
> g++ -c libbsc/bwt/libsais/libsais.c
> g++ -c libbsc/bwt/bwt.cpp
> g++ -c libbsc/coder/coder.cpp
> g++ -c libbsc/coder/qlfc/qlfc.cpp
> g++ -c libbsc/coder/qlfc/qlfc_model.cpp
> g++ -c libbsc/filters/detectors.cpp
> g++ -c libbsc/filters/preprocessing.cpp
> g++ -c libbsc/libbsc/libbsc.cpp
> g++ -c libbsc/lzp/lzp.cpp
> g++ -c libbsc/platform/platform.cpp
> g++ -c libbsc/st/st.cpp
> g++ -c bsc.cpp

Linking command:
> g++ -o bsc bsc.o adler32.o bwt.o coder.o detectors.o libbsc.o libsais.o lzp.o platform.o preprocessing.o qlfc.o qlfc_model.o st.o
==============================================================================================================================
