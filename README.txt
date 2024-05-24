				       ****************************** MtK-mer & RGCOK ****************************
							MtK-mer: Multi-Threaded Optimal k-mer Length
					RGCOK: Reference Genome Compression Algorithm using Optimal k-mer Length
							https://github.com/AnshuOishik/RGCOK
								Copyright (C) 2023 
==================================================================================================================================================================
If you use RGCOK then please cite the paper as: 
Subhankar Roy, Anirban Mukhopadhyay, A randomized optimal k-mer indexing approach for efficient parallel genome sequence compression,
Gene, Volume 907, 2024, ISSN 0378-1119, https://doi.org/10.1016/j.gene.2024.148235.
Link: https://www.sciencedirect.com/science/article/abs/pii/S0378111924001161
==================================================================================================================================================================
Introduction
To utilize the code, please use the Notepad++ editor.
Java has been utilized by us in the implementation.
Please use Linux or Windows as your operating system.
Please confirm that the physical memory on your computer is larger than 10 GB.
=============================================================================================================================
MtK-mer: Multi-Threaded Optimal k-mer Length
The MtK-mer algorithm used the randomization method to determine the ideal k-mer length. RGCOK uses the ideal k-mer length discovered using MtK-mer to compress the particular sequence.

# Compilation Command:
> javac -d . *.java

# Execution Command:
> java -Xms10240m mtkmer.Main chr.fa 9 30 8 

Notice:
# -Xms10240m is the initial allocation of memory (MB)
# The list of target file directories and the reference file path (the first line) are both found in chr.fa
# The k-mer length's lower and upper bounds are 9 and 30
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
# Please execute the procedure listed at the end of this file to create an executable file for the BSC compressor
# In the last phase, an alternative is to utilize a 7-zip compressor; the procedure is described at the end of this file
# The final compressed file created by the bsc compressor is called "BscC.bsc"
# "ZipC.7z" is the name of the compressed file that the 7-zip compressor produced at the end
# The number of threads is eight (4, by default, is the optional value)
# -Xms10240m is the initial allocation of memory (MB)
# Please place the executable "bsc" and "7za" in the main class file's directory
# For "bsc" and "7za" modes, kindly set "chmod 0777"
=============================================================================================================================
#Commands for "bsc" executable file generation from available code at https://github.com/IlyaGrebnov/libbsc
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
# Please change the platform.cpp file. In lines 51 and 66, change 'MEM_LARGE_PAGES' in Linux (Ubuntu) to 'MEM_4MB_PAGES' in Windows 10.
> g++ -c libbsc/st/st.cpp
> g++ -c bsc.cpp

#Linking command:
> g++ -o bsc bsc.o adler32.o bwt.o coder.o detectors.o libbsc.o libsais.o lzp.o platform.o preprocessing.o qlfc.o qlfc_model.o st.o

Notice:
# The created executable file name is bsc.
=============================================================================================================================
# Installing the 7-zip compressor for Windows can be done at https://www.7-zip.org. according to your operating system.
# Please set 7z path to environment variable.
# Please use the following command to install 7-zip on Linux.
> sudo apt-get update (If required)
> sudo apt-get install p7zip-full
=============================================================================================================================
### Contacts 
Please send an email to <subhankar.roy07@gmail.com> if you experience any issues.
