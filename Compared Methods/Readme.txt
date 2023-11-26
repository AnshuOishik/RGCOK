***********************Execution commands for compared methods*************************
General information
1. Executed on the Linux platform
2. Please install the "7za" compressor
You may use the following instructions to install "7za" on the AWS Ubuntu platform.
> sudo apt-get update
> sudo apt-get install p7zip-full
# Please update the "7za" file's permissions to "0777" after that.
=========================================================================================================================================
										:::HiRGC:::
:::Java:::
# Compilation command: 
> javac HiRGCModifiedUbuntu.java
# Execution command: 
> java HiRGCModifiedUbuntu

# Kindly provide the target sequence name and reference in the source files.
# There is no decompression code for Java because authors use it mainly to test HiRGC compression performance. 

:::C++:::
# Compilation commands:
> make hirgc
> make de_hirgc

# Execution commands:
# Compress:
> ./hirgc -m file -r r.fa -t t.fa
where,
# -m is the mode that the user want to used,
# -r is the reference (a FASTA file),
# -t is the target (a FASTA file),
# the mode "file" is to comrpess chromosomes,
# the target sequence is "t.fa" and the reference sequence is "r.fa",
# the compressed sequence is "t.fa_ref_r.fa.7z"

# Decompress:
> ./de_hirgc -m file -r r.fa -t t.fa_ref_r.fa.7z
where,
# compressed chromosome file "t.fa_ref_r.fa.7z" decompressed to "dec_t.fa_ref_r.fa".

Note: There are two more techniques to compress and decompress in the original HiRGC code.
=========================================================================================================================================
										:::SCCG:::
:::Compression of single sequence:::

# Compression commands:
# Compilation command: 
> javac SCCGC.java
# Please create a folder name "Result"
# Execution command: 
> java -Xms8192m SCCGC R.fa T1.fa Result
where,
# the target sequence is "T1.fa" and the reference sequence is "R.fa"
# The "T1.fa.7z" compressed sequence is kept in the "Result" folder

# Decompression commands:
# Compilation command: 
> javac SCCGD.java
# Execution command: 
> java -Xms8192m SCCGD R.fa Result/T1.fa.7z Result
where,
# "T1.fa.7z," a compressed chromosomal file, decompressed to "T1.fa," which is kept in the "Result" folder

:::Compression of muitiple sequence:::

# Compression commands:
# Compilation command: 
> javac SCCGC.java
# Please create a folder name "Result"
# Execution command: 
> java -Xms8192m SCCGC Result

where,
# The target sequences are "T1.fa," "T2.fa," and "T3.fa," respectively, and the reference sequence is "R.fa." The "Result" folder contains the compressed sequences "T1.fa.7z", "T2.fa.7z", and "T3.fa.7z", in that order
Note: Kindly provide the target sequence name and reference in the source files

# Decompression commands:
# Compilation command: 
> javac SCCGD.java
# Execution command: 
> java -Xms8192m SCCGD Result

where,
# the compressed chromosomal files "T1.fa.7z", "T2.fa.7z", and "T3.fa.7z" in the "Result" folder decompress to "T1.fa", "T2.fa", and "T3.fa" in that order
=========================================================================================================================================
										:::HRCM:::
# Compilation Commands:
> g++ -c  compress.cpp
> g++ -c  decompress.cpp

# Linking command:
> g++ -o main compress.o decompress.o

# Execution commands:
> ./main compress -r R.fa -f chr22.txt 20
> ./main decompress -r R.fa -f chr22.txt 20

where,
# percent = 20 (optional), parameters "compress" and "decompress" are used for compression and decompression, respectively, and R.fa is the reference file. All files to be compressed are contained in chr22.txt, with the exception of the reference file.
# -r is the reference (a FASTA file)
# -t is the target (a text file containing target genomes)
# "chr22.7z" is a compressed sequence that is decompressed to match the relevant target genomes.
=========================================================================================================================================
# The ERGC source code link is https://engineering.uconn.edu/~rajasek/ERGC.zip, and the NRGC source code link is https://engineering.uconn.edu/~rajasek/NRGC.zip. Please take note of these links. These two pages are currently unavailable. Therefore, using the same machine configuration as in the prior research, we used their result.
=========================================================================================================================================
### Contacts 
Please send an email to <subhankar.roy07@gmail.com> if you experience any issues.
