package mtkmer;
import java.io.*;
class ExtThread extends Thread {
    int seq_num;
	int k_mer_len;
	BufferedWriter bw;
    ExtThread(int seq_num1, int k_mer_len1, BufferedWriter bw1) {
        seq_num = seq_num1; 
		k_mer_len = k_mer_len1;
		bw = bw1;
    }
    public void run() {
        GenSequence sequence = new GenSequence(); 
		try {
			OptKmerSize.tarSeqExt(seq_num, sequence);
			OptKmerSize.codeMatch(seq_num, sequence, k_mer_len, bw);
		}
		catch(Exception e){
			System.out.println(e);
		}
    }
}
