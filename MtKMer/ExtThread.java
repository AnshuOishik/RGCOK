package mtkmer;
import java.io.*;
class ExtThread extends Thread {
    int seq_num;
	int k_mer_len;
    ExtThread(int k_mer_len1, int seq_num1) {
        k_mer_len = k_mer_len1;
		seq_num = seq_num1; 
    }
    public void run() {
        GenSequence sequence = new GenSequence(); 
		try {
			OptKmerSize.tarSeqExt(seq_num, sequence);
			OptKmerSize.codeMatch(seq_num, sequence, k_mer_len);
		}
		catch(Exception e){
			System.out.println(e);
		}
    }
}
