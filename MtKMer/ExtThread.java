package mtkmer;
import java.io.*;
class ExtThread extends Thread {
    int seq_num;
	int k_mer_len;
	BufferedWriter bw;
    ExtThread(int k_mer_len1, int seq_num1, BufferedWriter bw1) {
        k_mer_len = k_mer_len1;
		seq_num = seq_num1; 
		bw = bw1;
    }
    public void run() {
        GenSequence sequence = new GenSequence(); 
		try {
			OptKmerSize.kMerHash(k_mer_len , 1<<(2*k_mer_len));
			OptKmerSize.tarSeqExt(seq_num, sequence);
			OptKmerSize.codeMatch(seq_num, sequence, k_mer_len, bw);
		}
		catch(Exception e){
			System.out.println(e);
		}
    }
}
