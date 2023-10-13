import java.io.*;
class ExtThread extends Thread {
    int seq_num;
	int k_mer_len;
	BufferedWriter bw;
    ExtThread(int seq_num, int k_mer_len, BufferedWriter bw) {
        this.seq_num = seq_num; 
		this.k_mer_len = k_mer_len;
		this.bw = bw;
    }
    public void run() {
        GenSequence sequence = new GenSequence(); 
		try {
			OptKmerSize.targetSequenceExtraction(seq_num, sequence);
			OptKmerSize.codeFirstMatch(seq_num, sequence, k_mer_len, bw);
		}
		catch(Exception e){
			e.printStackTrace();
		}
    }
}
