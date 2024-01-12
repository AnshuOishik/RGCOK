//Main class of MtK-mer algorithm
package mtkmer;
import java.io.*;
class Main{
	public static void main(String[] a) throws IOException{ 
		String[] args = new String[4];
		String nt="4";
        args[0] = a[0];
		args[1] = a[1];
		args[2] = a[2];
        if (a.length == 3) {
            args[3] = nt;
        } else {
            args[3] = a[3];
        }
        String path = args[0];
        File path_file = new File(path);
		int k_mer_len1 = Integer.valueOf(args[1]);
		int k_mer_len2 = Integer.valueOf(args[2]);
        int pool_size = Integer.valueOf(args[3]);
       
		OptKmerSize.beginingSettings(path_file); 
		
		Long start_kmer_time = System.currentTimeMillis();
		
		OptKmerSize.findOptKmerSize(path_file.getAbsolutePath(), k_mer_len1, k_mer_len2, pool_size);
		
		Long end_kmer_time = System.currentTimeMillis();
        
		System.out.println("Time = " + 1.0*(end_kmer_time - start_kmer_time) / 1000 + "s.");	
	}
}