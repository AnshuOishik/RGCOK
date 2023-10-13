//Include multithreading concept, using executor thread
import java.io.*;
class Main{
	public static void main(String[] a) throws IOException,NewException{ 
		String[] args = new String[2];
		String nt="4";
        args[0] = a[0];
        if (a.length == 1) {
            args[1] = nt;
        } else {
            args[1] = a[1];
        }
        String path = args[0];
        File path_file = new File(path);
        int pool_size = Integer.valueOf(args[1]);
        
		OptKmerSize.begining(path_file); 
		
		Long start_kmer_time = System.currentTimeMillis();
		
		OptKmerSize.findOptKmerSize(path_file.getAbsolutePath(), pool_size);
		
		Long end_kmer_time = System.currentTimeMillis();
        
		System.out.println("Time = " + 1.0*(end_kmer_time - start_kmer_time) / 1000 + "s.");	
	}
}