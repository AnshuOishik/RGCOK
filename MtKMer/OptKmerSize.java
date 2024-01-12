//Finding optimal kmer length (Randomization approach)
package mtkmer;
import java.util.*; 
import java.io.*;
import java.util.concurrent.*;
class OptKmerSize{
	static List<String> seq_paths;
	static String[] id_vec;
	static int seq_num, max_seq_num = 120, max_cha_num = 1 << 28, ref_code_len; 
	static char[] ref_code;
	static int max_arr_num_bit = 30, max_arr_num = 1 << max_arr_num_bit; 
	static final int max_arr_num_bit_shift = max_arr_num_bit >> 1;
	static int[] ref_loc, ref_bucket;
	static File kmer_file;
	
	static void beginingSettings(File file_path){
		seq_num = readInputFile(file_path);
        id_vec = new String[seq_num];
    }
	
	static int readInputFile(File file_path){ 
        seq_paths = new ArrayList<>(max_seq_num);
        String f_path;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file_path));
            while ((f_path = br.readLine()) != null) {
                seq_paths.add(f_path);
            }
			br.close();
        } catch (Exception e) {
            System.out.println(e);
        } 
        return seq_paths.size();
    }
	
	static int twoBitIntCoding(char c) {
		int r;
        switch (c) {
            case 'A':	r = 0; break;
            case 'C':	r = 1; break;
            case 'G':	r = 2; break;
            case 'T':	r = 3; break;
            default:	r = -1;
        }
		return r;
    }
	
	static void kMerHash(int kmer_len) {  
		int min = 0, max = ref_code_len;
		int number_of_kmer = (int) Math.ceil((ref_code_len-kmer_len+1)/2); 
		int i, step_len, index=0, value;
		ref_loc = new int[number_of_kmer];
		ref_bucket = new int[max_arr_num];
        for (i = 0; i < max_arr_num; i++) { 
            ref_bucket[i] = -1;
        }
		while(true){
			step_len = (int)(Math.random()*(max-min+1)+min)-kmer_len;
			if(step_len < (kmer_len - 1)) 
				continue;
			long value1 = 0;
			int value2;
			for (i = step_len + kmer_len - 1; i >= step_len; i--) { 
				value1 <<= 2; 
				value1 += twoBitIntCoding(ref_code[i]);
			}
			
			if(kmer_len > 15)
				value2 = (int)(value1&(long)(max_arr_num - 1));
			else
				value2 = (int) value1;
			
			ref_loc[index] = ref_bucket[value2]; //s_arr
			ref_bucket[value2] = index++; //f_arr
			if(index >= number_of_kmer)
				break;
		}
    }
	
	static void refSeqExt(String ref_path){
		File file = new File(ref_path);
        ref_code = new char[max_cha_num];
        ref_code_len = 0;
        String info;
        char ch;
        Boolean flag = true;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            br.readLine();
            while ((info = br.readLine()) != null) {
                for (int i = 0; i < info.length(); i++) {
                    ch = info.charAt(i);
                    if (Character.isLowerCase(ch)) {  
                        ch = Character.toUpperCase(ch);
                    } 
                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        ref_code[ref_code_len++] = ch;
                    }
                }
            }
			br.close();
        } catch (Exception e) {
            System.out.println(e);
        } 
	}
	
	static void multiThreading(int pool_size, int k_mer_len) { 
        ExecutorService executorServiceLater = Executors.newFixedThreadPool(pool_size); 
        int i;
		for (i = 1; i < seq_num; i++) {  
			executorServiceLater.execute(new ExtThread(k_mer_len, i));
		}
        executorServiceLater.shutdown();
        while (!executorServiceLater.isTerminated()) {}
    }
	
	static void tarSeqExt(int seq_num, GenSequence seq) {
		String path = seq_paths.get(seq_num);
        File file = new File(path);
        int seq_code_len = 0;
        char[] seq_code = new char[max_cha_num];
        String info;
        char ch;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
			String id = br.readLine();
            while ((info = br.readLine()) != null) {
                for (int i = 0; i < info.length(); i++) {
                    ch = info.charAt(i);
                    if (Character.isLowerCase(ch)) {
                        ch = Character.toUpperCase(ch);
                    } 
                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        seq_code[seq_code_len++] = ch;
                    } 		
                }
            }
			br.close();
        } catch (Exception e) {
            System.out.println(e);
        } 
		seq.setCode(seq_code);
		seq.setCodeLen(seq_code_len);
    }
	
	static void codeMatch(int seq_num, GenSequence seq, int kmer_len)throws IOException { 
        char[] seq_code = seq.getCode();
        int seq_code_len = seq.getCodeLen(); 
		String tar_str=new String(seq_code,0,seq_code_len);
        int i , j , id;
		long count=0; 
        int step_len = seq_code_len - kmer_len + 1; 
		long tar_value;
        for (i = 0; i < step_len; i++) { 
            tar_value = 0;
            for (j = kmer_len - 1; j >= 0; j--) { 
                tar_value <<= 2; 
                tar_value += twoBitIntCoding(seq_code[i + j]); 
            }
			if(kmer_len > 15){
				id = ref_bucket[(int)(tar_value&(long)(max_arr_num - 1))];
				int ref_idx = id + max_arr_num_bit_shift;
                int tar_idx = i + max_arr_num_bit_shift;
				for(int k=0 ; k < kmer_len-15 ; k++){
					if(ref_code[ref_idx] == seq_code[tar_idx]) {
						ref_idx++;
						tar_idx++;
					}
					else{
						id = -1;
						break;
					}
				}
			}
			else
				id = ref_bucket[(int)tar_value];
            if (id != -1) 
              count++;
        }
		System.out.println("SeqNum, kmer_len, Count = "+seq_num+"  "+kmer_len+"  "+(kmer_len*count));
        seq.setCodeLen(0);
        seq.setCode(null);
    }
	
	static void findOptKmerSize(String path, int k_mer_len1, int k_mer_len2, int pool_size)throws IOException{ 
		try{
			if (seq_num > 1) {
				kmer_file = new File(path.replaceAll(".txt", "") + ".kmer");
				kmer_file.createNewFile();
			}
		}
		catch(Exception e){
			System.out.println(e);
		}
		
		refSeqExt(seq_paths.get(0));
		for(int k_mer_len = k_mer_len1; k_mer_len <= k_mer_len2; k_mer_len++){
			kMerHash(k_mer_len);
			multiThreading(pool_size, k_mer_len);
			ref_loc = ref_bucket = null;
			System.gc();
		}
	} 
}