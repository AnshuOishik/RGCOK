//Finding optimal kmer length (random approach)
import java.util.*; 
import java.io.*;
import java.util.concurrent.*;
class OptKmerSize{
	static List<String> seq_path;
	static String[] id_vec;
	static int seq_num,max_seq_num = 2000,max_cha_num = 1 << 28,ref_code_len;
	static char[] ref_code;
	static int kmer_len [] = new int[13]; //k-mer length = 9 to 21
	static int kmer_bit_num [] = new int[13]; 
    static int hash_table_len [] = new int[13]; 
	static int max_arr_num_bit = 30,max_arr_num = 1<<max_arr_num_bit; 
	static final int max_arr_num_bit_shift = max_arr_num_bit>>1;
	static int[] ref_loc;
    static int[] ref_bucket;
	static File kmer_file;
	static void begining(File path_file)throws NewException {
        if (!path_file.exists()) {
            throw new NewException("Fail to open the file input.fa!!!");
        }
		seq_num = readInputFile(path_file);
        id_vec = new String[seq_num];
    }
	static int readInputFile(File path_file)throws NewException { 
        if (!path_file.exists()) {
            throw new NewException("Fail to open the " + path_file.getAbsolutePath()+"file");
        }
        seq_path = new ArrayList<>(max_seq_num);
        String file_path;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(path_file));
            while ((file_path = br.readLine()) != null) {
                seq_path.add(file_path);
            }
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (Exception e) {
                System.out.println(e);
            }
        }
        return seq_path.size();
    }
	public static int toIintegerCoding(char c) {
        switch (c) {
            case 'A':return 0;
            case 'C':return 1;
            case 'G':return 2;
            case 'T':return 3;
            default:return -1;
        }
    }
	public static void kMerHashConstruction(int kmer_len, int hash_table_len) {  
		int min = 0, max = ref_code_len;
		int number_of_kmer=(int) Math.ceil((ref_code_len-kmer_len+1)/2); //Taking p=50% of reference code length
		int i, step_len, index=0, value;
		ref_loc = new int[max_cha_num];
		ref_bucket = new int[max_arr_num];
        for (i = 0; i < max_arr_num; i++) { 
            ref_bucket[i] = -1;
        }
		while(true){
			step_len = (int)(Math.random()*(max-min+1)+min)-kmer_len; //Maximum value of stepLen= 41-4=37
			if(step_len < (kmer_len - 1)) 
				continue;
			long value1 = 0;
			int value2;
			for (i = step_len+kmer_len - 1; i >= step_len; i--) { //3 to 0 for AGCT
				value1 <<= 2; 
				value1 += toIintegerCoding(ref_code[i]);//0+3 = 3, 12+1=13, 52+2=54, 216+0=216
			}
			if(kmer_len>15)
				value2 = (int)(value1&(long)(max_arr_num - 1));
			else
				value2 = (int) value1;
			ref_loc[index] = ref_bucket[value2]; 
			ref_bucket[value2] = index++; 
			if(index>number_of_kmer)
				break;
		}
    }
	public static void referenceSequenceExtraction(String referencePath, int kmer_len,int hash_table_len)throws NewException {
		File file = new File(referencePath);
        if (!file.exists()) {
            throw new NewException("fail to open the file and the name is " + file.getAbsolutePath());
        }
        ref_code = new char[max_cha_num];
        ref_code_len = 0;
        int lettersLen = 0;
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
					
                    if (Character.isLowerCase(ch)) {  //a
                        ch = Character.toUpperCase(ch);
                    } 
                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        ref_code[ref_code_len++] = ch;
                    }
                }
            }
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (Exception e) {
                System.out.println(e);
            }
        }
		kMerHashConstruction(kmer_len , hash_table_len);
	}
	private static void codeMatch(int pool_size, int kmer_len , BufferedWriter bw) { //Path of tmp, 4, Used for second level compression
        ExecutorService executorServiceLater = Executors.newFixedThreadPool(pool_size); //Thread pool is thread safe, Creates a thread pool that reuses a fixed number of threads operating off a shared unbounded queue. At any point, at most nThreads threads will be active processing tasks. If additional tasks are submitted when all threads are active, they will wait in the queue until a thread is available. If any thread terminates due to a failure during execution prior to shutdown, a new one will take its place if needed to execute subsequent tasks. The threads in the pool will exist until it is explicitly shutdown.
        for (int i = 1; i < seq_num; i++) {  
            executorServiceLater.execute(new ExtThread(i, kmer_len , bw));
        }
        executorServiceLater.shutdown();
        while (!executorServiceLater.isTerminated()) {}
    }
	public static void targetSequenceExtraction(int seq_num, GenSequence sequence) throws NewException{
		String path = seq_path.get(seq_num);
        File file = new File(path);
        if (!file.exists()) {
            throw new NewException("Fail to open the " + file.getAbsolutePath()+"file!!!");
        }
        int seq_code_len = 0;
        char[] seq_code = new char[max_cha_num];
        String info;
        char ch;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
			String id=br.readLine();
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
        } catch (Exception e) {
            System.out.println(e);
        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (Exception e) {
                System.out.println(e);
            }
        }
		sequence.setCode(seq_code);
		sequence.setCodeLen(seq_code_len);
    }
	public static void codeFirstMatch(int seq_num, GenSequence sequence, int kmer_len , BufferedWriter bw)throws IOException { //Target sequence
        char[] seq_code = sequence.getCode();
        int seq_code_len = sequence.getCodeLen(); 
		String tar_str=new String(seq_code,0,seq_code_len);
        int i , j , id;
		long count=0; 
        int step_len = seq_code_len - kmer_len + 1; 
		long tar_value;
        for (i = 0; i < step_len; i++) { 
            tar_value = 0;
            for (j = kmer_len - 1; j >= 0; j--) { 
                tar_value <<= 2; 
                tar_value += toIintegerCoding(seq_code[i + j]); 
            }
			if(kmer_len>15){
				id = ref_bucket[(int)(tar_value&(long)(max_arr_num - 1))];
				int ref_idx = id + max_arr_num_bit_shift;
                int tar_idx = i + max_arr_num_bit_shift;
				for(int k=0 ; k<kmer_len-15 ; k++){
					if(ref_code[ref_idx] == seq_code[tar_idx]) {
						ref_idx++;
						tar_idx++;
					}
					else{
						id=-1;
						break;
					}
				}
			}else
				id = ref_bucket[(int)tar_value];
            if (id != -1) 
              count++;
        }
		System.out.println("SeqNum, kmer_len, Count = "+seq_num+"  "+kmer_len+"  "+(kmer_len*count));
        sequence.setCodeLen(0);
        sequence.setCode(null);
    }
	public static void findOptKmerSize(String path, int pool_size)throws IOException,NewException{ 
		BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
		try{
			if (seq_num > 1) {
				kmer_file = new File(path.replaceAll(".txt", "") + ".kmer"); //chr22.kmer
				kmer_file.createNewFile();
			} else {
				throw new NewException("There is no target sequence!!!");
			}
		}
		catch(Exception e){
			System.out.println(e);
		}
		int K=9;
		for(int i=0;i<kmer_len.length;i++){
			kmer_len [i] = K++;
			kmer_bit_num [i] = 2*kmer_len [i];
			hash_table_len [i] = 1<<kmer_bit_num [i];
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(kmer_file));
		for(int i=0;i<kmer_len.length;i++){
			referenceSequenceExtraction(seq_path.get(0),kmer_len [i],hash_table_len [i]);
			codeMatch(pool_size,kmer_len [i],bw);
		}
	} 
}