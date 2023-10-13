package rgcok;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
class RGCOKCompress {
    static int max_seq_num = 2000, max_cha_num = 1 << 28, p = 25, k_mer_len = 11;  //max_cha_num = 1 << 28 and k_mer_len = 4(For Testing), Optimal k_mer_len value is selected using MtKmer algorithm                     
	static int kMer_bit_num = 2 * k_mer_len, hash_table_len = 1 << kMer_bit_num;     
	static final int max_arr_num_bit = 30, max_arr_num = 1<<max_arr_num_bit;   //max_arr_num_bit=20 (For Testing)
	static int vec_size = 1 << 20, min_rep_len = k_mer_len+1;
    static int seq_number, seq_bucket_len, sec_seq_num;                          
    static int ref_code_len, ref_low_vec_len, ref_n_vec_len;
    static char[] ref_code;
    static int[] ref_loc, ref_bucket, ref_low_vec_begin, ref_low_vec_length, ref_n_vec_begin, ref_n_vec_length, line_width_vec;
    static List<String> seq_paths;
    static String[] id_vec;
    static List<List<MatEntry>> match_result_vec;
    static List<int[]> seq_bucket_vec;
    static List<List<Integer>> seq_loc_vec; 
    
	static void secondLevelMatch(BufferedWriter bw, List<MatEntry> mat_res, int seq_num) {
        int hash_value, pre_seq_id = 1, id, pos, i, j, max_pos = 0, pre_pos = 0, delta_pos, len, max_len, delta_len, seq_id = 0, delta_seq_id;
        List<MatEntry> misMatchEntry = new ArrayList<>();
        try {
            misMatchEntry.add(mat_res.get(0));
			for (i = 1; i < mat_res.size() - 1; i++) {  
                hash_value = Math.abs(getTempHashValue(mat_res.get(i)) + getTempHashValue(mat_res.get(i + 1))) % seq_bucket_len; 		
				max_len = 0;
                for (j = 0; j < Math.min(seq_num - 1, sec_seq_num); j++) { 
                    id = seq_bucket_vec.get(j)[hash_value]; 		
                    if (id != -1) { 
                        for (pos = id; pos != -1; pos = seq_loc_vec.get(j).get(pos)) { 
                            len = getMatLen(match_result_vec.get(j), pos, mat_res, i);  
                            if (len > 1 && len > max_len) { 
                                seq_id = j + 1; 
                                max_pos = pos; 
                                max_len = len; 
                            }
                        }
                    }
                }
                if (max_len != 0) {
                    delta_seq_id = seq_id - pre_seq_id; 
                    delta_len = max_len - 2; 
                    delta_pos = max_pos - pre_pos;
                    pre_seq_id = seq_id; 
                    pre_pos = max_pos + max_len;
                    for (j = 0; j < misMatchEntry.size(); j++)
                        saveMatRes(bw, misMatchEntry.get(j));
                    bw.flush();
                    misMatchEntry = new ArrayList<>();
                    bw.write(delta_seq_id + " " + delta_pos + " " + delta_len + "\n");
                    i += max_len - 1; 
                } else
                    misMatchEntry.add(mat_res.get(i));
            }
            if (i == mat_res.size() - 1)
                misMatchEntry.add(mat_res.get(i));
            for (j = 0; j < misMatchEntry.size(); j++) 
               saveMatRes(bw, misMatchEntry.get(j));
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void saveFirstMatRes(BufferedWriter bw, List<MatEntry> mes) {
        int i;
		try {
            for (i = 0; i < mes.size(); i++)
                saveMatRes(bw, mes.get(i));
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void saveSplInfo(BufferedWriter bw, GenSequence seq) {
        int seq_low_vec_len = seq.getLowVecLen();
        int seq_n_vec_len = seq.getNVecLen();
        int seq_spe_cha_len = seq.getSpeChaLen();
        int flag = 0;
        try {
            if (seq_low_vec_len > 0 && ref_low_vec_len > 0) { 
                if ((2 * seq.getDiffLowVecLen()) < seq_low_vec_len) { 
                    flag = 1;
                    bw.write(flag + " "); 
                    eRLELowCase(bw, seq.getLowVecMatched(), seq_low_vec_len); 
                    //bw.write("\n"); Added
					savePosInfo(bw, seq.getDiffLowVecLen(), seq.getDiffLowVecBegin(), seq.getDiffLowVecLength());
				}
            }
            if (flag == 0) {
                bw.write(flag + " ");
                savePosInfo(bw, seq_low_vec_len, seq.getLowVecBegin(), seq.getLowVecLength());
            }
            seq.setLowVecLen(0);
            seq.setLowVecBegin(null);
            seq.setLowVecLength(null);
            seq.setLowVecMatched(null);
            seq.setDiffLowVecLen(0);
            seq.setDiffLowVecBegin(null);
			seq.setDiffLowVecLength(null);
			//bw.write("\n"); Added
			
			// RLE For N character matching:
			/*flag = 0;
			if (seq_n_vec_len > 0 && ref_low_vec_len > 0) { 
                if ((2 * seq.getDiffNVecLen()) < seq_n_vec_len) { 
                    flag = 1;
                    bw.write(flag + " ");
                    eRLENChar(bw, seq.getNVecMatched(), seq_n_vec_len);
                    bw.write("\n");
					savePosInfo(bw, seq.getDiffNVecLen(), seq.getDiffNVecBegin(), seq.getDiffNVecLength());
				}
            }
            if (flag == 0) {
                bw.write(flag + " ");
                savePosInfo(bw, seq_n_vec_len, seq.getNVecBegin(), seq.getNVecLength());
            }
            seq.setNVecLen(0);
            seq.setNVecBegin(null);
            seq.setNVecLength(null);
            seq.setNVecMatched(null);
            seq.setDiffNVecLen(0);
            seq.setDiffNVecBegin(null);
            seq.setDiffNVecLength(null);
			bw.write("\n");*/
            
			savePosInfo(bw, seq_n_vec_len, seq.getNVecBegin(), seq.getNVecLength()); //Original 
			seq.setNVecLen(0);
            seq.setNVecBegin(null);
            seq.setNVecLength(null);
			
			savePosInfo(bw, seq_spe_cha_len, seq.getSpeChaPos(), seq.getSpeChaCh()); 
            seq.setSpeChaLen(0);
            seq.setSpeChaPos(null);
            seq.setSpeChaCh(null);
			
            bw.write("\n");
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static int getTempHashValue(MatEntry me) { 
        int i, res = 0;
		String str = me.getMisMatStr();
        for (i = 0; i < str.length(); i++)
			res = res + 80021 * str.charAt(i);
		res = res + 40009 * me.getPos() + 50021 * me.getLen();
        res = res % seq_bucket_len;
        return res;
    }
	
    static void matResHash(List<MatEntry> mat_res) {
        int hash_value1, hash_value2, hash_value, i;
        List<Integer> seq_loc = new ArrayList<Integer>(vec_size);
        int[] seq_bucket = new int[seq_bucket_len];
        for (i = 0; i < seq_bucket_len; i++)
            seq_bucket[i] = -1;
        hash_value1 = getTempHashValue(mat_res.get(0)); 
        if (mat_res.size() < 2)
            hash_value2 = 0;
        else 
            hash_value2 = getTempHashValue(mat_res.get(1));
        hash_value = Math.abs(hash_value1 + hash_value2) % seq_bucket_len;
		seq_loc.add(seq_bucket[hash_value]);
        seq_bucket[hash_value] = 0;
        for (i = 1; i < mat_res.size() - 1; i++) {
            hash_value1 = hash_value2;
            hash_value2 = getTempHashValue(mat_res.get(i + 1));
            hash_value = Math.abs(hash_value1 + hash_value2) % seq_bucket_len;
            seq_loc.add(seq_bucket[hash_value]);
            seq_bucket[hash_value] = i;
        }
        seq_loc_vec.add(seq_loc);
		seq_bucket_vec.add(seq_bucket);
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
	
    static void firstLevelMatch(GenSequence seq) { 
        char[] seq_code = seq.getCode();
        int seq_code_len = seq.getCodeLen(), pre_pos = 0, step_len = seq_code_len - k_mer_len + 1; 
        int i, j, id, k, ref_id, tar_id, len, max_len, max_k;
		long tar_value; 
        MatEntry me = new MatEntry(); 
        List<MatEntry> mat_res = new ArrayList<>(vec_size);
        String mis_mat_str = String.valueOf(twoBitIntCoding(seq_code[0])); 
        for (i = 1; i < step_len; i++) { 
            tar_value = 0;
            for (j = k_mer_len - 1; j >= 0; j--) { 
                tar_value <<= 2; 
                tar_value += twoBitIntCoding(seq_code[i + j]);
            }
			id = ref_bucket[(int)(tar_value&(long)(max_arr_num - 1))];
			if(k_mer_len>15)
				id = ref_bucket[(int)(tar_value&(long)(max_arr_num - 1))];
			else
				id = ref_bucket[(int)tar_value];
            if (id > -1) {
                max_len = -1;
                max_k = -1;
                for (k = id; k != -1; k = ref_loc[k]) {
                    ref_id = k + k_mer_len;
                    tar_id = i + k_mer_len;
                    len = k_mer_len;
                    while (ref_id < ref_code_len && tar_id < seq_code_len && ref_code[ref_id++] == seq_code[tar_id++])
                        len++;   
                    min_rep_len = k_mer_len+1;
					if (len >= min_rep_len && len > max_len) {
                        max_len = len; 
                        max_k = k; 
                    }
                }
                if (max_len > -1) {
                    me.setMisMatStr(mis_mat_str); 
                    me.setPos(max_k - pre_pos); 
                    me.setLen(max_len - min_rep_len); 
                    mat_res.add(me);   
                    me = new MatEntry();
                    i += max_len;
                    pre_pos = max_k + max_len;
                    mis_mat_str = "";
                    if (i < seq_code_len) {
                        mis_mat_str += String.valueOf(twoBitIntCoding(seq_code[i]));
                    }
                    continue;
                }
            }
            mis_mat_str += String.valueOf(twoBitIntCoding(seq_code[i]));
        }
        if (i < seq_code_len) {
            for (; i < seq_code_len; i++) {
                mis_mat_str += String.valueOf(twoBitIntCoding(seq_code[i]));
            }
            me.setPos(0);
            me.setLen(-min_rep_len); 
            me.setMisMatStr(mis_mat_str); 
            mat_res.add(me);
        }
        seq.setMatchResult(mat_res);
        seq.setCodeLen(0);
        seq.setCode(null);
    }
	
	static void nMatch(GenSequence seq) {
        int seq_diff_n_vec_len = 0, start_pos = 1, i, j;
        int[] seq_n_vec_matched, seq_diff_n_vec_begin, seq_diff_n_vec_length, seq_n_vec_begin, seq_n_vec_length;
		seq_n_vec_matched = new int[vec_size];
		seq_diff_n_vec_begin = new int[vec_size]; 
        seq_diff_n_vec_length = new int[vec_size];
        seq_n_vec_begin = seq.getNVecBegin();
        seq_n_vec_length = seq.getNVecLength();
        for (i = 0; i < seq.getNVecLen(); i++) { 
            for (j = start_pos; j < ref_n_vec_len; j++) {
                if ((seq_n_vec_begin[i] == ref_n_vec_begin[j]) && (seq_n_vec_length[i] == ref_n_vec_length[j])) {
                    seq_n_vec_matched[i] = j;
                    start_pos = j + 1;
                    break;
                }
            }
            if (seq_n_vec_matched[i] == 0) {
                for (j = start_pos - 1; j > 0; j--) {
                    if ((seq_n_vec_begin[i] == ref_n_vec_begin[j]) && (seq_n_vec_length[i] == ref_n_vec_length[j])) {
                        seq_n_vec_matched[i] = j;
                        start_pos = j + 1;
                        break;
                    }
                }
            }
            if (seq_n_vec_matched[i] == 0) {
                seq_diff_n_vec_begin[seq_diff_n_vec_len] = seq.getNVecBegin()[i];
                seq_diff_n_vec_length[seq_diff_n_vec_len++] = seq.getNVecLength()[i];
            }
        }
        seq.setNVecMatched(seq_n_vec_matched);
        seq.setDiffNVecBegin(seq_diff_n_vec_begin);
        seq.setDiffNVecLength(seq_diff_n_vec_length);
        seq.setDiffNVecLen(seq_diff_n_vec_len);
    }
	
	static void lowerCaseMatch(GenSequence seq) {
        int seq_diff_low_vec_len = 0,start_pos = 1, i, j;
		int[] seq_low_vec_matched, seq_diff_low_vec_begin, seq_diff_low_vec_length, seq_low_vec_begin, seq_low_vec_length;
		seq_low_vec_matched = new int[vec_size];
        seq_diff_low_vec_begin = new int[vec_size]; 
        seq_diff_low_vec_length = new int[vec_size];
        seq_low_vec_begin = seq.getLowVecBegin();
        seq_low_vec_length = seq.getLowVecLength();
        for (i = 0; i < seq.getLowVecLen(); i++) { 
            for (j = start_pos; j < ref_low_vec_len; j++) { 
                if ((seq_low_vec_begin[i] == ref_low_vec_begin[j]) && (seq_low_vec_length[i] == ref_low_vec_length[j])) {
                    seq_low_vec_matched[i] = j; 
                    start_pos = j + 1; 
                    break;
                }
            }
            if (seq_low_vec_matched[i] == 0) {
                for (j = start_pos - 1; j > 0; j--) {
                    if ((seq_low_vec_begin[i] == ref_low_vec_begin[j]) && (seq_low_vec_length[i] == ref_low_vec_length[j])) {
                        seq_low_vec_matched[i] = j;
                        start_pos = j + 1;
                        break;
                    }
				}
            }
            if (seq_low_vec_matched[i] == 0) {
                seq_diff_low_vec_begin[seq_diff_low_vec_len] = seq.getLowVecBegin()[i];
                seq_diff_low_vec_length[seq_diff_low_vec_len++] = seq.getLowVecLength()[i];
            }
        }
        seq.setLowVecMatched(seq_low_vec_matched);
        seq.setDiffLowVecBegin(seq_diff_low_vec_begin);
        seq.setDiffLowVecLength(seq_diff_low_vec_length);
        seq.setDiffLowVecLen(seq_diff_low_vec_len);
    }
	
    static void tarSeqExt(int seq_num, GenSequence seq) { 
        String path = seq_paths.get(seq_num);
        File f = new File(path);
        int seq_code_len = 0,seq_low_vec_len = 0,seq_n_vec_len = 0, spe_cha_len = 0; 
        char[] seq_code = new char[max_cha_num];
		int[] seq_low_vec_begin, seq_low_vec_length, seq_n_vec_begin, seq_n_vec_length, spe_cha_pos, spe_cha_ch;
		seq_low_vec_begin = new int[vec_size];
		seq_low_vec_length = new int[vec_size];
		seq_n_vec_begin = new int[vec_size];
		seq_n_vec_length = new int[vec_size];
		spe_cha_pos =new int[vec_size]; 
		spe_cha_ch = new int[vec_size];
        String str;
        char ch;
        int lett_len = 0, n_lett_len = 0;
        Boolean flag = true, n_flag = false;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(f));
            id_vec[seq_num] = br.readLine(); 
            br.mark(1000);
            line_width_vec[seq_num] = br.readLine().length(); 
            br.reset();
            while ((str = br.readLine()) != null) {
                for (int i = 0; i < str.length(); i++) {
                    ch = str.charAt(i);
                    if (Character.isLowerCase(ch)) {
                        if (flag) {
                            flag = false;
                            seq_low_vec_begin[seq_low_vec_len] = lett_len;
                            lett_len = 0;
                        }
                        ch = Character.toUpperCase(ch);
                    } else {
                        if (!flag) {
                            flag = true;
                            seq_low_vec_length[seq_low_vec_len++] = lett_len;
                            lett_len = 0;
                        }
                    }
                    lett_len++;
                    if (ch == 'T' || ch == 'G' || ch == 'C' || ch == 'A') {
                        seq_code[seq_code_len++] = ch;
                    } else if (ch != 'N') {
                        spe_cha_pos[spe_cha_len] = seq_code_len;
                        spe_cha_ch[spe_cha_len++] = ch - 'A';
                    }
                    if (!n_flag) { 
                        if (ch == 'N') {
                            seq_n_vec_begin[seq_n_vec_len] = n_lett_len;
                            n_lett_len = 0;
                            n_flag = true;
                        }
                    } else {
                        if (ch != 'N') {
                            seq_n_vec_length[seq_n_vec_len++] = n_lett_len;
                            n_lett_len = 0;
                            n_flag = false;
                        }
                    }
                    n_lett_len++; 
                }
            }
            if (!flag)
                seq_low_vec_length[seq_low_vec_len++] = lett_len;
            if (n_flag)
                seq_n_vec_length[seq_n_vec_len++] = n_lett_len;
			br.close();
        } catch (Exception e) { 
            System.out.println(e);
        } 
        seq.setCode(seq_code);
        seq.setLowVecBegin(seq_low_vec_begin);
        seq.setLowVecLength(seq_low_vec_length);
        seq.setNVecBegin(seq_n_vec_begin);
        seq.setNVecLength(seq_n_vec_length);
        seq.setSpeChaPos(spe_cha_pos);
        seq.setSpeChaCh(spe_cha_ch);
        seq.setCodeLen(seq_code_len);
        seq.setLowVecLen(seq_low_vec_len);
        seq.setNVecLen(seq_n_vec_len);
        seq.setSpeChaLen(spe_cha_len);
    } 
	
	//Improved RLE
    static void eRLELowCase(BufferedWriter bw, int[] vec, int len) { 
        List<Integer> code = new ArrayList<>(vec_size);
        int i, c, code_len;
		if (len > 0) { 
            code.add(vec[0]); 
            c = 1;
            for (i = 1; i < len; i++) { 
                if (vec[i] - vec[i - 1] == 1)
                    c++; 
                else {
                    code.add(c); 
                    code.add(vec[i]); 
                    c = 1;
                }
            }
            code.add(c); 
        }
        code_len = code.size(); 
        try {
            bw.write(code_len + " "); 
            for (i = 0; i < code_len; i++)
                bw.write(code.get(i) + " "); 
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	//Modified RLE for N character, Included if result is improved.
	static void eRLENChar(BufferedWriter bw, int[] vec, int len) { 
        List<Integer> code = new ArrayList<>(vec_size);
        int i, c, code_len;
		if (len > 0) { 
            code.add(vec[0]); 
            c = 1;
            for (i = 1; i < len; i++) { 
                if (vec[i] - vec[i - 1] == 1)
                    c++; 
                else {
                    code.add(c); 
                    code.add(vec[i]); 
                    c = 1;
                }
            }
            code.add(c); 
        }
        code_len = code.size(); 
        try {
            bw.write(code_len + " "); 
            for (i = 0; i < code_len; i++)
                bw.write(code.get(i) + " "); 
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
    static void savePosInfo(BufferedWriter bw, int vec_len1, int[] vecBegin, int[] vec_len2) { //bw, 1, 12, 10; bw, 2, 19, 3
        int i;
		try {
            bw.write(vec_len1 + " "); 
            for (i = 0; i < vec_len1; i++)
                bw.write(vecBegin[i] + " " + vec_len2[i] + " "); //12 10; 19 3 12 3
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
    static void saveMatRes(BufferedWriter bw, MatEntry me) {
        try {
            bw.write(me.getMisMatStr() + "\n" + me.getPos() + " " + me.getLen() + "\n");
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
    static int getMatLen(List<MatEntry> ref_mat_res, int ref_idx, List<MatEntry> tar_mat_res, int tar_idx) {
        int len = 0;
        while (ref_idx < ref_mat_res.size() && tar_idx < tar_mat_res.size() && compMatEntry(ref_mat_res.get(ref_idx++), tar_mat_res.get(tar_idx++)))
            len++; 
        return len; 
    }
	
    static Boolean compMatEntry(MatEntry ref, MatEntry tar) {
        if (ref.getPos() == tar.getPos() && ref.getLen() == tar.getLen() && ref.getMisMatStr().equals(tar.getMisMatStr()))
            return true;
        else
            return false;
    }
	
    //BSC Compression
	public static void bscCompression() { 
		try {
            String tarCommand = "tar -cf " + "TarC.tar "+ "chr.id " + "chr.rgcok";
			Process p1 = Runtime.getRuntime().exec(tarCommand);
            p1.waitFor();
            String bscCommand = "./bsc e " + "TarC.tar " + "BscC.bsc -e2";
            Process p2 = Runtime.getRuntime().exec(bscCommand);
            p2.waitFor();

            deleteFile(new File("TarC.tar"));
            deleteFile(new File("chr.id"));
            deleteFile(new File("chr.rgcok"));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
	
	static boolean deleteFile(File f) { 
        if (!f.exists()) 
            return false;
        if (f.isDirectory()) {
            File[] fs = f.listFiles();
            for (File f1 : fs) 
                deleteFile(f1);
        }
        return f.delete();
    }
	
	//Identifiers stored as a raw string
    static void saveIdInfo(BufferedWriter bw, String[] id_vec) {
        int i;
		try {
            for (i = 1; i < id_vec.length; i++)
                bw.write(id_vec[i] + "\n");
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void eRLELineWidth(BufferedWriter bw, int[] vec, int len) { //desc file, number of sequence - 2, number of sequence - 2
        List<Integer> code = new ArrayList<>(vec_size);
        int i, c, code_len;
		if (len > 0) { 
            code.add(vec[1]); 
            c = 1;
            for (i = 2; i < len; i++) { //If there are more than one target seq
                if (vec[i] == vec[i - 1])
                    c++; 
                else {
                    code.add(c); 
                    code.add(vec[i]);
                    c = 1;
                }
            }
            code.add(c);
        }
        code_len = code.size(); 
        try {
            bw.write(code_len + " "); 
            for (i = 0; i < code_len; i++)
                bw.write(code.get(i) + " ");
            bw.write("\n");
            bw.flush();
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void multiThreading(String temp_dir_path, int pool_size) { 
        ExecutorService exe_service_first = Executors.newSingleThreadExecutor(); 
        int i;
		for (i = 1; i <= sec_seq_num; i++)
            exe_service_first.execute(new ExtThread(i, temp_dir_path));
        exe_service_first.shutdown(); 
        while (!exe_service_first.isTerminated()) {}
        ExecutorService exe_service_later = Executors.newFixedThreadPool(pool_size);
        for (i = sec_seq_num + 1; i < seq_number; i++)
            exe_service_later.execute(new ExtThread(i, temp_dir_path));
        exe_service_later.shutdown();
        while (!exe_service_later.isTerminated()) {}
    }
	
	static void kMerHash() { 
        ref_loc = new int[max_cha_num];
        ref_bucket = new int[hash_table_len]; 
		int i;
        for (i = 0; i < hash_table_len; i++)
            ref_bucket[i] = -1;
        long value1 = 0; 
        int value2, step_len = ref_code_len - k_mer_len + 1, shift_bit_num = 2 * (k_mer_len - 1);; 
        for (i = k_mer_len - 1; i >= 0; i--) { 
            value1 <<= 2; 
            value1 += twoBitIntCoding(ref_code[i]);
        }
		if(k_mer_len>15)
			value2 = (int)(value1&(long)(max_arr_num - 1));
		else
			value2 = (int) value1;
        ref_loc[0] = ref_bucket[value2]; 
        ref_bucket[value2] = 0; 
        for (i = 1; i < step_len; i++) { 
            value1 >>= 2; 
            value1 += (twoBitIntCoding(ref_code[i + k_mer_len - 1])) << shift_bit_num;
			if(k_mer_len>15)
				value2 = (int)(value1&(long)(max_arr_num - 1));
			else
				value2 = (int) value1;
            ref_loc[i] = ref_bucket[value2]; 
            ref_bucket[value2] = i; 
        }
    }
	
	static void refSeqExt(String ref_path) {
        File f = new File(ref_path);
        ref_code = new char[max_cha_num];
        ref_code_len = 0;
        ref_low_vec_begin = ref_low_vec_length = ref_n_vec_begin = ref_n_vec_length = new int[vec_size]; 
        ref_low_vec_len = ref_n_vec_len = 1; 
        int lett_len = 0, n_lett_len = 0, i; 
        String str;
        char ch;
        boolean flag = true, n_flag = false;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(f));
            br.readLine();
            while ((str = br.readLine()) != null) {
                for (i = 0; i < str.length(); i++) {
                    ch = str.charAt(i);
                    if (Character.isLowerCase(ch)) {
                        if (flag) {
                            flag = false;
                            ref_low_vec_begin[ref_low_vec_len] = lett_len;
                            lett_len = 0;
                        }
                        ch = Character.toUpperCase(ch);
                    } else {
                        if (!flag) {
                            flag = true;
                            ref_low_vec_length[ref_low_vec_len++] = lett_len;
                            lett_len = 0;
                        }
                    }
                    if (ch == 'T' || ch == 'G' || ch == 'C' || ch == 'A') {
                        ref_code[ref_code_len++] = ch;
                    }
                    lett_len++;		
					if (!n_flag) { 
                        if (ch == 'N') {
                            ref_n_vec_begin[ref_n_vec_len] = n_lett_len; 
                            n_lett_len = 0;
                            n_flag = true;
                        }
                    } else { 
                        if (ch != 'N') {
                            ref_n_vec_length[ref_n_vec_len++] = n_lett_len; 
                            n_lett_len = 0;
                            n_flag = false;
                        }
                    }
                    n_lett_len++; 
                }
            }
            if (!flag)
                ref_low_vec_length[ref_low_vec_len++] = lett_len;
			if (n_flag) 
                ref_n_vec_length[ref_n_vec_len++] = n_lett_len;
			br.close();
        } catch (Exception e) {
             System.out.println(e);
        } 
    }
	
	static String getTmpDirPath(String path) {
		String res = "";
		int i;
        String[] str = path.split("/");
        for (i = 0; i < str.length - 1; i++)
            res = res+str[i]+"/";
        File tmpDir = new File(res + "tempFile/");
        if (!tmpDir.exists())
            tmpDir.mkdir();
		return res+ "tempFile/";
    }
	
	static void seqCompress(String f_path, int pool_size) { 
        String temp_dir_path = getTmpDirPath(f_path);  
        File rgcok_file = new File(f_path.replaceAll(".fa", "") + ".rgcok"); 
        if (rgcok_file.exists()) rgcok_file.delete(); 
        File id_file = new File(f_path.replaceAll(".fa", "") + ".id"); 
        if (id_file.exists()) id_file.delete();
        refSeqExt(seq_paths.get(0)); 
        kMerHash();
        BufferedReader br = null;
        BufferedWriter bw1 = null,bw2 = null;
        try {
            multiThreading(temp_dir_path, pool_size); 
            File f;
            String str;
			int i;
            bw1 = new BufferedWriter(new FileWriter(rgcok_file, true));
            for (i = 1; i < seq_number; i++) { 
                f = new File(temp_dir_path + "compSeq-" + i);
                br = new BufferedReader(new FileReader(f));
                while ((str = br.readLine()) != null)
                    bw1.write(str + "\n"); 
                bw1.write("\n"); 
            }
            bw2 = new BufferedWriter(new FileWriter(id_file, true));
            eRLELineWidth(bw2, line_width_vec, seq_number);
            saveIdInfo(bw2, id_vec);
            deleteFile(new File(temp_dir_path));
			bw1.close();
			bw2.close();
			br.close();
        } catch (Exception e) { 
            System.out.println(e);
        }
    }
	
	static int getPrimeNumber(int n) { 
        int next_prime = n + 1,i;
        boolean p = false;
        while (!p) {
            p = true;
            for (i = 2; i < Math.sqrt(n) + 1; ++i) {
                if (next_prime % i == 0) {
                    p = false;
                    break;
                }
            }
            if (!p) {
                next_prime++;
            }
        }
        return next_prime;
    }  
	
	static int readInputFile(File file_path) { 
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
	
	static void initSettings(File file_path) { 
        seq_number = readInputFile(file_path);
        sec_seq_num = (int) Math.ceil((double) (p * seq_number) / 100); 
        seq_bucket_len = getPrimeNumber(vec_size);
        id_vec = new String[seq_number]; 
        line_width_vec = new int[seq_number]; 
        seq_bucket_vec = new ArrayList<>(seq_number); 
        seq_loc_vec = new ArrayList<>(seq_number); 
        match_result_vec = new ArrayList<>(sec_seq_num); 
    }
}
