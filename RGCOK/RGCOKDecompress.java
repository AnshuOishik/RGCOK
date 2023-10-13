package rgcok;
import java.io.*;
import java.util.*;
class RGCOKDecompress {
    static int max_seq_num = 2000, max_cha_num = 1 << 28, p = 25, k_mer_len = 11; //max_cha_num = 1 << 20 and k_mer_len = 4 (For Testing), Optimal k_mer_len value is selected using MtKmer algorithm           
    static int kMer_bit_num = 2 * k_mer_len, hash_table_len = 1 << kMer_bit_num; 
	static int max_arr_num_bit = 30, max_arr_num = 1<<max_arr_num_bit; //max_arr_num_bit = 20 (For Testing)
    static int vec_size = 1 << 20, min_rep_len = k_mer_len+1; 
    static int seq_number, sec_seq_num;
    static int ref_code_len, ref_low_vec_len;
    static char[] ref_code;
	static int[] ref_loc, ref_bucket, ref_low_vec_begin, ref_low_vec_length;
    static char[] int_code = new char[]{'A', 'C', 'G', 'T'};
    static List<String> seq_paths, id_vec;
    static List<Integer> line_width_vec;
    static List<List<MatEntry>> match_res_vec;
    static BufferedReader br1,br2;
    static BufferedWriter bw;
    
	static void saveSeqFile(GenSequence seq, int seq_num, String f_path) {
        String[] str = seq_paths.get(seq_num).split("/");
        File f = new File(getOutputFolderPath(f_path) + str[str.length - 1]);
        if (f.exists())
            f.delete();
        int code_len = seq.getCodeLen();
        int low_vec_len = seq.getLowVecLen();
        int n_vec_len = seq.getNVecLen();
        int spe_cha_len = seq.getSpeChaLen();
        int line_width = seq.getLineWidth();
        int i, j, k=0, temp, o_len = 0, n_len = 0, l = 0, tr = 0;;
        char[] temp_seq = new char[max_cha_num];
        char[] c_seq = seq.getCode();
        for (i = 0; i < spe_cha_len; i++) {
            while (o_len < seq.getSpeChaPos()[i] && o_len < code_len)
                temp_seq[n_len++] = c_seq[o_len++];
            temp_seq[n_len++] = (char) (seq.getSpeChaCh()[i] + 'A');
        }
        while (o_len < code_len) 
            temp_seq[n_len++] = c_seq[o_len++];
        char[] code_n = new char[max_cha_num];
        for (i = 0; i < n_vec_len; i++) {
            for (j = 0; j < seq.getNVecBegin()[i]; j++) 
                code_n[l++] = temp_seq[tr++];
            for (j = 0; j < seq.getNVecLength()[i]; j++) 
                code_n[l++] = 'N';
        }
        while (tr < n_len) 
            code_n[l++] = temp_seq[tr++];
        for (i = 0; i < low_vec_len; i++) {
            k += seq.getLowVecBegin()[i];
            temp = seq.getLowVecLength()[i];
            for (j = 0; j < temp; j++) {
                code_n[k] = Character.toLowerCase(code_n[k]);
                k++;
            }
        }
        try {
            bw = new BufferedWriter(new FileWriter(f, true));
            bw.write(seq.getId() + "\n");
            for (i = 0; i < l; i += line_width) {
                for (j = i; j < i + line_width; j++)
                    if (j < l) 
                        bw.write(code_n[j]);
                bw.write("\n");
                bw.flush();
            }
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void readCompTarSeq(GenSequence seq) {
        List<MatEntry> mat_res = seq.getMatchResult();
        int pos, pre_pos = 0, cur_pos, l;
        String mis_str;
        char[] code_arr = new char[max_cha_num];
        int i, j, k, n, code_len = 0;
        for (i = 0; i < mat_res.size(); i++) {
            pos = mat_res.get(i).getPos();
            cur_pos = pos + pre_pos;
            l = mat_res.get(i).getLen() + min_rep_len;
            pre_pos = cur_pos + l;
            mis_str = mat_res.get(i).getMisMatStr();
            for (j = 0; j < mis_str.length(); j++)
                code_arr[code_len++] = int_code[Integer.valueOf(mis_str.charAt(j)) - 48];
            for (k = cur_pos, n = 0; n < l; k++, n++)
                code_arr[code_len++] = ref_code[k];
        }
        seq.setCodeLen(code_len);
        seq.setCode(code_arr);
    }
	
	static void getMatRes(int id, int pos, int l, List<MatEntry> mat_res) {
        for (int i = 0; i < l; i++)
            mat_res.add(match_res_vec.get(id).get(pos++));
    }
	
	static void readMatRes(GenSequence seq) {
        List<MatEntry> mat_res = new ArrayList<>(vec_size);
        MatEntry mat_ent = new MatEntry();
        String s;
        String[] str;
        int id, pos, l, pre_seq_id = 0, pre_pos = 0;
        try {
            while (!(s = br1.readLine()).equals("")) {
                str = s.split("\\s+");
                if (str.length == 3) {
                    id = Integer.valueOf(str[0]);
                    pos = Integer.valueOf(str[1]);
                    l = Integer.valueOf(str[2]);
                    id += pre_seq_id;
                    pre_seq_id = id;
                    pos += pre_pos;
                    l += 2;
                    pre_pos = pos + l;
                    getMatRes(id, pos, l, mat_res);
                } else if (str.length == 2) {
                    pos = Integer.valueOf(str[0]);
                    l = Integer.valueOf(str[1]);
                    mat_ent.setPos(pos);
                    mat_ent.setLen(l);
                    mat_res.add(mat_ent);
                    mat_ent = new MatEntry();
                } else if (str.length == 1) 
                    mat_ent.setMisMatStr(s);
            }
            seq.setMatchResult(mat_res);
        } catch (IOException e) {
            System.out.println(e);
		}
    }
	
	static void readSplInfo(GenSequence seq) {
        try {
            String[] str = br1.readLine().split("\\s+");
            int flg = Integer.valueOf(str[0]);
            if (flg == 0) {
                int seq_low_vec_len = Integer.valueOf(str[1]);
                int i, k = 0;
                int[] seq_low_vec_begin = new int[seq_low_vec_len];
                int[] seq_low_vec_length = new int[seq_low_vec_len];
                for (i = 2; i < 2 + seq_low_vec_len * 2; i++) {
                    if (i % 2 == 0) 
                        seq_low_vec_begin[k] = Integer.valueOf(str[i]);
                    else 
                        seq_low_vec_length[k++] = Integer.valueOf(str[i]);
                }
                seq.setLowVecLen(seq_low_vec_len);
                seq.setLowVecBegin(seq_low_vec_begin);
                seq.setLowVecLength(seq_low_vec_length);
                int seq_n_vec_len = Integer.valueOf(str[2 + seq_low_vec_len * 2]);
                k = 0;
                int[] seq_n_vec_begin = new int[seq_n_vec_len];
                int[] seq_n_vec_length = new int[seq_n_vec_len];
                for (i = 3 + seq_low_vec_len * 2; i < 3 + seq_low_vec_len * 2 + seq_n_vec_len * 2; i++) {
                    if (i % 2 == 1) 
                        seq_n_vec_begin[k] = Integer.valueOf(str[i]);
                    else
                        seq_n_vec_length[k++] = Integer.valueOf(str[i]);
                }
                seq.setNVecLen(seq_n_vec_len);
                seq.setNVecBegin(seq_n_vec_begin);
                seq.setNVecLength(seq_n_vec_length);
                int seq_spe_cha_len = Integer.valueOf(str[3 + seq_low_vec_len * 2 + seq_n_vec_len * 2]);
                k = 0;
                if (seq_spe_cha_len > 0) {
                    int[] seq_spe_cha_pos = new int[seq_spe_cha_len];
                    int[] seq_spe_cha_ch = new int[seq_spe_cha_len];
                    for (i = 4 + seq_low_vec_len * 2 + seq_n_vec_len * 2; i < str.length; i++) {
                        if (i % 2 == 0)
                            seq_spe_cha_pos[k] = Integer.valueOf(str[i]);
						else 
                            seq_spe_cha_ch[k++] = Integer.valueOf(str[i]);
                    }
                    seq.setSpeChaPos(seq_spe_cha_pos);
                    seq.setSpeChaCh(seq_spe_cha_ch);
                }
                seq.setSpeChaLen(seq_spe_cha_len);
            } else {
                int i, j, l = Integer.valueOf(str[1]);
                List<Integer> v = new ArrayList<>(vec_size);
                for (i = 2; i < 2 + l; i++) 
                    v.add(Integer.valueOf(str[i]));
                int seq_low_vec_len = 0;
                for (i = 1; i < v.size(); i += 2) 
                    seq_low_vec_len += v.get(i);
                int[] seq_low_vec_matched = new int[seq_low_vec_len];
                int k = 0;
                for (i = 0; i < v.size(); i += 2)
                    for (j = 0; j < v.get(i + 1); j++)
                        seq_low_vec_matched[k++] = v.get(i) + j;
                int seq_diff_low_vec_len = Integer.valueOf(str[2 + l]);
                k = 0;
                int[] seq_diff_low_vec_begin = new int[seq_diff_low_vec_len];
                int[] seq_diff_low_vec_length = new int[seq_diff_low_vec_len];
                for (i = 3 + l; i < 3 + l + seq_diff_low_vec_len * 2; i++) {
                    if (i % 2 == 1) 
                        seq_diff_low_vec_begin[k] = Integer.valueOf(str[i]);
                    else
                        seq_diff_low_vec_length[k++] = Integer.valueOf(str[i]);
                }
                k = 0;
                int[] seq_low_vec_begin = new int[vec_size];
                int[] seq_low_vec_length = new int[vec_size];
                for (i = 0, j = 0; i < seq_low_vec_len; i++) {
                    k = seq_low_vec_matched[i];
                    if (k == 0) {
                        seq_low_vec_begin[i] = seq_diff_low_vec_begin[j];
                        seq_low_vec_length[i] = seq_diff_low_vec_length[j++];
                    } else {
                        seq_low_vec_begin[i] = ref_low_vec_begin[k];
                        seq_low_vec_length[i] = ref_low_vec_length[k];
                    }
                }
                seq.setLowVecLen(seq_low_vec_len);
                seq.setLowVecBegin(seq_low_vec_begin);
                seq.setLowVecLength(seq_low_vec_length);
                int seq_n_vec_len = Integer.valueOf(str[3 + l + seq_diff_low_vec_len * 2]);
                k = 0;
                int[] seq_n_vec_begin = new int[seq_n_vec_len];
                int[] seq_n_vec_length = new int[seq_n_vec_len];
                for (i = 4 + l + seq_diff_low_vec_len * 2; i < 4 + l + seq_diff_low_vec_len * 2 + seq_n_vec_len * 2; i++) {
                    if (i % 2 == 0)
                        seq_n_vec_begin[k] = Integer.valueOf(str[i]);
                    else
                        seq_n_vec_length[k++] = Integer.valueOf(str[i]);
                }
                seq.setNVecLen(seq_n_vec_len);
                seq.setNVecBegin(seq_n_vec_begin);
                seq.setNVecLength(seq_n_vec_length);
                int seq_spe_cha_len = Integer.valueOf(str[4 + l + seq_diff_low_vec_len * 2 + seq_n_vec_len * 2]);
                k = 0;
                if (seq_spe_cha_len > 0) {
                    int[] seq_spe_cha_pos = new int[seq_spe_cha_len];
                    int[] seq_spe_cha_ch = new int[seq_spe_cha_len];
                    for (i = 5 + l + seq_diff_low_vec_len * 2 + seq_n_vec_len * 2; i < str.length; i++) {
                        if (i % 2 == 1)
                            seq_spe_cha_pos[k] = Integer.valueOf(str[i]);
						else
                            seq_spe_cha_ch[k++] = Integer.valueOf(str[i]);
                    }
                    seq.setSpeChaPos(seq_spe_cha_pos);
                    seq.setSpeChaCh(seq_spe_cha_ch);
                }
                seq.setSpeChaLen(seq_spe_cha_len);
            }
        } catch (IOException e) {
            System.out.println(e);
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
	
	static String getOutputFolderPath(String f_path) {
        String out_folder_path = "";
        String s []= f_path.split("/");
        for (int i = 0; i < s.length - 1; i++) 
            out_folder_path += s[i]+"/";
        out_folder_path += "Output"+"/";
        return out_folder_path;
    }
	
    static void outputFolder(String f_path) {
        File out_folder = new File(getOutputFolderPath(f_path));
        if (out_folder.exists()) 
            deleteFile(out_folder);
        out_folder.mkdir();
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
	
	//construct hash table
    static void kMerHash() {
        int i;
		ref_loc = new int[max_cha_num];
        ref_bucket = new int[hash_table_len];
        for (i = 0; i < hash_table_len; i++)
            ref_bucket[i] = -1;
		long value1 = 0;
		int value2;
        int step_len = ref_code_len - k_mer_len + 1;
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
        int shift_bit_num = (k_mer_len * 2 - 2);
        for (i = 1; i < step_len; i++) {
            value1 >>= 2;
            value1 += (twoBitIntCoding(ref_code[i + k_mer_len - 1])) << shift_bit_num;
			value2 = (int)(value1&(long)(max_arr_num - 1));
            ref_loc[i] = ref_bucket[value2];
            ref_bucket[value2] = i;
        }
    }
	
	static void refSeqExt(String ref_path) {
        File f = new File(ref_path);
        ref_code = new char[max_cha_num];
        ref_code_len = 0;
        ref_low_vec_begin = new int[vec_size];
        ref_low_vec_length = new int[vec_size];
        ref_low_vec_len = 1;
        int lett_len = 0, i;
        String str;
        char c;
        Boolean flg = true;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(f));
            br.readLine();
            while ((str = br.readLine()) != null) {
                for (i = 0; i < str.length(); i++) {
                    c = str.charAt(i);
                    if (Character.isLowerCase(c)) {
                        if (flg) {
                            flg = false;
                            ref_low_vec_begin[ref_low_vec_len] = lett_len;
                            lett_len = 0;
                        }
                        c = Character.toUpperCase(c);
                    } else {
                        if (!flg) {
                            flg = true;
                            ref_low_vec_length[ref_low_vec_len++] = lett_len;
                            lett_len = 0;
                        }
                    }
                    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                        ref_code[ref_code_len++] = c;
                    lett_len++;
                }
            }
            if (!flg) 
                ref_low_vec_length[ref_low_vec_len++] = lett_len;
			br.close();
        } catch (IOException e) {
            System.out.println(e);
        } 
    }
	
	static void getIdentifierData(BufferedReader br) {
        String str;
        try {
            while ((str = br.readLine()) != null) 
                id_vec.add(str);
        } 
		catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void dRLELineWidth(BufferedReader br) {
        int l = 0, i, j;
        String s;
        try {
            s = br.readLine();
            String[] str = s.split("\\s+");
            int[] w = new int[Integer.valueOf(str[0])];
            for (i = 0; i < w.length; i++) 
                w[i] = Integer.valueOf(str[i + 1]);
            for (i = 1; i < Integer.valueOf(str[0]); i += 2) 
                l += w[i];
            if (l > 0) {
                line_width_vec = new ArrayList<>(l);
                for (i = 0; i < w.length; i += 2) 
                    for (j = 0; j < w[i + 1]; j++)
                        line_width_vec.add(w[i]);
            }
        } catch (IOException e) {
            System.out.println(e);
        }
    }
	
	static void seqDecompress(String f_path) {
        String rgcok_path = f_path.replaceAll("fa", "") + "rgcok";    //stored the base info of genomes
        String aux_path = f_path.replaceAll("fa", "") + "id";         //stored the line width and id
        File rgcok_file = new File(rgcok_path);
        File aux_file = new File(aux_path);
        try {
            br1 = new BufferedReader(new FileReader(rgcok_path));
            br2 = new BufferedReader(new FileReader(aux_path));
            dRLELineWidth(br2);
            getIdentifierData(br2);
            refSeqExt(seq_paths.get(0));
            kMerHash();
            GenSequence seq = new GenSequence();
            outputFolder(f_path);
            for (int i = 1; i < seq_number; i++) {
                seq.setId(id_vec.get(i - 1));
                seq.setLineWidth(line_width_vec.get(i - 1));
                readSplInfo(seq);
                readMatRes(seq);
                if (i <= sec_seq_num && i != seq_number - 1)
                    match_res_vec.add(seq.getMatchResult());
                readCompTarSeq(seq);
                saveSeqFile(seq, i, f_path);
                seq = new GenSequence();
            }
			br2.close();
			br1.close();
			bw.close();
        } catch (IOException e) { 
            System.out.println(e);
        } 
		deleteFile(rgcok_file);
		deleteFile(aux_file);
	}
		
	//BSC Decompression
	public static void bscDecompression() {
        try {
            String bscCommand = "./bsc d " + "BscC.bsc " + "TarC.tar";
            Process p1 = Runtime.getRuntime().exec(bscCommand);
            p1.waitFor();
            String tarCommand = "tar -xf " + "TarC.tar";
            Process p2 = Runtime.getRuntime().exec(tarCommand);
            p2.waitFor();
            String delCommand = "rm -f " + "TarC.tar "+ "BscC.bsc";
            Process p3 = Runtime.getRuntime().exec(delCommand);
            p3.waitFor();
        } catch (Exception e) {
            System.out.println(e);
        }
    }
	
	static int readInputFile(File file_path) {
        seq_paths = new ArrayList<>(max_seq_num);
        String f_path;
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file_path));
            while ((f_path = br.readLine()) != null) 
                seq_paths.add(f_path);
			br.close();
        } catch (IOException e) {
            System.out.println(e);
        } 
        return seq_paths.size();
    }
	
    static void initSettings(File file_path) {
        seq_number = readInputFile(file_path);
        sec_seq_num = (int) Math.ceil((double) (p * seq_number) / 100);
        id_vec = new ArrayList<>(seq_number);
        match_res_vec = new ArrayList<>(sec_seq_num);
    }
}