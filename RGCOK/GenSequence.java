package rgcok;
import java.util.*;
class GenSequence {
    String id;
    int line_width, code_len,low_vec_len,n_vec_len,spe_cha_len,diff_low_vec_len,diff_n_vec_len; 
    char[] code;
    int[] low_vec_begin,low_vec_length,n_vec_begin,n_vec_length,spe_cha_pos,spe_cha_ch,low_vec_matched,n_vec_matched,diff_low_vec_begin,diff_n_vec_begin,diff_low_vec_length,diff_n_vec_length;
    List<MatEntry> match_result;
    GenSequence() {
        int max_cha_num = 1 << 28, vec_size = 1 << 20; 
        id = "";   
        line_width = code_len = low_vec_len = diff_low_vec_len = diff_n_vec_len = n_vec_len = spe_cha_len = 0; 
        code = new char[max_cha_num]; 
        low_vec_matched = new int[vec_size];
		n_vec_matched = new int[vec_size];
        low_vec_begin = new int[vec_size]; 
        low_vec_length = new int[vec_size];
        diff_low_vec_begin = new int[vec_size];
		diff_n_vec_begin = new int[vec_size];
        diff_low_vec_length = new int[vec_size];
		diff_n_vec_length = new int[vec_size];
		n_vec_begin = new int[vec_size];
        n_vec_length = new int[vec_size];
        spe_cha_pos = new int[vec_size];
        spe_cha_ch = new int[vec_size];
        match_result = new ArrayList<>(vec_size);
    }
    String getId() {
        return id;
    }
    void setId(String id1) {
        id = id1;
    }
	int getCodeLen() {
        return code_len;
    }
    void setCodeLen(int code_len1) {
        code_len = code_len1;
    }
	int getLowVecLen() {
        return low_vec_len;
    }
    void setLowVecLen(int low_vec_len1) {
        low_vec_len = low_vec_len1;
    }
    int getDiffLowVecLen() {
        return diff_low_vec_len;
    }
    void setDiffLowVecLen(int diff_low_vec_len1) {
        diff_low_vec_len = diff_low_vec_len1;
    }
	int getDiffNVecLen() { 
        return diff_n_vec_len;
    }
    void setDiffNVecLen(int diff_n_vec_len1) { 
        diff_n_vec_len = diff_n_vec_len1;
    }
    int getNVecLen() {
        return n_vec_len;
    }
    void setNVecLen(int n_vec_len1) {
        n_vec_len = n_vec_len1;
    }
    int getSpeChaLen() {
        return spe_cha_len;
    }
    void setSpeChaLen(int spe_cha_len1) {
        spe_cha_len = spe_cha_len1;
    }
    int getLineWidth() {
        return line_width;
    }
    void setLineWidth(int line_width1) {
        line_width = line_width1;
    }
    char[] getCode() {
        return code;
    }
    void setCode(char[] code1) {
        code = code1;
    }
    int[] getLowVecMatched() {
        return low_vec_matched;
    }
    void setLowVecMatched(int[] low_vec_matched1) {
        low_vec_matched = low_vec_matched1;
    }
	int[] getNVecMatched() {
        return n_vec_matched;
    }
    void setNVecMatched(int[] n_vec_matched1) {
        n_vec_matched = n_vec_matched1;
    }
    int[] getLowVecBegin() {
        return low_vec_begin;
    }
    void setLowVecBegin(int[] low_vec_begin1) {
        low_vec_begin = low_vec_begin1;
    }
    int[] getLowVecLength() {
        return low_vec_length;
    }
    void setLowVecLength(int[] low_vec_length1) {
        low_vec_length = low_vec_length1;
    }
    int[] getNVecBegin() {
        return n_vec_begin;
    }
    void setNVecBegin(int[] n_vec_begin1) {
        n_vec_begin = n_vec_begin1;
    }
    int[] getNVecLength() {
        return n_vec_length;
    }
    void setNVecLength(int[] n_vec_length1) {
        n_vec_length = n_vec_length1;
    }
    int[] getSpeChaPos() {
        return spe_cha_pos;
    }
    void setSpeChaPos(int[] spe_cha_pos1) {
        spe_cha_pos = spe_cha_pos1;
    }
    int[] getSpeChaCh() {
        return spe_cha_ch;
    }
    void setSpeChaCh(int[] spe_cha_ch1) {
        spe_cha_ch = spe_cha_ch1;
    }
    int[] getDiffLowVecBegin() {
        return diff_low_vec_begin;
    }
    void setDiffLowVecBegin(int[] diff_low_vec_begin1) {
        diff_low_vec_begin = diff_low_vec_begin1;
    }
	int[] getDiffNVecBegin() {
        return diff_n_vec_begin;
    }
    void setDiffNVecBegin(int[] diff_n_vec_begin1) {
        diff_n_vec_begin = diff_n_vec_begin1;
    }
    int[] getDiffLowVecLength() {
        return diff_low_vec_length;
    }
    void setDiffLowVecLength(int[] diff_low_vec_length1) {
        diff_low_vec_length = diff_low_vec_length1;
    }
	int[] getDiffNVecLength() {
        return diff_n_vec_length;
    }
    void setDiffNVecLength(int[] diff_n_vec_length1) {
        diff_n_vec_length = diff_n_vec_length1;
    }
    List<MatEntry> getMatchResult() {
        return match_result;
    }
    void setMatchResult(List<MatEntry> match_result1) {
        match_result = match_result1;
    }
}
