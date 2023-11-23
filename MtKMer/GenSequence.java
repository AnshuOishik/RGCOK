package mtkmer;
class GenSequence {
    int code_len;
    char[] code;
    GenSequence() {
        int max_char_num = 1 << 28;
        int vec_size = 1 << 20;
        code_len = 0;   
        code = new char[max_char_num];  
    }
    int getCodeLen() {
        return code_len;
    }
    void setCodeLen(int code_len1) {
        code_len = code_len1;
    }
    char[] getCode() {
        return code;
    }
    void setCode(char[] code1) {
        code = code1;
    }
}
