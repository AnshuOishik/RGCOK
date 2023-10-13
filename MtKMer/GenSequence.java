class GenSequence {
    int code_len;
    char[] code;
    GenSequence() {
        int max_char_num = 1 << 28;
        int vec_size = 1 << 20;
        this.code_len = 0;   
        this.code = new char[max_char_num];  
    }
    int getCodeLen() {
        return code_len;
    }
    void setCodeLen(int code_len) {
        this.code_len = code_len;
    }
    char[] getCode() {
        return code;
    }
    void setCode(char[] code) {
        this.code = code;
    }
}
