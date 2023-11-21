package rgcok;
import java.io.*;
class ExtThread extends Thread {
    int seq_num;
    String temp_dir_path;
    ExtThread(int seq_num1, String temp_dir_path1) {
        seq_num = seq_num1;
        temp_dir_path = temp_dir_path1;
    }
    public void run() {
        GenSequence seq = new GenSequence();
        RGCOKCompress.tarSeqExt(seq_num, seq); 
        if (seq.getLowVecLen() > 0 && RGCOKCompress.ref_low_vec_len > 0)
            RGCOKCompress.lowerCaseMatch(seq);
		if (seq.getNVecLen() > 0 && RGCOKCompress.ref_n_vec_len > 0)
            RGCOKCompress.nMatch(seq);
        RGCOKCompress.firstLevelMatch(seq);
        if (seq_num <= RGCOKCompress.sec_seq_num) { 
            RGCOKCompress.match_result_vec.add(seq.getMatchResult());
            RGCOKCompress.matResHash(seq.getMatchResult());
        }
        BufferedWriter bw;
        try {
            bw = new BufferedWriter(new FileWriter(temp_dir_path + "compSeq-" + seq_num));
            RGCOKCompress.saveSplInfo(bw, seq); 
            if (seq_num != 1) { 
                if (seq.getMatchResult().size() < 3)
                    RGCOKCompress.saveFirstMatRes(bw, seq.getMatchResult());
                RGCOKCompress.secondLevelMatch(bw, seq.getMatchResult(), seq_num);
            } else {
                RGCOKCompress.saveFirstMatRes(bw, seq.getMatchResult());
            }
            seq.setMatchResult(null);
        } catch (IOException e) {
            System.out.println(e);
        }
    }
}
