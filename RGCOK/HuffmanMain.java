package huffman;
import java.io.*;
public class HuffmanMain {
	public static void main(String args1, String args2, String args3)throws IOException {											
		if(args1.equals("comp")){
			File inputFile1 = new File(args2);
			File outputFile1 = new File(args3);
			huffman.FrequencyTable freq1 = HuffmanCompress.getFrequencies(inputFile1);
			freq1.increment(256);  // EOF symbol gets a frequency of 1
			CodeTree code = freq1.buildCodeTree();
			CanonicalCode canonCode = new CanonicalCode(code, 257);
			code = canonCode.toCodeTree();  // Replace code tree with canonical one. For each symbol, the code value may change but the code length stays the same.
			// Read input file again, compress with Huffman coding, and write output file
			InputStream in4 = new BufferedInputStream(new FileInputStream(inputFile1));
			BitOutputStream out4 = new BitOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile1)));
			try {
				HuffmanCompress.writeCode(out4, canonCode);
				HuffmanCompress.compress(code, in4, out4);
				in4.close();
				delFile(inputFile1);
			} finally {
				out4.close();
				in4.close();
			}
		}  
		else if(args1.equals("decomp")){
			File inputFile2 = new File(args2);
			File outputFile2 = new File(args3);
			BitInputStream in5 = new BitInputStream(new BufferedInputStream(new FileInputStream(inputFile2)));
			OutputStream out5 = new BufferedOutputStream(new FileOutputStream(outputFile2));
			try {
				CanonicalCode canonCode1 = HuffmanDecompress.readCode(in5);
				CodeTree code1 = canonCode1.toCodeTree();
				HuffmanDecompress.decompress(code1, in5, out5);
			} finally {
				out5.close();
				in5.close();
			}
			//delFile(inputFile2);
		}			
	}//end of main
	public static boolean delFile(File file) {
        if (!file.exists()) {
            return false;
        }
        return file.delete();
    }
}//end of main class