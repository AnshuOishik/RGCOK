package huffman;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
// Uses static Huffman coding to compress an input file to an output file. Use HuffmanDecompress to decompress.
// Uses 257 symbols - 256 for byte values and 1 for EOF. The compressed file format contains the code length of each symbol under a canonical code, followed by the Huffman-coded data.
public final class HuffmanCompress {
	public  static FrequencyTable getFrequencies(File file) throws IOException {
		FrequencyTable freq = new FrequencyTable(new int[257]);
		InputStream input = new BufferedInputStream(new FileInputStream(file));
		try {
			while (true) {
				int b = input.read();
				if (b == -1)
					break;
				freq.increment(b);
			}
		} finally {
			input.close();
		}
		return freq;
	}
	public static void writeCode(BitOutputStream out, CanonicalCode canonCode) throws IOException {
		for (int i = 0; i < canonCode.getSymbolLimit(); i++) {
			int val = canonCode.getCodeLength(i);
			// For this file format, we only support codes up to 255 bits long
			if (val >= 256)
				throw new RuntimeException("The code for a symbol is too long");
			
			// Write value as 8 bits in big endian
			for (int j = 7; j >= 0; j--)
				out.write((val >>> j) & 1);
		}
	}
	public static void compress(CodeTree code, InputStream in, BitOutputStream out) throws IOException {
		HuffmanEncoder enc = new HuffmanEncoder(out);
		enc.codeTree = code;
		while (true) {
			int b = in.read();
			if (b == -1)
				break;
			enc.write(b);
		}
		enc.write(256);  // EOF
	}
}
