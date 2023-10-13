package huffman;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
// Decompresses an input file that was compressed with HuffmanCompress, to an output file.
public final class HuffmanDecompress {
	public static CanonicalCode readCode(BitInputStream in) throws IOException {
		int[] codeLengths = new int[257];
		for (int i = 0; i < codeLengths.length; i++) {
			// For this file format, we read 8 bits in big endian
			int val = 0;
			for (int j = 0; j < 8; j++) 
				val = val << 1 | in.readNoEof();
			codeLengths[i] = val;
		}
		return new CanonicalCode(codeLengths);
	}
	public static void decompress(CodeTree code, BitInputStream in, OutputStream out) throws IOException {
		HuffmanDecoder dec = new HuffmanDecoder(in);
		dec.codeTree = code;
		while (true) {
			int symbol = dec.read();
			if (symbol == 256)  // EOF symbol
				break;
			out.write(symbol);
		}
	}
}
