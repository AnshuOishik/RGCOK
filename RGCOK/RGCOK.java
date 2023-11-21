package rgcok;
import huffman.*;
import java.io.*;
class RGCOK {
    public static void main(String[] a) throws IOException{
        String[] args = new String[3];
        args[0] = a[0];	args[1] = a[1];
		String nt = "4";
        if (a.length == 2)	args[2] = nt;
        else	args[2] = a[2];
        int thread_pool_size = Integer.parseInt(args[2]);
        if (args[1].equals("comp")) {
			double beforeUsedMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			Long start_time = System.currentTimeMillis();
			System.out.println("Compression Started...");
			
            RGCOKCompress.initSettings(new File(args[0]));
			RGCOKCompress.seqCompress(args[0], thread_pool_size);
			
			//Compression using bsc encoder
            RGCOKCompress.bscCompression(); //Output of BSC compressor is "BscC.bsc"
			
			//Compression using Huffman encoder
			//HuffmanMain.main(args[1],"chr.id","HCchr.id");	
			//HuffmanMain.main(args[1],"chr.rgcok","HCchr.rgcok");
            
			Long end_time = System.currentTimeMillis();
            System.out.println("Total compression time = " + (1.0*end_time - start_time)/ 1000 + " S");
			double afterUsedMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			double actualMemUsed = afterUsedMem-beforeUsedMem;
			System.out.printf("Actual memory used = %.3f GB\n", (actualMemUsed/1024/1024/1024));
			System.out.println("Compression Completed...");
        } else if (args[1].equals("decomp")) {
			double beforeUsedMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			Long start_time = System.currentTimeMillis();
			System.out.println("Decompression Started...");

			RGCOKDecompress.initSettings(new File(args[0]));
			
			//Decompression using Huffman decoder
			//HuffmanMain.main(args[1],"HCchr.id","chr.id");
			//HuffmanMain.main(args[1],"HCchr.rgcok","chr.rgcok");
			
			//Decompression using bsc decoder
			RGCOKDecompress.bscDecompression();
			
			RGCOKDecompress.seqDecompress(args[0]);
			
            Long end_time = System.currentTimeMillis();
            System.out.println("Total decompression time = " + (1.0*end_time - start_time)/ 1000 + " S");
			double afterUsedMem = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
			double actualMemUsed = afterUsedMem-beforeUsedMem;
			System.out.printf("Actual memory used = %.3f GB\n", (actualMemUsed/1024/1024/1024));
			System.out.println("Decompression Completed...");
        } 
    }
}
