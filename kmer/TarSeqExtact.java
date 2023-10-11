import java.io.*;
class TarSeqExtact{
	public static int MAX_CHA_NUM = 1 << 24;   //maximum length of a chromosome in HRCM 1<<28
	public static int VEC_SIZE = 1 << 20;      //length for other character arrays
	//public static int refCodeLen, refLowVecLen;
	//public static int[] refLowVecBegin;
    //public static int[] refLowVecLength;
	//public static char[] refCode;
	//static int[] speChaPos = new int[VEC_SIZE];
	public static String targetSequenceExtraction(String targetPath) {
        //String path = seqPaths.get(seqNum);
        File file = new File(targetPath);
        if (!file.exists()) {
            throw new MyException("fail to open the file and the name is " + file.getAbsolutePath());
        }

        int seqCodeLen = 0;
        int seqLowVecLen = 0;
        int seqNVecLen = 0;
        int speChaLen = 0;

        char[] seqCode = new char[MAX_CHA_NUM];
        int[] seqLowVecBegin = new int[VEC_SIZE];
        int[] seqLowVecLength = new int[VEC_SIZE];
        int[] seqNVecBegin = new int[VEC_SIZE];
        int[] seqNVecLength = new int[VEC_SIZE];
        int[] speChaPos = new int[VEC_SIZE];
        int[] speChaCh = new int[VEC_SIZE];

        String info;
        char ch;
        int lettersLen = 0, nLettersLen = 0, oLettersLen = 0;
        Boolean flag = true, nFlag = false;

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            //identifierVec[seqNum] = br.readLine();

            //br.mark(1000);
            //lineWidthVec[seqNum] = br.readLine().length();
            //br.reset();
			String identifier=br.readLine();
			System.out.println(identifier);
			
            while ((info = br.readLine()) != null) {
                for (int i = 0; i < info.length(); i++) {
                    ch = info.charAt(i);

                    if (Character.isLowerCase(ch)) {
                        if (flag) {
                            flag = false;
                            seqLowVecBegin[seqLowVecLen] = lettersLen;
                            lettersLen = 0;
                        }
                        ch = Character.toUpperCase(ch);
                    } else {
                        if (!flag) {
                            flag = true;
                            seqLowVecLength[seqLowVecLen++] = lettersLen;
                            lettersLen = 0;
                        }
                    }
                    lettersLen++;

                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        seqCode[seqCodeLen++] = ch;
                    } 
					/*else if (ch != 'N') {
                        speChaPos[speChaLen] = seqCodeLen+nLettersLen; //Problem
                        speChaCh[speChaLen++] = ch - 'A';
                    }*/
					else if (ch != 'N') {
                        speChaPos[speChaLen] = oLettersLen; //Problem
                        speChaCh[speChaLen++] = ch - 'A';
                    }


                    if (!nFlag) {
                        if (ch == 'N') {
                            seqNVecBegin[seqNVecLen] = nLettersLen;
                            nLettersLen = 0;
                            nFlag = true;
                        }
                    } else {
                        if (ch != 'N') {
                            seqNVecLength[seqNVecLen++] = nLettersLen;
                            nLettersLen = 0;
                            nFlag = false;
                        }
                    }
                    nLettersLen++;
					oLettersLen++;
                }
            }

            if (!flag) {
                seqLowVecLength[seqLowVecLen++] = lettersLen;
            }

            if (nFlag) {
                seqNVecLength[seqNVecLen++] = nLettersLen;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
		
		System.out.println(seqCodeLen);
		
		//for(int i=0;i<seqCodeLen;i++)
			//System.out.print(seqCode[i]);
		//System.out.println();
		
		String tarStr=new String(seqCode,0,seqCodeLen);
		
		for(int i=0;i<seqLowVecLen;i++)
			System.out.print(seqLowVecBegin[i]+" ");
		System.out.println();
		
		for(int i=0;i<seqLowVecLen;i++)
			System.out.print(seqLowVecLength[i]+" ");
		System.out.println();
		
		for(int i=0;i<seqNVecLen;i++)
			System.out.print(seqNVecBegin[i]+" ");
		System.out.println();
		
		for(int i=0;i<seqNVecLen;i++)
			System.out.print(seqNVecLength[i]+" ");
		System.out.println();
		
		//for(int i=0;i<speChaLen;i++)
			//System.out.print(speChaPos[i]+" "); //Use delta coding here
		//System.out.println();
		
		//Using delta coding here
		for(int i=speChaLen-1;i>0;i--)
			speChaPos[i]=speChaPos[i]-speChaPos[i-1]-1;
		for(int i=0;i<speChaLen;i++)
			System.out.print(speChaPos[i]+" "); 
		System.out.println();
		
		for(int i=0;i<speChaLen;i++)
			System.out.print(speChaCh[i]+" "); //Actually display speChaCh[i]-'A'
		System.out.println();
		
		return tarStr;
		
    }
	/*public static void main(String ...a){
		targetSequenceExtraction("TestTar.txt");
	}*/
}