import java.io.*;
class RefSeqExtact{
	public static int MAX_CHA_NUM = 1 << 24;   //maximum length of a chromosome in HRCM 1<<28
	public static int VEC_SIZE = 1 << 20;      //length for other character arrays
	public static int refCodeLen, refLowVecLen;
	public static int[] refLowVecBegin;
    public static int[] refLowVecLength;
	public static char[] refCode;
	public static String referenceSequenceExtraction(String referencePath) {
        File file = new File(referencePath);
        if (!file.exists()) {
            throw new MyException("fail to open the file and the name is " + file.getAbsolutePath());
        }

        refCode = new char[MAX_CHA_NUM];
        refCodeLen = 0;
        refLowVecBegin = new int[VEC_SIZE];
        refLowVecLength = new int[VEC_SIZE];
        refLowVecLen = 1;

        int lettersLen = 0;
        String info;
        char ch;
        Boolean flag = true;

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            //br.readLine();
			System.out.println(br.readLine());
            while ((info = br.readLine()) != null) {
				//System.out.println(info); //GCGTggNNnnAA
                for (int i = 0; i < info.length(); i++) {
                    ch = info.charAt(i);
					
                    if (Character.isLowerCase(ch)) {  //a
                        if (flag) {
                            flag = false;
                            refLowVecBegin[refLowVecLen] = lettersLen;
                            lettersLen = 0;
                        }
                        ch = Character.toUpperCase(ch);
                    } else {
                        if (!flag) {
                            flag = true;
                            refLowVecLength[refLowVecLen++] = lettersLen;
                            lettersLen = 0;
                        }
                    }

                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        refCode[refCodeLen++] = ch; //G
                    }

                    lettersLen++; //1
                }
            }

            if (!flag) {
                refLowVecLength[refLowVecLen++] = lettersLen;
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
    //}
	//public static void main(String ...a){
		//referenceSequenceExtraction("TestRef.txt");
		System.out.println(refCodeLen);
		
		//for(int i=0;i<refCodeLen;i++)
			//System.out.print(refCode[i]);
		//System.out.println();
		
		String refStr=new String(refCode,0,refCodeLen);
		//System.out.println(refStr);
		
		for(int i=1;i<refLowVecLen;i++)
			System.out.print(refLowVecBegin[i]+" ");
		System.out.println();
		
		for(int i=1;i<refLowVecLen;i++)
			System.out.print(refLowVecLength[i]+" ");
		System.out.println();
		
		return refStr;
	}
}