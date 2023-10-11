class Main{
	public static void main(String[] args) { 
	
		String refStr = RefSeqExtact.referenceSequenceExtraction("BuEb.txt");
		//System.out.println(refStr);
		
		String tarStr = TarSeqExtact.targetSequenceExtraction("AgPh.txt");
		//System.out.println(tarStr);
		
		new OptKmerSize().findOptKmerSize(refStr,tarStr);

	}
}