//Finding optimal kmer length (random approach)
import java.nio.file.*; 
import java.util.*; 
class OptKmerSize{
	String refStr,tarStr;
	public void findOptKmerSize(String refStr,String tarStr){ 
		System.out.println("Optimal k-mer calculation");
		this.refStr=refStr;
		this.tarStr=tarStr;

		Set<String> kmer=new HashSet<>();
		Map<String,Integer> kmerc=new HashMap<>();
		
		//Random Selection of 'number_of_kmer' of size 'kmerlen' take at runtime from reference sequence
		int i=0, kmerlen=10, min = 0, max = refStr.length()-kmerlen,count=0,number_of_kmer=20000 ;//refStr.length()/kmerlen;	
		for(;i<refStr.length();){
			i = (int)(Math.random()*(max-min+1)+min); //In the above formula, the min value is inclusive while the max value is exclusive.   
			kmer.add(refStr.substring(i,i+kmerlen));
			if(kmer.size()>=number_of_kmer)
				break;
		}
		
		//Count the occurance of each random kmer in target sequence
		Iterator<String> itr = kmer.iterator(); 
        while (itr.hasNext()){ 
			String str=itr.next();
			//count=0;
			//System.out.println(str+" "+count); 
			for(i=0;i<tarStr.length();){
				if((tarStr.length()-i>=kmerlen)&&str.equals(tarStr.substring(i,i+kmerlen))){
					if(kmerc.containsKey(str)){
						count=kmerc.get(str);
						kmerc.put(str,++count);
						i+=kmerlen;
					}
					else if(!kmerc.containsKey(str)){
						count=1;
						kmerc.put(str,count);
						i+=kmerlen;
					}
				} 
				else {
					i+=kmerlen;
				}
			}  
		}
		//System.out.println(kmerc);
		
		//Difference b/w Original size & mapped character size
		int original_length=tarStr.length(),mapped_length=0;
		
		List<Integer> mapValues = new ArrayList<>(kmerc.values());
		Collections.sort(mapValues);
		Collections.reverse(mapValues);
		//System.out.println(mapValues);
		
		Iterator<Integer> valueIt = mapValues.iterator();
		while (valueIt.hasNext()){
			int value=valueIt.next();
			mapped_length+=value*kmerlen;
			if(mapped_length>original_length){
				mapped_length-=value*kmerlen;
				break;
			}
		}
		System.out.println(original_length);
		System.out.println(mapped_length);
		System.out.println(original_length-mapped_length);
		
		//Original size vs. mapped character size
		/*Set<Map.Entry<String,Integer>> set = kmerc.entrySet();
		for(Map.Entry<String,Integer> me:set){
			temp=me.getValue();
			mapped_length+=temp*kmerlen;
		}
		System.out.println(original_length-mapped_length);*/
	} 
}