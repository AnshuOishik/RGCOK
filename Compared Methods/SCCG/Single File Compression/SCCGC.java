

import java.io.*;
import java.util.*;

import java.lang.management.*;

public class SCCGC {

	public static String reference, target;
	public static String meta_data = "";
	public static String text = "";
	public static int kmer_length = 21;//21;
	public static int sub_length = 30000, limit = 100;
	public static double T1 = 0.5;
	public static int T2 = 4;
	public static int sot, eot, sor, eor;// start of target,end of target,start of reference,end of reference
	public static int length = 0, mismatch, endref = sub_length - 1;
	public static final int maxchar = 26843545; // 2^28, 268435456
	public static final int maxseq = 10737418; // 2^30, 1073741824
	public static int[] next_kmer = new int[maxchar];
	public static int[] kmer_location = new int[maxseq];// global hash table
	public static boolean local = true;
	public static HashMap<Integer, List<newkmer>> hashmap;


	public SCCGC() throws IOException {
		text = "";
		hashmap = new HashMap<Integer, List<newkmer>>();
	}

	public static long getCPUTime() {
		ThreadMXBean bean = ManagementFactory.getThreadMXBean();
		return bean.isCurrentThreadCpuTimeSupported() ? bean.getCurrentThreadCpuTime() : 0L;
	}
	
	// read sequence for local matching
	public static String LreadSeq(String sequenceFileName) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(sequenceFileName);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));
		StringBuilder stringbuilder = new StringBuilder();
		String line;
		boolean base = false;
		
		meta_data = bufferedreader.readLine();// Chromosome label
		while ((line = bufferedreader.readLine()) != null) {
			stringbuilder.append(line);
			if (!base) {
				length = line.length();
				base = true;
			}
		}
		bufferedreader.close();
		return stringbuilder.toString();
	}
	
	// read reference sequence for global matching
	public static String GreadrefSeq(String sequenceFileName) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(sequenceFileName);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));
		String line;
		int line_length;
		char temp_ch;

		StringBuilder stringbuilder = new StringBuilder();
		line = bufferedreader.readLine();
		while ((line = bufferedreader.readLine()) != null) {
			line_length = line.length();
			for (int i = 0; i < line_length; i++) {
				temp_ch = line.charAt(i);
				if (!Character.isUpperCase(temp_ch)) {
					temp_ch = Character.toUpperCase(temp_ch);
				}
				if (temp_ch != 'N' ) {
					stringbuilder.append(temp_ch);
				}
			}
		}
		bufferedreader.close();
		return stringbuilder.toString();
	}
	
	// read target sequence for global matching
	public static String GreadtarSeq(String sequenceFileName, String fileName) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(sequenceFileName);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));

		FileWriter filewriter = new FileWriter(fileName, true);
		BufferedWriter out = new BufferedWriter(filewriter);

		StringBuilder Llist = new StringBuilder();
		StringBuilder Nlist = new StringBuilder();

		String meta_data, line;

		int line_length, totallength = 0, Llen = 0, Nlen = 0, last_L = 0, last_N = 0, start_L = 0,
				start_N = 0;
		char temp_ch;
		boolean is_low = false;
		boolean is_N = false;
		boolean linelen = true;

		StringBuilder stringbuilder = new StringBuilder();

		meta_data = bufferedreader.readLine();
		out.write(meta_data);
		out.write('\n');

		while ((line = bufferedreader.readLine()) != null) {
			line_length = line.length();

			if (linelen) {
				out.write(Integer.toString(line_length));
				out.write('\n');
				linelen = false;
			}

			for (int i = 0; i < line_length; i++) {
				temp_ch = line.charAt(i);

				if (Character.isLowerCase(temp_ch)) {
					if (is_low) {
						Llen++;
					} else {
						is_low = true;
						start_L = i + totallength;
						Llist.append(start_L - last_L).append(' ');
						Llen++;
					}
					temp_ch = Character.toUpperCase(temp_ch);
				} else {
					if (is_low) {
						Llen--;
						Llist.append(Llen).append(' ');
						last_L = start_L + Llen;
					}
					is_low = false;
					Llen = 0;
				}

				if (temp_ch == 'N') {
					if (is_N) {
						Nlen++;
					} else {
						is_N = true;
						start_N = i + totallength;
						Nlist.append(start_N - last_N).append(' ');
						Nlen++;
					}
				} else {
					if (is_N) {
						Nlist.append(Nlen).append(' ');
						last_N = start_N + Nlen;
					}
					is_N = false;
					Nlen = 0;
					
					stringbuilder.append(temp_ch);
				}
			}
			totallength += line_length;
		}

		if (Llist.length() > 0) {
			if (is_low) {
				Llist.append(Llen);
			} else {
				Llist.deleteCharAt(Llist.length() - 1);
			}
		}

		if (Nlist.length() > 0) {
			if (is_N) {
				Nlist.append(Nlen);
			} else {
				Nlist.deleteCharAt(Nlist.length() - 1);
			}

		}

		bufferedreader.close();

		out.write(Llist.toString());
		out.write('\n');
		out.write(Nlist.toString());
		out.write('\n');
		out.flush();
		out.close();
		return stringbuilder.toString();
	}
	
	// build local hash table
	public static void buildLhashtable(String read, int kmer_length) {
		int length = read.length(), i = kmer_length;
		String Nkmer = "";
		while (i > 0) {
			Nkmer += "N";
			i--;
		}

		int Nkey = Nkmer.hashCode();

		while (i < length - kmer_length + 1) {
			String kmer = read.substring(i, i + kmer_length);
			newkmer newKMer = new newkmer();
			newKMer.setkmerstart(i);
			newKMer.setkmer(kmer);
			int key = kmer.hashCode();
			if (hashmap.containsKey(key)) {
				{
					List<newkmer> list = hashmap.get(key);
					list.add(newKMer);
				}
			} else {
				List<newkmer> list = new ArrayList<newkmer>();
				list.add(newKMer);
				hashmap.put(key, list);
			}
			i++;
			if (key == Nkey) {
				while (read.substring(i, i + 1).equalsIgnoreCase("n") && i < length - kmer_length + 1) {
					i++;
				}
			}
		}
	}

	// build golbal hash table
	public static void buildGhashtable(String read, int kmer_length) {
		int iteration = read.length() - kmer_length + 1;

		for (int i = 0; i < maxseq; i++) {
			kmer_location[i] = -1;
		}
		for (int i = 0; i < iteration; i++) {
			String kmer = read.substring(i, i + kmer_length);
			long key = Math.abs(kmer.hashCode());
			
			if (key == -2147483648) {
				key = 0;
			}
			
			while (key > maxseq - 1) {
				key = key / 2;
			}
			next_kmer[i] = kmer_location[(int) key];
			kmer_location[(int) key] = i;
		}
	}

	public static List<Position> lowercase_position(String sequence) {
		List<Position> list = new ArrayList<Position>();
		boolean successive = false;
		int start = 0, end = 0;

		for (int i = 0; i < sequence.length(); i++) {
			if (Character.isLowerCase(sequence.charAt(i))) {
				if (successive) {
					end += 1;
				} else {
					start = i;
					end += 1;
					successive = true;
				}
			} else {
				if (successive) {
					Position position = new Position();
					position.setstartinTar(start);
					position.setendinTar(end - 1);
					list.add(position);
				}
				successive = false;
				start = 0;
				end = i + 1;
			}
		}

		if (successive) {
			Position position = new Position();
			position.setstartinTar(start);
			position.setendinTar(end);
			list.add(position);
		}
		return list;
	}

	public static int get_incre(int endinRef, int endinTar, int Refendindex, int Tarendindex) {
		int position = 0;
		int endIndex;

		if (Refendindex - endinRef <= Tarendindex - endinTar) {
			endIndex = Refendindex - endinRef + 1;
		} else {
			endIndex = Tarendindex - endinTar + 1;
		}

		for (int i = 1; i < endIndex; i++) {
			if (target.charAt(endinTar + i) == reference.charAt(endinRef + i)) {
				position++;
			} else {
				break;
			}
		}
		return position;
	}
	
	// local match
	public static List<Position> Lmatch(String reference, String target, int kmer_length) {
		List<Position> list = new ArrayList<Position>();
		int index = 0, startIndex;
		int increment, most_incre, key;
		int kmerstart, endinRef, endinTar, Refendindex, Tarendindex;
		String kmer;

		buildLhashtable(reference, kmer_length);

		while (true) {
			increment = 0;
			most_incre = 0;
			kmer = target.substring(index, index + kmer_length);
			key = kmer.hashCode();

			//not match
			if (!hashmap.containsKey(key)) {
				startIndex = Integer.MAX_VALUE;
				index = index + 1;
				if (index + kmer_length > target.length()) {
					break;
				}
				continue;
			}

			List<newkmer> klist = hashmap.get(key);//all matched positions
			startIndex = Integer.MAX_VALUE;
			most_incre = 0;
			for (newkmer newKMer : klist) {
				if (newKMer.getkmer().equals(kmer)) {
					kmerstart = newKMer.getkmerstart();
					endinRef = kmerstart + kmer_length - 1;
					endinTar = index + kmer_length - 1;
					Refendindex = reference.length() - 1;
					Tarendindex = target.length() - 1;
					increment = get_incre(endinRef, endinTar, Refendindex, Tarendindex);
					//search the longest match
					if (klist.size() > 1) {
						if (increment == most_incre) {// choose nearest position
							if (list.size() > 1) {
								int lastEIR = list.get(list.size() - 1).getendinRef();
								if (kmerstart - lastEIR < startIndex - lastEIR)
									startIndex = kmerstart;
							}
						} else if (increment > most_incre) {
							most_incre = increment;
							startIndex = kmerstart;
						}
					} else {
						most_incre = increment;
						startIndex = kmerstart;
						break;
					}
				}
			}

			if (startIndex == Integer.MAX_VALUE) {
				index = index + 1;
				if (index + kmer_length > target.length()) {
					break;
				}
				continue;
			}

			Position position = new Position();
			position.setstartinTar(index);
			position.setendinTar(index + kmer_length + most_incre - 1);
			position.setstartinRef(startIndex);
			position.setendinRef(startIndex + kmer_length + most_incre - 1);
			list.add(position);
			index = index + kmer_length + most_incre + 1;
			if (index + kmer_length > target.length()) {
				break;
			}
		}
		return list;
	}
	
	// global match
	public static List<Position> Gmatch(String reference, String target, int kmer_length) {
		List<Position> list = new ArrayList<Position>();
		int index = 0, startIndex, lastEIR = 0;

		String kmer;

		buildGhashtable(reference, kmer_length);

		while (true) {
			int increment = 0, most_incre = 0;
			kmer = target.substring(index, index + kmer_length);
			int key = Math.abs(kmer.hashCode());

			if (key == -2147483648) {
				key = 0;
			}
			while (key > maxseq - 1) {
				key = key / 2;
			}
			// not match
			if (kmer_location[key] == -1) {
				startIndex = Integer.MAX_VALUE;
				index = index + 1;
				if (index + kmer_length > target.length()) {
					break;
				}
				continue;
			}

			startIndex = Integer.MAX_VALUE;
			most_incre = 0;
			boolean match = false;

			for (int k = kmer_location[key]; k != -1; k = next_kmer[k]) {
				increment = 0;
				String Rkmer = reference.substring(k, k + kmer_length);
				if (!kmer.equals(Rkmer)) {
					continue;
				}
				
				try {
					lastEIR = list.get(list.size() - 1).getendinRef();
				} catch (Exception x) {
					lastEIR = 0;
				}
				//try to find matches in a limit range
				if (k - lastEIR > limit || k - lastEIR < -limit) {
					continue;
				}

				match = true;
				int ref_idx = k + kmer_length;
				int tar_idx = index + kmer_length;
				while (ref_idx < reference.length() && tar_idx < target.length()
						&& (reference.substring(ref_idx, ref_idx + 1).equals(target.substring(tar_idx, tar_idx + 1)))) {
					ref_idx++;
					tar_idx++;
					increment++;
				}
				if (increment == most_incre) {
					if (list.size() > 1) {
						if (k - lastEIR < startIndex - lastEIR)
							startIndex = k;
					}
				} else if (increment > most_incre) {
					most_incre = increment;
					startIndex = k;
				}
			}
			// search matches in whole sequence 
			if (!match) {
				for (int k = kmer_location[key]; k != -1; k = next_kmer[k]) {
					increment = 0;
					String Rkmer = reference.substring(k, k + kmer_length);
					if (!kmer.equals(Rkmer)) {
						continue;
					}

					int ref_idx = k + kmer_length;
					int tar_idx = index + kmer_length;

					while (ref_idx < reference.length() && tar_idx < target.length() && (reference
							.substring(ref_idx, ref_idx + 1).equals(target.substring(tar_idx, tar_idx + 1)))) {
						ref_idx++;
						tar_idx++;
						increment++;
					}
					if (increment == most_incre) {
						if (list.size() > 1) {
							if (k - lastEIR < startIndex - lastEIR)
								startIndex = k;
						}
					} else if (increment > most_incre) {
						most_incre = increment;
						startIndex = k;
					}
				}
			}

			if (startIndex == Integer.MAX_VALUE) {
				index = index + 1;
				if (index + kmer_length > target.length()) {
					break;
				}
				continue;
			}

			Position position = new Position();
			position.setstartinTar(index);
			position.setendinTar(index + kmer_length + most_incre - 1);
			position.setstartinRef(startIndex);
			position.setendinRef(startIndex + kmer_length + most_incre - 1);
			list.add(position);
			index = index + kmer_length + most_incre + 1;
			if (index + kmer_length > target.length()) {
				break;
			}
		}
		return list;
	}

	//list to matching pairs
	public static Position format_matches(List<Position> list) {
		int startinTar, startinRef, endinRef;
		int trouble = 0;
		for (int i = 0; i < list.size(); i++) {

			if (i == 0) {
				startinTar = list.get(i).getstartinTar();
				startinRef = list.get(i).getstartinRef();
				endinRef = list.get(i).getendinRef();
				if (endinRef >= endref) {
					endref = endinRef;
				}

				if (startinTar > 0) {
					String preamble = target.substring(0, startinTar);
					text += preamble + "\n";
					trouble += preamble.length();
				}
				text += "" + (startinRef + sor) + "," + (endinRef + sor) + "\n";
				continue;
			}

			startinTar = list.get(i).getstartinTar();
			startinRef = list.get(i).getstartinRef();
			endinRef = list.get(i).getendinRef();
			if (endinRef >= endref) {
				endref = endinRef;
			}

			int endinTarPrev = list.get(i - 1).getendinTar();
			String mismatch = target.substring(endinTarPrev + 1, startinTar);

			if (mismatch.length() > 0) {
				text += mismatch + "\n";
				trouble += mismatch.length();
			}

			text += (startinRef + sor) + "," + (endinRef + sor) + "\n";
		}
		if (trouble > (sub_length * T1))// T1
		{
			mismatch++;
		}

		return list.get(list.size() - 1);
	}

	public static void format_matches(List<Position> list, String fileName) throws IOException {

		StringBuilder stringbuilder = new StringBuilder();

		int startinTar, startinRef, endinRef, endinTar = 0;
		for (int i = 0; i < list.size(); i++) {
			if (i == 0) {
				startinTar = list.get(i).getstartinTar();
				endinTar = list.get(i).getendinTar();
				startinRef = list.get(i).getstartinRef();
				endinRef = list.get(i).getendinRef();
				if (startinTar > 0) {
					String preamble = target.substring(0, startinTar);
					stringbuilder.append(preamble).append("\n");
				}
				stringbuilder.append(startinRef).append(",").append(endinRef).append("\n");
				continue;
			}

			startinTar = list.get(i).getstartinTar();
			startinRef = list.get(i).getstartinRef();
			endinRef = list.get(i).getendinRef();
			endinTar = list.get(i).getendinTar();
			int endinTarPrev = list.get(i - 1).getendinTar();
			String mismatch = target.substring(endinTarPrev + 1, startinTar);
			if (mismatch.length() > 0) {
				stringbuilder.append(mismatch).append("\n");
			}
			stringbuilder.append(startinRef).append(",").append(endinRef).append("\n");
		}
		if (endinTar < (target.length() - 1)) {
			stringbuilder.append(target.substring(endinTar + 1, (target.length())));
		}
		write(fileName, stringbuilder.toString(), true);
	}

	public static void write(String filename, String text, boolean append) throws IOException {
		//System.out.println("filename "+filename);
		FileWriter filewriter = new FileWriter(filename, append);
		BufferedWriter output = new BufferedWriter(filewriter);
		output.write(text);
		//output.flush();
		output.close();
	}

	public static void write(String filename, List<Position> list, boolean append, String auxiliary)
			throws IOException {
		//System.out.println("filename..."+filename);
		FileWriter filewriter = new FileWriter(filename, append);
		BufferedWriter output = new BufferedWriter(filewriter);
		StringBuilder text = new StringBuilder();
		int temp_end = 0;
		output.write(auxiliary);

		for (Position position : list) {
			int start = position.getstartinTar();
			int end = position.getendinTar();
			text.append(start - temp_end).append(' ').append(end - start).append(' ');
			temp_end = end;
		}
		output.write(text.toString());
		output.write("\n");
		//output.flush();
		output.close();
	}

	public static void postprocess(String filename, String final_file) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(filename);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));
		String line;
		StringBuilder stringbuilder = new StringBuilder();
		List<Integer> list = new ArrayList<Integer>();
		// merge continuous matches
		while ((line = bufferedreader.readLine()) != null) {
			if (line.contains(",")) {
				String[] strings = line.split(",");
				int begin = Integer.parseInt(strings[0]);
				int end = Integer.parseInt(strings[1]);
				if (list.size() > 0) {
					int prevEnd = list.get(list.size() - 1);
					if (begin != prevEnd + 1) {
						stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");
						list = new ArrayList<Integer>();
					}
				}
				list.add(begin);
				list.add(end);
			} else if (line.length() > 0) {
				try {
					stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");
					if (!line.contains("^")) {
						stringbuilder.append(line).append("\n");
					}
				} catch (Exception x) {
					if (!line.contains("^")) {
						stringbuilder.append(line).append("\n");
					}
				}
				list = new ArrayList<Integer>();
			}
		}
		if (list.size() > 0) {
			stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");
		}
		bufferedreader.close();
		
		//delta encoding
		InputStream inputStream = new ByteArrayInputStream(stringbuilder.toString().getBytes());
		datainputstream = new DataInputStream(inputStream);
		bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));

		stringbuilder = new StringBuilder();
		list = new ArrayList<Integer>();
		List<String> stringList = new ArrayList<String>();

		while ((line = bufferedreader.readLine()) != null) {
			stringList.add(line);
		}

		int prev = 0;
		boolean successive = false;
		for (int i = 0; i < stringList.size(); i++) {
			String str = stringList.get(i);

			if (str.contains(",")) {
				if (!successive) {
					String[] strings = str.split(",");
					int begin = Integer.parseInt(strings[0]);
					int end = Integer.parseInt(strings[1]);
					list.add(begin);
					list.add(end - begin);
					prev = end;
					stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");
					successive = true;
				} else {
					String[] strings = str.split(",");
					int begin = Integer.parseInt(strings[0]);
					int end = Integer.parseInt(strings[1]);
					list.add(begin - prev);
					list.add(end - begin);
					prev = end;

					stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");

				}
				list = new ArrayList<Integer>();
			} else if (str.length() > 0) {
				stringbuilder.append(str).append("\n");
			}
		}

		if (list.size() > 0) {
			stringbuilder.append(list.get(0)).append(",").append(list.get(list.size() - 1)).append("\n");
		}
		bufferedreader.close();

		write(final_file, stringbuilder.toString(), true);

	}

	public static class Position {

		private int startinRef;
		private int endinRef;
		private int startinTar;
		private int endinTar;

		public int getstartinRef() {
			return startinRef;
		}

		public void setstartinRef(int startinRef) {
			this.startinRef = startinRef;
		}

		public int getendinRef() {
			return endinRef;
		}

		public void setendinRef(int endinRef) {
			this.endinRef = endinRef;
		}

		public int getstartinTar() {
			return startinTar;
		}

		public void setstartinTar(int startinTar) {
			this.startinTar = startinTar;
		}

		public int getendinTar() {
			return endinTar;
		}

		public void setendinTar(int endinTar) {
			this.endinTar = endinTar;
		}
	}

	public static class newkmer {
	private String kmer;
		private int kmerstart;

		public String getkmer() {
			return kmer;
		}

		public void setkmer(String kmer) {
			this.kmer = kmer;
		}

		public int getkmerstart() {
			return kmerstart;
		}

		public void setkmerstart(int kmerstart) {
			this.kmerstart = kmerstart;
		}

	}

	// 7zip
	public static void use7zip(String filename) throws IOException {
		System.out.println("Hello");
		File zipFile = new File(filename);
		if (!zipFile.exists()) {
			return;
		}
		//String zipname = filename.replaceAll(".txt", "");
		String exec = "7za a " + filename + ".7z " + filename + " -m0=PPMD";
		try{Process process = Runtime.getRuntime().exec(exec);
			process.waitFor();}
		catch(Exception e){
			e.printStackTrace();
		}
	}

	// main
	public static void main(String[] args) throws IOException, InterruptedException {
		kmer_length = 21; //21
		int  controuble = 0;
		boolean is_con = false;
		
		if (args.length != 3) {
			System.out.println("Make sure you have inputted 3 arguments.");
			System.exit(0);
		}
		
		String reference_genome = args[0]; // reference file .fa path
		String target_genome = args[1]; // target file .fa path
		String final_folder = args[2];// output file folder
		//String final_folder = "Result";// output file folder
		String[] chrs = target_genome.split("/");
		String final_file = final_folder + "/" + chrs[(chrs.length-1)];
		System.out.println(final_file);
		Runtime currRuntime = Runtime.getRuntime();
		int nMAX = (int) (currRuntime.maxMemory() / 1024 / 1024);
		System.out.println("MAX RAM :" + nMAX + "M");

		Date startDate = new Date();
		long startCpuTimeNano = getCPUTime();
		System.out.println("Start time: " + startDate);


		while (local) {
			mismatch = 0;
			System.out.println(target_genome + " is compressing...\n");

			File file = new File(final_file);
			if (file.exists()) {
				file.delete();
			}

			//load data
			
			String greference = reference_genome;
			String gtarget = target_genome;
			String tempfile = final_folder + "//interim.txt";

			String reference_seq = LreadSeq(greference);

			String target_seq = LreadSeq(gtarget);
			int target_length = target_seq.length();
			
			if(target_length< sub_length*5){ 
				local=false;
				break;
			}else if(target_length< sub_length*1333){
				T1=0.1;T2=0;
			}
			
			//pre-processing
			
			String auxiliary = meta_data + "\n" + length + "\n";

			List<Position> L_list = lowercase_position(target_seq);

			write(final_file, L_list, false, auxiliary);

			write(final_file, "\n", true);
			reference_seq = reference_seq.toUpperCase();
			target_seq = target_seq.toUpperCase();

			file = new File(tempfile);
			if (file.exists()) {
				file.delete();
			}
			
			sot = 0;
			eot = sub_length;
			sor = 0;
			eor = sub_length;
			Position position = new Position();

			while (true) {

				SCCGC utilities = new SCCGC();
				int kmerlength = SCCGC.kmer_length;

				if (eor > reference_seq.length()) {
					text = (target_seq.substring(sot));
					if (text.length() <= 0) {
						break;
					} else {
						write(tempfile, text, true);
						break;
					}
				}
				if (eot > target_seq.length()) {
					text = (target_seq.substring(sot));
					if (text.length() <= 0) {
						break;
					} else {
						write(tempfile, text, true);
						break;
					}
				}

				reference = (reference_seq.substring(sor, eor));
				target = (target_seq.substring(sot, eot));

				//segmentation-based Local Matching Phase
				
				List<Position> list = Lmatch(reference, target, kmerlength);

				if (list.size() <= 0) {
					kmerlength = 11;   // k`
					list = Lmatch(reference, target, kmerlength);
				}

				if (list.size() <= 0) {
					mismatch++;
					
					if (eot >= target_seq.length() - 1) {
						text = (target_seq.substring(sot));
						write(tempfile, text, true);
						break;
					}
					
					if (is_con) {
						controuble++;
					}
					is_con = true;
					
					text += target + "\n";
					write(tempfile, text, true);
					sot += sub_length;
					eot = sot + sub_length;
					eor += sub_length;
					int difference = target_seq.length() - sot;

					if (difference <= kmer_length) {
						text = (target_seq.substring(sot));
						write(tempfile, text, true);
						break;
					} else if (difference < sub_length) {
						eot = target_seq.length() - 1;
					}

					int difference_ref = reference_seq.length() - sor;

					if (difference_ref < sub_length) {
						eor = reference_seq.length() - 1;
					}
					if (eot >= target_seq.length()) {
						break;
					}
					if (difference_ref <= kmer_length) {
						text = (target_seq.substring(sot));
						if (text.length() <= 0) {
							break;
						} else {
							if (text.length() > (sub_length * T1))// T1
							{
								mismatch++;
							}
							if (mismatch > T2) {// T2
								local = false;
								break;
							}
							write(tempfile, text, true);
							break;
						}
					}
					continue;
				}

				is_con = false;
				if (controuble > 2)
					mismatch-=controuble;				
				controuble = 1;
				
				position = format_matches(list);

				if (mismatch > T2) {// T2
					local = false;
					break;
				}

				sot += position.getendinTar() + 1;
				eot = sot + sub_length;
				sor += endref + 1;
				endref = 0;
				eor = sor + sub_length;

				write(tempfile, text, true);

				int difference = target_seq.length() - sot;

				if (difference <= kmer_length) {
					text = (target_seq.substring(sot));
					if (text.length() <= 0) {
						break;
					} else {
						write(tempfile, text, true);
						break;
					}

				} else if (difference < sub_length) {
					eot = target_seq.length() - 1;
				}

				int difference_ref = reference_seq.length() - sor;

				if (difference_ref < sub_length) {
					eor = reference_seq.length() - 1;
				}
				if (difference_ref <= kmer_length) {
					text = (target_seq.substring(sot));
					if (text.length() <= 0) {
						break;
					} else {
						if (text.length() > (sub_length * T1))// T1
						{
							mismatch++;
						}
						if (mismatch > T2) {// T2
							local = false;
							break;
						}
						write(tempfile, text, true);
						break;
					}
				}

			}

			if (!local) {
//				long localtime = getCPUTime() - startCpuTimeNano;
//				System.out.println( "Local attempt time: " + (double) localtime / 1000000000.0 + " seconds.");
				break;
			}
			
			//post-processing
			
			postprocess(tempfile, final_file);

			use7zip(final_file);

			file = new File(tempfile);
			file.delete();

			file = new File(final_file);
			file.delete();
			
			long taskCPUTimeNano = getCPUTime() - startCpuTimeNano;
			System.out.println( "Compressed time: " + (double) taskCPUTimeNano / 1000000000.0 + " seconds.");
			break;
		}

		if (!local) {
			//System.out.println("Hello");
			
			File file = new File(final_file);
			if (file.exists()) {
				file.delete();
			}
			
			//load data
			
			String greference = reference_genome;
			String gtarget = target_genome;
			String tempfile = final_folder +"//interim.txt";

			String reference_seq = GreadrefSeq(greference);
			String target_seq = GreadtarSeq(gtarget, final_file);

			file = new File(tempfile);
			if (file.exists()) {
				file.delete();
			}

			//pre-processing
			
			int kmer_length = SCCGC.kmer_length;
			
			reference = (reference_seq.substring(0, reference_seq.length()));
			target = (target_seq.substring(0, target_seq.length()));
			
			//global Matching Phase
			
			List<Position> list = Gmatch(reference, target, kmer_length);
			
			//post-processing
			
			format_matches(list, tempfile);
			postprocess(tempfile, final_file);

			use7zip(final_file);

			file = new File(tempfile);
			file.delete();
			
			file = new File(final_file);
			file.delete();

			long taskCPUTimeNano = getCPUTime() - startCpuTimeNano;
			System.out.println( "Compressed time: " + (double) taskCPUTimeNano / 1000000000.0 + " seconds.");
		}
		System.out.println( "Compressed file is: "+final_folder + "/" +chrs[(chrs.length-1)]+ ".7z");
		System.out.println("Done\n" +"-----------------------------------------------------");
	}
}
