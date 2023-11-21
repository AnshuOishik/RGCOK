import java.io.*;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class SCCGD {

	public static List<Position> L_list;
	public static List<Position> N_list;
	public static String meta_data;
	public static int length;
	public static int line_length;

	public static String readSeq(String sequenceFileName) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(sequenceFileName);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));
		String line;

		StringBuilder stringbuilder = new StringBuilder();
		line = bufferedreader.readLine();
		while ((line = bufferedreader.readLine()) != null) {
			stringbuilder.append(line);
		}
		bufferedreader.close();
		return stringbuilder.toString();
	}

	public static String readrefSeq(String sequenceFileName) throws IOException {
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
				if (temp_ch != 'N') {
					stringbuilder.append(temp_ch);
				}
			}
		}
		bufferedreader.close();
		return stringbuilder.toString();
	}

	public static void use7zip(String filename, String Dfilename) throws IOException {

		File zipFile = new File(filename);
		if (!zipFile.exists()) {
			return;
		}
		String exec = "7za e " + filename + " -o" + Dfilename + " -aos";
		try {
			Process process = Runtime.getRuntime().exec(exec);
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static class Position {
		int startinRef;
		int endinRef;
	}

	public static StringBuilder reconstruct(String inFileName, String reference) throws IOException {
		FileInputStream fileinputstream = new FileInputStream(inFileName);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));

		String line;
		StringBuilder stringbuilder = new StringBuilder();

		int prev_end = 0, index = 0;

		while ((line = bufferedreader.readLine()) != null) {

			if (index == 4) {

				if (line.contains(",")) {
					String[] strings = line.split(",");
					int begin = Integer.parseInt(strings[0]);
					int end = Integer.parseInt(strings[1]);
					begin += prev_end;
					end += begin;
					try {
						String text = reference.substring(begin, end + 1);
						stringbuilder.append(text);
						prev_end = end;
					} catch (Exception x) {
						System.out.println("come" + prev_end + begin + end);
					}
				} else if (line.length() > 0) {
					stringbuilder.append(line);
				}
				continue;

			} else if (index < 4) {
				index++;
				continue;
			}
		}
		bufferedreader.close();
		return stringbuilder;
	}

	public static void write(String fileName, String text) throws IOException {
		FileWriter filewriter = new FileWriter(fileName, false);
		BufferedWriter out = new BufferedWriter(filewriter);
		out.write(text);
		out.flush();
		out.close();
	}

	public static long getCpuTime() {
		ThreadMXBean bean = ManagementFactory.getThreadMXBean();
		return bean.isCurrentThreadCpuTimeSupported() ? bean.getCurrentThreadCpuTime() : 0L;
	}

	public static void test(String ref, String tar) {
		for (int i = 0; i < ref.length(); i++) {
			if (ref.charAt(i) != tar.charAt(i)) {
				System.out.println("error:" + i);
			}
		}
	}

	public static void main(String[] args) throws IOException {

		if (args.length != 3) {
			System.out.println("Make sure you have inputted 3 arguments.");
			System.exit(0);
		}

		String reference_file = args[0];
		String in_file_7z = args[1];
		String[] chrs = in_file_7z.split("/");
		//String final_file = args[2] + "/" + chrs[(chrs.length-1)];
		String final_file = args[2] + "/" + chrs[(chrs.length-1)].replaceAll(".7z", "");
		
		//String in_file_name = final_file.replaceAll(".fa", ".txt");
				
		Date startDate = new Date();
		long startCpuTimeNano = getCpuTime();
		System.out.println("Start time: " + startDate);

		System.out.println(args[1] + " is decompressing...");

		File file = new File(final_file);
		file.delete();

		use7zip(in_file_7z, args[2]+"/");

		FileInputStream fileinputstream = new FileInputStream(final_file);
		DataInputStream datainputstream = new DataInputStream(fileinputstream);
		BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(datainputstream));
		String line;
		int Pindex = 0;

		L_list = new ArrayList<Position>();
		N_list = new ArrayList<Position>();

		// load pre-processing data
		while ((line = bufferedreader.readLine()) != null) {
			if (Pindex == 1) {
				line_length = Integer.parseInt(line);
				Pindex++;
				continue;
			} else if (Pindex == 2) {
				if (line == null || line.length() <= 0) {
					Pindex++;
					continue;
				}
				String[] strings = line.split(" ");
				int j = 0, prev = 0;
				while (j < strings.length) {
					Position position = new Position();
					position.startinRef = prev + Integer.parseInt(strings[j]);
					j++;
					position.endinRef = position.startinRef + Integer.parseInt(strings[j]);
					j++;
					L_list.add(position);
					prev = position.endinRef;
				}
				Pindex++;
				continue;
			} else if (Pindex == 3) {
				if (line == null || line.length() <= 0) {
					Pindex++;
					continue;
				}
				String[] strings = line.split(" ");
				int j = 0, prev = 0;
				while (j < strings.length) {
					Position position = new Position();
					position.startinRef = prev + Integer.parseInt(strings[j]);
					j++;
					position.endinRef = position.startinRef + Integer.parseInt(strings[j]);
					j++;
					N_list.add(position);
					prev = position.endinRef;
				}
				Pindex++;
				continue;
			} else if (Pindex == 0) {
				meta_data = line + "\n";
				Pindex++;
				continue;
			}
		}
		bufferedreader.close();

		String reference = "";

		if (N_list.size() > 0) {
			reference = readrefSeq(reference_file);
		} else if (N_list.size() <= 0) {
			reference = readSeq(reference_file);
			reference = reference.toUpperCase();
		}

		// reconstruct target sequence
		StringBuilder target_string = reconstruct(final_file, reference);

		// complete target
		StringBuilder interim = new StringBuilder();
		int index = 0, iterator = 0;
		boolean accept = false;
		for (Position position : N_list) {

			while (true) {

				if (iterator >= position.startinRef && iterator < position.endinRef) {
					interim.append('N');
				} else if (iterator < position.startinRef && !accept) {
					interim.append(target_string.charAt(index));
					index = index + 1;
					if (index >= target_string.length()) {
						accept = true;
					}
				}
				iterator++;

				if (iterator >= position.endinRef) {
					break;
				}
			}
		}
		if (N_list.size() > 0) {
			while (index < target_string.length()) {
				interim.append(target_string.charAt(index));
				index++;
				iterator++;
			}
		}

		if (interim.length() > 0) {
			target_string = new StringBuilder(interim);
			// System.out.println("FINAL SIZE: " + target_string.length());
		} else {
			// System.out.println("FINAL SIZE: " + target_string.length());
		}

		for (Position position : L_list) {
			for (int j = position.startinRef; j <= position.endinRef; j++) {
				System.out.println(position.startinRef+" "+position.endinRef);
				if (j == target_string.length()) {
					break;
				}
				target_string.setCharAt(j, Character.toLowerCase(target_string.charAt(j)));
			}
		}

		String final_string = target_string.toString();

		target_string = new StringBuilder();
		target_string.append(final_string.charAt(0));

		for (int t = 1; t < final_string.length(); t++) {
			if (t % line_length == 0) {
				target_string.append("\n");
			}
			target_string.append(final_string.charAt(t));
		}

		final_string = meta_data + (target_string.append("\n")).toString();

		// test(res,tar);//debug

		write(final_file, final_string);

//		file = new File(in_file_name);
//		if (file.exists()) {
//			file.delete();
//		}

		long taskCpuTimeNano = getCpuTime() - startCpuTimeNano;
		System.out.println("Decompressed time: " + (double) taskCpuTimeNano / 1000000000.0 + " seconds.");
		System.out.println("Decompressed file is " + args[2] + "/" + chrs[(chrs.length-1)]);
		System.out.println("Done\n" + "----------------------------------------------------------------------------");
	}
}
