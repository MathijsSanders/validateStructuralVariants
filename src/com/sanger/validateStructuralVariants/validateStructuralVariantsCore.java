package com.sanger.validateStructuralVariants;

import java.util.*;
import java.util.concurrent.*;
import java.io.*;
import java.util.zip.*;

import htsjdk.samtools.*;


public class validateStructuralVariantsCore {
	public validateStructuralVariantsCore(File pacbioBam, File input_bed, String output_bed, int extract_width, int search_width, int threads) {
		System.out.println("Retrieving SVs...");
		var structuralVariants = retrieveStructuralVariants(input_bed);
		var svList = retrieveStatistics(structuralVariants.getStructuralVariants(), pacbioBam, extract_width, search_width, threads);
		Collections.sort(svList, Comparator.comparing(svInformation::getChromLeft).thenComparing(svInformation::getStartLeft));
		writeResults(structuralVariants.getHeader(), svList, output_bed);
	}
	
	private List<svInformation> retrieveStatistics(List<svInformation> svList, File pbBam, int extract_width, int search_width, int threads) {
		var forkJoinPool = new ForkJoinPool(threads);
		try {
			forkJoinPool.submit(() -> svList.parallelStream().forEach(i -> addStatistics(i, pbBam, extract_width, search_width))).get();
		} catch(InterruptedException | ExecutionException e) {
			e.printStackTrace();
			System.exit(-3);
		}
		return svList;
	}
	
	private void writeResults(String header, List<svInformation> svList, String output_bed_file) {
		try {
			BufferedWriter output = new BufferedWriter(new FileWriter(output_bed_file));
			StringBuffer buffer = new StringBuffer();
			int counter = 0;
			output.write(String.join("\t", header, "Total read count", "Total supported reads (ZMW)"));
			output.flush();
			for(svInformation info : svList) {
				buffer.append("\n" + String.join("\t", info.getLine(), Integer.toString(info.size()), Integer.toString(info.getCount())));
				counter++;
				if(counter > 100) {
					output.write(buffer.toString());
					output.flush();
					counter = 0;
					buffer = new StringBuffer();
				}
			}
			if(counter > 0) {
				output.write(buffer.toString());
				output.flush();
			}
			output.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.out.println("Could not write the results to the output BED file.");
			System.exit(-6);
		}
	}
	
	private void addStatistics(svInformation info, File bam, int ew, int sw) {
		var inputSam = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).samRecordFactory(DefaultSAMRecordFactory.getInstance()).open(bam);
		var posSet = new HashSet<String>(1000, 0.99999f);
		SAMRecord currentRecord = null;
		var intervals = QueryInterval.optimizeIntervals(new QueryInterval[]{new QueryInterval(inputSam.getFileHeader().getSequenceIndex(info.getChromRight()), info.getStartRight() - sw, info.getEndRight() + sw)});
		var it = inputSam.query(intervals, false);
		while(it.hasNext()) {
			currentRecord = it.next();
			info.addRead(currentRecord.getReadName());
			if(!posSet.contains(currentRecord.getReadName())) {
				if(fitSV(currentRecord, info, ew)) {
					posSet.add(currentRecord.getReadName());
					info.increaseCount();
				} else if(info.isDeletion() && reviewCigar(currentRecord, info, 50)) {
					posSet.add(currentRecord.getReadName());
					info.increaseCount();
				}
			}
		}
	}
	private boolean reviewCigar(SAMRecord current, svInformation info, int width) {
		var trackerRef = current.getAlignmentStart();
		var it = current.getCigar().getCigarElements().iterator();
		CigarElement tel = null;
		while(it.hasNext()) {
			tel = it.next();
			if(tel.getOperator().toString().equals("H"))
				continue;
			else if(tel.getOperator().toString().equals("S"))
				continue;
			else if(tel.getOperator().toString().equals("M"))
				trackerRef += tel.getLength();
			else if(tel.getOperator().toString().equals("I"))
				continue;
			else if(tel.getOperator().toString().equals("D")) {
				if(Math.abs(info.getStartLeft() - trackerRef) <= width && Math.abs(info.getStartRight() - (trackerRef + tel.getLength())) <= width)
					return true;
				trackerRef += tel.getLength();
			} else {
				System.out.println("You are missing the following operator: " + tel.getOperator().toString());
				System.exit(0);
			}
		}
		return false;
	}
	private boolean fitSV(SAMRecord current, svInformation info, int ew) {
		var searchList = new ArrayList<search>(100);
		var sa = current.getStringAttribute("SA");
		boolean fitLeft = false, fitRight = false;		
		if(sa == null || sa.equals(""))
			return false;
		searchList.add(new search(current.getReferenceName(), current.getAlignmentStart(), current.getAlignmentEnd()));
		var saTokens = sa.split(";");
		for(final String token : saTokens) {
			var elements = token.split(",");
			searchList.add(new search(elements[0], Integer.parseInt(elements[1]), Integer.parseInt(elements[1]) + TextCigarCodec.decode(elements[3]).getReferenceLength()));
		}
		for(final search el : searchList) {
			if(info.getChromLeft().equals(el.chr) && (Math.abs(info.getStartLeft() - el.start) <= ew || Math.abs(info.getStartLeft() - el.end) <= ew))
				fitLeft = true;
			else if(info.getChromRight().equals(el.chr) && (Math.abs(info.getStartRight() - el.start) <= ew || Math.abs(info.getStartRight() - el.end) <= ew))
				fitRight = true;
			if(fitLeft && fitRight)
				return true;
		}
		return false;
	}
	private retrieveResults retrieveStructuralVariants(File brass_bed_file) {
		List<String> header = new ArrayList<String>(1000);
		List<svInformation> svList = new ArrayList<svInformation>(1000);
		try {
			BufferedReader br = null;
			if(isGzipped(brass_bed_file))
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(brass_bed_file))));
			else
				br = new BufferedReader(new FileReader(brass_bed_file));
			String line = null;
			String[] tokens = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")) {
					header.add(line);
					continue;
				}
				tokens = line.split("\t");
				svList.add(new svInformation(line,tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), tokens[3], Integer.parseInt(tokens[4]), Integer.parseInt(tokens[5]), tokens[11].contains("deletion")));
			}
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.out.println("Error reading BRASS BED file.");
			System.exit(-1);
		}
		return new retrieveResults(String.join("\n", header), svList);
	}
	private boolean isGzipped(File db) {
		int magic = 0;
		try {
			var raf = new RandomAccessFile(db, "r");
			magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
			raf.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-9);
		}
		return magic == GZIPInputStream.GZIP_MAGIC;
	}
	private boolean intersect(svInformation info, int search) {
		if(!info.getChromLeft().equals(info.getChromRight()))
			return false;
		else if(info.getEndLeft() + search + 1 < info.getStartRight()-search)
			return false;
		return true;
	}
}
class retrieveResults {
	private String header = null;
	private List<svInformation> svList = new ArrayList<svInformation>(1000);
	public retrieveResults(String header, List<svInformation> svList) {
		this.header = header;
		this.svList = svList;
	}
	public String getHeader() {
		return header;
	}
	public List<svInformation> getStructuralVariants() {
		return svList;
	}
}

class search {
	public String chr = null;
	public Integer start = -1;
	public Integer end = -1;
	public search(String chr, Integer start, Integer end) {
		this.chr = chr;
		this.start = start;
		this.end = end;
	}
}