package com.sanger.validateStructuralVariants;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.validators.PositiveInteger;

public class validateStructuralVariants {
	private static String versionNumber = "0.1";
	@Parameter
	private List<String> parameters = new ArrayList<String>();
	
	@Parameter(names = "--input-pacbio-bam-file", description = "Input PacBio BAM file.", required = true, converter = FileConverter.class, validateWith = FileValidator.class, order=0)
	public File input_pacbio_bam_file = null;
	
	@Parameter(names = "--input-filtered-bed", description = "Input filtered BRASS BED file", required = true, converter = FileConverter.class, validateWith=FileValidator.class, order=1)
	public File input_bed_file = null;
	
	@Parameter(names = "--output-bed-file", description = "Output BED file to store results.", required = true, order=2)
	public String output_bed_file = null;
	
	@Parameter(names = "--width-extract", description = "Window to extract PacBio reads (sv_position +- extract_width).", validateWith = PositiveInteger.class, order=3)
	public Integer extract_width = 250;
	
	@Parameter(names = "--width-search", description = "Window for searching additional PacBio reads (sv +- search_width)", validateWith = PositiveInteger.class, order=4)
	public Integer search_width = 5000;
		
	@Parameter(names = "--threads", description = "Number of threads.", validateWith = PositiveInteger.class, order=5)
	public Integer threads = 1;
	
	@Parameter(names = {"--help","-help"}, help = true, description = "Get usage information", order=6)
	private boolean help;
	
	@Parameter(names = {"--version","-version"}, description = "Get current version", order=7)
	private Boolean version = null;
	
	public static void main(String[] args) {
		var vsv  = new validateStructuralVariants();
		JCommander jCommander = new JCommander(vsv);
		jCommander.setProgramName("validateStructuralVariants.jar");
		JCommander.newBuilder().addObject(vsv).build().parse(args);
		if(vsv.version != null && vsv.version) {
			System.out.println("Validate Illumina SVs in matched PacBio data: " + versionNumber);
			System.exit(0);
		}
		else if(vsv.help) {
			jCommander.usage();
			System.exit(0);
		} else {
			var nThreads = Runtime.getRuntime().availableProcessors();
			if(vsv.threads > nThreads)
				System.out.println("Warning: Number of threads exceeds number of available cores");
			new validateStructuralVariantsCore(vsv.input_pacbio_bam_file, vsv.input_bed_file, vsv.output_bed_file, vsv.extract_width, vsv.search_width, vsv.threads);
		}
		
	}
}
