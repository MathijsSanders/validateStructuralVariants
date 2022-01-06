# validateStructuralVariants
Redetect structural variants determined from Illumina WGS data by BRASS in matched PacBio WGS data

## How do I run it?

First, BRASS needs to be run on your sample of interest and matched control to determine somatic structural variants (SV). Second, it is preferred, but not necessary, that the SVs are annotated and filtered by AnnotateBRASS. Finally, the resultant BED files is used as input in combination with the matched PacBio BAM file to determine whether the SVs are detectable in the PacBio data.

### The recommended way

The pre-compiled JAR file is included with the repository, but in case the package needs to be recompiled, please run:

```bash
mvn package clean
```

The following command adds a single column to the BRASS input BED file indicating whether the SV is detectable in the PacBio data (0: absent, 1: present).

```bash
java -Xmx2G -jar validateStructuralVariants.jar --input-pacbio-bam-file input_pacbio_bam_file --input-filtered-bed BRASS_bed_file --output-bed-file output_bed_file --width-extract extract_window --width-search search_window --threads number_of_threads
```

- --input-pacbio-bam-file*: Input PacBio BAM file.
- --input-filtered-bed*: Input BRASS BED file.
- --output-bed-file*: Output BRASS BED file.
- --width-extract: Window for extracting reads around the SV breakpoints (default: 250nt).
- --width-search: Window for searching for reads supporting the SV around the breakpoints (default: 5000nt).
- --threads: Comma-separated list of decoy names (default: 1).    
- --help, -help: Get usage information
- --version, -version: Get current version
- \* Required.

*Dependencies*
- Maven version 3+ (For compiling only).
- Java JDK 11+
