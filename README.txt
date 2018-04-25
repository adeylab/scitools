scitools Version: 0.1.1
Adey Lab (www.adeylab.org)

scitools is a set of scripts designed for working with single-cell
combinatorial indexing data. It includes tools to go from fastq files
off of the sequencer (after bcl2fastq) to a processed dataset. It
includes a number of external tools and R packages that are called
by scitools. If scitools is used, be sure to cite those tools!

scitools commands are in the form of [class]-[operation]. Most can
be specified in reverse order, or if the operation is unique to the
class of files or the class type can be determined by the files in
the arguments then just the operation name can be specified.

Dependencies: (command-line callable, can be specified in options)
Run scitools dependencies to check software and R requirements.
The default executables can be altered at the start of the scitools
code in the # GLOBAL DEFAULTS section.

###Executables:
   gzip         For gzipped fastq files. Default: gzip & zcat
                (these are hardoded at the beginning of the scitools
                 code and not available as options)
   samtools     Bam-related commands. Default: samtools
   bedtools     Bed-related commands. Default: bedtools
   Rscript      Numerous R operations. Default: Rscript
   bwa          For alignment only. Default: bwa
   macs2        For atac-callpeak only. Default: macs2
   scitools     Can call itself. Default: scitools

###R packages:
   ggplot2         For plotting commands
   svd             For Latent Semantic Indexing (LSI)
   Rtsne           For tSNE visualization
   methods         For PCA
   dbscan          For density-based clustering

Some default locations / shortcuts:
   Fastq directory (where bcl2fastq outputs fastq files)
      fastq_input_directory
   Output fastq directory (for processed SCI fastq files)
      SCI_fastq_directory
   SCI index file (should comtain all barcodes in the proper format)
      SCI_index_file
   hg19 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      hg19_ref
   hg38 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      hg38_ref
   mm10 ([reference.fa], [reference.fa].fai, [reference.fa].[bwa_index], and [reference].refGene.txt)
      mm10_ref

File Types used by scitools:

   fastq    Standard fastq files, input can be gzipped or not, output fastq
               files will always be gzipped.
   
   bam      Bam files. For scitools, it is expected the barcode sequence
               or cellID is in the read name ([barcode]:[read_number])
               Sam files are not supported to encourage (force) space saving
			   
   annot    Annotation file: tab delimited, 2 columns
               Col 1 = cellID (barcode), could also be feature name
               Col 2 = annotation information (e.g. experimental condition)
            values file A special case of an annotaiton file where the
               annotation is a continuous variable and not discrete.
			   
   matrix   Matrix file, tab delimited
               Row 1 = CellIDs, has 1 less column than all other rows
               Rows = field 1 is name, then 1 field for each cellID
			   
   dims     Dimensions file, tab delimited
               Col 1 = cellID
               Col 2-N = dimensions for each cell
            range    Range format is not a file type but a way to input a
               range of values (e.g. which dimensions to use from a dims file)
               This is a comma separated list and can include dashes for a
               range of values (e.g. 1,3-6,8-10,13)

Typical scitools analysis for sci-ATAC-seq:

   1) Make an annotation file for demultiplexing a sequencing run:
      scitools annot-make -p   # this will output an example annot csv file, edit in excel
      scitools annot-make -P [myExperimentDescriptorFile.csv] > [mySamples.annot]

   2) Demultiplex and reformat raw fastq files:
      scitools fastq-dump -A [mySamples.annot] -R [RUN_NAME, must be in fastq directory]

   3) Align fastq files (if reference defaults are set up use hg38, hg19, or mm10)
      scitools align -t [n threads] [reference_prefix or shortcut] [reads.sampleA.1.fq.gz] [reads.sampleA.2.fq.gz] [sampleA_prefix]

   4) Remove duplicates and create complexity file:
      scitools rmdup [sampleA].bam

   5) Plot complexity file to assess library performance and determine bam filters:
      scitools plot-complexity [sampleA].complexity.txt
         Note: if multiple conditions are present in the sample, create a condition annot file as in step 1
               you can then add '-A [myConditions.annot]' to the above command to plot by conditions.
               It is also possible to specify colors for plotting conditions using -C or -c

   6) Determine performance and filter bam to remove noise reads:
      scitools bam-filter -N [read threshold, ~1000] -C [sampleA].complexity.txt -c [min compl, 0],[max compl, 60]

   7) Examine index-based perfoemance on filtered bam:
      scitools index-performance [sampleA].bbrd.q10.filt.bam

   8) Call ATAC peaks (macs2):
      scitools callpeak -f [reference.fai file] [sampleA].bbrd.q10.filt.bam

   9) Make a counts matrix:
      scitools atac-counts [sampleA].bbrd.q10.filt.bam [sampleA].bbrd.q10.filt.bed

   10) Filter the counts matrix:
      scitools matrix-filter [sampleA counts matrix]

   11) Perform term-frequency inverse-document-frequency tarnsformation:
      scitools tfidf [sampleA filtered matrix]

   12) Latent semantic indexing:
      scitools lsi [sampleA tfidf matrix]

   13) Visualize via tSNE:
      scitools tsne [sampleA LSI matrix]

   14) Plot tSNE:
      scitools plot-dims -A [myConditions.annot] [sampleA LSI].tsne.dims
