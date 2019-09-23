package sci_commands::da_enrichment;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("da_enrichment");

sub da_enrichment {

@ARGV = @_;
# Defaults

getopts("O:R:Xg:l:LP:p:n:B:r:q:", \%opt);

$die2 = "
scitools da_enrichment [options] [Input differential accessibility File]
   or    enrichment_da

This script will take the full differential accessibility table output 
(Differential_acc_[tested samples]_as_ref_shrunk_q02_wald.txt) from the matrix-da scitools function.
It will then take user-defined top hits of differential accessibility and look for
known transcription factor motif enrichment compared to all peaks within the table output.

If -L flagged, will perform LOLA analysis for enrichment of accessibility across various
annotations. (e.g. cpgislands,H3K24me3 peaks, promoters etc.)

If -B and -r are flagged, will perform analysis for enrichment of accessibility over a set of provided bed files (such as topic bed files)


Options:
   -O   [STR]   Output Prefix (default file name output is [Input differential accessibility File].homer.motifs)
   -R   [STR]   Rscript call (def = $Rscript)
   -g   [STR]   Genome call [hg38,mm10, or a reference genome fasta]. Default: hg38
   -l   [FLT]   Log2 Fold Change of Accessibility filter to be used for motif discovery. Def: NULL
                Sites have to be equal to or above the given float.
   -p   [FLT]   adjusted p-value filter to be usef for motif discovery. Def: NULL
                Sites have to be equal to or below the given float.
   -q   [FLT]   q-value filter to be used for motif discovery. Def: NULL
                Sites have to be equal to or below the given float.
   -P   [FLT]   Top percentage of differential accessibility peaks to be used for motif discovery. Def: NULL
                This filter is run after -l and -p filtering, if they are specified.
                Sites have to be greater or equal to the top percentage (sorted by lowest to highest q values) 
   -B		[STR]	  If a text file with a list of bedfile names is provided analysis will be done on these with opt r as the ref 
   -r		[STR]	  All peaks bed file
   -L   [FLAG]  If flagged, will perform LOLA analysis on data sets from Sheffield lab. 
   -n   [INT]   Number of cores to be used for LOLA analysis. (Default = 1)
   -X   [FLAG]  Retain intermediate files (Default = delete)

";


if (!defined $opt{'g'}) {$opt{'g'} = "hg38"};
if (!defined $opt{'n'}) {$opt{'n'} = 1};
if ((!defined $opt{'B'} ) && (!defined $opt{'r'} )){
if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
#if (!defined $opt{'P'} && !defined $opt{'p'} && !defined $opt{'l'}) {$opt{'P'} = 5};



open R, ">$opt{'O'}.da.enrichment.r";
print R "
dat<-read.delim(\"$ARGV[0]\")
full_peaks<-as.character(dat\$annotation)
dat\$chr<-sapply(strsplit(full_peaks,\"_\"),\"[\",1)
dat\$start<-sapply(strsplit(full_peaks,\"_\"),\"[\",2)
dat\$end<-sapply(strsplit(full_peaks,\"_\"),\"[\",3)
dat_full_peaks<-cbind(dat\$chr,dat\$start,dat\$end)
write.table(dat_full_peaks,\"$opt{'O'}.full_peaks.bed\",quote=F,sep=\"\\t\",col.names=F,row.names=F)

";

if (defined $opt{'p'}) {
print R "
message(\"Filtering Peaks to those with p-values less than or equal to $opt{'p'}\")
dat<-dat[dat\$padj<=$opt{'p'},]


";
};

if (defined $opt{'q'}) {
print R "
if (\"qval\" %in% colnames(dat)){
message(\"Filtering Peaks to those with q-values less than or equal to $opt{'q'}\")
dat<-dat[dat\$qval<=$opt{'q'},]
}else{
message(\"q Values not found in given file. Skipping q-value filter.\")

}

";
};
if (defined $opt{'l'}) {
print R "
message(\"Filtering Peaks to those with Log2 Fold Change greater than or equal to $opt{'l'}\")
dat<-dat[dat\$log2fold>=$opt{'l'},]
";
};

if (defined $opt{'P'}) {
print R "
message(\"Filtering Peaks to $opt{'P'}\% most significant.\")
dat<-head(dat[order(dat\$qval),],n=floor(($opt{'P'}/100)*nrow(dat)))
dat_sig_peaks<-cbind(dat\$chr,dat\$start,dat\$end)
write.table(dat_sig_peaks,\"$opt{'O'}.top.$opt{'P'}.percent.significant_peaks.bed\",quote=F,sep=\"\\t\",col.names=F,row.names=F)
write.table(dat_sig_peaks,\"$opt{'O'}.significant_peaks.bed\",quote=F,sep=\"\\t\",col.names=F,row.names=F)
";
} else  {
print R "

dat_sig_peaks<-cbind(dat\$chr,dat\$start,dat\$end)
write.table(dat_sig_peaks,\"$opt{'O'}.significant_peaks.bed\",quote=F,sep=\"\\t\",col.names=F,row.names=F)
"; 
};

if (defined $opt{'L'}) {
print R "

library(LOLA)
library(GenomicRanges)
message(\"Conducting LOLA Analysis.\")
#/home/groups/oroaklab/src/R/R-3.5.1/library2/LOLA/nm/t1/resources/regions/LOLACore/$opt{'g'}

dbPath = system.file(\"nm/t1/resources/regions/LOLACore\", \"$opt{'g'}\", package=\"LOLA\")
if (dbPath == \"\"){
  message(\"LOLA Database not found. Skipping LOLA analysis.\")
  }else{
regionDB = loadRegionDB(dbLocation=dbPath)
dat_sig_peaks<-as.data.frame(dat_sig_peaks)
colnames(dat_sig_peaks)<-c(\"chr\",\"start\",\"end\")
sig_peaks<-makeGRangesFromDataFrame(dat_sig_peaks,ignore.strand=TRUE)

dat_full_peaks<-as.data.frame(dat_full_peaks)
colnames(dat_full_peaks)<-c(\"chr\",\"start\",\"end\")
full_peaks<-makeGRangesFromDataFrame(dat_full_peaks,ignore.strand=TRUE)

locResults = runLOLA(userSets=sig_peaks,userUniverse=full_peaks,regionDB=regionDB,cores=$opt{'n'})

writeCombinedEnrichment(locResults, outFolder= \"$opt{'O'}.Enrichment_Results\", includeSplits=TRUE)}

";
};

close R;

system("$Rscript $opt{'O'}.da.enrichment.r");

system("/home/groups/oroaklab/src/homer/bin/findMotifsGenome.pl $opt{'O'}.significant_peaks.bed $opt{'g'} $opt{'O'}.Enrichment_Results -size 500 -bg $opt{'O'}.full_peaks.bed");

if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}.da.enrichment.r");
    system("rm -f $opt{'O'}.full_peaks.bed");
    system("rm -f $opt{'O'}.significant_peaks.bed");
} 


}else {
open IN, "$opt{'B'}";
   while ($l = <IN>) {
      chomp $l;
      $bed_include{$l}++;
  } close IN;

  foreach $bed (sort keys %bed_include) {
  $opt{'O'}=$bed;
  $opt{'O'} =~ s/\.bed//;
  
  if (defined $opt{'L'}) {
	open R, ">$opt{'O'}.da.enrichment.r";
	print R "

	library(LOLA)
	library(GenomicRanges)
	message(\"Conducting LOLA Analysis.\")
	#/home/groups/oroaklab/src/R/R-3.5.1/library2/LOLA/nm/t1/resources/regions/LOLACore/$opt{'g'}

	dbPath = system.file(\"nm/t1/resources/regions/LOLACore\", \"$opt{'g'}\", package=\"LOLA\")
	regionDB = loadRegionDB(dbLocation=dbPath)
	dat_sig_peaks<-read.table(file=\"$bed\",header=F,sep=\"\\t\")
	colnames(dat_sig_peaks)<-c(\"chr\",\"start\",\"end\")
	sig_peaks<-makeGRangesFromDataFrame(dat_sig_peaks,ignore.strand=TRUE)

	dat_full_peaks<-read.table(file=\"$opt{'r'}\",header=F,sep=\"\\t\")
	colnames(dat_full_peaks)<-c(\"chr\",\"start\",\"end\")
	full_peaks<-makeGRangesFromDataFrame(dat_full_peaks,ignore.strand=TRUE)

	locResults = runLOLA(userSets=sig_peaks,userUniverse=full_peaks,regionDB=regionDB,cores=$opt{'n'})

	writeCombinedEnrichment(locResults, outFolder= \"$opt{'O'}.LOLA.Enrichment_Results\", includeSplits=TRUE)

";
	close(R);
	system("$Rscript $opt{'O'}.da.enrichment.r");
  };
   system("/home/groups/oroaklab/src/homer/bin/findMotifsGenome.pl $bed $opt{'g'} $opt{'O'}.homer.Enrichment_Results -size 500 -bg $opt{'r'}"); 
  
  if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}.da.enrichment.r");
	}
  
  }

}





}
1;

