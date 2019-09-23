package sci_commands::met_metfilter;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("met_metfilter");

sub met_metfilter {

@ARGV = @_;

getopts("O:n:p:", \%opt);
$die2 = "
scitools met_metfilter [options] [readID.meth File] [any_C File]

sciMET Filtering of Reads which may be incorrectly bisulfite converted.
[readID.meth File] 			= 		Tab-separated File Generated during scitools met_metextract.
									Format: <cellID><readID><cgmet><cgnonmet><chmet><chnonmet>
[any_C File] 				= 		File generated from scitools met_metextract call.

1. Takes a readID.meth File and calculates number of CH measurements.
2. Filters out reads exceeding [-p] CH methylation with [-n] defined minimum CH counts.
3. Recreates a new readID filtered any_C file.
4. Regenerates readID percentile bin read methylation plots.

Options:

   -O 	[STR]	Output Prefix. 
   				(Default: all text within last \"/\" and first \"\.\"\ of [Input File])
   -n 	[INT]	Minimum number of CH measurements in a readID to apply percentile filter.
   				(Default: 5)
   -p 	[INT]	Maximum CH percent methylation. Reads above will be excluded.
   				(Default: 70 [cutoff based on Luo et al. 2017])	
	
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'n'}) {$opt{'n'} = 5};
if (!defined $opt{'p'}) {$opt{'p'} = 70};
if (!defined $opt{'O'}) {$opt{'O'}=$ARGV[1]; my @O = split(/\//,$opt{'O'}); $opt{'O'}=$O[-1]; @O = split(/\./,$opt{'O'}); $opt{'O'}=$O[0]};

if ($ARGV[1] =~ /gz$/) {
$arg1_dgz=~ s/\.gz\$//;
$command="gzip -d $ARGV[1]; awk '{if(\$5+\$6>=$opt{'n'}&&(\$5/(\$5+\$6)*100>=$opt{'p'})) print \$2}' $ARGV[0] | grep -vf - $arg1_dgz > $opt{'O'}.readID.filt.txt; gzip $arg1_dgz";
} else {
$command="awk '{if(\$5+\$6>=$opt{'n'}&&(\$5/(\$5+\$6)*100>=$opt{'p'})) print \$2}' $ARGV[0] | grep -vf - $ARGV[1] > $opt{'O'}.readID.filt.txt";
}

open LOG, ">readID.meth.filter.log";
$ts = localtime(time);
print LOG "$ts scitools anyC-readfilter
readID.meth File = $ARGV[0]
any_C File = $ARGV[1]

Output Prefix = $opt{'O'}

 Options:
 ";
 foreach $option (keys %opt) {
     print LOG "   $option   $opt{$option}\n";
 }
print LOG "$command\n"; 
close LOG;

system($command);

open R, ">any_C_context_$opt{'O'}.readIDfilt.percentile_bin_meth.r";
print R "
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(plyr)

percent_met_row <- function(x,y){
	working_row<-c(paste(x),nrow(y[y\$cellID==x & y\$value<=0.01,]),
	nrow(y[y\$cellID==x & y\$value>0.01 & y\$value<=0.02,]),
	nrow(y[y\$cellID==x & y\$value>0.02 & y\$value<=0.03,]),
	nrow(y[y\$cellID==x & y\$value>0.03 & y\$value<=0.04,]),
	nrow(y[y\$cellID==x & y\$value>0.04 & y\$value<=0.05,]),
	nrow(y[y\$cellID==x & y\$value>0.05 & y\$value<=0.06,]),
	nrow(y[y\$cellID==x & y\$value>0.06 & y\$value<=0.07,]),
	nrow(y[y\$cellID==x & y\$value>0.07 & y\$value<=0.08,]),
	nrow(y[y\$cellID==x & y\$value>0.08 & y\$value<=0.09,]),
	nrow(y[y\$cellID==x & y\$value>0.09 & y\$value<=0.1,]),
	nrow(y[y\$cellID==x & y\$value>0.1 & y\$value<=0.11,]),
	nrow(y[y\$cellID==x & y\$value>0.11 & y\$value<=0.12,]),
	nrow(y[y\$cellID==x & y\$value>0.12 & y\$value<=0.13,]),
	nrow(y[y\$cellID==x & y\$value>0.13 & y\$value<=0.14,]),
	nrow(y[y\$cellID==x & y\$value>0.14 & y\$value<=0.15,]),
	nrow(y[y\$cellID==x & y\$value>0.15 & y\$value<=0.16,]),
	nrow(y[y\$cellID==x & y\$value>0.16 & y\$value<=0.17,]),
	nrow(y[y\$cellID==x & y\$value>0.17 & y\$value<=0.18,]),
	nrow(y[y\$cellID==x & y\$value>0.18 & y\$value<=0.19,]),
	nrow(y[y\$cellID==x & y\$value>0.19 & y\$value<=0.2,]),
	nrow(y[y\$cellID==x & y\$value>0.2 & y\$value<=0.21,]),
	nrow(y[y\$cellID==x & y\$value>0.21 & y\$value<=0.22,]),
	nrow(y[y\$cellID==x & y\$value>0.22 & y\$value<=0.23,]),
	nrow(y[y\$cellID==x & y\$value>0.23 & y\$value<=0.24,]),
	nrow(y[y\$cellID==x & y\$value>0.24 & y\$value<=0.25,]),
	nrow(y[y\$cellID==x & y\$value>0.25 & y\$value<=0.26,]),
	nrow(y[y\$cellID==x & y\$value>0.26 & y\$value<=0.27,]),
	nrow(y[y\$cellID==x & y\$value>0.27 & y\$value<=0.28,]),
	nrow(y[y\$cellID==x & y\$value>0.28 & y\$value<=0.29,]),
	nrow(y[y\$cellID==x & y\$value>0.29 & y\$value<=0.3,]),
	nrow(y[y\$cellID==x & y\$value>0.3 & y\$value<=0.31,]),
	nrow(y[y\$cellID==x & y\$value>0.31 & y\$value<=0.32,]),
	nrow(y[y\$cellID==x & y\$value>0.32 & y\$value<=0.33,]),
	nrow(y[y\$cellID==x & y\$value>0.33 & y\$value<=0.34,]),
	nrow(y[y\$cellID==x & y\$value>0.34 & y\$value<=0.35,]),
	nrow(y[y\$cellID==x & y\$value>0.35 & y\$value<=0.36,]),
	nrow(y[y\$cellID==x & y\$value>0.36 & y\$value<=0.37,]),
	nrow(y[y\$cellID==x & y\$value>0.37 & y\$value<=0.38,]),
	nrow(y[y\$cellID==x & y\$value>0.38 & y\$value<=0.39,]),
	nrow(y[y\$cellID==x & y\$value>0.39 & y\$value<=0.4,]),
	nrow(y[y\$cellID==x & y\$value>0.4 & y\$value<=0.41,]),
	nrow(y[y\$cellID==x & y\$value>0.41 & y\$value<=0.42,]),
	nrow(y[y\$cellID==x & y\$value>0.42 & y\$value<=0.43,]),
	nrow(y[y\$cellID==x & y\$value>0.43 & y\$value<=0.44,]),
	nrow(y[y\$cellID==x & y\$value>0.44 & y\$value<=0.45,]),
	nrow(y[y\$cellID==x & y\$value>0.45 & y\$value<=0.46,]),
	nrow(y[y\$cellID==x & y\$value>0.46 & y\$value<=0.47,]),
	nrow(y[y\$cellID==x & y\$value>0.47 & y\$value<=0.48,]),
	nrow(y[y\$cellID==x & y\$value>0.48 & y\$value<=0.49,]),
	nrow(y[y\$cellID==x & y\$value>0.49 & y\$value<=0.5,]),

	nrow(y[y\$cellID==x & y\$value>0.5 & y\$value<=0.51,]),
	nrow(y[y\$cellID==x & y\$value>0.51 & y\$value<=0.52,]),
	nrow(y[y\$cellID==x & y\$value>0.52 & y\$value<=0.53,]),
	nrow(y[y\$cellID==x & y\$value>0.53 & y\$value<=0.54,]),
	nrow(y[y\$cellID==x & y\$value>0.54 & y\$value<=0.55,]),
	nrow(y[y\$cellID==x & y\$value>0.55 & y\$value<=0.56,]),
	nrow(y[y\$cellID==x & y\$value>0.56 & y\$value<=0.57,]),
	nrow(y[y\$cellID==x & y\$value>0.57 & y\$value<=0.58,]),
	nrow(y[y\$cellID==x & y\$value>0.58 & y\$value<=0.59,]),
	nrow(y[y\$cellID==x & y\$value>0.59 & y\$value<=0.6,]),
	nrow(y[y\$cellID==x & y\$value>0.6 & y\$value<=0.61,]),
	nrow(y[y\$cellID==x & y\$value>0.61 & y\$value<=0.62,]),
	nrow(y[y\$cellID==x & y\$value>0.62 & y\$value<=0.63,]),
	nrow(y[y\$cellID==x & y\$value>0.63 & y\$value<=0.64,]),
	nrow(y[y\$cellID==x & y\$value>0.64 & y\$value<=0.65,]),
	nrow(y[y\$cellID==x & y\$value>0.65 & y\$value<=0.66,]),
	nrow(y[y\$cellID==x & y\$value>0.66 & y\$value<=0.67,]),
	nrow(y[y\$cellID==x & y\$value>0.67 & y\$value<=0.68,]),
	nrow(y[y\$cellID==x & y\$value>0.68 & y\$value<=0.69,]),
	nrow(y[y\$cellID==x & y\$value>0.69 & y\$value<=0.7,]),
	nrow(y[y\$cellID==x & y\$value>0.7 & y\$value<=0.71,]),
	nrow(y[y\$cellID==x & y\$value>0.71 & y\$value<=0.72,]),
	nrow(y[y\$cellID==x & y\$value>0.72 & y\$value<=0.73,]),
	nrow(y[y\$cellID==x & y\$value>0.73 & y\$value<=0.74,]),
	nrow(y[y\$cellID==x & y\$value>0.74 & y\$value<=0.75,]),
	nrow(y[y\$cellID==x & y\$value>0.75 & y\$value<=0.76,]),
	nrow(y[y\$cellID==x & y\$value>0.76 & y\$value<=0.77,]),
	nrow(y[y\$cellID==x & y\$value>0.77 & y\$value<=0.78,]),
	nrow(y[y\$cellID==x & y\$value>0.78 & y\$value<=0.79,]),
	nrow(y[y\$cellID==x & y\$value>0.79 & y\$value<=0.8,]),
	nrow(y[y\$cellID==x & y\$value>0.8 & y\$value<=0.81,]),
	nrow(y[y\$cellID==x & y\$value>0.81 & y\$value<=0.82,]),
	nrow(y[y\$cellID==x & y\$value>0.82 & y\$value<=0.83,]),
	nrow(y[y\$cellID==x & y\$value>0.83 & y\$value<=0.84,]),
	nrow(y[y\$cellID==x & y\$value>0.84 & y\$value<=0.85,]),
	nrow(y[y\$cellID==x & y\$value>0.85 & y\$value<=0.86,]),
	nrow(y[y\$cellID==x & y\$value>0.86 & y\$value<=0.87,]),
	nrow(y[y\$cellID==x & y\$value>0.87 & y\$value<=0.88,]),
	nrow(y[y\$cellID==x & y\$value>0.88 & y\$value<=0.89,]),
	nrow(y[y\$cellID==x & y\$value>0.89 & y\$value<=0.9,]),
	nrow(y[y\$cellID==x & y\$value>0.9 & y\$value<=0.91,]),
	nrow(y[y\$cellID==x & y\$value>0.91 & y\$value<=0.92,]),
	nrow(y[y\$cellID==x & y\$value>0.92 & y\$value<=0.93,]),
	nrow(y[y\$cellID==x & y\$value>0.93 & y\$value<=0.94,]),
	nrow(y[y\$cellID==x & y\$value>0.94 & y\$value<=0.95,]),
	nrow(y[y\$cellID==x & y\$value>0.95 & y\$value<=0.96,]),
	nrow(y[y\$cellID==x & y\$value>0.96 & y\$value<=0.97,]),
	nrow(y[y\$cellID==x & y\$value>0.97 & y\$value<=0.98,]),
	nrow(y[y\$cellID==x & y\$value>0.98 & y\$value<=0.99,]),
	nrow(y[y\$cellID==x & y\$value>0.99,]))
	return(working_row)
}

IN<-read.table(file=\"$opt{'O'}.readID.filt.txt\",comment.char=\"\")
colnames(IN)<-c(\"cellID\",\"readID\",\"cgmet\",\"cgnonmet\",\"chmet\",\"chnonmet\")
IN\$cg.met<-IN\$cgmet/(IN\$cgmet+IN\$cgnonmet)
IN\$ch.met<-IN\$chmet/(IN\$chmet+IN\$chnonmet)
IN\$tot.met<-(IN\$chmet+IN\$cgmet)/(IN\$chmet+IN\$chnonmet+IN\$cgmet+IN\$cgnonmet)

IN<-na.omit(melt(IN[c(\"cellID\",\"cg.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(file=\"$opt{'O'}.readID.percentile_bin_meth.CG.readIDfilt.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

IN<-na.omit(melt(IN[c(\"cellID\",\"ch.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(file=\"$opt{'O'}.readID.percentile_bin_meth.CH.readIDfilt.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

IN<-na.omit(melt(IN[c(\"cellID\",\"tot.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(file=\"$opt{'O'}.readID.percentile_bin_meth.allC.readIDfilt.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

"; close R;

system("$Rscript any_C_context_$opt{'O'}.readIDfilt.percentile_bin_meth.r");

}
1;