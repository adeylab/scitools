package sci_commands::met_metextract;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("met_metextract");

sub met_metextract {

@ARGV = @_;

$threads = 1;

getopts("O:x:t:", \%opt);
$die2 = "
scitools met_metextract [options] [Bismark Reference Genome] [De-duplicated Filtered Bismark Aligned Bam] [Output Directory]

sciMET methylation extraction per cell.
Takes a Bismark alignment file, following removal of duplicates and filtering via bam-filter.
Generates CG and CH methylation cov files through Bismark.

Bismark Reference Genomes are Bismark pre-processed bisulfite converted genomes. Will accept hg19b, hg38b and mm10b.

Output files are modified bed files:
<CHROM> <START> <END> <methylated C count> <non-methylated C count> <Percent Methylation at C> <Barcode>

Options:
   -x   [STR]   Bismark call (Default: $bismark)
   -t   [INT]   Number of threads to use for extraction (Default: 1)
   -O   [STR]   Output Prefix (default is bam file prefix)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $ARGV[2]) {die $die2};

if (!defined $REF{$ARGV[0]}) {die "\n\nERROR: Genome $ARGV[0] is not a proper genome reference! Exiting!\n"} else {$genome = $REF{$ARGV[0]}};

if (defined $opt{'x'}) {$bismark = $opt{'x'}};
if (!defined $opt{'t'}) {$opt{'t'} = $threads};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bam$//};

$command="${bismark}_methylation_extractor --merge_non_CpG --comprehensive --mbias_off --no_header --gzip --yacht --parallel  $opt{'t'} -s --genome_folder $genome -o $ARGV[2] $ARGV[1] 2>>$ARGV[2]/methylation_extraction.log";

$command_2="awk 'NR==FNR{split(\$1,a,\":\"); barc=a[1];chnonmet[barc]=0;chmet[barc]=0; cgnonmet[barc]=0;cgmet[barc]=0 ;next} {split(\$1,a,\":\"); barc=a[1];if(\$5==\"h\" || \$5==\"x\") chnonmet[barc]++; else if (\$5==\"H\" || \$5==\"X\") chmet[barc]++; else if (\$5==\"Z\") cgmet[barc]++; else if (\$5==\"z\") cgnonmet[barc]++} END {for (i in cgmet) print i,cgmet[i],cgnonmet[i],chmet[i],chnonmet[i]}' <(zcat any_C_context_$opt{'O'}.txt.gz) <(zcat any_C_context_$opt{'O'}.txt.gz) > any_C_context_$opt{'O'}.methylation_extraction.cellID.meth";

$command_3="awk 'NR==FNR{split(\$1,a,\":\"); barc=a[1];chnonmet[\$1]=0;chmet[\$1]=0; cgnonmet[\$1]=0;cgmet[\$1]=0 ;next} {if(\$5==\"h\" || \$5==\"x\") chnonmet[\$1]++; else if (\$5==\"H\" || \$5==\"X\") chmet[\$1]++; else if (\$5==\"Z\") cgmet[\$1]++; else if (\$5==\"z\") cgnonmet[\$1]++} END {for (i in cgmet) print i,cgmet[i],cgnonmet[i],chmet[i],chnonmet[i]}' <(zcat any_C_context_$opt{'O'}.txt.gz) <(zcat any_C_context_$opt{'O'}.txt.gz) | awk 'OFS=\"\\t\"{split(\$1,a,\":\") print a[1],\$1,\$2,\$3,\$4,\$5} -'> any_C_context_$opt{'O'}.methylation_extraction.readID.meth";

open LOG, ">$ARGV[2]/methylation_extraction.log";
$ts = localtime(time);
print LOG "$ts scitools methylationextraction-bam
Alignment File = $ARGV[1]
Reference Genome = $genome
Output Directory = $ARGV[2]

 Options:
 ";
 foreach $option (keys %opt) {
     print LOG "   $option   $opt{$option}\n";
 }
print LOG "$command\n";
print LOG "$command_2\n";
print LOG "$command_3\n";
close LOG;

#system($command);

#Additional Plotting of Global Methylation

open R, ">any_C_context_$opt{'O'}.cellID_globalmethylation_counts.r";
print R "
library(ggplot2)
library(reshape2)
IN<-read.table(file=\"any_C_context_$opt{'O'}.methylation_extraction.cellID.meth\")
colnames(IN)<-c(\"cellID\",\"cgmet\",\"cgnonmet\",\"chmet\",\"chnonmet\")
IN\$count_cg<-IN\$cgmet+IN\$cgnonmet
IN\$count_ch<-IN\$chmet+IN\$chnonmet
IN\$count_c<-IN\$count_cg+IN\$count_ch
IN\$perc_ch<-IN\$chmet/(IN\$count_ch)*100
IN\$perc_cg<-IN\$cgmet/(IN\$count_cg)*100

IN.perc<-melt(IN[c(\"cellID\",\"perc_cg\",\"perc_ch\")])
IN.count<-melt(IN[c(\"cellID\",\"count_cg\",\"count_ch\",\"count_c\")])

PLT<-ggplot(IN.perc,aes(x=variable,y=value))+
theme_bw()+
geom_violin(aes(color=as.factor(variable)))+
geom_jitter(height=0,width=0.5,aes(color=as.factor(variable)),size=0.5,alpha=0.5)+
coord_flip()+
xlab(\"Methylation Mark\")+
scale_x_discrete(labels=c(\"CG\",\"CH\"))+
ylab(\"Global Percent Methylation Per Barcode\")+
scale_color_discrete(guide=F)+ylim(c(0,100))
ggsave(plot=PLT,filename=\"$opt{'O'}.cellID.globalmethylation.pdf\")
ggsave(plot=PLT,filename=\"$opt{'O'}.cellID.globalmethylation.png\")


PLT<-ggplot(IN.count,aes(x=variable,y=log10(value)))+
theme_bw()+
geom_violin(aes(color=as.factor(variable)))+
geom_jitter(height=0,width=0.5,aes(color=as.factor(variable)),size=0.5,alpha=0.5)+
coord_flip()+
xlab(\"Methylation Mark\")+
scale_x_discrete(labels=c(\"CG\",\"CH\",\"C\"))+
ylab(\"Measured Cytosines Per Barcode (log10)\")+
scale_color_discrete(guide=F)
ggsave(plot=PLT,filename=\"$opt{'O'}.cellID.Ccounts.pdf\")
ggsave(plot=PLT,filename=\"$opt{'O'}.cellID.Ccounts.png\")

"; close R;

#system($command_2);

system("$Rscript any_C_context_$opt{'O'}.cellID_globalmethylation_counts.r");

#Additional Plotting of Per Read Methylation



open R, ">any_C_context_$opt{'O'}.readID.percentile_bin_meth.r";
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

IN<-read.table(file=\"any_C_context_$opt{'O'}.methylation_extraction.readID.meth\",comment.char=\"\")
colnames(IN)<-c(\"cellID\",\"readID\",\"cgmet\",\"cgnonmet\",\"chmet\",\"chnonmet\")
IN\$cg.met<-IN\$cgmet/(IN\$cgmet+IN\$cgnonmet)
IN\$ch.met<-IN\$chmet/(IN\$chmet+IN\$chnonmet)
IN\$tot.met<-(IN\$chmet+IN\$cgmet)/(IN\$chmet+IN\$chnonmet+IN\$cgmet+IN\$cgnonmet)

IN<-na.omit(melt(IN[c(\"cellID\",\"cg.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(\"$opt{'O'}.readID.percentile_bin_meth.CG.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

IN<-na.omit(melt(IN[c(\"cellID\",\"ch.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(\"$opt{'O'}.readID.percentile_bin_meth.CH.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

IN<-na.omit(melt(IN[c(\"cellID\",\"tot.met\")]))
percentile_met_plot<-data.frame()
percentile_met_plot<-ldply(unique(IN\$cellID)[1:100], function(x) percent_met_row(x,IN))
row.names(percentile_met_plot)<-percentile_met_plot\$V1
percentile_met_plot<-percentile_met_plot[-c(1)]
pdf(\"$opt{'O'}.readID.percentile_bin_meth.allC.pdf\")
Heatmap(t(apply(data.matrix(percentile_met_plot),1,scale)),cluster_columns=F,row_names_gp = gpar(fontsize = 3),name=\"cellID Z-scored Read Count in Binned Percentile Methylation\")
dev.off()

"; close R;

system($command_3);

system("$Rscript any_C_context_$opt{'O'}.readID.percentile_bin_meth.r");

}
1;