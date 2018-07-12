package sci_commands::met_probesummary;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("met_probesummary");

sub met_probesummary {

@ARGV = @_;

use Getopt::Std; %opt = ();
getopts("B:R:F:O:", \%opt);
$die2 = "
scitools met_probesummary [options] [FeatureSet] [NOMe.sorted.bed.gz File]

sciNOMe generation of methylation summaries per cellID per feature.
[FeatureSet]                            =               Tab-separated Bed-format file of <CHR><START><END><ANNOTATION(optional)>
[NOMe.sorted.bed.gz File]               =               File generated from scitools met-NOMeextract call.

1. Takes a [NOMe.sorted.bed.gz File] either HCH, HCG, or GCH file.
                HCG, HCH, GCH sorted Bed Files are modified bed files:
                <CHROM> <START> <END> <methylated C count> <non-methylated C count> <Percent Methylation at C> <Barcode>

2. Then intersects with a feature set bed window and collapses into two FeatureSet by CellID WindowMatrix.
                One contains C's measured per FeatureSet x CellID [CountWindowMatrix]
                One contains Percent Methylation per FeatureSet x CellID [PercWindowMatrix]

                [Perc/Count] WindowMatrix Files (tab-separated):
                FeatureSet              cellID_1                cellID_2                cellID_3...
                Feature1                Perc/Count              Perc/Count              Perc/Count...
                Feature2                Perc/Count              Perc/Count              Perc/Count...
                Feature3                Perc/Count              Perc/Count              Perc/Count...
                ...

Options:
   -B   [STR]   Bedtools call (Default: $bedtools)
   -R   [STR]   Rscript call (Default: $Rscript)
   -F   [STR]   Feature Set Name for output file appending.
                                (Default: all text within last \"/\" and first \"\.\"\ of [FeatureSet])
   -O   [STR]   Output Prefix.
                                (Default: all text within last \"/\" and first \"\.\"\ of [Input File])
";

if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};

if (defined $opt{'B'}) {$bedtools = $opt{'B'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'F'}) {$opt{'F'}=$ARGV[0]; my @F = split(/\//,$opt{'F'}); $opt{'F'}=$F[-1]; @F = split(/\./,$opt{'F'}); $opt{'F'}=$F[0]};
if (!defined $opt{'O'}) {$opt{'O'}=$ARGV[1]; my @O = split(/\//,$opt{'O'}); $opt{'O'}=$O[-1]; @O = split(/\./,$opt{'O'}); $opt{'O'}=$O[0] .".". $O[1]};

if ($ARGV[1] =~ /sorted.bed$/){
print "Found a Sorted Bed File Which is Non-gzipped.\n";
$command="Sorted Bed File supplied. Skipping Generation.\n";
$input_bed=$ARGV[1];
} elsif ($ARGV[1] =~ /sorted.bed.gz$/) {
print "Found a Sorted Bed File Which is gzipped.\n";
print "Unzipping.\n";
$command="Sorted Bed File supplied. Skipping Generation.\n";
system("gzip -d $ARGV[1]");
my @input_bed = split(/\.gz/,$ARGV[1]);
$input_bed = @input_bed[0];
} 

system($command);
print "Using $input_bed as sorted bed file.\n";

# DETECT IF FEATURE SET FILE IS ANNOTATED (GENE NAMES ETC)
open IN, "head -n 1 $ARGV[0] |";
while ($l = <IN>) {
chomp $l;
@col_count = split(/\t/, $l);
}

@out_bed = split(/\./,$input_bed);
$out_bed = join(/\./,$out_bed[0,1]);

if (scalar @col_count ==4){
        print "Generating tab.bed file with annotation column.\n";
        $command_2="bedtools intersect -wa -wb -sorted -a $ARGV[0] -b $input_bed | bedtools groupby -g 1,2,3,4,10 -c 8,9 -o sum,sum> $out_bed.$opt{'F'}.tab.bed";
} elsif (scalar @col_count ==3){
        print "Generating tab.bed file without annotation column.\n";
        $command_2="bedtools intersect -wa -wb -sorted -a $ARGV[0] -b $input_bed | bedtools groupby -g 1,2,3,9 -c 7,8 -o sum,sum > $out_bed.$opt{'F'}.tab.bed";
} else {die "Please Check that Feature Set File is properly formatted Bed with 1 or 0 annotation columns."}

open LOG, ">>window_generation.log";
$ts = localtime(time);
print LOG "$ts scitools NOMeBed-windowsummary
Feature Set File = $ARGV[0]
Input File = $ARGV[1]

Output Prefix = $out_bed.$opt{'F'}

 Options:
 ";
 foreach $option (keys %opt) {
     print LOG "   $option   $opt{$option}\n";
 }
print LOG "$command\n";
print LOG "$command_2\n";
close LOG;

system($command_2);

open R, ">$out_bed.$opt{'F'}.windowsummarize.r";
print R "
#TAKE OUTPUT TAB.BED FILE AND GENERATE THE COUNT AND MEAN WINDOWMATRIX

args = commandArgs(trailingOnly=TRUE)
library(reshape2)

tab<-read.table(file=args[1],header=F)
write(paste(\"Read in Long Format File:\",as.character(args[1])), stderr())

#check if the windows tab file are annotated
if(ncol(tab) > 6){
colnames(tab)<-c(\"chr\",\"start\",\"end\",\"anno\",\"barc\",\"met\",\"unmet\")
tab\$window<-paste(tab\$chr,tab\$start,tab\$end,tab\$anno,sep=\"_\")
}else{
colnames(tab)<-c(\"chr\",\"start\",\"end\",\"barc\",\"met\",\"unmet\")
tab\$window<-paste(tab\$chr,tab\$start,tab\$end,sep=\"_\")
}

tab\$count<-as.numeric(tab\$met)+as.numeric(tab\$unmet)
tab\$perc<-(as.numeric(tab\$met)/(as.numeric(tab\$met)+as.numeric(tab\$unmet)))*100

write(paste(\"Casting Count Matrix.\", stderr()))
mat_count<-dcast(tab,window~barc,fun.aggregate=sum,value.var=\"count\")
rownames(mat_count)<-mat_count[,1]
mat_count<-mat_count[,-1]

write(paste(\"Casting Percent Matrix.\", stderr()))
mat_perc<-dcast(tab,window~barc,fun.aggregate=mean,value.var=\"perc\")
rownames(mat_perc)<-mat_perc[,1]
mat_perc<-mat_perc[,-1]

write(paste(\"Writing out Matrices to tab-separated data frames.\", stderr()))
write.table(x=mat_count,file=paste0(args[1],\".CountWindowMatrix.txt\"),quote=F,col.names=T,row.names=T,sep=\"\t\")
write.table(x=mat_perc,file=paste0(args[1],\".PercWindowMatrix.txt\"),quote=F,col.names=T,row.names=T,sep=\"\t\")

"; close R;

print "Generating WindowMatrix Count and Percentage Summaries.";
system("$Rscript $out_bed.$opt{'F'}.windowsummarize.r $out_bed.$opt{'F'}.tab.bed");
}
1;