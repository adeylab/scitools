package sci_commands::atac_chromvar;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_chromvar");

sub atac_chromvar {

@ARGV = @_;

# Defaults
$count_cutoff = 1000;
$frac_in_peaks = 0.15;
$motif_set = "human_pwms_v2";
%MOTIF_SETS = (
   "human_pwms_v1" => 1,
   "mouse_pwms_v1" => 1,
   "human_pwms_v2" => 1,
   "mouse_pwms_v2" => 1,
   "homer_pwms" => 1,
   "encode_pwms" => 1
);
%GENOMES = (
   "hg19" => "BSgenome.Hsapiens.UCSC.hg19",
   "hg38" => "BSgenome.Hsapiens.UCSC.hg38",
   "mm10" => "BSgenome.Mmusculus.UCSC.mm10"
);

getopts("O:R:Xg:p:M:c:", \%opt);

$die2 = "
scitools atac-chromvar [options] [bam file] [peaks bed file]
   or    chromvar

Bam file MUST have RG lines present. If they are not, add them using
'scitools bam-addrg'. This wrapper wille xecute chromVAR and output
deviation and deviation z-score matrix files, along with a variaiton
file for each motif or peak set assessed.

Options:
   -O   [STR]   Output prefix (default is bam prefix)
                (creates prefix.chromVAR directory)
   -g   [STR]   Genome (hg38, hg19, or mm10; def = hg38)
   -M   [STR]   Motif set (def = $motif_set)
   -c   [INT,FRAC]   count_cutoff,frac_in_peaks (def = $count_cutoff,$frac_in_peaks)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Rewrite existing directory

Note: Requires the chromVAR R package to be installed and loadable.
      It can be found here: github.com/GreenleafLab/chromVAR
      It is also recommended to install chromVARmotifs vaialble here:
	  github.com/GreenleafLab/chromVARmotifs

Motif Sets: human_pwms_v1, mouse_pwms_v1, human_pwms_v2, mouse_pwms_v2,
            homer_pwms, encode_pwms

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bam//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'g'}) {$opt{'g'} = "hg38"};
if (!defined $GENOMES{$opt{'g'}}) {die "\n\nERROR: Genome $opt{'g'} is not a prper genome reference! Exiting!\n"} else {$genome = $GENOMES{$opt{'g'}}};
if (!defined $opt{'M'}) {$opt{'M'} = "human_pwms_v2"};
if (!defined $MOTIF_SETS{$opt{'M'}}) {die "\n\nERROR: The motif set provided ($opt{'M'}) does not exist in chromVARmotifs\n"};
if (defined $opt{'c'}) {($count_cutoff,$frac_in_peaks) = split(/,/, $opt{'c'})};

if (-e "$opt{'O'}.chromVAR") {
	if (defined $opt{'X'}) {
		print STDERR "\n\nWARNING: $opt{'O'}.chromVAR exists, will rewrite contents!\n";
	} else {
		die "\n\nERROR: $opt{'O'}.chromVAR exists! Exiting!\n";
	}
}

system("mkdir $opt{'O'}.chromVAR");

$peakfile = $ARGV[1];
@BAMS = split(/,/, $ARGV[0]);
$bam_list = "c(";
foreach $bam (@BAMS) {$bam_list .= "$bam,"};
$bam_list =~ s/,$/\)/;

open R, ">$opt{'O'}.chromVAR/chromVAR.r";
$ts = localtime(time);
print R "# chromVAR script generated: $ts
# bam: $ARGV[0]
# peak bed: $ARGV[1]
# genome: $opt{'g'}
# motif set: $opt{'M'}
# count_cutoff,frac_cutoff = $count_cutoff,$frac_in_peaks

# load libraries
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
register(MulticoreParam(8))
library(JASPAR2016)
library(ggplot2)

# load motifs and genome
library($genome)
data($opt{'M'})

# read in peaks and filter them
peaks <- getPeaks(\"$ARGV[1]\",sort=TRUE)
peaks <- resize(peaks, width = 500, fix = \"center\")

# make counst matrix and normalize / filter
counts<-getCounts(\"$ARGV[0]\", peaks, paired = TRUE, by_rg = TRUE, format = \"bam\", colData = DataFrame(celltype = \"cells\"))
counts <- addGCBias(counts, genome = $genome)
counts <- filterSamples(counts, min_depth = $count_cutoff, min_in_peaks = $frac_in_peaks, shiny = FALSE)
counts <- filterPeaks(counts, non_overlapping = TRUE)

# make motif index
motif_ix <- matchMotifs($opt{'M'}, counts, genome = $genome)

# calculate & print deviations
dev <- computeDeviations(object = counts, annotations = motif_ix)
write.table(as.matrix(deviations(dev)),file = \"$opt{'O'}.chromVAR/deviations.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
write.table(as.matrix(deviationScores(dev)),file = \"$opt{'O'}.chromVAR/deviation_scores.matrix\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

# calculate & print variabilities
var <- computeVariability(dev)
write.table(as.matrix(var),file = \"$opt{'O'}.chromVAR/variability.txt\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)
plot<-plotVariability(var,use_plotly=FALSE)
ggsave(plot,file = \"$opt{'O'}.chromVAR/variability.png\")

# generate tSNE on the deviations
#############tsne causes segmentation faults############
###################needs to be corrected###############
#tsne <- deviationsTsne(dev, threshold = 1.5, perplexity = 10, shiny = FALSE)
#write.table(as.matrix(tsne),file = \"$opt{'O'}.chromVAR/tsne.dims\", col.names = TRUE, row.names = TRUE, sep = \"\\t\", quote = FALSE)

"; close R;

system("$Rscript $opt{'O'}.chromVAR/chromVAR.r >> $opt{'O'}.chromVAR/chromVAR.log 2>> $opt{'O'}.chromVAR/chromVAR.log");


}
1;
