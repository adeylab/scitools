package sci_commands::bam_project;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_project");

sub bam_project {

@ARGV = @_;

# Defaults
$training_increment = 0.01;
$complexity_increment = 0.05;
$min_unique = 1000;
$read_increments = 10000000;
$gradient_def = "BuYlRd";

getopts("s:O:n:t:c:R:r:XPfG:", \%opt);

$die2 = "
scitools bam-project [options] [non rmdup bam OR rand_reads.txt]
   or    project-bam
   or    project

Will use the current sequencing to build a projection curve
and project out the expected unique reads at additional levels
of read depth. (use only a NON RMDUP bam). Projects to 0.05
complexity (median of cells).

If run with a bam, will create a [output].rand_reads.txt file
which can be used for subsequent runs for a faster runtime.

This model typically predicts within 1-2% of the actual unique
read counts; however, its intended use is just to get a rough
estimate of the target reads desired to continue sequencing a
library to the desired complexity point.

Note: It is also important to consider the bam input file has
already been barcode matched - therefore to predict raw read
counts for additional sequencing, you need to account for the
unmatching reads. (ie divide by the fraction matching)

Options:
   -O   [STR]   Output dir prefix (default is bam file prefix)
   -n   [INT]   Minimum unique read count to include cellID
                (def = $min_unique)
   -t   [FLT]   Fraction of reads to use as each input data point
                (def = $training_increment)
   -c   [FLT]   Complexity increments (median) to print at
                (def = $complexity_increment)
   -r   [INT]   Read counts to increment projections
                (def = $read_increments)
   -G   [GRD]   Color gradient in plot (def = $gradient_def)
                For all available gradients, run 'scitools gradient'
   -s   [STR]   Samtools call (def = $samtools)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Keep intermediate files (def = delete)
   -f           Filter out models
                (agressive; may fail, but if it works, it may
                 be a more accurate projection)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//; $opt{'O'} =~ s/\.rand_reads.txt$//;
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'t'}) {$training_increment = $opt{'t'}};
if (defined $opt{'n'}) {$min_unique = $opt{'n'}};
if (defined $opt{'c'}) {$complexity_increment = $opt{'c'}};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});

if (-e "$opt{'O'}.read_projections") {
	print STDERR "\n\nWARNING: $opt{'O'}.read_projections Directory Exists! Will Rewrite Contents!\n";
} else {
	system("mkdir $opt{'O'}.read_projections");
}

$readCT_for_proj = 0;
if ($ARGV[0] =~ /\.bam$/) {
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		$readCT_for_proj++;
		@P = split(/\t/, $l);
		($cellID,$null) = split(/:/, $P[0]);
		$seed = 1;
		while (defined $SEED_info{$seed}) {$seed = rand(1e20)};
		if ($P[1] & 4) {$P[4] = 0};
		if ($P[2] =~ /(M|Y|L|K|G|Un|Random|Alt)/i) {$P[4] = 0};
		$P[0] =~ s/#.+$//;
		$SEED_info{$seed} = "$P[0]\t$P[4]\t$P[2]\t$P[3]";
	} close IN;
	open OUT, ">$opt{'O'}.rand_reads.txt";
	foreach $seed (sort keys %SEED_info) {
		print OUT "$SEED_info{$seed}\n";
	} close OUT;
	%SEED_info = ();
	open IN, "$opt{'O'}.rand_reads.txt";
} elsif ($ARGV[0] =~ /\.rand_reads\.txt$/) {
	open IN, "$ARGV[0]";
	while ($l = <IN>) {$readCT_for_proj++};
	close IN;
	open IN, "$ARGV[0]";
} else {
	die "\nCannot determine file type. Must be .bam or .rand_reads.txt\n$die2";
}

$totalCT = 0; $allKept = 0;
$increment = int(($readCT_for_proj*$training_increment)+1);
$report = $increment; $report_num = 0;
%CELLID_total = (); %CELLID_kept = ();
open OUT, ">$opt{'O'}.read_projections/model_all.txt";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$null) = split(/:/, $P[0]);
	$totalCT++;
	if (defined $KEEP{$P[0]}) {
		$CELLID_total{$cellID}++;
		$CELLID_kept{$cellID}++;
		$allKept++;
	} elsif ($P[1] < 10) {} else {
		$CELLID_total{$cellID}++;
		if (!defined $CELLID_POS_ISIZE{$cellID}{"$P[2]:$P[3]"} && !defined $OBSERVED{$P[0]}) {
			$CELLID_POS_ISIZE{$cellID}{"$P[2]:$P[3]"} = 1;
			$KEEP{$P[0]} = 1;
			$CELLID_kept{$cellID}++;
			$allKept++;
		}
		$OBSERVED{$P[0]} = 1;
	}
	if ($totalCT>=$report) {
		$report_num++;
		foreach $cellID (keys %CELLID_total) {
			$frac = sprintf("%.3f", $CELLID_kept{$cellID}/$CELLID_total{$cellID});
			print OUT "$cellID\t$report_num\t$report\t$CELLID_total{$cellID}\t$CELLID_kept{$cellID}\t$frac\t".(1/$CELLID_total{$cellID})."\t".(1/$CELLID_kept{$cellID})."\n";
		}
		$report += $increment;
	}
} close IN; close OUT;

open IN, "$opt{'O'}.read_projections/model_all.txt";
open OUT, ">$opt{'O'}.read_projections/model_increments.txt";
while ($l = <IN>) {
	chomp $l; @P = split(/\t/, $l);
	$cellID = $P[0];
	if ($CELLID_kept{$cellID} >= $min_unique) {
		print OUT "$l\n";
	}
} close IN; close OUT;
if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.read_projections/model_all.txt");
}

open R, ">$opt{'O'}.read_projections/model.r";
print R "
INCREMENTS<-read.table(\"$opt{'O'}.read_projections/model_increments.txt\")
colnames(INCREMENTS) <- c(\"cellID\",\"ReportID\",\"ReportCT\",\"Total\",\"Uniq\",\"Complexity\",\"X\",\"Y\")\n";

$cellIDs_included = 0;
foreach $cellID (keys %CELLID_kept) {
	if ($CELLID_kept{$cellID} >= $min_unique) {
		$CELLID_fraction{$cellID} = $CELLID_total{$cellID}/$totalCT;
		if ($cellIDs_included<1) {
			print R "MODEL<-rbind(c(\"$cellID\",t(as.matrix(lm(subset(INCREMENTS,cellID==\"$cellID\")\$Y ~ subset(INCREMENTS,cellID==\"$cellID\")\$X)\$coefficients))))\n";
		} else {
			print R "MODEL<-rbind(MODEL,c(\"$cellID\",t(as.matrix(lm(subset(INCREMENTS,cellID==\"$cellID\")\$Y ~ subset(INCREMENTS,cellID==\"$cellID\")\$X)\$coefficients))))\n";
		}
		$cellIDs_included++;
	}
}
print R "write.table(MODEL,file=\"$opt{'O'}.read_projections/model.txt\",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=\"\\t\")\n";
close R;

system("$Rscript $opt{'O'}.read_projections/model.r");

$test_count = 10000; $failed_models = 0; $passed_models = 0;
open MODEL, "$opt{'O'}.read_projections/model.txt";
while ($l = <MODEL>) {
	chomp $l;
	($cellID,$int,$slope) = split(/\t/, $l);
	$Vmax = (1/$int); $Km = ($slope/$int);
	$test_uniq = (($Vmax*$test_count)/($Km+$test_count));
	if ($test_count<$test_uniq || !defined $opt{'f'}) {
		$CELLID_Vmax{$cellID} = $Vmax;
		$CELLID_Km{$cellID} = $Km;
		$passed_models++;
	} else {
		$failed_models++;
	}
} close MODEL;

if ($passed_models>0) {
	$fail_frac = sprintf("%.3f", $failed_models/$passed_models);
} else {$fail_frac=1};
if ($fail_frac > 0.1) {
	print STDERR "\nWARNING: $failed_models cell models failed, $passed_models passed ($fail_frac)\n";
	if ($fail_frac > 0.5) {
		die "ERROR: failed model fraction less than 0.5, not enough data to proceed. Exiting.\n";
	}
}

open CELL_PROJ, ">$opt{'O'}.read_projections/cell_projections.txt";
open MED_PROJ, ">$opt{'O'}.read_projections/summary_projections.txt";

$projected_complexity = 1;
$previous_complexity = 1.1;
$projected_read_total = $read_increments;
while ($projected_complexity > 0.05) {
	@COMPLEXITY = (); @UNIQ_RDS = (); $uniq_read_sum = 0;
	foreach $cellID (keys %CELLID_Vmax) {
		$cell_total_read_count = int($projected_read_total*$CELLID_fraction{$cellID});
		$cell_unique_read_count = int(($CELLID_Vmax{$cellID}*$cell_total_read_count)/($CELLID_Km{$cellID}+$cell_total_read_count));
		$cell_unique_read_frac = sprintf("%.4f", $cell_unique_read_count/$cell_total_read_count);
		print CELL_PROJ "$cellID\t$projected_read_total\t$cell_total_read_count\t$cell_unique_read_count\t$cell_unique_read_frac\n";
		push @COMPLEXITY, $cell_unique_read_frac;
		push @UNIQ_RDS, $cell_unique_read_count;
		$uniq_read_sum += $cell_unique_read_count;
	}
	@COMPLEXITYs = sort {$a<=>$b} @COMPLEXITY; $len = @COMPLEXITYs;
	$projected_complexity = sprintf("%.2f", $COMPLEXITYs[int($len/2)]);
	@UNIQ_RDSs = sort {$a<=>$b} @UNIQ_RDS;
	$projected_unique = $UNIQ_RDSs[int($len/2)];
	$projected_mean = int($uniq_read_sum/$len);
	if ($projected_complexity != $previous_complexity) {
		print MED_PROJ "$projected_complexity\t$projected_read_total\t$projected_unique\t$projected_mean\n";
		$previous_complexity = $projected_complexity;
	}
	$projected_read_total += $read_increments;
}
close CELL_PROJ; close MED_PROJ;

open R, ">$opt{'O'}.read_projections/plot_projections.r";
print R "
library(ggplot2)
$gradient_function
CELLS<-read.table(\"$opt{'O'}.read_projections/cell_projections.txt\")
SUMMARY<-read.table(\"$opt{'O'}.read_projections/summary_projections.txt\")
PLT<-ggplot() + theme_bw() +
	geom_line(aes((CELLS\$V5*100),log10(CELLS\$V4),group=CELLS\$V1),alpha=0.15,size=0.15,color=\"lightsteelblue4\") +
	geom_line(aes((SUMMARY\$V1*100),log10(SUMMARY\$V4)),color=\"black\",size=1,linetype=\"dashed\") +
	geom_line(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"black\",size=1) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3)),color=\"black\",size=3) +
	geom_point(aes((SUMMARY\$V1*100),log10(SUMMARY\$V3),color=log10(SUMMARY\$V2)),size=2) +
	scale_color_gradientn(colours=gradient_funct(99)) +
	scale_x_continuous(limits=c(0,100)) +
	scale_y_continuous(limits=c(2,6)) +
	xlab(\"Complexity\") +
	ylab(\"log10 Unique Reads\") +
	labs(color=\"Log10 Total\\nReads\")
ggsave(plot=PLT,filename=\"$opt{'O'}.read_projections/projected_complexity.png\",width=6,height=5)
ggsave(plot=PLT,filename=\"$opt{'O'}.read_projections/projected_complexity.pdf\",width=6,height=5)
"; close R;

system("$Rscript $opt{'O'}.read_projections/plot_projections.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.read_projections/plot_projections.r $opt{'O'}.read_projections/model.r $opt{'O'}.read_projections/model.txt $opt{'O'}.read_projections/model_increments.txt");
}

}
1;
