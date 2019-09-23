package sci_commands::atac_callpeak;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_callpeak");

sub atac_callpeak {

@ARGV = @_;

# Defaults:
$min_feature_size = 500;

getopts("s:m:b:O:l:f:Xe", \%opt);

$die2 = "
scitools atac-callpeak [options] [duplicate removed and filtered bam file]
   or    callpeak(s)

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -l   [INT]   Minimum feature size (def = $min_feature_size);
   -b   [STR]   Bedtools call (def = $bedtools)
   -m   [STR]   Macs2 call (def = $macs2)
   -s   [STR]   Samtools call (def = $samtools)
   -f   [STR]   Fai file for chr lengths (shorcut examples: hg19, hg38, and mm10 if in .cfg)
                If toggled will ensure n peaks extend beyond
   -e 			If toggled, will not filter chromosomes by standard filter pattern (i.e.(M|Y|L|K|G|Un|Random|Alt))
   -X           Retain intermediate files (def = remove)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'l'}) {$min_feature_size = $opt{'l'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

if (defined $opt{'f'}) {
	if (defined $REF{$opt{'f'}}) {
		open FAI, "$REF{$opt{'f'}}.fai";
	} else {
		open FAI, "$opt{'f'}";
	}
	while ($l = <FAI>) {
		chomp $l;
		@P = split(/\t/, $l);
		$CHR_length{$P[0]} = $P[1];
	} close FAI;
}

system("$macs2 callpeak --keep-dup all -t $ARGV[0] -n $opt{'O'} >> $opt{'O'}.macs2.log 2>> $opt{'O'}.macs2.log");

open IN, "$bedtools merge -i $opt{'O'}_peaks.narrowPeak 2>/dev/null |";
open OUT, "| $bedtools sort -i - 2>/dev/null | $bedtools merge -i - > $opt{'O'}.$min_feature_size.tmp 2>/dev/null";

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $opt{'e'} || $P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
		if (($P[2]-$P[1])<$min_feature_size) {
			$mid = ($P[2]+$P[1])/2;
			$start = int($mid-($min_feature_size/2));
			$end = int($mid+($min_feature_size/2));
			print OUT "$P[0]\t$start\t$end\n";
		} else {
			print OUT "$P[0]\t$P[1]\t$P[2]\n";
		}
	}
} close IN; close OUT;

open IN, "$opt{'O'}.$min_feature_size.tmp";
open OUT, ">$opt{'O'}.$min_feature_size.bed";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $opt{'f'}) {
		if ($P[2] < $CHR_length{$P[0]}) {
			print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
		} else {
			print OUT "$P[0]\t$P[1]\t$CHR_length{$P[0]}\t$P[0]_$P[1]_$CHR_length{$P[0]}\n";
		}
	} else {
		print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	}
} close IN; close OUT;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.$min_feature_size.tmp $opt{'O'}_peaks.narrowPeak $opt{'O'}_model.r $opt{'O'}_peaks.xls $opt{'O'}_summits.bed");
}

}
1;
