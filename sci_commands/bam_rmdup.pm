package sci_commands::bam_rmdup;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_rmdup");

sub bam_rmdup {

@ARGV = @_;
getopts("s:O:", \%opt);

$die2 = "
scitools bam-rmdup [options] [sorted bam file]
   or    rmdup

Will produce a barcode-based duplicate removed bam
and associated complexity file.

Will exclude /(M|Y|L|K|G|Un|Random|Alt)/i chroms

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

open OUT, "| $samtools view -bS - > $opt{'O'}.bbrd.q10.bam 2>/dev/null";

open H, "$samtools view -H $ARGV[0] |";
while ($l = <H>) {print OUT $l};
close H;

open IN, "$samtools view -q 10 $ARGV[0] |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($barc,$null) = split(/:/, $P[0]);
	if (defined $KEEP{$P[0]}) {
		print OUT "$l\n";
		$BARC_total{$barc}++;
		$BARC_kept{$barc}++;
	} elsif ($P[1] & 4) {} else {
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
			$BARC_total{$barc}++;
			if (!defined $BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]"} && !defined $OBSERVED{$P[0]}) {
				$BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]"} = 1;
				$KEEP{$P[0]} = 1;
				print OUT "$l\n";
				$BARC_kept{$barc}++;
			}
			$OBSERVED{$P[0]} = 1;
		}
	}
} close IN; close OUT;


open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;

}
1;
