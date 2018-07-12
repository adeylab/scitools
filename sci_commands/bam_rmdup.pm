package sci_commands::bam_rmdup;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_rmdup");

sub bam_rmdup {

@ARGV = @_;
getopts("s:O:xm:H:", \%opt);

$memory = "2G";

$die2 = "
scitools bam-rmdup [options] [sorted bam file] (additonal sorted bam file(s)) ...
   or    rmdup

Will produce a barcode-based duplicate removed bam and associated
complexity file. Will merge bams if multipe are specified.

Uses header from the first bam provided. If -O is not specified, the name
of the first bam provided will be used as the output prefix.

Will exclude /(M|Y|L|K|G|Un|Random|Alt)/i chroms

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
                 It is highly recommended to specify -O if
                 multiple bams are used as the input.
   -x           If multiple bams, do not include BAMID field
   -m   [MEM]   Samtools sort mex memory K/M/G (def = $memory)
                 (only for multiple bams)
   -H   [BAM]   Use header from this bam instead.
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'m'}) {$memory = $opt{'m'}};
if (!defined $ARGV[1]) {$opt{'x'} = 1};

if (!defined $ARGV[1]) {
	open OUT, "| $samtools view -bS - > $opt{'O'}.bbrd.q10.bam 2>/dev/null";
} else {
	$out_prefix = "$opt{'O'}.bbrd.q10";
	open OUT, "| $samtools view -bSu - | $samtools sort -m $memory -T $out_prefix.TMP - > $out_prefix.bam";
}

if (!defined $opt{'H'}) {$opt{'H'} = $ARGV[0]};
open H, "$samtools view -H $opt{'H'} |";
while ($l = <H>) {print OUT $l};
close H;
$reads_q10_to_other_chr = 0;

for ($bamID = 0; $bamID < @ARGV; $bamID++) {
	open IN, "$samtools view -q 10 $ARGV[$bamID] |";
	while ($l = <IN>) {
		$q10_reads++;
		chomp $l;
		@P = split(/\t/, $l);
		if (!defined $opt{'x'}) {$P[0] .= ":BAMID=$bamID"; $l = join("\t", @P)};
		$barc = $P[0]; $barc =~ s/:.+$//;
		if (defined $KEEP{$P[0]}) {
			print OUT "$l\n";
			$BARC_total{$barc}++;
			$BARC_kept{$barc}++;
			$total_kept++;
		} elsif ($P[1] & 4) {} else {
			if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt)/i) {
				$BARC_total{$barc}++;
				if (!defined $BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]:$P[8]"} && !defined $OBSERVED{$P[0]}) {
					$BARC_POS_ISIZE{$barc}{"$P[2]:$P[3]:$P[8]"} = 1;
					$KEEP{$P[0]} = 1;
					print OUT "$l\n";
					$BARC_kept{$barc}++;
					$total_kept++;
				}
				$OBSERVED{$P[0]} = 1;
			} else {
				$reads_q10_to_other_chr++;
			}
		}
	} close IN;
	%KEEP = ();
}

close OUT;


open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;

open OUT, ">$opt{'O'}.complexity.log";
$args = join(" ", @ARGV);
$ts = localtime(time);
$pct_off = sprintf("%.2f", ($reads_q10_to_other_chr/$q10_reads)*100);
$pct_kept = sprintf("%.2f", ($total_kept/$q10_reads)*100);
print OUT "$ts bam-rmdup completed. Args:\n\t$args\nQ10 reads = $q10_reads\nAligning to other chroms: $reads_q10_to_other_chr ($pct_off %)\nTotal reads kept: $total_kept ($pct_kept %)\n";
close OUT;

}
1;
