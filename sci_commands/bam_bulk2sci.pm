package sci_commands::bam_bulk2sci;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_bulk2sci");

sub bam_bulk2sci {

# Defaults
$memory = "2G";

@ARGV = @_;
getopts("s:O:Gm:", \%opt);

$die2 = "
scitools bam-bulk2sci [options] [output bam] [name1]=[bam1] [name2]=[bam2] ...
   or    bulk2sci

Will produce a SCI bam file where each input bam is treated
as a separate cell / barcode combo.

Will exclude /(M|Y|L|K|G|Un|Random|Alt|_)/i chroms

Options:
   -G           Do not add RG lines (def = add them)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)
   -s   [STR]   Samtools call (def = $samtools)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};

$ARGV[0] =~ s/\.bam$//;

$RG_header = "";
for ($i = 1; $i < @ARGV; $i++) {
	($name,$bam) = split(/=/, $ARGV[$i]);
	$NAME_bam{$name} = $bam;
	$RG_header .= "\@RG\tID:$name\tSM:$name\tLB:$name\tPL:bulk2sci\n";
}

$header = "";
open H, "$samtools view -H $bam 2>/dev/null |";
while ($l = <H>) {
	chomp $l;
	$header .= "$l\n";
} close H;

open OUT, "| $samtools view -bSu - 2>/dev/null | $samtools sort -m $memory -T $ARGV[0].TMP - > $ARGV[0].bam 2>/dev/null";
print OUT "$header";
if (!defined $opt{'G'}) {
	print OUT "$RG_header";
}
foreach $name (keys %NAME_bam) {
	$num = 0;
	%READ_num = ();
	open IN, "$samtools view -q 10 $NAME_bam{$name} 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt|_)/i) {
			$tag = shift(@P); $tag =~ s/#.+$//;
			$line = join("\t", @P);
			if (defined $READ_num{$tag}) {
				$read_num = $READ_num{$tag};
			} else {
				$READ_num{$tag} = $num;
				$read_num = $num;
				$num++;
			}
			print OUT "$name:$read_num\t$line\n";
		}
	} close IN;
} close OUT;

}
1;
