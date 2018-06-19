package sci_commands::bam_addrg;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_addrg");

sub bam_addrg {

@ARGV = @_;
getopts("s:A:L:O:", \%opt);

$die2 = "
scitools bam-addrg [options] [input bam]
   or    addrg

Produces a sorted bam file with read name and RG barcodes.

Options:
   -O   [STR]   Output (default is bam file prefix, adds .RG.bam)
   -A   [STR]   Annotation file (optional, for faster RG pre-processing)
                (Note: will exclude barcodes not in the annot file)
   -L   [STR]   List of CellIDs (optional, for faster RG pre-processing)
                (Note: will exclude barodes not in the list file)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.bam$//; $opt{'O'} =~ s/\.RG//;

if (defined $opt{'L'}) {
	%CELLID_annot = ();
	open IN, "$opt{'L'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$CELLID_annot{$cellID} = "Cell";
	} close IN;
}
if (defined $opt{'A'}) {read_annot($opt{'A'})};

$RG_header = "";
if (defined $opt{'A'} || defined $opt{'L'}) {
	foreach $cellID (keys %CELLID_annot) {
		$RG_header .= "\@RG\tID:$cellID\tSM:$cellID\tLB:$cellID\tPL:SCI\n";
		$CELLID_include{$cellID} = 1;
	}
} else {
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.*$//;
		$RG_header .= "\@RG\tID:$cellID\tSM:$cellID\tLB:$cellID\tPL:SCI\n";
		$CELLID_include{$cellID} = 1;
	} close IN;
}

$out_header = "";
open IN, "$samtools view -h $ARGV[0] 2>/dev/null |";
open OUT, "| $samtools view -bS - 2>/dev/null > $opt{'O'}.RG.bam";
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^\@/) {
		$out_header .= "$l\n";
	} else {
		if ($out_header ne "done") {
			print OUT $out_header.$RG_header."\@PG\tID:scitools_bam-addrg\tVN:$version\n";
			$out_header = "done";
		} else {
			@P = split(/\t/, $l);
			$cellID = $P[0]; $cellID =~ s/:.*$//;
			print OUT "$l\tRG:Z:$cellID\n";
		}
	}
} close IN; close OUT;

}
1;
