package sci_commands::bam_split;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_split");

sub bam_split {

@ARGV = @_;
getopts("s:O:A:a:", \%opt);

$die2 = "
scitools bam-split [options] [input bam]
   or    split-bam

Options:
   -A   [STR]   Annotation file (required)
   -O   [STR]   Output prefix (default is bam file prefix)
   -a   [STR]   Comma separated list of annotations to include (requires -A)
                (Default = a seaprate bam for each annotation)
   -s   [STR]   Samtools call (def = $samtools)

";

if (!defined $ARGV[0] || !defined $opt{'A'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

read_annot($opt{'A'});

$header = "";
open IN, "$samtools view -H $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	$header .= "$l\n";
} close IN;

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {$ANNOT_include{$annot} = 1};
} else {
	foreach $annot (keys %ANNOT_count) {$ANNOT_include{$annot} = 1};
}

foreach $annot (keys %ANNOT_include) { 
	$file = "$annot\_out";
	open $file, "| $samtools view -bS - 2>/dev/null > $opt{'O'}.$annot.bam";
	print $file "$header";
}

open IN, "$samtools view $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$null) = split(/:/, $P[0]);
	$annot = $CELLID_annot{$cellID};
	if (defined $ANNOT_include{$annot}) {
		$file = "$annot\_out";
		print $file "$l\n";
	}
} close IN;

foreach $annot (keys %ANNOT_include) {
	$file = $ANNOT_file{$annot};
	close $file;
}

}
1;
