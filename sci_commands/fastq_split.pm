package sci_commands::fastq_split;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_split");

sub fastq_split {

@ARGV = @_;

getopts("A:O:X", \%opt);

$die2 = "
scitools fastq-split [options] read1.fq(.gz) read2.fq(.gz)
   or    split-fastq

Will split your fastq files by annotation. Reads that do not
match the annotastion will go to [out_prefix].unassigned.
Matching will go to [out_prefix].[annot]...

Options:
   -A   [STR]   Annotation file (comma separated for more than one)
                (If multiple, must be non-conflicting)
				(Required)
   -O   [STR]   Output prefix (def = annot file prefix)
   -X           Do not output unassigned reads

";

if (!defined $ARGV[0] || !defined $ARGV[1] || !defined $opt{'A'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $opt{'A'}; $opt{'O'} =~ s/\.annot$//};

$r1 = $ARGV[0]; $r2 = $ARGV[1];

read_annot($opt{'A'});
print STDERR "$annot_count total annotations found.\n";

# setup output files
foreach $annot (keys %ANNOT_count) {
	$out1_handle = "$annot.1";
	open $out1_handle, "| $gzip > $opt{'O'}.$annot.1.fq.gz";
	$HANDLE1{$annot} = $out1_handle;
	$out2_handle = "$annot.2";
	open $out2_handle, "| $gzip > $opt{'O'}.$annot.2.fq.gz";
	$HANDLE2{$annot} = $out2_handle;
}

if (!defined $opt{'X'}) {
	open O1, "| $gzip > $opt{'O'}.unassigned.1.fq.gz";
	open O2, "| $gzip > $opt{'O'}.unassigned.2.fq.gz";
}

if ($r1 =~ /\.gz$/) {open R1, "$zcat $r1 |"} else {open R1, "$r1"};
if ($r2 =~ /\.gz$/) {open R2, "$zcat $r2 |"} else {open R2, "$r2"};

while ($r1tag = <R1>) {
	$r2tag = <R2>; chomp $r1tag; chomp $r2tag;
	$r1seq = <R1>; chomp $r1seq; $null = <R1>; $r1qual = <R1>; chomp $r1qual;
	$r2seq = <R2>; chomp $r2seq; $null = <R2>; $r2qual = <R2>; chomp $r2qual;
	$cellID = $r1tag; $cellID =~ s/^\@//; $cellID =~ s/:.+$//;
	if (defined $CELLID_annot{$cellID}) {
		$annot = $CELLID_annot{$cellID};
		$ANNOT_count{$annot}++;
		$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
		print $out1_handle "$r1tag\n$r1seq\n\+\n$r1qual\n";
		print $out2_handle "$r2tag\n$r2seq\n\+\n$r2qual\n";
	} else {
		$non_annot_count++;
		if (!defined $opt{'X'}) {
			print O1 "$r1tag\n$r1seq\n\+\n$r1qual\n";
			print O2 "$r2tag\n$r2seq\n\+\n$r2qual\n";
		}
	}
}

foreach $annot (keys %ANNOT_count) {
	$out1_handle = $HANDLE1{$annot}; $out2_handle = $HANDLE2{$annot};
	close $out1_handle; close $out2_handle;
	print STDERR "Annot: $annot, count = $ANNOT_count{$annot}\n";
}
if (!defined $opt{'X'}) {close O1; close O2};
print STDERR "Not in annotaiton file: $non_annot_count\n";

}
1;
