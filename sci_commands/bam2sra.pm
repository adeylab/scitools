package sci_commands::bam2sra;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam2sra");

sub bam2sra {

@ARGV = @_;

# Defaults
$memory = "2G";

getopts("s:Mm:A:a:", \%opt);

$die2 = "
scitools bam2sra [output prefix] [unfiltered bam] (unfiltered bam2, optional) (etc...)
   or    bam2geo

This tool will take in a bam file with barcodes as read names and convert it to a
set of fastq files for the paired reads and then a third for the barcode.

Adds .(1|2|ix).fq.gz to the output prefix.

It is designed to prepare reads for GEO / SRA upload.

Options:
   -A   [STR]   Only include barcodes present in provided annotation file
   -a   [STR]   Comma separated list of annotations within -A to include
   -s   [STR]   samtools call (def = $samtools)
   -M           also run md5sum on each fastq (def = no)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)
   
";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'m'}) {$memory = $opt{'m'}};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

open O1, "| $gzip > $ARGV[0].1.fq.gz";
open O2, "| $gzip > $ARGV[0].2.fq.gz";
open O3, "| $gzip > $ARGV[0].ix.fq.gz";

$readNumber = 0; $orphans = 0;
for ($bamID = 1; $bamID < @ARGV; $bamID++) {
	$prevID = "null"; $prevSeq = ""; $prevQual = ""; $match = 0; # 1 means loaded first & looking for mate, 0 means looking for new first read
	%ORPHAN = ();
	open IN, "$samtools sort -m $memory -T $ARGV[0].TMP $ARGV[$bamID] -n -O SAM |";
	while ($l = <IN>) {
		if ($l !~ /^\@/) {
			chomp $l;
			@P = split(/\t/, $l);
			$readID = $P[0]; $readID =~ s/#.+$//;
			$barc = $readID; $barc =~ s/:.+$//;
			if (!defined $opt{'A'} ||
			   (defined $opt{'A'} && !defined $opt{'a'} && defined $CELLID_annot{$barc}) ||
			   (defined $opt{'A'} && defined $opt{'a'} && defined $ANNOT_include{$CELLID_annot{$barc}})) { # annot matching
				$barcQual = $barc; $barcQual =~ s/./#/g;
				if ($readID eq $prevID) {
					$match = 0;
					$readNumber++;
					$outID = sprintf("%012d", $readNumber);
					if ($P[1] & 64) { # first
						print O1 "\@$outID#0/1\n$P[9]\n\+\n$P[10]\n";
						print O2 "\@$outID#0/2\n$prevSeq\n\+\n$prevQual\n";
					} else {
						print O2 "\@$outID#0/1\n$P[9]\n\+\n$P[10]\n";
						print O1 "\@$outID#0/2\n$prevSeq\n\+\n$prevQual\n";
					}
					print O3 "\@$outID#0/3\n$barc\n\+\n$barcQual\n";
				} elsif ($match == 0) { # signal to load first read
					$match = 1;
				} elsif ($match == 1) { # store prev as orphan and signal storing of current read
					if (defined $ORPHAN{$readID}) {
						$orphans--;
						$match = 0; # signal found match and looking for new first read
						($orphanSeq,$orphanQual) = split(/\t/, $ORPHAN{$readID});
						$readNumber++;
						$outID = sprintf("%012d", $readNumber);
						if ($P[1] & 64) { # first
							print O1 "\@$outID#0/1\n$P[9]\n\+\n$P[10]\n";
							print O2 "\@$outID#0/2\n$orphanSeq\n\+\n$orphanQual\n";
						} else {
							print O2 "\@$outID#0/1\n$P[9]\n\+\n$P[10]\n";
							print O1 "\@$outID#0/2\n$orphanSeq\n\+\n$orphanQual\n";
						}
						print O3 "\@$outID#0/3\n$barc\n\+\n$barcQual\n";
					} else {
						$orphans++;
						$ORPHAN{$readID} = "$prevSeq\t$prevQual";
					}
				}
				if ($match == 1) { # load read
					$prevID = $readID; $prevSeq = $P[9]; $prevQual = $P[10];
				}
			}
		}
	}
	close IN;
}

close O1; close O2; close O3;

print STDERR "INFO: $orphans reads did not have a mate.\n";

if (defined $opt{'M'}) {
	system("md5sum $ARGV[0].1.fq.gz > $ARGV[0].1.fq.gz.md5sum");
	system("md5sum $ARGV[0].2.fq.gz > $ARGV[0].2.fq.gz.md5sum");
	system("md5sum $ARGV[0].ix.fq.gz > $ARGV[0].ix.fq.gz.md5sum");
}

}

1;