package sci_commands::bam_aggregate;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_aggregate");

sub bam_aggregate {

@ARGV = @_;
getopts("O:s:r", \%opt);

$die2 = "
scitools bam-aggregate [options] [bam file] [annotation file]
   or    aggregate-bam

Note: if a non counts matrix is provided, it will still sum the values.
Support for comma-separated annotations (i.e. cell aggregation with
oversampling) has been added.

Options:
   -O   [STR]   Output prefix (default is [input annot].aggregate.bam)
   -s   [STR]   Samtools call (def = $samtools)
   -r           Add RG fields (def = no)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.annot$//};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
read_annot($ARGV[1]);

open OUT, "| $samtools view -bS - 2>/dev/null >$opt{'O'}.aggregate.bam";
open HEAD, "$samtools view -H $ARGV[0] 2>/dev/null |";
while ($l = <HEAD>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] !~ /RG/) {
		print OUT "$l\n";
	}
} close HEAD;

if (defined $opt{'r'}) {
	foreach $annot_field (keys %ANNOT_count) {
		@ANNOTS = split(/,/, $annot_field);
		foreach $annot (@ANNOTS) {
			if (!defined $ANNOT_inRG{$annot}) {
				print OUT "\@RG\tID:$annot\tSM:$annot\tLB:$annot\tPL:SCI_aggregate\n";
				$ANNOT_inRG{$annot}++;
			}
		}
	}
}

open IN, "$samtools view $ARGV[0] 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	($cellID,$other) = split(/:/, $P[0]);
	if (defined $CELLID_annot{$cellID}) {
		@ANNOTS = split(/,/, $CELLID_annot{$cellID});
		foreach $annot (@ANNOTS) {
			print OUT "$annot:$cellID.$other";
			for ($i = 1; $i < @P; $i++) {
				if ($P[$i] !~ /^RG:Z:/) {
					print OUT "\t$P[$i]";
				}
			}
			if (defined $opt{'r'}) {
				print OUT "\tRG:Z:$annot";
			}
			print OUT "\n";
		}
	}
} close IN; close OUT;

}
1;
