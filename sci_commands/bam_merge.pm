package sci_commands::bam_merge;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_merge");

sub bam_merge {

@ARGV = @_;
getopts("s:H:m:x", \%opt);

$memory = "2G";

$die2 = "
scitools bam-merge [options] [output bam] [input bam1] [input bam2] ...
   or    merge-bam

Simple wrapper for samtools merge, can add BAMID field to reads.
Uses header from first bam listed.

Options:
   -s   [STR]   Samtools call (def = $samtools)
   -H   [BAM]   Use header from this bam instead.
   -m   [MEM]   Samtools sort mex memory K/M/G (def = $memory)
   -x           Do not add bam id.

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'m'}) {$memory = $opt{'m'}};

$ARGV[0] =~ s/\.bam$//;

open OUT, "| $samtools view -bSu - 2>/dev/null | $samtools sort -m $memory -T $ARGV[0].TMP - > $ARGV[0].bam";

if (!defined $opt{'H'}) {$opt{'H'} = $ARGV[1]};
open H, "$samtools view -H $opt{'H'} |";
while ($l = <H>) {print OUT $l};
close H;

for ($bamID = 1; $bamID < @ARGV; $bamID++) {
        open IN, "$samtools view $ARGV[$bamID] |";
        while ($l = <IN>) {
                chomp $l;
                @P = split(/\t/, $l);
                if (!defined $opt{'x'}) {$P[0] .= ":BAMID=$bamID"; $l = join("\t", @P)};
                print OUT "$l\n";
        } close IN;
}

close OUT;

}
1;
