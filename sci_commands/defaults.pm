package sci_commands::defaults;

use sci_utils::general;
use Getopt::Std; %opt = ();
use FindBin '$RealBin';
use Exporter "import";
@EXPORT = ("defaults");

sub defaults {

$ts = localtime(time);
print "
scitools defaults (or check-defaults, defaults-check)
Run at $ts
Scitools location: $RealBin

Default command executable calls:
   gzip:        $gzip
   zcat:        $zcat
   bwa:         $bwa
   samtools:    $samtools
   scitools:    $scitools
   macs2:       $macs2
   bedtools:    $bedtools
   R scripts:   $Rscript
   Py scripts:  $Pscript
   bowtie2:     $bowtie2
   bismark:     $bismark

Reference genomes:\n";
foreach $refID (sort keys %REF) {
	print "   $refID = $REF{$refID}\n";
}
print "
Other defaults:\n";
foreach $var (sort keys %VAR) {
	print "   $var = $VAR{$var}\n";
}
print "\n";

}
1;
