package sci_commands::annot_rename;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_rename");

sub annot_rename {

$joiner = "_";

@ARGV = @_;
getopts("O:x", \%opt);

$die2 = "
scitools annot-rename [options] [annot file to rename] [annot mapping]
   or    rename-annot

Will rename annotaitons using a mapping file with:
   current_annot (tab) new_annot

Options:
   -O   [STR]   Output file name / prefix (def = annotPfx.mapPfx.annot)
   -x           Exclude annotaitons not mapped (default is to include them unchanged)
   
";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/annot$//;
	$opt{'O'} .= $ARGV[1];
	$opt{'O'} =~ s/\.annot$//;
	$opt{'O'} =~ s/\.txt$//;
}

open IN, $ARGV[1];
while ($l = <IN>) {
	chomp $l;
	($old,$new) = split(/\t/, $l);
	$RENAME{$old} = $new;
} close IN;


open OUT, ">$opt{'O'}.annot";
open IN, $ARGV[0];
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	if (defined $RENAME{$annot}) {
		print OUT "$cellID\t$RENAME{$annot}\n";
	} elsif (!defined $opt{'x'}) {
		print OUT "$cellID\t$annot\n";
	}
} close IN; close OUT;

}
1;