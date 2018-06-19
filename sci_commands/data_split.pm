package sci_commands::data_split;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("data_split");

sub data_split {

@ARGV = @_;
getopts("O:", \%opt);

$die2 = "
scitools split-data [options] [input_data_file]

Will take a combined data file and split into scitools format files.
Can be gzipped or not.
   
For cells without information in a file, NA will be reported.

Options:
   -O   [STR]   Output prefix (def = [input - .data].[suffix])

";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.gz//; $opt{'O'} =~ s/\.data//;

if ($ARGV[0] =~ /\.gz$/) {
	open IN, "$zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}

$info = <IN>;
$cellID_data = <IN>; chomp $cellID_data;
@CELLIDS = split(/\t/, $cellID_data);

$out_type = "null";
while ($l = <IN>) {
	chomp $l; @P = split(/\t/, $l);
	if ($l =~ /^#/) {
		if ($out_type eq "dims") {
			for ($i = 2; $i < @CELLIDS; $i++) {
				print OUT "$CELLIDS[$i]";
				foreach $dim (sort {$a<=>$b} @DIMS) {
					print OUT "\t$CELLID_DIMS{$CELLIDS[$i]}{$dim}";
				} print OUT "\n";
			} close OUT;
		} elsif ($out_type eq "matrix") {close OUT};
		$name = $P[1]; $name =~ s/NAME=//;
		if ($P[0] =~ /(ANNOTATION|LAMBDA|VALUES)/) {
			if ($P[0] eq "#ANNOTATION_DATA") {
				open OUT, ">$opt{'O'}.$name.annot";
			} elsif ($P[0] eq "#VALUES_DATA") {
				open OUT, ">$opt{'O'}.$name.values";
			} elsif ($P[0] eq "#LAMBDA_DATA") {
				open OUT, ">$opt{'O'}.$name.lambda";
			}
			$l = <IN>; chomp $l; @P = split(/\t/, $l);
			for ($i = 2; $i < @P; $i++) {
				print OUT "$CELLIDS[$i]\t$P[$i]\n";
			} close OUT;
		} elsif ($P[0] eq "#COMPLEXITY_DATA") {
			%COMP_DATA = ();
			for ($comp_line = 1; $comp_line <= 4; $comp_line++) {
				$l = <IN>; chomp $l; @P = split(/\t/, $l);
				for ($i = 2; $i < @P; $i++) {
					$COMP_DATA{$P[1]}{$CELLIDS[$i]} = $P[$i];
				}
			}
			open OUT, ">$opt{'O'}.$name.complexity.txt";
			for ($i = 2; $i < @P; $i++) {
				print OUT "$COMP_DATA{'RANK'}{$CELLIDS[$i]}\t$CELLIDS[$i]\t$COMP_DATA{'RAW_READS'}{$CELLIDS[$i]}\t$COMP_DATA{'UNIQUE_READS'}{$CELLIDS[$i]}\t$COMP_DATA{'COMPLEXITY'}{$CELLIDS[$i]}\n";
			} close OUT;
		} elsif ($P[0] eq "#DIMENSIONS_DATA") {
			open OUT, ">$opt{'O'}.$name.dims";
			$out_type = "dims";
			%CELLID_DIMS = (); @DIMS = ();
		} elsif ($P[0] eq "#MATRIX_DATA") {
			open OUT, ">$opt{'O'}.$name.matrix";
			print OUT "$CELLIDS[2]";
			for ($i = 3; $i < @CELLIDS; $i++) {
				print OUT "\t$CELLIDS[$i]";
			} print OUT "\n";
			$out_type = "matrix";
		}
	} else {
		if ($out_type eq "dims") {
			($null,$dim) = split(/_/, $P[1]);
			push @DIMS, $dim;
			for ($i = 2; $i < @P; $i++) {
				$CELLID_DIMS{$CELLIDS[$i]}{$dim} = $P[$i];
			}
		} elsif ($out_type eq "matrix") {
			print OUT "$P[1]";
			for ($i = 2; $i < @P; $i++) {
				print OUT "\t$P[$i]";
			} print OUT "\n";
		}
	}
} close IN;

}
1;
