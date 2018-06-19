package sci_commands::annot_make;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_make");

sub annot_make {

@ARGV = @_;
getopts("O:I:P:ph", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");

$die2 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...
   or    make-annot

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -P   [STR]   Plate descriptor file (instead of written descriptors)
   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)
   -h           More detailed description of plate / combo specification

";

$die3 = "
scitools annot-make [options] [annotation_description_1] [annotaiton_description_2] ...

Options:
   -O   [STR]   Output annotation file (default = STDOUT)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -P   [STR]   Plate descriptor file (instead of written descriptors)
   -p           Print out sample plate file for modificaiton and exit
                (ExamplePlateDescriptor.csv)
   -h           More detailed description of plate / combo specification

Annotation descriptors are provided as:

[Annotation_Name]+[Transposase or PCR descriptor]+[Transposase or PCR descriptor]+[etc...]
  must provide at least 1 transposase descriptor and at least 1 PCR descriptor for
  each annotation.

Transposase descriptor:
[TN5],[TN5 index set combination]=[row],[row]:[columns],[row],etc...
 or [NEX]

PCR Descriptor:
[PCR],[PCR index set combination]=[row]:[columns],[row],[row]:[columns],etc...

Before the \"=\" are two comma separated fields. The first is TN5/NEX or PCR to
define what stage of indexing is specified. The second is the index set IDs in
the order of i5 and then i7, (e.g. AA).

After the \"=\" are a series of comma separated fields, where each field
corresponds to a column of a plate, where columns are numbered 1-12. The column
field can be further specified with a subset of rows after a colon, where rows
are A-H.

Note: each index stage can be specified multiple times. The result is an all by
all of the transposase and PCR index sets that are specified.

Example:
My_Sample_1+NEX,AA=ALL+PCR,AC=1-8+PCR,AD=1:A-D,2-5,6:ACDFH,7-12 My_Sample_2+NEX,AA=ALL+NEX,BB=ALL+PCR,GE=1-8

   row specifications can be listed as a range OR letters with no spacing (e.g. ABCGH)
   commas in the column specification are ONLY for separating out columns NOT rows
      (i.e. 1,2,3:A,B would NOT be OK because of the comma between row letters)

";

if (defined $opt{'h'}) {die $die3};

if (defined $opt{'p'}) {
open OUT, ">ExamplePlateDescriptor.csv";
	print OUT "#NEX,MySampleID1,AA,Partial
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
1,1,1,1,1,1,0,0,0,0,0,0
#NEX,MySampleID1,BB,All
#PCR,MySampleID1,CE,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
#NEX,MySampleID2,AA,Partial
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
0,0,0,0,0,0,1,1,1,1,1,1
#PCR,MySampleID2,EE,All
#PCR,MySampleID2,DF,Partial
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
1,1,1,1,1,1,1,1,1,1,1,1
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0\n";
close OUT;
exit;
}

if (!defined $ARGV[0] && !defined $opt{'P'}) {die $die2};

# Read in index file
open IN, $VAR{'SCI_index_file'};
while ($l = <IN>) {
	chomp $l;
	($ID,$pos,$seq) = split(/\t/, $l);
	($tier,$set,$side,$wells) = split(/_/, $ID);
	if ($tier =~ /(Tn5|Nex)/i) {
		if ($side =~ /i5/) {
			$TN5SET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$TN5SET_i7WELLS_seq{$set}{$wells} = $seq;
		}
		$tier = "Tn5";
	} else {
		if ($side =~ /i5/) {
			$PCRSET_i5WELLS_seq{$set}{$wells} = $seq;
		} else {
			$PCRSET_i7WELLS_seq{$set}{$wells} = $seq;
		}
	}
} close IN;

if (defined $opt{'O'}) {open OUT, ">$opt{'O'}"};

if (defined $opt{'P'}) {
	%NEX_ID_i5_i7_pair = ();
	%PCR_ID_i5_i7_pair = ();
	open IN, "$opt{'P'}";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^#/) {
			($class,$annot,$combo,$subset) = split(/,/, $l);;
			$class =~ s/^#//; $annot =~ s/ /_/g;
			if (!defined $ANNOT_flag{$annot}) {
				$ANNOT_flag{$annot} = $class;
			} else {
				$ANNOT_flag{$annot} .= ",$class";
			}
			($i5_set,$i7_set) = split(//, $combo);
			if ($subset =~ /all/i) {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$rowLetter = $LETTERS[$rowNum];
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
							$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
							$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
						}
					}
				}
			} else {
				for ($rowNum = 1; $rowNum <= 8; $rowNum++) {
					$row = <IN>; chomp $row; $rowLetter = $LETTERS[$rowNum];
					@ROW_COLS = split(/,/, $row); unshift @ROW_COLS, "0";
					for ($colNum = 1; $colNum <= 12; $colNum++) {
						if ($ROW_COLS[$colNum]>0) {
							if ($class =~ /(Tn5|Nex)/i) {
								$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$rowLetter},$TN5SET_i7WELLS_seq{$i7_set}{$colNum}";
								$NEX_ID_i5_i7_pair{$annot}{$pair} = 1;
							} else {
								$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$rowLetter},$PCRSET_i7WELLS_seq{$i7_set}{$colNum}";
								$PCR_ID_i5_i7_pair{$annot}{$pair} = 1;
							}
						}
					}
				}
			}
		}
	} close IN;
	
	foreach $annot (keys %ANNOT_flag) {
		if ($ANNOT_flag{$annot} =~ /(Tn5|Nex)/i && $ANNOT_flag{$annot} =~ /pcr/i) {
			print STDERR "Printing $annot index combinations.\n";
			foreach $NEX_pair (keys %{$NEX_ID_i5_i7_pair{$annot}}) {
				($ix3,$ix1) = split(/,/, $NEX_pair);
				foreach $PCR_pair (keys %{$PCR_ID_i5_i7_pair{$annot}}) {
					($ix4,$ix2) = split(/,/, $PCR_pair);
					if (defined $opt{'O'}) {
						print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
					} else {
						print "$ix1$ix2$ix3$ix4\t$annot\n";
					}
				}
			}
		} else {
			print STDERR "\nWARNING: Transposase (Nex/Tn5) AND PCR specifications must both be included for each annotation!\nBoth were not found for $annot! - SKIPPING!\n";
		}
	}
	
} else {
	foreach $annot_descriptor (@ARGV) {

		@DESCRIPTORS = split(/\+/, $annot_descriptor);
		$annot = shift(@DESCRIPTORS);
		%NEX_i5_i7_pair = ();
		%PCR_i5_i7_pair = ();
		
		print STDERR "\n##### PARSING ANNOT $annot #####\n";
		
		for ($i = 0; $i < @DESCRIPTORS; $i++) {
		
			print STDERR "\n##### Input Parse Descriptor: $DESCRIPTORS[$i] #####\n";
			($class,$columns) = split(/=/, $DESCRIPTORS[$i]);
			print STDERR "\tClass=$class, columns=$columns\n";
			($tier,$combo) = split(/,/, $class);
			print STDERR "\tTier=$tier, combo=$combo\n";
			($i5_set,$i7_set) = split(//, $combo);
			print STDERR "\ti5=$i5_set, i7=$i7_set\n";
			@COLUMN_DESCRIPTORS = split(/,/, $columns);
			for ($j = 0; $j < @COLUMN_DESCRIPTORS; $j++) {
				
				if ($COLUMN_DESCRIPTORS[$j] =~ /ALL/i) {$COLUMN_DESCRIPTORS[$j] = "1-12"};
				
				# determine if a subset of rows are specified (after :)
				if ($COLUMN_DESCRIPTORS[$j] =~ /:/) {
					($col,$rows) = split(/:/, $COLUMN_DESCRIPTORS[$j]);
				} else {
					$col = $COLUMN_DESCRIPTORS[$j];
					$rows = "A-H";
				}
				
				print STDERR "\tRows=$rows, which is:";
				# Add rows to the row list for the column(s)
				@ROW_LIST = split(//, $rows);
				@ROW_INCLUDE = ();
				for ($k = 0; $k < @ROW_LIST; $k++) {
					if ($ROW_LIST[($k+1)] eq "-") {
						$rowStart = $LETTER_NUM{$ROW_LIST[$k]};
						$k++; $k++;
						$rowEnd = $LETTER_NUM{$ROW_LIST[$k]};
						for ($l = $rowStart; $l <= $rowEnd; $l++) {
							push @ROW_INCLUDE, $LETTERS[$l];
							print STDERR " $LETTERS[$l]";
						}
					} else {
						push @ROW_INCLUDE, $ROW_LIST[$k];
						print STDERR " $ROW_LIST[$k]";
					}
				}
				
				# Go through the column(s)
				print STDERR "\n\tCols=$col, which is:";
				@COL_INCLUDE = ();
				if ($col =~ /-/) {
					($startCol,$endCol) = split(/-/, $col);
					for ($l = $startCol; $l <= $endCol; $l++) {
						push @COL_INCLUDE, $l;
						print STDERR " $l";
					}
				} else {
					push @COL_INCLUDE, $col;
					print STDERR " $col";
				}
				
				print STDERR "\n\t  --> Adding all pairs\n";
				# Now add all relevant barcodes to their respective groupings
				foreach $add_col (@COL_INCLUDE) {
					foreach $add_row (@ROW_INCLUDE) {
						if ($class =~ /(Tn5|Nex)/i) {
							$pair = "$TN5SET_i5WELLS_seq{$i5_set}{$add_row},$TN5SET_i7WELLS_seq{$i7_set}{$add_col}";
							$NEX_i5_i7_pair{$pair} = 1;
						} else {
							$pair = "$PCRSET_i5WELLS_seq{$i5_set}{$add_row},$PCRSET_i7WELLS_seq{$i7_set}{$add_col}";
							$PCR_i5_i7_pair{$pair} = 1;
						}
					}
				}
			}
		}

		foreach $NEX_pair (keys %NEX_i5_i7_pair) {
			($ix3,$ix1) = split(/,/, $NEX_pair);
			foreach $PCR_pair (keys %PCR_i5_i7_pair) {
				($ix4,$ix2) = split(/,/, $PCR_pair);
				if (defined $opt{'O'}) {
					print OUT "$ix1$ix2$ix3$ix4\t$annot\n";
				} else {
					print "$ix1$ix2$ix3$ix4\t$annot\n";
				}
			}
		}
	}
}

if (defined $opt{'O'}) {close OUT};

}
1;
