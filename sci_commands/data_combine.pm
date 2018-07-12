package sci_commands::data_combine;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("data_combine");

sub data_combine {

@ARGV = @_;
# Defaults
$maxDim = 15;

getopts("O:D:A:R:zF:", \%opt);

$die2 = "
scitools combine-data [output_file] [data1_name]=[data1_file] [data2_name]=[data2_file] etc...
   or    data-combine

Will report all data from files into the specified output file. If the data name is not
specified, the file name will be used instead.

Will auto detect the following filetypes from names:
   annot, dims, matrix, lambda, values, complexity
If it is a different file type it will include as class 'other', to force specifying 'other'
add another argument after the file that is the specified class - will force 'other'
   [data_name]=[data_file]=[other_type]
This can be used for other types of matrix that do not have cellID columns.
   e.g.   Mean_DevZ=myMeanDeviationByCluster.matrix=matrix
Which will output with a suffix '.matrix', but will be of class 'other'

For cells without information in a file, NA will be reported.

Options:
   -D   [INT]   Maximum dimension to include for dims files (def = $maxDim)
   -F   [STR]   File specifying files to combine - will read each line as an input argument.
   -A   [STR]   Only include cells within specified annot file.
   -a   [STR]   Include only specified annotations (in -A), comma sep.
   -R   [STR]   Rename cells using rename.annot file
                (if absent will create as the provided file name)
   -z           Gzip the output

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if ($ARGV[0] =~ /[,=]/) {die "ERROR: First argument must be the output file, it looks like the first argument specified is an input file.\n$die2"};
if (defined $opt{'D'}) {$maxDim = $opt{'D'}};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};

if (defined $opt{'A'}) {read_annot($opt{'A'})};

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (keys %ANNOT_count) {
		$ANNOT_include{$annot} = 1;
	}
}

if (defined $opt{'F'}) {
	open IN, "$opt{'F'}";
	while ($l = <IN>) {
		chomp $l;
		push @ARGV, $l;
	} close IN;
}

# get full file type info and cells present
$included_cell_ct = 0; $cellID_out = ""; @INCLUDED_CELLIDS = ();
for ($i = 1; $i < @ARGV; $i++) {
	if ($ARGV[$i] =~ /[,=]/) {
		$type = "";
		($name,$file,$type) = split(/[,=]/, $ARGV[$i]);
		if ($type ne "") {
			$CLASS{$i} = "other";
			$TYPE{$i} = $type;
		}
	} else {
		$name = $ARGV[$i]; $file = $ARGV[$i];
	}
	if (-e "$file") {
		if ($file =~ /\.annot$/ && !defined $CLASS{$i}) {$CLASS{$i} = "annot"}
		elsif ($file =~ /\.dims$/ && !defined $CLASS{$i}) {$CLASS{$i} = "dims"}
		elsif ($file =~ /\.complexity\.txt$/ && !defined $CLASS{$i}) {$CLASS{$i} = "complexity"}
		elsif ($file =~ /\.lambda$/ && !defined $CLASS{$i}) {$CLASS{$i} = "lambda"}
		elsif ($file =~ /\.values$/ && !defined $CLASS{$i}) {$CLASS{$i} = "values"}
		elsif ($file =~ /\.matrix$/ && !defined $CLASS{$i}) {$CLASS{$i} = "matrix"}
		else {
			$CLASS{$i} = "other";
		}
		if (!defined $SKIP{$i}) {
			$FILES{$i} = $file;
			$NAMES{$i} = $name;
			if ($CLASS{$i} =~ /(annot|dims|lambda|values|compelxity)/) {
				open IN, "$file";
				while ($l = <IN>) {
					chomp $l;
					@P = split(/\t/, $l);
					$cellID = $P[0];
					if (!defined $opt{'A'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) {
						if (!defined $CELLID_include{$cellID}) {
							$CELLID_include{$cellID} = 1;
							$included_cell_ct++; $cellID_out .= "$cellID\t";
							push @INCLUDED_CELLIDS, $cellID;
						} else {
							$CELLID_include{$cellID}++;
						}
					}
				}
				close IN;
			} else {
				open IN, "$file"; $h = <IN>; close IN; chomp $h;
				@H = split(/\t/, $h);
				foreach $cellID (@H) {
					if (!defined $opt{'A'} || defined $ANNOT_include{$CELLID_annot{$cellID}}) {
						if (!defined $CELLID_include{$cellID}) {
							$CELLID_include{$cellID} = 1;
							$included_cell_ct++; $cellID_out .= "$cellID\t";
						push @INCLUDED_CELLIDS, $cellID;
						} else {
							$CELLID_include{$cellID}++;
						}
					}
				}
			}
		}
	} else {
		print STDERR "\nWARNING: Cannot find file: $file, skipping!\n";
	}
}
if ($included_cell_ct==0) {
	die "\nERROR: Included cell count is 0! Check your file cellID compatability!\n";
} elsif ($included_cell_ct<100) {
	print STDERR "\nWARNING: Included cell count is < 100, if this is less than expected, re-check your files and cellIDs!\n";
}

# rename cellIDs
%CELLID_cellName = ();
if (defined $opt{'R'}) {
	if (-e "$opt{'R'}") {
		open RENAME, "$opt{'R'}";
		while ($l = <RENAME>) {
			chomp $l;
			($cellID,$cellName) = split(/\t/, $l);
			$CELLID_cellName{$cellID} = $cellName;
		} close RENAME;
	} else {
		$opt{'R'} =~ s/\.annot$//; $opt{'R'} =~ s/\.rename$//;
		$cellNum = 1;
		open RENAME, ">$opt{'R'}.rename.annot";
		for ($i = 0; $i < @INCLUDED_CELLIDS; $i++) {
			$cellName = "CellID_".sprintf("%09s", $cellNum);
			$cellNum++;
			$cellID = $INCLUDED_CELLIDS[$i];
			$CELLID_cellName{$cellID} = $cellName;
			print RENAME "$cellID\t$cellName\n";
		} close RENAME;
	}
	$cellID_out = "";
	for ($i = 0; $i < @INCLUDED_CELLIDS; $i++) {
		$cellID = $INCLUDED_CELLIDS[$i];
		if (defined $CELLID_cellName{$cellID}) {
			$cellID_out .= "$CELLID_cellName{$cellID}\t";
		} else {
			die "\nERROR: Renaming file $opt{'R'} does not contain all cellIDs! Cannot proceed!\n";
		}
	}
}

$cellID_out =~ s/\t$//;

# open OUT and print cellIDs and some basic info
if (!defined $opt{'z'}) {
	open OUT, ">$ARGV[0]";
} else {
	$ARGV[0] =~ s/\.gz$//;
	open OUT, "| $gzip > $ARGV[0].gz";
}
if (defined $opt{'A'}) {
	print OUT "#CELL_IDS\tCELL_INCLUSION_ANNOTATION_FILE=$opt{'A'}";
	if (defined $opt{'a'}) {
		print OUT "\tINCLUDED_ANNOTATIONS=$opt{'a'}";
	}
} else {
	print OUT "#CELL_IDS";
}

for ($i = 1; $i < @ARGV; $i++) {
	print OUT "\tDATA=$ARGV[$i]";
} print OUT "\n";

print OUT "CELL_IDS\tUniqueIdentifier\t$cellID_out\n";
for ($i = 1; $i < @ARGV; $i++) {
	if (!defined $SKIP{$i}) {
		if ($CLASS{$i} eq "annot") {
			print OUT "#ANNOTATION_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_annot($FILES{$i});
			print OUT "$NAMES{$i}\tANNOTATION";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_annot{$cellID}) {
					print OUT "\t$CELLID_annot{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "values") {
			print OUT "#VALUES_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_values($FILES{$i});
			print OUT "$NAMES{$i}\tVALUES";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_value{$cellID}) {
					print OUT "\t$CELLID_value{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "lambda") {
			print OUT "#LAMBDA_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_values($FILES{$i});
			print OUT "$NAMES{$i}\tLAMBDA";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_value{$cellID}) {
					print OUT "\t$CELLID_value{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "complexity") {
			print OUT "#COMPLEXITY_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_complexity($FILES{$i});
			print OUT "$NAMES{$i}\tRAW_READS";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_raw_reads{$cellID}) {
					print OUT "\t$CELLID_raw_reads{$cellID}";
				} else {
					print OUT "\tNA";
				}
			}
			print OUT "\n$NAMES{$i}\tUNIQUE_READS";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_uniq_reads{$cellID}) {
					print OUT "\t$CELLID_uniq_reads{$cellID}";
				} else {
					print OUT "\tNA";
				}
			}
			print OUT "\n$NAMES{$i}\tCOMPLEXITY";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_complexity{$cellID}) {
					print OUT "\t$CELLID_complexity{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n$NAMES{$i}\nRANK";
			for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
				$cellID = $INCLUDED_CELLIDS[$j];
				if (defined $CELLID_complexity_rank{$cellID}) {
					print OUT "\t$CELLID_complexity_rank{$cellID}";
				} else {
					print OUT "\tNA";
				}
			} print OUT "\n";
		} elsif ($CLASS{$i} eq "dims") {
			print OUT "#DIMENSIONS_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			read_dims($FILES{$i});
			if ($Ndims<$maxDim) {$lastDim = $Ndims} else {$lastDim = $maxDim};
			for ($dim = 1; $dim < $lastDim; $dim++) {
				print OUT "$NAMES{$i}\tDIMENSION_$dim";
				for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
					$cellID = $INCLUDED_CELLIDS[$j];
					if (defined $CELLID_DIMS{$cellID}[$dim]) {
						print OUT "\t$CELLID_DIMS{$cellID}[$dim]";
					} else {
						print OUT "\tNA";
					}
				} print OUT "\n";
			}
		} elsif ($CLASS{$i} eq "matrix") {
			print OUT "#MATRIX_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\n";
			open MAT, "$FILES{$i}";
			$mat_head = <MAT>; chomp $mat_head; @MATRIX_CELLS = split(/\t/, $mat_head);
			%MATRIX_CELLID_pos = ();
			for ($cellPos = 0; $cellPos < @MATRIX_CELLS; $cellPos++) {
				$MATRIX_CELLID_pos{$MATRIX_CELLS[$cellPos]} = $cellPos;
			}
			while ($mat_l = <MAT>) {
				chomp $mat_l;
				@MAT_FIELDS = split(/\t/, $mat_l);
				$feature = shift(@MAT_FIELDS);
				print OUT "$NAMES{$i}\t$feature";
				for ($j = 0; $j < @INCLUDED_CELLIDS; $j++) {
					$cellID = $INCLUDED_CELLIDS[$j];
					if (defined $MATRIX_CELLID_pos{$cellID}) {
						$value = $MAT_FIELDS[$MATRIX_CELLID_pos{$cellID}];
					} else {
						$value = "NA";
					}
					print OUT "\t$value";
				}
				print OUT "\n";
			}
			close MAT;
		} elsif ($CLASS{$i} eq "other") {
			if (!defined $TYPE{$i}) {$TYPE{$i} = "txt"};
			print OUT "#OTHER_DATA\tNAME=$NAMES{$i}\tFILE=$FILES{$i}\tTYPE=$TYPE{$i}\n";
			open IN, "$FILES{$i}";
			while ($in_l = <IN>) {
				chomp $in_l;
				print OUT "$NAMES{$i}\t$in_l\n";
			} close IN;
		}
	}
}
close OUT;

}
1;
