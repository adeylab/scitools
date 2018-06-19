package sci_utils::general;

use Getopt::Std; %opt = ();
use Exporter "import";

# Export all of the subroutines and variables that should be passed to scitools
@EXPORT = (
	"load_defaults",
		qw($color_mapping),qw($ref_shortcuts),qw(@BASES),qw(%REF),qw(%VAR),qw($gzip),qw($zcat),
		qw($bwa),qw($samtools),qw($scitools),qw($macs2),qw($bedtools),qw($Rscript),qw($Pscript),
	"read_annot",
		qw(%CELLID_annot),qw(%ANNOT_count),qw($annot_count),qw(@ANNOT_FILES),
	"read_complexity",
		qw(%CELLID_uniq_reads),qw(%CELLID_raw_reads),qw(%CELLID_complexity),qw(%CELLID_complexity_rank),
	"read_matrix",
		qw(%CELLID_FEATURE_value),
	"read_matrix_stats",
		qw(@MATRIX_COLNAMES),qw(@MATRIX_ROWNAMES),qw(%MATRIX_CellID_nonZero),qw(%MATRIX_feature_nonZero),
		qw(%MATRIX_CellID_signal),qw(%MATRIX_feature_signal),qw($matrix_colNum),qw($matrix_rowNum),
	"read_color_string","read_color_file",
		qw(%ANNOT_color),
	"read_dims",
		qw(%CELLID_DIMS),qw($Ndims),
	"read_pcurve_dims",
		qw(%CELLID_PCURVE_DIMS),
	"read_values",
		qw(%CELLID_value),qw(@VALUES),qw($value_min),qw($value_max),qw($value_mean),qw($value_sum),qw($value_median),
	"read_ranges",
		qw(@RANGE_VALUES),qw($range_R_set),
	"read_mode",
		qw(%MODE_GROUP_CLASS_PARTS),qw(%ALIAS_mode),
	"read_indexes",
		qw(%INDEX_POS_SEQ_id),qw(%INDEX_POS_SEQ_well),
	"read_indexdir",
		qw(%INDEX_TYPE_SEQ_seq),qw(%INDEX_TYPE_SEQ_id),qw(%INDEX_TYPE_length),qw($indexes_loaded),
	"read_refgene",
		qw(%GENENAME_coords),qw(%GENEID_coords),qw(%GENECOORDS_refGene),
	"get_gradient",
	"load_gradient_defaults",
		qw(%COLOR_GRADIENT),
	"print_gradients"
);


sub load_defaults {
	# Some global ones that are not configured
	$color_mapping = "none";
	$ref_shortcuts = "";
	@BASES = ("A", "C", "G", "T", "N");
	%REF; %VAR;
	if (-e "$_[0]") {
		open DEF, "$_[0]";
		while ($def = <DEF>) {
			if ($def !~ /^#/) {
				chomp $def;
				($var,$val) = split(/=/, $def);
				if ($var =~ /_ref$/) {
					$refID = $var; $refID =~ s/_ref$//;
					$REF{$refID} = $val;
					$ref_shortcuts .= "\t$refID = $val\n";
				} else {
					if ($var eq "gzip") {$gzip = $val}
					elsif ($var eq "zcat") {$zcat = $val}
					elsif ($var eq "bwa") {$bwa = $val}
					elsif ($var eq "samtools") {$samtools = $val}
					elsif ($var eq "scitools") {$scitools = $val}
					elsif ($var eq "macs2") {$macs2 = $val}
					elsif ($var eq "bedtools") {$bedtools = $val}
					elsif ($var eq "Rscript") {$Rscript = $val}
					elsif ($var eq "Pscript") {$Pscript = $val}
					else {$VAR{$var} = $val};
				}
			}
		} close DEF;
	} else {
		# No config file = load default executables for convenience
		$gzip = "gzip"; #DEFAULT=gzip
		$zcat = "zcat"; #DEFAULT=zcat
		$bwa = "bwa"; #DEFAULT=bwa
		$samtools = "samtools"; #DEFAULT=samtools
		$scitools = "scitools"; #DEFAULT=scitools
		$macs2 = "macs2"; #DEFAULT=macs2
		$bedtools = "bedtools"; #DEFAULT=bedtools
		$Rscript = "Rscript"; #DEFAULT=Rscript
		$Pscript = "python"; #DEFAULT=Pscript
	}
}

sub read_annot {
	%CELLID_annot = (); %ANNOT_count = (); $annot_count = 0;
	@ANNOT_FILES = split(/,/, $_[0]);
	foreach $annot_file (@ANNOT_FILES) {
		open ANNOT, "$annot_file";
		while ($annot_line = <ANNOT>) {
			chomp $annot_line;
			($annot_cellID,$annot) = split(/\t/, $annot_line);
			$CELLID_annot{$annot_cellID} = $annot;
			if (!defined $ANNOT_count{$annot}) {
				$annot_count++;
				$ANNOT_count{$annot}=0;
			}
		} close ANNOT;
	}
}

sub read_complexity {
	%CELLID_uniq_reads = ();
	%CELLID_raw_reads = ();
	%CELLID_complexity = ();
	%CELLID_complexity_rank = ();
	@COMPLEXITY_FILES = split(/,/, $_[0]);
	foreach $complexity_file (@COMPLEXITY_FILES) {
		open COMPL, "$complexity_file";
		while ($comp_line = <COMPL>) {
			chomp $comp_line;
			($num,$cellID,$raw,$uniq,$pct) = split(/\t/, $comp_line);
			$CELLID_complexity_rank{$cellID} = $num;
			$CELLID_uniq_reads{$cellID} = $uniq;
			$CELLID_raw_reads{$cellID} = $raw;
			$CELLID_complexity{$cellID} = $pct;
		} close COMPL;
	}
}

sub read_matrix {
	%CELLID_FEATURE_value = (); @MATRIX_COLNAMES = (); @MATRIX_ROWNAMES = ();
	%MATRIX_CellID_nonZero = (); %MATRIX_feature_nonZero = ();
	%MATRIX_CellID_signal = (); %MATRIX_feature_signal = ();
	$matrix_colNum = 0; $matrix_rowNum = 0;
	open MAT, "$_[0]";
	$column_line = <MAT>; chomp $column_line; @MATRIX_COLNAMES = split(/\t/, $column_line);
	foreach $column (@MATRIX_COLNAMES) {$matrix_colNum++};
	while ($matrix_line = <MAT>) {
		chomp $matrix_line;
		@MATRIX_ROW = split(/\t/, $matrix_line);
		$featureName = shift(@MATRIX_ROW);
		push @MATRIX_ROWNAMES, $featureName;
		for ($colNum = 0; $colNum < @MATRIX_ROW; $colNum++) {
			$CELLID_FEATURE_value{$MATRIX_COLNAMES[$colNum]}{$featureName} = $MATRIX_ROW[$colNum];
			if ($MATRIX_ROW[$colNum]>0) {
				$MATRIX_CellID_nonZero{$MATRIX_COLNAMES[$colNum]}++;
				$MATRIX_feature_nonZero{$featureName}++;
				$MATRIX_CellID_signal{$MATRIX_COLNAMES[$colNum]}+=$MATRIX_ROW[$colNum];
				$MATRIX_feature_signal{$featureName}+=$MATRIX_ROW[$colNum];
			}
		}
		$matrix_rowNum++;
	} close MAT;
}

sub read_matrix_stats {
	@MATRIX_COLNAMES = (); @MATRIX_ROWNAMES = ();
	%MATRIX_CellID_nonZero = (); %MATRIX_feature_nonZero = ();
	%MATRIX_CellID_signal = (); %MATRIX_feature_signal = ();
	$matrix_colNum = 0; $matrix_rowNum = 0;
	open MAT, "$_[0]";
	$column_line = <MAT>; chomp $column_line; @MATRIX_COLNAMES = split(/\t/, $column_line);
	foreach $column (@MATRIX_COLNAMES) {$matrix_colNum++};
	while ($matrix_line = <MAT>) {
		chomp $matrix_line;
		@MATRIX_ROW = split(/\t/, $matrix_line);
		$featureName = shift(@MATRIX_ROW);
		push @MATRIX_ROWNAMES, $featureName;
		for ($colNum = 0; $colNum < @MATRIX_ROW; $colNum++) {
			if ($MATRIX_ROW[$colNum]>0) {
				$MATRIX_CellID_nonZero{$MATRIX_COLNAMES[$colNum]}++;
				$MATRIX_feature_nonZero{$featureName}++;
				$MATRIX_CellID_signal{$MATRIX_COLNAMES[$colNum]}+=$MATRIX_ROW[$colNum];
				$MATRIX_feature_signal{$featureName}+=$MATRIX_ROW[$colNum];
			}
		}
		$matrix_rowNum++;
	} close MAT;
}

sub read_color_string {
	%ANNOT_color = ();
	@COL_STRING = split(/,/, $_[0]);
	$color_mapping = "\"Cell\" = \"lightsteelblue4\",";
	foreach $color_assignment (@COL_STRING) {
		($annot,$color) = split(/=/, $color_assignment);
		$ANNOT_color{$annot} = $color;
		$color_mapping .= "\"$annot\" = \"$color\",";
	} $color_mapping =~ s/,$//;
}

sub read_color_file {
	%ANNOT_color = ();
	open COL_FILE, "$_[0]";
	$color_mapping = "\"Cell\" = \"lightsteelblue4\",";
	while ($color_assignment = <COL_FILE>) {
		chomp $color_assignment;
		($annot,$color) = split(/\t/, $color_assignment);
		$ANNOT_color{$annot} = $color;
		$color_mapping .= "\"$annot\" = \"$color\",";
	} $color_mapping =~ s/,$//; close COL_FILE;
}

sub read_dims {
	%CELLID_DIMS = ();
	@DIM_FILES = split(/,/, $_[0]);
	foreach $dim_file (@DIM_FILES) {
		open DIMS, "$dim_file";
		while ($dim_line = <DIMS>) {
			chomp $dim_line;
			@DIM_SET = split(/\t/, $dim_line);
			$cellID = $DIM_SET[0];
			@{$CELLID_DIMS{$cellID}} = @DIM_SET;
			$Ndims = @DIM_SET;
		} close DIMS;
	}
}

sub read_pcurve_dims {
	%CELLID_PCURVE_DIMS = ();
	@PC_DIM_FILES = split(/,/, $_[0]);
	foreach $pc_dim_file (@PC_DIM_FILES) {
		open DIMS, "$pc_dim_file";
		while ($dim_line = <DIMS>) {
			chomp $dim_line;
			@DIM_SET = split(/\t/, $dim_line);
			$cellID = $DIM_SET[0];
			@{$CELLID_PCURVE_DIMS{$cellID}} = @DIM_SET;
		} close DIMS;
	}
}

sub read_values {
	%CELLID_value = ();
	@VALUES = ();
	$value_min = "na"; $value_max = "na";
	$value_mean = 0; $value_sum = 0; $value_median = 0;
	@VAL_FILES = split(/,/, $_[0]);
	foreach $val_file (@VAL_FILES) {
		open VALS, "$val_file";
		while ($val_line = <VALS>) {
			chomp $val_line;
			($cellID,$value) = split(/\t/, $val_line);
			$CELLID_value{$cellID} = $value;
			push @VALUES, $value;
			$value_sum += $value;
			if ($value > $value_max || $value_max eq "na") {$value_max = $value};
			if ($value < $value_min || $value_min eq "na") {$value_min = $value};
		} close VALS;
	}
	$value_range = $value_max - $value_min;
	$value_mean = $value_sum/@VALUES;
	@SORTED_VALUES = sort {$a<=>$b} @VALUES;
	$value_median = $SORTED_VALUES[int(@VALUES/2)];
}

sub read_ranges {
	@RANGE_VALUES = ();
	$range_R_set = "c(";
	@RANGE_SETS = split(/,/, $_[0]);
	for ($range_pos = 0; $range_pos < @RANGE_SETS; $range_pos++) {
		$range_set = $RANGE_SETS[$range_pos];
		if ($range_set =~ /-/) {
			($range_start,$range_end) = split(/-/, $range_set);
			for ($range_value = $range_start; $range_value <= $range_end; $range_value++) {
				push @RANGE_VALUES, $range_value;
				$range_R_set .= "$range_value,";
			}
		} else {
			push @RANGE_VALUES, $range_set;
			$range_R_set .= "$range_set,";
		}
	}
	$range_R_set =~ s/,$//;
	$range_R_set .= ")";
}

sub read_indexdir {
	%INDEX_TYPE_SEQ_seq = ();
	%INDEX_TYPE_SEQ_id = ();
	%INDEX_TYPE_length = ();
	$indexes_loaded = 0;
	@INDEX_FILES = split(/,/, $_[0]);
	foreach $index_specification (@INDEX_FILES) {
		($index_name,$target) = split(/=/, $index_specification);
		if ($index_name =~ /DIR/i) {
			$target =~ s/\/$//;
			open INDEX, "cat $target/*.txt |";
		} else {
			open INDEX, "$target";
		}
		while ($index_line = <INDEX>) {
			chomp $index_line;
			($index_id,$index_type_raw,$index_seq_raw) = split(/\t/, $index_line);
			$index_type = lc($index_type_raw);
			$index_seq = uc($index_seq_raw);
			($id_tier,$id_set,$id_side,$id_well) = split(/_/, $index_id);
			$INDEX_TYPE_SEQ_seq{$index_type}{$index_seq} = $index_seq;
			$INDEX_TYPE_length{$index_type} = length($index_seq);
			$INDEX_TYPE_SEQ_id{$index_type}{$index_seq} = $index_id;
			$indexes_loaded++;
		} close INDEX;
	}
}

sub read_mode { # FOR MODE IMPLEMENTATION
	$sci_modes_file = $_[0];
	%MODE_GROUP_CLASS_PARTS = ();
	%ALIAS_mode = ();
	open MODES, "$sci_modes_file" || die "ERROR: Cannot open sci_modes file: $sci_modes_file.\n";
	while ($mode_l = <MODES>) {
		if ($mode_l !~ /^#/) {
			chomp $mode_l;
			if ($mode_l =~ /^>/) {
				$mode_l =~ s/^>//;
				@MODE_ALIASES = split(/,/, $mode_l);
				$mode_name = $MODE_ALIASES[0];
				for ($aliasID = 0; $aliasID < @MODE_ALIASES; $aliasID++) {
					$ALIAS_mode{$MODE_ALIASES[$aliasID]} = $mode_name;
					#print STDERR "Adding alias $MODE_ALIASES[$aliasID] to $mode_name\n";
				}
				$mode_group = 0;
			} elsif ($mode_l =~ /\&/) {
				$mode_group++;
				#print STDERR " >>> Jumping to a second mode group <<<\n";
			} else {
				($item,$spec) = split(/=/, $mode_l);
				#print STDERR " item = $item, spec = $spec\n";
				if ($item =~ /name/) {
					$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'name'} = $spec;
					#print STDERR "   $item is a name, adding to names.\n";
				} else {
					@READ_PARTS = split(/,/, $spec);
					$current_offset = 0;
					#print STDERR "   Parsing read parts:\n";
					for ($read_partID = 0; $read_partID < @READ_PARTS; $read_partID++) {
						($read_item,$item_length) = split(/:/, $READ_PARTS[$read_partID]);
						#print STDERR "      Part = $read_item, length = $item_length, current offset = $current_offset\n";
						if ($read_item =~ /null/) {
							$current_offset += $item_length; # do not load anything for null bases, just add offset
						} else {
							if ($read_item =~ /^read_/) {
								$read_name = $read_item;
								$read_name =~ s/^read_//;
								push @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{'read_outputs'}}, $read_name;
								#print STDERR "         Part is a read - adding to read outputs - named: $read_name\n";
							}
							push @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'name'}}, $read_item;
							push @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'length'}}, $item_length;
							push @{$MODE_GROUP_CLASS_PARTS{$mode_name}[$mode_group]{$item}{'offset'}}, $current_offset;
							if ($item_length =~ /^\d+$/) {$current_offset += $item_length};
						}
					}
				}
			}
		}
	} close MODES;
}

sub read_indexes {
	%INDEX_POS_SEQ_id = ();
	%INDEX_POS_SEQ_well = ();
	open INDEX, "$_[0]";
	while ($index_line = <INDEX>) {
		chomp $index_line;
		($index_id,$index_pos,$index_seq) = split(/\t/, $index_line);
		($id_tier,$id_set,$id_side,$id_well) = split(/_/, $index_id);
		$INDEX_POS_SEQ_id{$index_pos}{$index_seq} = $id_set;
		$INDEX_POS_SEQ_well{$index_pos}{$index_seq} = $id_well;
	} close INDEX;
}

sub read_refgene {
	%GENENAME_coords = ();
	%GENEID_coords = ();
	%GENECOORDS_refGene = ();
	open REFGENE, "$_[0]";
	while ($refgene_line = <REFGENE>) {
		chomp $refgene_line;
		@REFGENE = split(/\t/, $refgene_line);
		$gene_coords = "$REFGENE[2]:$REFGENE[4]-$REFGENE[5]";
		$GENEID_coords{$REFGENE[1]} = $gene_coords;
		$GENENAME_coords{$REFGENE[12]} = $gene_coords;
		$GENECOORDS_refGene{$gene_coords} = $refgene_line;
	} close REFGENE;
}

sub get_gradient {
	if ($_[0] =~ /,/) {
		@GRADIENT_COLORS = split(/,/, $_[0]);
		$gradient_specification = "gradient_funct<-colorRampPalette(c(";
		for ($grad_pos = 0; $grad_pos < @GRADIENT_COLORS; $grad_pos++) {
			$gradient_specification .= "\"$GRADIENT_COLORS[$grad_pos]\",";
		} $gradient_specification =~ s/,$//;
		$gradient_specification .= "))";
	} elsif (defined $COLOR_GRADIENT{$_[0]}) {
		$gradient_specification = $COLOR_GRADIENT{$_[0]};
	} else {
		die "
ERROR: Color gradients must either be one of the pre-set scitools color gradients,
       or defined as a comma separated list of two or more colors. For more information
       run 'scitools gradient'.
";
	}
	return $gradient_specification;
}

sub load_gradient_defaults {
# COLOR GRADIENT DEFAULTS:
%COLOR_GRADIENT = ();

# diverging
$COLOR_GRADIENT{'PuOr'} = "gradient_funct<-colorRampPalette(c(\"#542788\",\"#8073ac\",\"#b2abd2\",\"#d8daeb\",\"#f7f7f7\",\"#fee0b6\",\"#fdb863\",\"#e08214\",\"#b35806\"))";
$COLOR_GRADIENT{'BuRd'} = "gradient_funct<-colorRampPalette(c(\"#2166ac\",\"#4393c3\",\"#92c5de\",\"#d1e5f0\",\"#f7f7f7\",\"#fddbc7\",\"#f4a582\",\"#d6604d\",\"#b2182b\"))";
$COLOR_GRADIENT{'PuRd'} = "gradient_funct<-colorRampPalette(c(\"#542788\",\"#8073ac\",\"#b2abd2\",\"#d8daeb\",\"#f7f7f7\",\"#fddbc7\",\"#f4a582\",\"#d6604d\",\"#b2182b\"))";
$COLOR_GRADIENT{'RdYlGn'} = "gradient_funct<-colorRampPalette(c(\"#d73027\",\"#f46d43\",\"#fdae61\",\"#fee08b\",\"#ffffbf\",\"#d9ef8b\",\"#a6d96a\",\"#66bd63\",\"#1a9850\"))";
$COLOR_GRADIENT{'BrBG'} = "gradient_funct<-colorRampPalette(c(\"#8c510a\",\"#bf812d\",\"#dfc27d\",\"#f6e8c3\",\"#f5f5f5\",\"#c7eae5\",\"#80cdc1\",\"#35978f\",\"#01665e\"))";
$COLOR_GRADIENT{'PiYG'} = "gradient_funct<-colorRampPalette(c(\"#c51b7d\",\"#de77ae\",\"#f1b6da\",\"#fde0ef\",\"#f7f7f7\",\"#e6f5d0\",\"#b8e186\",\"#7fbc41\",\"#4d9221\"))";
$COLOR_GRADIENT{'PuRdGn'} = "gradient_funct<-colorRampPalette(c(\"#762a83\",\"#9970ab\",\"#c2a5cf\",\"#e7d4e8\",\"#f7f7f7\",\"#d9f0d3\",\"#a6dba0\",\"#5aae61\",\"#1b7837\"))";
$COLOR_GRADIENT{'BuYlRd'} = "gradient_funct<-colorRampPalette(c(\"#2c7bb6\",\"#abd9e9\",\"#ffffbf\",\"#fdae61\",\"#d7191c\"))";

# sequential multi-hue
$COLOR_GRADIENT{'YlOrRd'} = "gradient_funct<-colorRampPalette(c(\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#bd0026\",\"#800026\"))";
$COLOR_GRADIENT{'WtYlOrRd'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#ffffcc\",\"#ffeda0\",\"#fed976\",\"#feb24c\",\"#fd8d3c\",\"#fc4e2a\",\"#e31a1c\",\"#bd0026\",\"#800026\"))";
$COLOR_GRADIENT{'RdPu'} = "gradient_funct<-colorRampPalette(c(\"#fff7f3\",\"#fde0dd\",\"#fcc5c0\",\"#fa9fb5\",\"#f768a1\",\"#dd3497\",\"#ae017e\",\"#7a0177\",\"#49006a\"))";
$COLOR_GRADIENT{'YlGnBu'} = "gradient_funct<-colorRampPalette(c(\"#ffffd9\",\"#edf8b1\",\"#c7e9b4\",\"#7fcdbb\",\"#41b6c4\",\"#1d91c0\",\"#225ea8\",\"#253494\",\"#081d58\"))";
$COLOR_GRADIENT{'BuGnYl'} = "gradient_funct<-colorRampPalette(c(\"#081d58\",\"#253494\",\"#225ea8\",\"#1d91c0\",\"#41b6c4\",\"#7fcdbb\",\"#c7e9b4\",\"#edf8b1\",\"#ffffd9\"))";
$COLOR_GRADIENT{'YlGn'} = "gradient_funct<-colorRampPalette(c(\"#ffffe5\",\"#f7fcb9\",\"#d9f0a3\",\"#addd8e\",\"#78c679\",\"#41ab5d\",\"#238443\",\"#006837\",\"#004529\"))";
$COLOR_GRADIENT{'Spect'} = "gradient_funct<-colorRampPalette(c(\"#3288bd\",\"#66c2a5\",\"#abdda4\",\"#e6f598\",\"#ffffbf\",\"#fee08b\",\"#fdae61\",\"#f46d43\",\"#d53e4f\"))";
$COLOR_GRADIENT{'WtSpect'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#3288bd\",\"#66c2a5\",\"#abdda4\",\"#e6f598\",\"#ffffbf\",\"#fee08b\",\"#fdae61\",\"#f46d43\",\"#d53e4f\"))";


# sequential single-hue
$COLOR_GRADIENT{'WtPu'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#fcfbfd\",\"#efedf5\",\"#dadaeb\",\"#bcbddc\",\"#9e9ac8\",\"#807dba\",\"#6a51a3\",\"#54278f\",\"#3f007d\"))";
$COLOR_GRADIENT{'WtRd'} = "gradient_funct<-colorRampPalette(c(\"white\",\"#fff5f0\",\"#fee0d2\",\"#fcbba1\",\"#fc9272\",\"#fb6a4a\",\"#ef3b2c\",\"#cb181d\",\"#a50f15\",\"#67000d\"))";

# non colorbrewer
$COLOR_GRADIENT{'BuGo'} = "gradient_funct<-colorRampPalette(c(\"lightblue\",\"blue\",\"red\",\"orange\",\"gold\"))";
$COLOR_GRADIENT{'BuG90Rd'} = "gradient_funct<-colorRampPalette(c(\"blue\",\"gray90\",\"red\"))";

}

sub print_gradients {
print "
scitools available gradients and gradient specification. For a number 
of plotting options, a gradient can be specified. Scitools has several
default gradients shown below, or a custom gradient can be specified.
These commands will add in a color ramp function in R.

Default Gradients: (name, R color function)
";

foreach $gradient (sort keys %COLOR_GRADIENT) {
	print "   $gradient\t($COLOR_GRADIENT{$gradient})\n";
}

print "
To specify a custom gradient, specify a comma-separated list of R
colors, or hex codes:
   eg.   red,yellow,blue   would be a gradient of those three colors.
   
";
}

1;
