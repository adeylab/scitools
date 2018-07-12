package sci_commands::check_vars;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("check_vars");

sub check_vars {

@ARGV = @_;
if (!defined $ARGV[0]) {print "DEBUG: No variables were passed to this command!\n"};
for ($i = 0; $i < @ARGV; $i++) {
	print "DEBUG: Var ARGV $i passed is $ARGV[$i]\n";
}

}
1;
