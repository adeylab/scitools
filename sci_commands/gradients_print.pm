package sci_commands::gradients_print;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("gradients_print");

sub gradients_print {

print_gradients();

}
1;
