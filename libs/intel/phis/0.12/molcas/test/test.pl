# run from the command line as
# > perl test.pl
# 

use 5.6.0;

use File::Path;
use File::Copy;
use File::Temp qw/tempfile tmpnam/;
use FileHandle;

@Results = (  -0.1879625707, -0.1856695536, -0.1758996984, -0.1840491744 );

$USER = $ENV{USER};

$thisdir = qx(pwd); chomp $thisdir;
$exedir = "$thisdir/../../examples";

$workdir = "/tmp/$USER/phis-test/";

foreach $test_no (0..3) {

  $inp = new FileHandle;
  $inp->open("$thisdir/test.$test_no") 
    or die "Cannot open test file test.$test_no";
  &molcas($inp, 'seward', 'scf', 'motra');

  print "test case no. $test_no\n";

  print "dump_all, energy-first ";
  pass_or_fail(system("$exedir/dump_all -m >/dev/null") == 0);
 
  print "dump_all, symmetry-first ";
  pass_or_fail(system("$exedir/dump_all -ms >/dev/null") == 0);

  @a = `$exedir/sample -m`;
  (undef,undef,$a) = split(/\s+/, pop @a);
  printf "MP2\t%f\t%f ", $Results[$test_no], $a;
  pass_or_fail(abs($Results[$test_no] - $a) < 1.0e-5);
}

sub molcas {
#
# Expects a file handle to the input file and a list of modules to execute.
# Returns a file handle to the output file.
#
  my $input = shift;

  my $MOLCAS = $ENV{MOLCAS};

  my ($name, $module, $exe);
  my ($rd, $wt, $out);

  # this is for Molcas version 5
  $DataDir="data"; $ext=".exe";

  mkpath($workdir);
  chdir $workdir or die "Cannot change to $workdir: $!";

  unlink("BASLIB");
  symlink("$MOLCAS/basis_library",'BASLIB');

  foreach $name ('ABDATA', 'RYSRW', 'SYSDB') {
    symlink("$MOLCAS/$DataDir/$name", "$name");
  }

  symlink ('SCFORB', 'INPORB');

# execution
  $out = tmpnam( DIR => $workdir );
  while (@_) {

    $module = shift;
    $exe="$MOLCAS/bin/$module$ext";

    # this makes error messages work
    symlink("$MOLCAS/$DataDir/$module.MsgDB", 'ERRDB');

    $wt = new FileHandle;
    $wt->open("| $exe >>$out") or die "Cannot execute $exe: $!";
    $input->seek(0,0);
    while (<$input>) { print $wt $_; }
    undef $wt;

  }

  $rd = new FileHandle;
  $rd->open($out) or die "Cannot open MOLCAS output: $!\n";
  return $rd;
}

sub pass_or_fail {

    my $flag = shift;

    if ($flag) {
	printf "PASS\n";
    } else {
	printf "FAIL\n";
	exit(1);
    }
}
