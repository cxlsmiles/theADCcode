# run from the command line as 
# > perl test.pl [path_to_guk]

use File::Path;

# the path to the GUK executable may be given on the command line.
$gukexe = ($ARGV[0] || "guk");

@Results = (  -0.1879625707, -0.1856695536, -0.1758996984, -0.1840491744 );

$USER = $ENV{USER};

$thisdir = `pwd`; chomp $thisdir;
$exedir = "$thisdir/../../examples";
$scratch_dir = "/tmp/$USER/phis-test/"; mkpath $scratch_dir;

$ENV{ed3} = "$scratch_dir/dfile";
$ENV{ed6} = "$scratch_dir/mfile";

chdir $scratch_dir or die "Cannot change to $scratch_dir: $!";

foreach $test_no (0..3) {

  system("$gukexe <$thisdir/test.$test_no >gukout");

  print "test case no. $test_no\n";

  print "dump_all, energy-first ";
  pass_or_fail(system("$exedir/dump_all -g >/dev/null") == 0);
 
  print "dump_all, symmetry-first ";
  pass_or_fail(system("$exedir/dump_all -gs >/dev/null") == 0);

  @a = `$exedir/sample -g`;
  (undef,undef,$a) = split(/\s+/, pop @a);

  printf "MP2\t%f\t%f ", $Results[$test_no], $a;
  pass_or_fail(abs($Results[$test_no] - $a) < 1.0e-5);
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
