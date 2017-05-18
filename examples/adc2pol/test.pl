#!/usr/bin/perl

#Compare the six first (largest transitions) eigenvalues

# Initialize a test structure
# Use  Kirill's excitation propagator 
my %orig = ();
open(TESTFILE, '<', "soeren_kirill_pol.out") or die $!;
my $roots = 0;
while (<TESTFILE>) {
    if (/^\s*\d+\)\s*(\d+\.\d{4})\s*(\d{1,2}\.\d)/)  {
        last if ($roots++ > 5);
        my $val = $1;
	$orig{"$val"}++;
    }
}
close(TESTFILE);

# Loop over folders
foreach my $dir ( "guk", "gus", "molcas") {
    my %test =(); 
    $roots = 0;
    chdir "$ENV{PWD}/$dir";
    my $output =  `./pol.in 2>&1`;
    while( $output =~ /\s*\d+\:\s*(\d+.\d{6})\,\s*(\d+.\d{2})\,\s*(\d\.\d{6})/g) {
        last if ($roots++ > 5);
	my $val = sprintf("%1.4f",$1/27.21139);
	$test{"$val"}++;
    }
    if( &compare(\%orig, \%test) ) {
	print uc $dir, " test passed.\n";
    } else {
	print uc $dir ," test failed.\n";
    }  
}

exit;

##### Subs
sub compare  {
  my %hash1 = %{$_[0]};
  my %hash2 = %{$_[1]};
  my $size1 = scalar keys %hash1;
  my $size2 = scalar keys %hash2;
  return 0 if $size2 != $size1;
  my @keys1 = sort keys %hash1;
  my @keys2 = sort keys %hash2;
  for my $i (0 .. ($size1 - 1)) {
    if (abs($keys1[$i] - $keys2[$i]) > 1.1e-4) {
      print $i, " ", $keys1[$i] , " ", $keys2[$i], " ", abs($keys1[$i] - $keys2[$i]), "\n";
      return 0;
    }
  }  
  return 1;
}	


