#!/usr/bin/perl

# Initialize a test structure
# Use J Chem Phys 123 (2005) 144115, TABLE V, \Sigma(4+), N2
my %orig = ();
$orig{"15.62"}=1;
$orig{"16.79"}=2;
$orig{"18.95"}=1;


# Loop over folders
foreach my $dir ("guk", "gus", "molcas") {
    my %test =(); 
    chdir "$ENV{PWD}/$dir";
    my $output =  `./ip.in 2>&1`;
    while( $output =~ /\s*\d+\:\s*(\d+.\d{6})\,\s*(\d+.\d{2})\,\s*(\d\.\d{6})/g) {
	my $val = sprintf("%1.2f",$1);
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
    if (abs($keys1[$i] - $keys2[$i]) > 1.1e-2) {
      print $i, " ", $keys1[$i] , " ", $keys2[$i], " ", abs($keys1[$i] - $keys2[$i]), "\n";
      return 0;
    }
  }  
  return 1;
}	


