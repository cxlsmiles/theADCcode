#!/usr/bin/perl

# Initialize a test structure
# Use Francesco's popana output
my %orig = ();
open(TESTFILE, '<', "adcpop.out") or die $!;
my $line = 0;
my %pops;
my @labels;
while (<TESTFILE>) {
        if (!$line) {
                @labels = split /\s+/;
		$line ++;
                next;
        }
        next if (($_ !~ /82\.0872/) && ($_ !~ /29\.5627/));
	my @lin = split  /\s+/;
	foreach my $i (2 .. 16) {
                my $lb =  join ('/',sort {$a cmp $b} split(/\//,$labels[$i-1]));
		$pops{sprintf("%1.4f %s",$lin[1],$lb)} = $lin[$i];
	}
}
close(TESTFILE);
# Loop over folders
foreach my $dir ("guk", "gus", "molcas") {
    my %test =(); 
    chdir "$ENV{PWD}/$dir";
    my $output =  `./popana.in 2>&1`;
    $output =~ s/.*?Eigenvalue\s+Orbital contributions\s+\-+\s//s;
    $output =~ s/Computing the spectrum.*?Eigenvalue\s+Orbital contributions\s+\-+\s//s;
    $line = 0;
    my $num_lbls = 0;
    for (split /^/, $output) {
       if (/^\s*$/){$line = 0; next;}
       if (!$line) {
                @labels = split /\s+/;
                $line ++;
                $num_lbls = scalar @labels;
                next;
        }
        next if (($_ !~ /82\.087/) && ($_ !~ /29\.5627/));
        my @lin = split  /\s+/;
        foreach my $i (2 .. $num_lbls) {
		my $lb =  join ('/',sort {$a cmp $b} split(/\//,$labels[$i-1]));
                $test{sprintf("%1.4f %s",$lin[1],$lb)} = $lin[$i];
        }

    }
    if( &compare(\%pops, \%test) ) {
	print uc $dir, " test passed.\n";
    } else {
	print uc $dir ," test failed.\n";
    } 
}



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
    my @vals1 = split /\s/,$keys1[$i]; @vals2 = split /\s/,$keys2[$i];
    return 0 if (($vals1[0] - $vals2[0]) > 1.1e-4);
    return 0 if ($vals1[1] ne $vals2[1]);
    return 0 if ($hash1{$keys1[$i]} != $hash2{$keys2[$i]});
  }  
  return 1;
}	


