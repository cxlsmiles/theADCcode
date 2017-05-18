#!/usr/bin/perl

# Extracts the spectrum as a contribution from atomic orbitals
use lib "/home/yasen/work/extract_spectra/";
use CreateXmgr;

#parse the adc output file
my $popana_flag = 0; # a flag raised when the population data is to be read
my %pop;
my @eigvals;
my $parts_read;      # counts the parts of the population data set
my $not_first_reading = 0; # a flag raised after the first population data set has been read

while(<>) {

    if (/^\s*ADC two\-hole population analysis$/)
    {
	$popana_flag = 1 ;
	$parts_read = -1;
    }
 
    if (/^\s*Computing spectrum for symmetry/)
    {
	$not_first_reading = 1 if $popana_flag;
	$popana_flag = 0 ;
    }
    next unless $popana_flag;
 
    
    next if /^\s*ADC two\-hole population analysis$/;
    next if /^\s+Eigenvalue\s+Orbital\scontributions$/;
    next if /^\s{4}-+$/;
    next if /^\s$/;
    
    
    @line = split /\s+/, $_;
    shift @line;
    
    if ($line[0] =~ /\d+\.\d+/)
    {
	my $eigval = shift @line;
	# all eigenvalues are saved in the first part of the population data set
	push @eigvals, $eigval unless $parts_read; 
	
	foreach (0..$#line) {
	    push @{$pop{$groups[$parts_read * 8 + $_]}}, $line[$_];
	}
	
    } else 
    {
	$parts_read++;
	next if $not_first_reading; 
	push @groups, @line;
    }
    
}
#done parsing


#write the data sets into a temporary file
$tempfile = ".temp_" . time . ".dat";
open(OUTFILE, ">$tempfile");

#get the total intensity 
my @ps = (0) x scalar @eigvals;

for (my $i = 0; $i < scalar @eigvals; $i++){
    foreach my $key ( keys %pop) {
	$ps[$i] += $pop{$key}[$i];
    }
}

foreach  my $key (sort keys %pop) {
    print OUTFILE "\#$key\n";
    for(my $i = 0; $i < scalar @eigvals; $i++) {
	
	my $num = 0;
	$num =  $pop{$key}[$i] if $pop{$key}[$i] > 0;
	print OUTFILE "$eigvals[$i] $ps[$i]\n";
	$ps[$i] -= $num;
	
    }
    print OUTFILE "&\n";
}

close(OUTFILE);
#done writing

#plot using grace
&create_xmgr(\%pop, $tempfile);


