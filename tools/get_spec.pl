#!/usr/bin/perl

# Extracts the spectrum from an ADC calculation output
use lib "/home/yasen/work/extract_spectra/";
use CreateXmgr;

#parse the adc output file
my @eigvals;
my %eigcomp;
my $mainspace_flag = 0;
while(<>) {
    # get the eigenvalues and polestrengths throughout the file
    if (/^\s*\d+\:\s(\d+\.\d{6})\,\s(\d{1,2}\.\d{2})/) {
	push(@eigvals,$1);
    }
    
    if (/^\s*Overlaps with main\-space configurations:/) {
	$mainspace_flag = 1;
	next;
    }
    
    if (!/(\<[\d\,]+\|):([-\s]\d\.\d{6})/) {
	$mainspace_flag = 0;
    }

    #save the contributions for each eigenvalue
    while ($mainspace_flag && /(\<[\d\,]+\|):([-\s]\d\.\d{6})/g) 
    {
	$eigcomp{$1}[$#eigvals] = $2;
    }
}
#done parsing

#write the data sets into a temporary file
$tempfile = ".temp_" . time . ".dat";

open(OUTFILE, ">$tempfile");

#get the total intensity 
my @ps = (0) x scalar @eigvals;
for(my $i = 0; $i < @eigvals; $i++) {
    foreach  my $key (sort keys %eigcomp) {
	$ps[$i] += $eigcomp{$key}[$i]**2 if exists $eigcomp{$key}[$i];
    }
}

#save to the temporary file
foreach  my $key (sort keys %eigcomp) {
    print OUTFILE "\#$key\n";
    for(my $i = 0; $i < @eigvals; $i++) {	
	my $num = 0;
	$num =  $eigcomp{$key}[$i]**2 if exists $eigcomp{$key}[$i];
	print OUTFILE "$eigvals[$i] $ps[$i]\n" unless $num == 0;
	$ps[$i] -= $num;

    }
    print OUTFILE "&\n";
}

close(OUTFILE);
#done writing

#plot the temp file using grace
&create_xmgr(\%eigcomp, $tempfile);




