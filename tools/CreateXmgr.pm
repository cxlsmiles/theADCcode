
our @EXPORT = qw(create_xmgr);


#creates a xmgr file from a previously created temp file 
#uses the hash just to get the name of the data sets
sub create_xmgr 
{
    my %hash = %{shift @_};
    my $temp_file = shift @_;

    # get the title of the file from the current path
    my @title = split(/\/|\s/ , `pwd`);
    shift @title; 
    shift @title;
    shift @title;
    
    # set the xmgrace parameters
    my $par_string = "-pexec 'title \"". join(',', @title) ."\"' ";
    $par_string .= "-pexec 'yaxis label \"Pole strength\"' ";
    $par_string .= "-pexec 'xaxis label \"Ionization energy (eV)\"' ";
    my @key = sort keys %hash;

    for (0 .. $#key) {
	$par_string .= "-pexec \"s$_ type bar\" ";
	$par_string .= "-pexec \"s$_ symbol size 0\" ";
	$par_string .= "-pexec \"s$_ line type 0\" ";
	$par_string .= "-pexec 's$_ legend  \"$key[$_]\"' ";
    }
    
    my $command = "gracebat $par_string -settype bar $temp_file -timestamp -saveall " . join('-', @title) . 
	".agr -noprint";
    `$command; rm -f $temp_file`;
}



