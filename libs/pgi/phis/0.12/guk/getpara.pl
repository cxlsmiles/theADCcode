#!/usr/bin/perl 

use Getopt::Std;
use FileHandle;
use POSIX;

getopts( 'hd' );
if( $opt_h or not @ARGV ) {

  print STDERR << "EOF";
getpara.pl [-hdv] <dumpfile>
tries to determine the values of various constants compiled into GAMESS-UK
from its dumpfile.
-h\t prints this text
-d\t debug output
EOF

  exit;
}

# use a typeglob to map debugging output to stdout if verbose operation
# is requested. Send to nirvana otherwise.
if( $opt_d ) {
  *DEBUG = *STDERR; 
} else { open( DEBUG, "> /dev/null" ); }

# determine the endianess
if( unpack( "n", pack( "s", 10 )) != 10 ) {
  $little_endian = 1;
  print DEBUG "System is little endian.\n";
} else {
  print DEBUG "System is big endian.\n";
}

$path = $ARGV[0];
# see if path points to dumpfile or GUK executable
if( -x $path ) {
  # $path is the executable: create temporary dumpfile
  $df = POSIX::tmpnam();
  # delete temporary file unless in debug mode
  $deltmp = 1 unless $opt_d; 

  print DEBUG "Creating temporary dumpfile $df\n";
  $ENV{ed3} = $df;

  open( GUKIN, "| $path >/dev/null" );
  print GUKIN << "EOF";
zmat
he
end
basis sto3g
enter 1
EOF
  close GUKIN;
} else {
  $df = $path;  # $path points to existing dumpfile
}

open( FH, "$df" ) or die "Cannot open dumpfile $df: $!\n";

sysread FH, $buf, 4088;
@buf = unpack "S2040", $buf;
while( @buf ) {
  ++$sec;
  (undef, $blk, $typ, $len) = splice( @buf, 0, 4 );
  next if( $typ == 0 );
  warn "Section type $typ multiply defined.\n" if( defined $SecByType{typ} );
  $SecByType{$typ} = [ $sec, $blk, $len ];
}

print DEBUG join( "\t", "typ", "sec", "blk", "len" ), "\n";
foreach $typ ( sort { $a <=> $b } keys %SecByType ) {
  print DEBUG join( "\t", $typ, @{$SecByType{$typ}} ), "\n";
}

# *** process section type 1 ***
$blk = ${$SecByType{1}}[1];

# determine MXPRIM
$mxprim = CountWords( $blk );
#print VERBOSE "MXPRIM seems to be $mxprim\n";

# $tmp = MAXAQM*MAXPRIM - MAXAT
$tmp = CountWords();
print DEBUG "Going on...\n";

# determine MAXAT
$maxat = CountWords();
$maxat = ($maxat-10);
#print VERBOSE "MAXAT seems to be $maxat\n";

# calculate MAXAQM
$maxaqm = ($tmp-$maxat)/$mxprim;
#print VERBOSE "MAXAQM seems to be $maxaqm\n";

# determine MAXSHELL
$mxshel = CountWords();
$mxshel = ($mxshel - 2)/7*2; # FIXME: *2 because 2 4byte integers are one word
#print VERBOSE "MXSHEL seems to be $mxshel\n";

# *** process section type 3 ***
$blk = ${$SecByType{3}}[1];
++$blk;  # skip job specification

# determine MAXORB
$maxorb = CountWords( $blk );
$maxorb = ($maxorb - 4)/2;
#print VERBOSE "MAXORB seems to be $maxorb\n";

# print output suitable for shell execution and cleanup
print "maxorb=$maxorb maxat=$maxat mxshel=$mxshel mxprim=$mxprim maxaqm=$maxaqm";
unlink $df if( $deltmp );

sub CountWords {

  my $tot = 0;
  my( $buf, @buf );
  my( $pos, $cnt );
  
  if( $_[0] ) { 
    print DEBUG "Seeking block $_[0]\n";
    $pos = $_[0]*4096;
    seek FH, $pos, 0;
  }

  do {
    sysread FH, $buf, 4096;
    @buf = unpack "d511L2", $buf;
    @tmp = splice( @buf, 511, 4);
    if( $little_endian ) { pop @tmp; }
    $cnt = pop @tmp;
    print DEBUG "Block: $blk\t Count: $cnt\t Tot: $tot\n";
    ++$blk;
    $tot += $cnt;
  } while( $cnt == 511 );

return $tot;
}

