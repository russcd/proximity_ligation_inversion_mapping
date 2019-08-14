use strict ; 
use warnings ; 

### will replace GATCGATC with GATC and cut the remainder of the read and quality score
### usage :
#### gzip -dc fastq.gz | perl prune_GATC.pl | gzip - > fastq.pruned.gz 

while (<STDIN>) { 

	my $read = <STDIN> ; 
	my $mid = <STDIN> ; 
	my $qual = <STDIN> ; 

	chomp $qual ; 
	chomp $read ; 
	
	$read =~ s/GATCGATC.+/GATC/ ; 
	while ( length( $read ) < length( $qual ) ) { 
		chop $qual ; 
	}		

	print $_ ; 
	print $read, "\n" ;
	print $mid ; 
	print $qual, "\n" ;
}
