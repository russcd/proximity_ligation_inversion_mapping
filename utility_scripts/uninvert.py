### this is a simple script to "univert" proximity ligation data generated for a sample with a large chromosomeal inversion relative to the reference genome
### it is provided as a utility script in hopes it will be useful but is not considered a core component of the pacakge

### usage:

##      python2 univert.py file.csv.gz bp1 bp2 | gzip - > uninverted.csv.gz

### this will convert the coordinates of each read pair that overlaps a predicted inversion to the "univerted" coordinates
### we have foudn it useful for reversing inversions in multiply inverted arrangements


     

#!/usr/bin/python
from sys import argv
import csv
import gzip

### initial prediction
breakpoint = [ float( argv[2] ), float( argv[3] ) ]
pos1, pos2 = breakpoint

### parameters
distance = 100e6

### data objects
p1 = []
p2 = []

### read data into paired-lists
with gzip.open( argv[1] ) as tsv :

	### split on tab
    for line in csv.reader(tsv, delimiter="\t") :

    	### check to make sure we're close enough to consider the read
    	if ( ( abs( float(line[1]) - breakpoint[0] ) < distance or abs( float(line[1]) - breakpoint[1] ) < distance ) and ( abs( float(line[2]) - breakpoint[0] ) < distance or abs( float(line[2]) - breakpoint[1] ) < distance ) ):

    		if ( float( line[1] ) < float( line[2] ) ) :
	    		p1.append( float( line[1] ) )
    			p2.append( float( line[2] ) )

    		### or if alternative, add in other order
       		elif ( float( line[1] ) > float( line[2] ) ) :
	    		p1.append( float ( line[2] ) )
    			p2.append( float ( line[1] ) )

### iterate through positions
for i in range( len(p1)-1 ) :
		
	## now go through and update our lines and compute the distance for the points
	if ( ( p1[i] < pos1 and p2[i] < pos1 ) or ( p1[i] < pos1 and p2[i] > pos2 ) or ( pos2 < p1[i] and pos2 < p2[i] ) ) :
		print '\t'.join(map(str,['2R',p1[i],p2[i]]))
	elif ( p1[i] < pos1 and pos1 < p2[i] < pos2 ) :
		p2[i] = pos1 + pos2 - p2[i] 
		print '\t'.join(map(str,['2R',p1[i],p2[i]]))
	elif ( pos1 < p1[i] < pos2 and pos2 < p2[i] ) : 
		p1[i] = pos1 + pos2 - p1[i] 
		print '\t'.join(map(str,['2R',p1[i],p2[i]]))
	elif ( pos1 < p1[i] < pos2 and pos1 < p2[i] < pos2 ) :
		p1[i] = pos1 + pos2 - p1[i] 
		p2[i] = pos1 + pos2 - p2[i] 
		print '\t'.join(map(str,['2R',p2[i],p1[i]]))

