#!/usr/bin/python
from scipy.optimize import minimize 
import argparse
from sys import argv
import csv
import gzip
import random
import numpy

### command line arg parser
parser = argparse.ArgumentParser(description='Command Line Options')

### set required
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument('--bp1', type=int, required=True, help='Breakpoint One [int]')
required.add_argument('--bp2', type=int, required=True, help='Breakpoint Two [int]')
required.add_argument('--csv', type=str, required=True, help='Read Pair Coordinates [gzipped csv]')

### set optionals
optional = parser.add_argument_group('optional arguments')
optional.add_argument('--max_distance', type=int, default = 5e6, help='Maximum Distance To Breakpoint [int]' )
optional.add_argument('--min_distance', type=int, default = 0, help='Minimum Distance Between Links [int]' )
optional.add_argument('--bootstrap', type=int, default = 0, help='Number of Bootstraps [int]' )

## parse
args = parser.parse_args()

### data objects to store read positions
p1 = []
p2 = []

### compute distance function
def compute_distance( breakpoint ) : 

    ### take in the breakpoint estimates
    pos1, pos2 = breakpoint

    ### store distance here 
    dist = 0 

    ### iterate through positions
    for i in range( len(p1)-1 ) :

        ## now go through and update our lines and compute the distance for the points
        if ( ( p1[i] < pos1 and p2[i] > pos2 ) or ( pos1 < p1[i] < pos2 and pos1 < p2[i] < pos2 ) or ( p1[i] <= p2[i] < pos1 ) or ( p1[i] > pos2 and p2[i] > pos2 ) ) : 
            dist += ( p2[i] - p1[i] ) 
        elif ( p1[i] < pos1 and pos1 < p2[i] < pos2 ) :
            dist += ( pos1 + ( pos2 - p2[i] ) - p1[i] ) 
        elif ( pos1 < p1[i] < pos2 and pos2 < p2[i] ) : 
            dist += ( p2[i] - ( pos2 - ( p1[i] - pos1 ) ) )

    ### return the total distance spanned by the reads 
    return dist 

### read data into paired-lists 
with gzip.open( args.csv ) as tsv :

    ### split on tab
    for line in csv.reader(tsv, delimiter="\t") :

    	### remove read pairs that are too close
        if ( abs( float(line[1]) - float(line[2]) ) < args.min_distance ):
            continue

        ### check to make sure we're close enough to consider the read
        if ( ( abs( float(line[1]) - args.bp1 ) < args.max_distance or abs( float(line[1]) - args.bp2 ) < args.max_distance ) and ( abs( float(line[2]) - args.bp1 ) < args.max_distance or abs( float(line[2]) - args.bp2 ) < args.max_distance ) ):

        ### append read positions to list 
            if ( float( line[1] ) < float( line[2] ) ) :
                p1.append( float( line[1] ) )
                p2.append( float( line[2] ) )

                ### or if alternative, add in other order
            elif ( float( line[1] ) > float( line[2] ) ) :
                p1.append( float( line[2] ) )
                p2.append( float( line[1] ) )              

### now do the optimization
estimate = minimize(compute_distance, [args.bp1,args.bp2], method="Nelder-Mead", options={'maxiter':5000,'maxfev':5000,'ftol':100000} )

### print the output of the point estimate
print "Estimated Breakpoint Positions:\t", estimate.x[0], "\t", estimate.x[1]

### bootstrap confidence if desired
if ( args.bootstrap == 0 ) :
	exit()

### store bootstraps
p1_boot = p1
p2_boot = p2

### store positions of bootstrap estimates
boot1 = []
boot2 = []

### perform number of requested bootstraps
for b in range(args.bootstrap) :

	### clear lists
	p1 = []
	p2 = []

	### create draws
	for d in range(len(p1_boot)) :
		draw = random.randint( 0, len(p1_boot)-1 )
		p1.append(p1_boot[draw])
		p2.append(p2_boot[draw])

	### now run optimization
	estimate = minimize(compute_distance, [args.bp1,args.bp2], method="Nelder-Mead", options={'maxiter':5000,'maxfev':5000} )
	boot1.append( int(estimate.x[0]) )
	boot2.append( int(estimate.x[1]) )

	print b, estimate.x[0], estimate.x[1]

## sort each and output 95% CI position estimates
boot1.sort()
boot2.sort()

print boot1
print boot2

### print output
print "Boostrapped Confidence Intervals:\t", numpy.percentile(boot1, [2.5, 97.5]), "\t", numpy.percentile(boot2,[2.5,97.5])
