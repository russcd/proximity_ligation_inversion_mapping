#!/usr/bin/python
from scipy.optimize import minimize 
import argparse
from sys import argv
from sys import stderr
import csv
import gzip
import random
import numpy

### command line arg parser
parser = argparse.ArgumentParser(description='Command Line Options')

### set required
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument('--csv', type=str, required=True, help='Read Pair Coordinates [gzipped csv]')

### set optionals
optional = parser.add_argument_group('optional arguments')
optional.add_argument('--bp1', type=int, default = 0, help='Breakpoint One [int]')
optional.add_argument('--bp2', type=int, default = 0, help='Breakpoint Two [int]')
optional.add_argument('--max_distance', type=int, default = 5e12, help='Maximum Distance To Breakpoint [int]')
optional.add_argument('--min_distance', type=int, default = 0, help='Minimum Distance Between Links [int]')
optional.add_argument('--bootstrap', type=int, default = 0, help='Number of Bootstraps [int]')
optional.add_argument('--grid', type=int, default = 0, help='Distance between Gridpoints [int]')
optional.add_argument('--log', type=int, default = 1, help='Use log scaled distances [0,false | 1,true]')

## parse
args = parser.parse_args()

### check bp1 and bp2 are set if bootstrap or esimation
if args.grid == 0 and ( args.bp1 == 0 or args.bp2 == 0 ) :
	exit("must specify starting position estimate for optimization or bootstrapping using --bp1 and --bp2")

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
            if ( args.log == 0 ) :
                dist += ( p2[i] - p1[i] )
            else :
                dist += numpy.log( p2[i] - p1[i] )
        elif ( p1[i] < pos1 and pos1 < p2[i] < pos2 ) :
            if ( args.log == 0 ) :
                dist += ( pos1 + ( pos2 - p2[i] ) - p1[i] )
            else :
                numpy.log( pos1 + ( pos2 - p2[i] ) - p1[i] )
        elif ( pos1 < p1[i] < pos2 and pos2 < p2[i] ) :
            if ( args.log == 0 ) :
                dist += ( p2[i] - ( pos2 - ( p1[i] - pos1 ) ) )
            else :
                dist += numpy.log( p2[i] - ( pos2 - ( p1[i] - pos1 ) ) )

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

### now do the optimization if this is neither a bootstrap or grid search
if ( args.bootstrap == 0 and args.grid == 0 ) :

	### estimate the inversion breakpoint position
	estimate = minimize(compute_distance, [args.bp1,args.bp2], method="Nelder-Mead", options={'maxiter':5000,'maxfev':5000,'ftol':100000} )

	### print the output of the point estimate
	print "Estimated Breakpoint Positions:\t", estimate.x[0], "\t", estimate.x[1]
	print "Total Distance:\t", estimate.fun

### bootstrap confidence if desired
if args.bootstrap > 0 :

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

##		print b, estimate.x[0], estimate.x[1]

	## sort each and output 95% CI position estimates
	boot1.sort()
	boot2.sort()

	### print output
	print "95% Confidence Interval for Lower Breakpoint:\t", numpy.percentile(boot1, [2.5, 97.5]), "\n"
	print "95% Confidence Interval for Upper Breakpoint:\t", numpy.percentile(boot2,[2.5,97.5]), "\n"

### grid search if requested
if args.grid > 0 :

	### get total distance spanned without changes
	unmodified_distance = compute_distance( [0,0] )

	### x,y,z data objects for 3-d plot
	x_list = []
	y_list = []
	z_list = []

	### mins
	xmin = 0
	ymin = 0
	zmin = 1e1000

	### iterate through x,y points
	for x in range(int(min(p1)),int(max(p2)),args.grid) :
		for y in range(int(min(p1)),int(max(p2)),args.grid) :

			### skip duplicates
			if ( y < x ) :
				continue

			### compute point and update lists
			x_list.append(x)
			y_list.append(y)
			z = compute_distance( [x, y] )
			z_list.append(z)

			print x, "\t", y, "\t", z, "\t", z/unmodified_distance

			### find mins
			if zmin > z :
				zmin = z
				xmin = x
				ymin = y

	### provide users with possible starting point
	print >> stderr, "Suggesting Starting Point Is:\t", xmin, "\t", ymin
