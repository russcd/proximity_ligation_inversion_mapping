# proximity_ligation_inversion_mapping

# Overview:
This repository contains scripts for mapping inversion breakpoints using proximity ligation (Hi-C) sequencing data and a simple two-parameter optimization algorithm. An evaluation of the properties of this approach and application to fine-mapping inversion breakpoints of the An. gambiae species complex is availble from our preprint here: https://www.biorxiv.org/content/10.1101/662114v1.abstract. Here, we focus on practical details associated with using our implementation. 

# Dependencies:
The primary script "map_breakpoints.py" is a python 2 implementation of our algorithm with the following python module dependencies:
1. scipy.optimize
2. argparse
3. numpy

Each can be installed using pip. E.g., "pip install numpy"

# Recommended Usage:
The primary use for this script is fine-mapping inversion breakpoints using proximity ligation sequencing data.

python2 map_breakpoints.py --bp1 [int] --bp2 [int] --csv [gzipped csv file] --max_distance [int] --min_distance [int] --bootstrap [int] --grid [int] --log [int, 0|false, 1|true]

--bp1, the starting position for the optimization algorithm for the first breakpoint in basepairs 

--bp2, the starting positionat as above, but for the second breakpoint

In general, the starting position has not substantially affected breakpoint estimates when we run this program. However, it will generally be wise to vary the initial position provided in subsequent runs to make sure that optimized breakpoint position is estimated robustly. 

--csv, the gzipped paired read mapping position (see below)

--max_distance, the largest distance from the predicted starting breakpoints to consider, generally this parameter has little effect on the resulting breakpoint estimates, but it can sometimes speed up the program significantly

--min_distance, the minimum distance between two read pairs to consider, we have used 1,000 bp for this parameter. This will tend to remove self-ligated read pairs and those for which there is no cut site.  

--bootstrap, number of bootstraps to use to estimate breakpoint confidence intervals, the default is no bootstrapping.

--grid, distance between grid points for a coarse grid search covering the entire chromosome. If run, this will also produce a starting point for fine-mapping via the programs default behavior and output this to standard error. 

--log, optimize using the log-scaled distance spanned by read pairs

# Inversion Detection:
We have found that this program can also be used for inversion breakpoint detection by simply running it across an entire chromosome as a coarse grid search. This may be necessary if it is not known whether an inversion is present in the sample. If run in grid mode, the program will produce for all pairs of points across the grid the total distance spanned by all read pairs in the CSV file assuming an inversions breakpoints map to the coordinates of the point, and the ratio of the total distance spanned relative to their mapping positions with no assumed inversion. In general, if there is no inversion there should be no decrease in the total distance spanned by read pairs and the ratio will rarely be substantially lower than one. Therefore, a ratio less than one might be evidence of an inversion, and the minimum ratio across the grid can be used as a starting point for fine-mapping the breakpoint positions. However, in practice, many other features of the genome, e.g. translocations or missassemblies, might also result in false positives. We therefore caution that the results of this detection method should be interrogated graphically if possible. 

To run a grid search, and example command line that includes points spaced every 250Kb is:

python2 map_breakpoints.py --csv [CSV_FILE] --grid 250000 > [OUTPUT_GRID_fILE] 2> [SUGGESTED OPTIMIZATION STARTING POINT]

Please be aware that if grid points are spaced very densly (e.g. every 10 Kb), this may take a very long time to run. 

# Recommended Preparation of Input Data:
Generally, proximity ligated DNA is fragmented using a restriction enzyme and then end-repaired and blunt-end ligated. Therefore, we expect junctions between ligated DNA molecules to be demarcated with two tandem copies of the recognition sequence. Other approaches are possible so it is important to understand the specifics of each method. In the utility scripts, users will find a short perl script that will simply truncate all reads to the first tandem copy of the recognition sequence in our study, Mbol.

Once reads have been processed to remove junctions, map them to the appropriate reference genome using the default bwa mem alignment function. Remove PCR duplicates and from these, remove the subset of read pairs for which one or both do not map above a specified mapping quality threshhold and those for which read pairs map to different chromosomes. Finally, from the resulting bam file, extract the mapping coordinates for each read using cut and compress with gzip. E.g.,

samtools view BAMFILE | awk '$9 > 0' | cut -f3,4,8 | gzip - > inputfile.txt.gz 


