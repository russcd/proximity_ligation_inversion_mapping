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

python2 map_breakpoints.py --bp1 [int] --bp2 [int] --csv [gzipped csv file] --max_distance [int] --min_distance [int] --bootstrap [int]

--bp1, the starting position for the optimization algorithm for the first breakpoint in basepairs 

--bp2, the starting positionat as above, but for the second breakpoint

In general, the starting position has not substantially affected breakpoint estimates when we run this program. However, it will generally be wise to vary the initial position provided in subsequent runs to make sure that optimized breakpoint position is estimated robustly. 

--csv, the gzipped paired read mapping position (see below)

--max_distance, the largest distance from the predicted starting breakpoints to consider, generally this parameter has little effect on the resulting breakpoint estimates, but it can sometimes speed up the program significantly

--min_distance, the minimum distance between two read pairs to consider, we have used 1,000 bp for this parameter. This will tend to remove self-ligated read pairs and those for which there is no cut site.  

--bootstrap, number of bootstraps to use to estimate breakpoint confidence intervals, the default is no bootstrapping. 

# Inversion Detection:
We have found that this program can also be used for inversion breakpoint detection by simply running it across an entire chromosome and providing evenly spaced point (at 1/3 and 2/3rds across the chromosome) as the intiial estimates. When no inversion is present, i.e., in a genome that is completely colinear with the refnerence, breakpoint estimates are unrealistically close together (less than 10 Kb), or span across the entire length of the chromosome. Visual inspect of the contact map associated with any predicted breakpoint should be performed. 

# Recommended Preparation of Input Data:
Generally, proximity ligated DNA is fragmented using a restriction enzyme and then end-repaired and blunt-end ligated. Therefore, we expect junctions between ligated DNA molecules to be demarcated with two tandem copies of the recognition sequence. Other approaches are possible so it is important to understand the specifics of each method. In the utility scripts, users will find a short perl script that will simply truncate all reads to the first tandem copy of the recognition sequence in our study, Mbol.

Once reads have been processed to remove junctions, map them to the appropriate reference genome using the default bwa mem alignment function. Remove PCR duplicates and from these, remove the subset of read pairs for which one or both do not map above a specified mapping quality threshhold and those for which read pairs map to different chromosomes. Finally, from the resulting bam file, extract the mapping coordinates for each read using cut and compress with gzip. E.g.,

samtools view BAMFILE | awk '$9 > 0' | cut -f3,4,9 | gzip - > inputfile.txt.gz 


