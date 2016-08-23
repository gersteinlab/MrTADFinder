# MrTADFinder
MrTADFinder aims to identify topologically associating domains (TADs) in multiple resolutions.

##INPUT FILES:

MrTADFinder takes an whole-genome-to-whole genome contact map as an input. The contact map should be a sparse matrix stored in a tab-delimited file as follow

1	1464	39  
1	1768	19  
1	4455	20  
1	4458	70  
1	7368	29  
1	10413	25  
1	10563	30    
1	10687	37  
1	10690	24    
1	11123	28  
...  

The first 2 columns are indices of genomic bins, and the third column represents the contact frequency (corresponding matrix elements). This simple format has been widely used by mapping tools such as HiCPro.

The mapping between chromosome and genomics bins is specified by two annotation files. 

To bin the human genome (hg19) in 40kb bin size, an annotation file in the following format is provided (see ./data/bins_file1):

1	chr1	0	6231  
2	chr2	6232	12311  
3	chr3	12312	17262  
4	chr4	17263	22041  
5	chr5	22042	26564  
6	chr6	26565	30842  
7	chr7	30843	34821  
8	chr8	34822	38481  
9	chr9	38482	42012  
10	chr10	42013	45401  
11	chr11	45402	48777  
12	chr12	48778	52124  
13	chr13	52125	55004  
14	chr14	55005	57688  
15	chr15	57689	60252  
16	chr16	60253	62511  
17	chr17	62512	64541  
18	chr18	64542	66493  
19	chr19	66494	67972  
20	chr20	67973	69548  
21	chr21	69549	70752  
22	chr22	70753	72035  
23	chrX	72036	75917  
24	chrY	75918	77402  
25	chrM	77403	77403  

The first and second columns show the indices and names of various chromosome. Based on a bin size of 40kb, chromosome 1 is divided into 6232 bins (from bin 0 to bin 6231). Apart from the last bin, all the bins are 40kb. In this example, the whole human genome is divided int 77404 bins, and therefore the corresponding contact map is a square matrix of size 77404. Numbers running from 1 to 77404 are used to indexing the contact map file.

The 2nd file has the form

0	1	40000  
0	40001	80000  
0	80001	120000  
0	120001	160000  
0	160001	200000  
0	200001	240000  
0	240001	280000  
0	280001	320000  
0	320001	360000  
0	360001	400000  

IN this case, there are 77404 lines, representing all genomic bins (chromosome number (0=chr1 etc), start and end points).  

Annotation files based on binning human genome (hg19) in 40kb are provided. Users can use different annotation files for different bin sizes or different organisms.

##OUTPUT FILE:

Output file is simply a bed file that stores a list of TADs (chromosome number, start and end coordinates).  


##USAGE:

MrTADFinder is written in Julia. It has been tested in Julia v0.4.3. If Julia and the required packages are installed (see the first few lines in MrTADFinder.jl), one could simply run in the command prompt

> julia run_MrTADFinder.jl contact_map ./data/bins_file1 ./data/bins_file2 res=2.5 10 TAD_chr10.bed

The 1st agrument: contact map.  
The 2nd and 3rd agruments are the 2 annotation files.  
The 4th argument is the resolution parameter.  
The 5th argument is the chromosome of interest.  
The 6th argument is the path and name of the output file.


##REFERENCE:



