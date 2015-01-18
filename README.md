## About GEO

Gene Expression Omnibus(GEO) is a publicly available repository of community
submitted high throughput functional genomic datasets.

The entire list is as follows:

- Gene expression profiling by microarray or next-generation sequencing
- Non-coding RNA profiling by microarray or next-generation sequencing
- Chromatin immunoprecipitation (ChIP) profiling by microarray or next-generation sequencing
- Genome methylation profiling by microarray or next-generation sequencing
- Genome variation profiling by array
- SNP arrays
- Serial Analysis of Gene Expression (SAGE)
- Protein arrays

An experiment is summarized in the form of a 'GEO Series' which is a 
collection of 'GEO Samples' ([example](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM81022)).
Any 'GEO Sample' has an associated metadata such as the platform used for sequencing, the type of 
library preparation, protocol & technology used and other related information.

The datasets are available for public download via NCBI’s ftp server.
GEO datasets can also be programmatically retrieved using the 
[Entrez programming utilities](http://www.ncbi.nlm.nih.gov/books/NBK25499/)

Using Entrez utilities, it is possible to search for datasets belonging to 
a particular category of experiments, such as 
‘all experiments that involve whole genome bisulfite sequencing’.

Entrez also supports querying all the associated metadata with a particular 
record as is displayed on the NCBI website. 
It is thus possible to programmatically download all data
along with the metadata for any study.


## GEO Profiles

[GEO Profiles database](http://www.ncbi.nlm.nih.gov/geoprofiles/) is  a *curated* dataset of gene
expression profiles as derived from GEO datasets. They can be queried programmatically for various keywords,
[gene symbols, gene names etc.](http://www.ncbi.nlm.nih.gov/geo/info/profiles.html).

GEO profiles can be used to query nearest neighbors based on 
profile expression([example](http://www.ncbi.nlm.nih.gov/geoprofiles?LinkName=geoprofiles_geoprofiles_prof&from_uid=112040738) and can be used to find ‘related data’.

 Any query returns individual gene expressions. A
 ‘“smoking cancer”’ query will return all those gene expression studies
 that have “smoking cancer” annotation appearing somewhere. The ‘nearest profile neighbors’
 returns nearest neighbors based on precalculated pearson correlations. These correlations 
 are calculated for samples coming from the SAME dataset. Thus two independent datasets will
 never show up as nearest neighbors, even if the nature of experiment was similar.
 
 In short, ‘Profile neighbors’ is “within” the dataset.
 ‘Sequence neighbors’ are “across” the dataset


## MG-RAST
MG-RAST now has an API available: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004008


The raw and processed sequences are stored on their servers.. Can be retrieved easily too:
https://github.com/MG-RAST/MG-RAST-Tools/blob/master/examples/python/download_metagenome_sequences.py
(The script is broken right now, I submitted a PR for that)


## KBase.us

Looks promising, though I have not played around with it; 
especcialy interesting: ‘Test microbial ecological hypotheses through taxonomic and functional
analysis of quality-assessed metagenomic data’

They have an API too: http://kbase.us/developer-zone/services


## Periodically querying new data

In order to periodically updated the metadata we rely on the [XML feed](http://www.ncbi.nlm.nih.gov/geo/feed/series/) provided by NCBI.
A cron job downloads the xml every 24 hours and parses it. The parsed XML provides us new GEO datasets than have been uploaded
, a simple query in the local database can determine if there are new datasets available for download.

The metadata also has a [field](http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html) indicating the date of update and
it should thus be possible to capture the changes in already downloaded datasets too.

