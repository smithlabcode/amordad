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


## Periodically querying new data

In order to periodically updated the metadata we rely on the [XML feed](http://www.ncbi.nlm.nih.gov/geo/feed/series/) provided by NCBI.
A cron job downloads the xml every 24 hours and parses it. The parsed XML provides us new GEO datasets than have been uploaded
, a simple query in the local database can determine if there are new datasets available for download.

The metadata also has a [field](http://www.ncbi.nlm.nih.gov/geo/info/qqtutorial.html) indicating the date of update and
it should thus be possible to capture the changes in already downloaded datasets too.


## More about GEO Profiles

For obtaining profile neighbors, they have pre-calculated pairwise
pearson correlations. So any query like this:
http://www.ncbi.nlm.nih.gov/geoprofiles?LinkName=geoprofiles_geoprofiles_prof&from_uid=33741191
randomly throws up top 200 hits(arbitrarily). The thresholds are arbitrary too.

There is something called ‘Chromosome neighbors’ too which outputs the
20 ‘physically’ close genes for the query gene in the same dataset.
E.g. http://www.ncbi.nlm.nih.gov/geoprofiles?LinkName=geoprofiles_geoprofiles_chr&from_uid=33741191

Then there are ‘Sequence Neighbors’ that simply gives link to all
similar datasets, taking sequence simiairty into account.
