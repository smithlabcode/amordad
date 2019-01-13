# Amordad

## Introduction
Amordad is an efficient database engine for alignment-free, content-based indexing of metagenomic datasets.
Amordad places the metagenome comparison problem in a geometric
context, and uses an indexing strategy that combines random hashing
with a regular nearest neighbor graph. This framework allows refinement
of the database over time by continual application of random
hash functions, with the effect of each hash function encoded in the
nearest neighbor graph. This eliminates the need to explicitly maintain
the hash functions in order for query efficiency to benefit from the
accumulated randomness. Results on real and simulated data show
that Amordad can support logarithmic query time for identifying similar
metagenomes even as the database size reaches into the millions.

## Contacts and bug reports
Andrew D. Smith
andrewds@usc.edu

Ehsan Behnam
behnamgh@usc.edu

Wenzheng Li
wenzhenl@usc.edu

If you found a bug or mistake in this project, we would like to know about it. Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been fixed.
2. Check that your input is in the correct format and you have selected the correct options.
3. Please reduce your input to the smallest possible size that still produces the bug; we will need your input data to reproduce the problem, and the smaller you can make it, the easier it will be.

## LICENSE
The Amordad database engine for metagenomics Copyright (C) 2014 Andrew D Smith, Ehsan Behnam, Wenzheng Li, and the University of Southern California

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
