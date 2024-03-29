# mapmaker

## Notes on the 2019 update
This is an update of the MapMaker from Eric Lander's lab. The branches are
* release-3.1: Updates to compile with a C17 compiler.
* release-3.2: Removed unused code and reformatted code. 
* release-3.3: Added flip and 'all lods' commands. Bugs fixes are mostly applied in this and subsequent branches.
* release-3.4: Added multithreading.
* master: The most recent release

It's been a pleasure to bring MapMaker up to date. This a historic piece of software that is very well cited. Compilation should be straightforward using CMake. Speedup of the threaded code is fairly linear with core count. Additional performance appears available through vectorization, but datasets available to me were small and so I didn't implement vectorization. Please contact me if you are willing to provide larger data sets. To verify correctness of changes, I compared results from before each change, but I likely have introduced bugs. If you discover one, please contact me, supply a dataset and instructions on how to reproduce the problem and I'll look at creating a fix. Finally, I mainly focused on MapMaker/EXP, and I don't know about the state of MapMaker/QTL.

Kent A. Vander Velden, PhD
kent.vandervelden@gmail.com


## Original MapMaker abstract

MAPMAKER/EXP 3.0 AND MAPMAKER/QTL 1.1

From https://www.ncbi.nlm.nih.gov/pubmed/3692487 :
Also https://www.semanticscholar.org/paper/MAPMAKER%3A-an-interactive-computer-package-for-maps-Lander-Green/b941c4617ef321e4d7b42c598ad773fd55bf9c9d

1. Genomics. 1987 Oct;1(2):174-81.

MAPMAKER: an interactive computer package for constructing primary genetic
linkage maps of experimental and natural populations.

Lander ES(1), Green P, Abrahamson J, Barlow A, Daly MJ, Lincoln SE, Newberg LA.

Author information: 
(1)Whitehead Institute for Biomedical Research, Cambridge, Massachusetts 02142.

Erratum in
    Genomics. 2009 Apr;93(4):398. Newburg, L [corrected to Newberg, L A].

With the advent of RFLPs, genetic linkage maps are now being assembled for a
number of organisms including both inbred experimental populations such as maize 
and outbred natural populations such as humans. Accurate construction of such
genetic maps requires multipoint linkage analysis of particular types of
pedigrees. We describe here a computer package, called MAPMAKER, designed
specifically for this purpose. The program uses an efficient algorithm that
allows simultaneous multipoint analysis of any number of loci. MAPMAKER also
includes an interactive command language that makes it easy for a geneticist to
explore linkage data. MAPMAKER has been applied to the construction of linkage
maps in a number of organisms, including the human and several plants, and we
outline the mapping strategies that have been used.

DOI: 10.1016/0888-7543(87)90010-3 
PMID: 3692487  [Indexed for MEDLINE]


