# ComMeta
ComMeta is a program usable for analyzing metagenome-derived genomes with a comparative approach. It provides statistics such as the Average Nucleotide Identity (ANI) and the Average Amino acid Identity (AAI) as well as a Reciprocal Best Hit analysis and core- and variable genome analysis. ComMeta also uses Circos to create plots showing either the overlapping region between genomes or information about single genomes (GC%, GC skew, features and contigs). To do the analyses various published tools have been used (listed below).

## Instalation
ComMeta is python 2.7 written and makes use of the Biopython and jinja2 python packages. To use ComMeta the following 3rd party dependencies should be on your system path.

- [Prokka](https://github.com/tseemann/prokka) >= v1.12: Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069. 
- [OrthoANIu](https://www.ezbiocloud.net/tools/orthoani) >= v1.2: Lee, I., Kim, Y. O., Park, S. C., & Chun, J. (2016). OrthoANI: an improved algorithm and software for calculating average nucleotide identity. International journal of systematic and evolutionary microbiology, 66(2), 1100-1103.
- [CompareM](https://github.com/dparks1134/CompareM) >= v0.0.23: Parks, D. (2016). CompareM A Toolbox for Comparative Genomics. Retrieved from https://github.com/dparks1134/CompareM
- seblastian: Frank, J. not publicly available
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/): Camacho, C. N. C. for B. I. (US). (2008). BLAST (r) Command Line Applications User Manual. National Center for Biotechnology Information (US).
- [Circos](http://circos.ca/software/download/circos/) >= v0.69-6: Krzywinski, M., Schein, J., Birol, İ., Connors, J., Gascoyne, R., Horsman, D., … Marra, M. A. (2009). Circos: An information aesthetic for comparative genomics. Genome Research , 19(9), 1639–1645. https://doi.org/10.1101/gr.092759.109
- [Mummer 3](https://sourceforge.net/projects/mummer/files/): Kurtz, S., Phillippy, A., Delcher, A. L., Smoot, M., Shumway, M., Antonescu, C., & Salzberg, S. L. (2004). Versatile and open software for comparing large genomes. Genome Biology, 5(2), R12. https://doi.org/10.1186/gb-2004-5-2-r12
- [R package](https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf) >= v1.6.19
: Chen, H. (2017). VennDiagram: Generate High-Resolution Venn and Euler Plots. Retrieved from https://cran.r-project.org/package=VennDiagram

The above stated 3rd party dependencies may also have dependencies that need to be downloaded on your system.

## Workflow

The basic command for ComMeta is running the entire pipeline including ANI, AAI, RBH, Dot plots, Circos plots and Venn diagrams. This can be performed using the following command:

```> ComMeta/console.py -i <input_dir> -o <output_dir> -s <subject_file>```

The -i <input_dir> requires the directory with the set of genomes to compare, which includes reference files in GBFF format and one subject file in FASTA format. The -o <output_dir> is the desired directory in which the output will be placed. The -s <subject_file> indicates the FASTA file of the genome you want to compare to the references.

A number of optional arguments can also be specified. These arguments are listed below:

- __-t__: The number of threads to use (default: 8)
- __--usearch__: The path to where usearch9 (for OrthoANIu) is installed
- __--window__: The window used for calculating the GC skew and GC percentage (default: 500)
- __--kmer__: The minimum alignment length to be shown in the overlapping region Circos plot (default: 1000)
- __--treshold__: The treshold for e-value for showing in the overlapping region Circos plot (default: 1e-6)

If you want to skip certain analyses, some can be turned off using the following arguments:

- __--no_dp__: Skip dot plots
- __--no_venn__: Skip Venn diagram
- __--no_scirc__: Skip single genome Circos plot
- __--no_ccirc__: Skip overlapping region Circos plot


