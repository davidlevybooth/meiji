# Meiji 0.1 Snakemake 

 Welcome the Meiji Snakemake pipeline. This suite of scripts
 is designed around the idea that by aligning short reads to a 
 comprehensive set of species-level genomes, we can quickly and 
 accurately quantify the microbiome. Species-level (ANI > 95%) 
 genomes and their revised taxonomic classifications were sourced from
 the Genome Taxonomy Database (GTDB) https://gtdb.ecogenomic.org/
 For ultra-fast and accurate read mapping we use the bbmap suite of 
 tools developed by the Joint Genome Institute. 
 https://jgi.doe.gov/data-and-tools/bbtools/ 

 Mapping to sequenced genomes also gives us the advantage of having 
 annotated functional data at our fingertips. Here, we wrap the functional
 profiling module from SHOGUN: https://github.com/knights-lab/SHOGUN
 (hence the name Meiji). We pulled functional annotations for each genome 
 from ANNOTREE http://annotree.uwaterloo.ca/, which are created using 
 UniRef100 clusters  https://www.uniprot.org/, KEGG https://www.kegg.jp/ 
 KOs and PFAM http://pfam.xfam.org/ annotations. 

 This work is in its earliest developmental iteration. It is built for use with 
 relatively powerful (0.5TB RAM) servers. An implementation for the cloud is 
 in development. Should you have any issues or concerns, please contact 
 David Levy-Booth, dlevyboo@mail.ubc.ca
 And thank you for choosing Meiji.
