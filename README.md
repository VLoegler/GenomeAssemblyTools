# GenomeAssemblyTools

This repository contains usefull python tools to deal with de novo genome assemblies using MUMmer4. 

### reorderContigs.py
This tools reorder the contigs of an assembly according to a reference assembly. The script changes order and orientation of contigs if needed. 

### NbContigsToX.py
Get the number of an assembly contigs that cover X% of a reference assembly. 

### maskTelomeres.py
Mask extremities of each chromosome (with a given size, default is 20kb). 

### findMergedChromosomes
Find contigs corresponding to the merge of several chromosomes of a reference. To be considered as merged, more than 90% of several reference chromosomes must align on the contig. 

### Rename contigs
Rename contigs with a reference genome. Contigs will be renamed with reference chromosome name if more than 20kb align on the contig. If several chromosomes align, name of all chromosomes will be in the new contig name.  
