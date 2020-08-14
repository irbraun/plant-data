## plant-data

### Description

This repository organizes information about plant genes, phenotypes, and annotations from a variety of data sources, and contains notebooks and scripts for preprocessing and combining this information. The emphasis is on obtaining datasets for looking into how well text mining approaches work with this data.

* All files are described in `file_descriptions.csv`, with links to the original data where applicable.
* The `pipeline.sh` script runs all of the other scripts for preprocessing and and merging the information present in those files. This pipeline generates files in the reshaped data directory. Full output files are not including here because of file sizes, but the first few lines of each are included in the reshaped samples directory. Reproducing the results of the pipeline in full requires all of the input files mentioned above, some of which are available through database subscriptions or requests. 
* The primary dataset of interest produced is available here as `genes_texts_annots.csv` and `genes_texts_annots.tsv`.
* A truncated sample from this dataset that is viewable on GitHub is `genes_texts_annots_sample.tsv`.

### Publication

A publication that uses this dataset to look at how computational methods can handle phenotype descriptions is in progress.

### Feedback
Send any feedback, questions, or suggestions to irbraun at iastate dot edu.
