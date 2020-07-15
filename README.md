## plant-data

This repository organizes information about plant phenotypes and annotations from a variety of data sources, and contains notebooks and scripts for preprocessing and merging this information. The emphasis is on obtaining datasets for looking into how well text mining approaches work with this data.

* All files are described in `file_descriptions.tsv`, with links to the original data where applicable.
* Python notebooks that preprocess and reshape the data from these files are in `preprocessing/`.
* Python scripts that generate both pickle and TSV versions of a combined dataset are in `scripts/`.
* The combined preprocessed dataset is available here as `genes_texts_annots.tsv`.
* A truncated sample from this dataset that is viewable on GitHub is `genes_texts_annots_sample.tsv`.

