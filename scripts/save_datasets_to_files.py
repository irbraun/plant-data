import sys
import pandas as pd
import numpy as np
import glob
import os
import warnings
warnings.simplefilter('ignore')
sys.path.append("../../oats")
from oats.biology.dataset import Dataset






def truncate_string(text, char_limit):
	truncated_text = text[:char_limit]
	if len(text)>char_limit:
		truncated_text = "{}...".format(truncated_text)
	return(truncated_text)



def save_combined_dataset_to_files(csv_path, tsv_path, tsv_sample_path, input_dir, *input_filenames):

	print("\nworking on creating", os.path.basename(csv_path))
	dataset = Dataset()
	for filename in input_filenames:
		filepath = os.path.join(input_dir, filename)
		dataset.add_data(pd.read_csv(filepath, lineterminator="\n"))
		print("finished adding data from {}".format(filepath))

	# Saving the version of the dataset after merging based on gene names.
	print("merging rows based on gene names...")
	dataset.filter_has_description()

	# Write the dataset to both tsv and csv files that can be used for loading it later.
	df = dataset.to_pandas()
	df = df.sort_values(by="id",inplace=False)
	df.to_csv(tsv_path, sep="\t", index=False)
	df.to_csv(csv_path, index=False)

	# Prepare a sample file from that whole dataset that makes it easy to understand what the context is.
	df = df.head(100)
	df["unique_gene_identifiers"] = df["unique_gene_identifiers"].map(lambda x: truncate_string(x, 30))
	df["other_gene_identifiers"] = df["other_gene_identifiers"].map(lambda x: truncate_string(x, 20))
	df["gene_models"] = df["gene_models"].map(lambda x: truncate_string(x, 30))
	df["descriptions"] = df["descriptions"].map(lambda x: truncate_string(x, 100))
	df["annotations"] = df["annotations"].map(lambda x: truncate_string(x, 60))
	df.to_csv(tsv_sample_path, sep="\t", index=False)
	print("done")


reshaped_dir = "../reshaped_data"

save_combined_dataset_to_files("../genes_texts_annots.csv", "../genes_texts_annots.tsv",  "../genes_texts_annots_sample.tsv", reshaped_dir, 
	"oellrich_walls_phene_descriptions.csv", 
	"oellrich_walls_phenotype_descriptions.csv", 
	#"oellrich_walls_annotations.csv",
	"sgn_phenotype_descriptions.csv", 
	"maizegdb_phenotype_descriptions.csv", 
	"maizegdb_curated_go_annotations.csv",
	"tair_phenotype_descriptions.csv",
	"tair_curated_go_annotations.csv", 
	"tair_curated_po_annotations.csv",
	"planteome_curated_annotations.csv"
	)

