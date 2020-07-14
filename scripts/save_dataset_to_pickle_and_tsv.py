import sys
import pandas as pd
import numpy as np
import glob
import os
import warnings
warnings.simplefilter('ignore')
sys.path.append("../../oats")
from oats.biology.dataset import Dataset
from oats.utils.utils import save_to_pickle






def truncate_string(text, char_limit):
	truncated_text = text[:char_limit]
	if len(text)>char_limit:
		truncated_text = "{}...".format(truncated_text)
	return(truncated_text)



def save_combined_dataset_to_files(pickle_path, tsv_path, tsv_sample_path, input_dir, *input_filenames):

	print("\nworking on creating", os.path.basename(pickle_path))
	dataset = Dataset()
	for filename in input_filenames:
		filepath = os.path.join(input_dir, filename)
		dataset.add_data(pd.read_csv(filepath, lineterminator="\n"))
		print("finished adding data from {}".format(filepath))

	# Saving a version of the dataset prior to merging based on gene names.
	split_basename = os.path.basename(pickle_path).split(".")
	unmerged_pickle_path = os.path.join(os.path.dirname(pickle_path), "{}_unmerged.{}".format(split_basename[0],split_basename[1]))
	save_to_pickle(obj=dataset, path=unmerged_pickle_path)	

	# Saving the version of the dataset after merging based on gene names.
	print("merging rows based on gene names...")
	dataset.collapse_by_all_gene_names(case_sensitive=False)
	dataset.filter_has_description()

	# Write the dataset both to a pickle that can be loaded as an object and a tsv file for checking.
	save_to_pickle(obj=dataset, path=pickle_path)
	df = dataset.to_pandas()
	df = df.sort_values(by="id",inplace=False)
	df.to_csv(tsv_path, sep="\t", index=False)

	# Prepare a sample file from that whole dataset that makes it easy to understand what the context is.
	df = df.head(100)
	df["gene_names"] = df["gene_names"].map(lambda x: truncate_string(x, 30))
	df["gene_synonyms"] = df["gene_synonyms"].map(lambda x: truncate_string(x, 20))
	df["description"] = df["description"].map(lambda x: truncate_string(x, 100))
	df["term_ids"] = df["term_ids"].map(lambda x: truncate_string(x, 60))
	df.to_csv(tsv_sample_path, sep="\t", index=False)
	print("done")


reshaped_dir = "../reshaped_data"


save_combined_dataset_to_files("../pickles/genes_texts_annots.pickle", "../genes_texts_annots.tsv",  "../genes_texts_annots_sample.tsv", reshaped_dir, 
	"oellrich_walls_phene_descriptions.csv", 
	"oellrich_walls_phenotype_descriptions.csv", 
	"oellrich_walls_annotations.csv",
	"sgn_phenotype_descriptions.csv", 
	"maizegdb_phenotype_descriptions.csv", 
	"maizegdb_curated_go_annotations.csv",
	"tair_phenotype_descriptions.csv",
	"tair_curated_go_annotations.csv", 
	"tair_curated_po_annotations.csv"
	)

