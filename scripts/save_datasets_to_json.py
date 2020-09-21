import sys
import pandas as pd
import numpy as np
import glob
import os
import warnings
import json
import re
warnings.simplefilter('ignore')
sys.path.append("../../oats")
from oats.biology.dataset import Dataset










def truncate_string(text, char_limit):
	truncated_text = text[:char_limit]
	if len(text)>char_limit:
		truncated_text = "{}...".format(truncated_text)
	return(truncated_text)


def truncate_list(list_, item_limit):
	truncated_list = list_[:item_limit]
	if len(list_)>item_limit:
		truncated_list.append("...")
	return(truncated_list)




# Load the dataset that was created and saved in the previous step in the pipeline.
path = "../genes_texts_annots.csv"
dataset = Dataset(path)




# Output this dataset as a new json file.
json_data = dataset.to_json()
json_path = "../genes_texts_annots.json"
with open(json_path, "w") as f:
	json.dump(json_data, f, indent=4)





# Create a sample version of the file by truncating some of the strings and lists.
json_data = dataset.to_json()
json_path = "../genes_texts_annots_sample.json"


list_limit = 4
description_char_limit = 100
for gene in json_data:
	gene["unique_gene_identifiers"] = truncate_list(gene["unique_gene_identifiers"], list_limit)
	gene["other_gene_identifiers"] = truncate_list(gene["other_gene_identifiers"], list_limit)
	gene["gene_models"] = truncate_list(gene["gene_models"], list_limit)
	gene["descriptions"] = truncate_string(gene["descriptions"], description_char_limit)
	gene["annotations"] = truncate_list(gene["annotations"], list_limit)
	gene["sources"] = truncate_list(gene["sources"], list_limit)


s = json.dumps(json_data, indent=4)
s = re.sub(r'": \[\s+', '": [', s)
s = re.sub(r'",\s+', '", ', s)
s = re.sub(r'"\s+\]', '"]', s)
with open(json_path, "w") as f:
	f.write(s)






