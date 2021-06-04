#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pandas as pd
import numpy as np
import glob
import os
import json
import re
from itertools import chain
from collections import defaultdict
import pandas as pd
import itertools
import networkx as nx
from gensim.parsing.preprocessing import strip_non_alphanum, stem_text, preprocess_string
from gensim.parsing.preprocessing import remove_stopwords, strip_punctuation
from nltk.corpus import brown, stopwords
from nltk.tokenize import sent_tokenize, word_tokenize
import warnings
warnings.simplefilter('ignore')


# In[ ]:


def to_json(df):
    infinite_defaultdict = lambda: defaultdict(infinite_defaultdict)
    split_on_bar_without_empty_strings = lambda x: [y.strip() for y in x.split("|") if y.strip() != ""]
    json_data = []
    for row in df.itertuples():
        d = infinite_defaultdict()
        d["_gene_id"] = row._1
        d["species_code"] = row.species_code
        d["species_name"] = row.species_name
        d["unique_gene_identifiers"] = split_on_bar_without_empty_strings(row.unique_gene_identifiers)
        d["other_gene_identifiers"] = split_on_bar_without_empty_strings(row.other_gene_identifiers)
        d["gene_models"] = split_on_bar_without_empty_strings(row.gene_models)
        d["text_unprocessed"] = row.text_unprocessed
        d["text_tokenized_sents"] = row.text_tokenized_sents
        d["text_tokenized_words"] = row.text_tokenized_words
        d["text_tokenized_stems"] = row.text_tokenized_stems
        d["annotations"] = split_on_bar_without_empty_strings(row.annotations)
        d["annotations_nc"] = split_on_bar_without_empty_strings(row.annotations_nc)
        d["reference_name"] = row.reference_name
        d["reference_file"] = row.reference_file
        d["reference_link"] = row.reference_link
        json_data.append(d)
    return(json_data)


# In[2]:


def truncate_string(text, char_limit):
    truncated_text = text[:char_limit]
    if len(text)>char_limit:
        truncated_text = "{}...".format(truncated_text)
    return(truncated_text)


def truncate_list(list_, item_limit):
    truncated_list = list_[:item_limit]
    return(truncated_list)


# In[3]:


# Make the full size json object and save to file.
path = "../final/data/genes_texts_annotations.csv"
df = pd.read_csv(path)
df = df.fillna("")
json_data = to_json(df)
json_path = "../final/data/genes_texts_annotations.json"
with open(json_path, "w") as f:
    json.dump(json_data, f, indent=4)


    
 
# Create a sample version of the file by truncating some of the strings and lists.
json_data = to_json(df)
json_path = "../final/samples/genes_texts_annotations.json"

# Subset both the number of entries in the dataset and truncate information in each field.
list_limit = 4
char_limit = 100
num_genes = 100
json_data = json_data[:num_genes]
for gene in json_data:
    gene["unique_gene_identifiers"] = truncate_list(gene["unique_gene_identifiers"], list_limit)
    gene["other_gene_identifiers"] = truncate_list(gene["other_gene_identifiers"], list_limit)
    gene["gene_models"] = truncate_list(gene["gene_models"], list_limit)
    gene["text_unprocessed"] = truncate_string(gene["text_unprocessed"], char_limit)
    gene["text_tokenized_sents"] = truncate_string(gene["text_tokenized_sents"], char_limit)
    gene["text_tokenized_words"] = truncate_string(gene["text_tokenized_words"], char_limit)
    gene["text_tokenized_stems"] = truncate_string(gene["text_tokenized_stems"], char_limit)
    gene["annotations"] = truncate_list(gene["annotations"], list_limit)
    gene["annotations_nc"] = truncate_list(gene["annotations_nc"], list_limit)

# This is an inelegant solution to formatting the json string the way we want to for the small sample file.
# Highly dependent on what structure of the dictinoary is, will break if that is changed.
indent_size = 4
s = json.dumps(json_data, indent=indent_size)
s = re.sub(r'": \[\s+', '": [', s)
s = re.sub(r'",\s+', '", ', s)
s = re.sub(r'"\s+\]', '"]', s)


# These are necessary because the above inelegant part misses newlines following str:str relationships in the json file.
s = s.replace(' "species_name":', '\n{}"species_name":'.format(" "*(indent_size*2)))
s = s.replace(' "annotations":', '\n{}"annotations":'.format(" "*(indent_size*2)))
s = s.replace(' "text_tokenized_sents":', '\n{}"text_tokenized_sents":'.format(" "*(indent_size*2)))
s = s.replace(' "text_tokenized_words":', '\n{}"text_tokenized_words":'.format(" "*(indent_size*2)))
s = s.replace(' "text_tokenized_stems":', '\n{}"text_tokenized_stems":'.format(" "*(indent_size*2)))
s = s.replace(' "reference_file":', '\n{}"reference_file":'.format(" "*(indent_size*2)))
s = s.replace(' "reference_link":', '\n{}"reference_link":'.format(" "*(indent_size*2)))
s = s.replace(' "unique_gene_identifiers":', '\n{}"unique_gene_identifiers":'.format(" "*(indent_size*2)))

with open(json_path, "w") as f:
    f.write(s)

print("done")

