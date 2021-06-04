#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available through the Planteome database
# The purpose of this notebook is to read in and do a preliminary survey of the data related to text and ontology annotations of text that is available through TAIR. The datasets need to be organized and also restructured into a standard format that will allow it be combined with datasets from other resources. This notebook takes the following input files that were downloaded from the TAIR and produces a set of files that have standard columns, which are listed and described below.
# 
# ### Files read
# ```
# plant-data/databases/planteome/biological_process.txt
# plant-data/databases/planteome/cellular_component.txt
# plant-data/databases/planteome/molecular_function.txt
# plant-data/databases/planteome/plant_anatomical_entity.txt
# plant-data/databases/planteome/plant_structure_development_stage.txt
# plant-data/databases/planteome/quality.txt
# ```
# 
# ### Files created
# ```
# plant-data/reshaped/data/planteome_curated_annotations.csv
# ```
# 
# ### Columns in the created files
# * **species_name**: String is the name of the species.
# * **species_code**: String identifier for the species, uses the 3-letter codes from KEGG.
# * **unique_gene_identifiers**: Pipe delimited list of gene identifers, names, models, etc that uniquely refer to this gene.
# * **other_gene_identifiers**: Same as the previous, but may not uniquely refer to a given gene.
# * **gene_models**: Pipe delimited list of gene model names, subset of unique_gene_identifiers.
# * **text_unprocessed**: A free text field for any descriptions of phenotyes associated with this gene.
# * **annotations**: Pipe delimited list of gene ontology term identifiers.
# * **reference_name**: String naming the database or paper that was the source for this data.
# * **reference_link**: The link to the reference resource if applicable.
# * **reference_file**: The specific name of the file from which this data comes if applicable.

# In[64]:


import sys
import os
import re
import warnings
import pandas as pd
import numpy as np
import itertools
import re

sys.path.append("../utils")
from constants import EVIDENCE_CODES

OUTPUT_DIR = "../reshaped/data"
warnings.simplefilter('ignore')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


# In[65]:


# Creating a list of lambdas for finding gene model strings.
gene_model_patterns = []
gene_model_patterns.append(re.compile("grmzm.+"))
gene_model_patterns.append(re.compile("zm[0-9]+d[0-9]+"))
gene_model_patterns.append(re.compile("solyc[0-9]+g[0-9]+"))
gene_model_patterns.append(re.compile("at[0-9]{1}g[0-9]+"))
gene_model_patterns.append(re.compile("medtr[0-9]{1}g[0-9]+"))
gene_model_patterns.append(re.compile("os[0-9]+g[0-9]+"))
gene_model_patterns.append(re.compile("loc_os[0-9]+g[0-9]+"))
is_gene_model = lambda s: any([bool(pattern.match(s.lower())) for pattern in gene_model_patterns])


# In[80]:


# Columns that should be in the final reshaped files.
reshaped_columns = [
 "species_name",
 "species_code",
 "unique_gene_identifiers", 
 "other_gene_identifiers", 
 "gene_models", 
 "text_unprocessed", 
 "annotations", 
 "reference_name",
 "reference_link",
 "reference_file"]

# The files that were obtained through the Planteome browser queries.
planteome_annotation_filepaths = [
    "../databases/planteome/biological_process.txt",
    "../databases/planteome/cellular_component.txt",
    "../databases/planteome/molecular_function.txt",
    "../databases/planteome/plant_anatomical_entity.txt",
    "../databases/planteome/plant_structure_development_stage.txt",
    "../databases/planteome/quality.txt"
]
# The fields that were included for each of those downloaded files.
columns = [
    "bioentity",
    "bioentity_name",
    "type",
    "bioentity_label",
    "annotation_class",
    "aspect",
    "annotation_extension_json",
    "taxon",
    "evidence_type",
    "evidence_with",
    "reference",
    "assigned_by"
]
# Read in all the files and stack rows, because they all have the same formatting and fields.
dfs = [pd.read_csv(path, sep="\t", names=columns) for path in planteome_annotation_filepaths]
df = pd.concat(dfs)
df.head(20)


# In[67]:


df.tail(20)


# In[68]:


# How many of each evidence type are there in this combined?
code_quantities = {c:len([x for x in df["evidence_type"] if EVIDENCE_CODES[x] in c]) for c in list(set(EVIDENCE_CODES.values()))}
for k,v in code_quantities.items():
    print("{:25}{:8}".format(k,v))


# In[69]:


df.shape


# In[70]:


# Retain just the annotations that we're considering high confidence.
high_confidence_categories = ["experimental","author_statement","curator_statement"]
df["high_confidence"] = df["evidence_type"].apply(lambda x: EVIDENCE_CODES[x] in high_confidence_categories)
df = df[df["high_confidence"]]
df.shape


# In[71]:


# Taxon IDs are given in the file but we need to use the naming scheme used across files to make it compatible.
ncbi_taxon_ids = ["3702","3847","3880","4081","4530","4577"]
species_code_strings = ["ath","gmx","osa","mtr","sly","zma"]
species_name_strings = ["Arabidopsis","soybean","rice","Medicago","tomato","maize"]
mapping = dict(zip(ncbi_taxon_ids,species_code_strings))
df["species_code"] = df["taxon"].map(lambda x: mapping.get(x[-4:],None))
mapping = dict(zip(ncbi_taxon_ids,species_name_strings))
df["species_name"] = df["taxon"].map(lambda x: mapping.get(x[-4:],None))
df.dropna(subset=["species_code"], axis=0, inplace=True)
print(pd.unique(df["species_name"]))
print(pd.unique(df["species_code"]))


# In[72]:


df.shape


# In[73]:


df.head(20)


# In[74]:


# Formatting the gene identifier columns to be the same as the other files.
df["unique_gene_identifiers"] = df.apply(lambda row: "{}|{}".format(row["bioentity_name"],row["bioentity_label"]), axis=1)
df["other_gene_identifiers"] = ""
df.head(20)


# In[75]:


# Other columns that are needed.
df["gene_models"] = df["unique_gene_identifiers"].map(lambda x: "".join([s for s in x.split("|") if is_gene_model(s)]))
df["text_unprocessed"] = ""
df["annotations"] = df["annotation_class"]
df["reference_name"] = "Planteome"
df["reference_link"] = "https://planteome.org/"


aspects = ["P", "C", "F", "A", "T", "G"]
filenames = [
    "biological_process.txt",
    "cellular_component.txt",
    "molecular_function.txt",
    "plant_anatomical_entity.txt",
    "quality.txt",
    "plant_structural_development_stage.txt"]

df["reference_file"] = df["aspect"].map(dict(zip(aspects,filenames)))


# In[76]:


df["ontology"] = df["annotation_class"].map(lambda x: x.split(":")[0])
pd.unique(df["ontology"])
df = df[df["ontology"].isin(["GO","PO","PATO"])]
df.shape


# In[77]:


df = df[reshaped_columns]
df.head(30)


# In[78]:


# Outputting the dataset of annotations to a csv file.
path = os.path.join(OUTPUT_DIR,"planteome_curated_annotations.csv")
df.to_csv(path, index=False)
df.head(30)


# In[79]:


df.sample(30)

