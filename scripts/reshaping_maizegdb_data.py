#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available through Maize GDB (Maize Genetics and Genomics Database)
# The purpose of this notebook is to read in and do a preliminary analysis of the data related to text descriptions that are available through Maize GDB. The data was provided in the form of the input file by a request through Maize GDB curators, rather than obtained through an already available file from the database. The data needs to be organized and also restructured into a standard format that will allow it to be easily combined with datasets from other resources. This notebook takes the following input files that were obtained from MaizeGDB and produces a set of files that have standard columns that are listed and described below.
# 
# ### Files read
# ```
# plant-data/databases/maizegdb/pheno_genes.txt
# plant-data/databases/maizegdb/maize_v3.gold.gaf
# ```
# 
# 
# ### Files created
# ```
# plant-data/reshaped/data/maizegdb_phenotype_descriptions.csv
# plant-data/reshaped/data/maizegdb_go_annotations.csv
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

# In[1]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import sys
import os
import warnings
import pandas as pd
import numpy as np
import itertools
import re
import matplotlib.pyplot as plt
import nltk
nltk.download('punkt')
from nltk.tokenize import word_tokenize
from nltk.tokenize import sent_tokenize

sys.path.append("../utils")
from constants import NCBI_TAG, UNIPROT_TAG, EVIDENCE_CODES

sys.path.append("../../oats")
from oats.nlp.preprocess import concatenate_with_delim, subtract_string_lists, replace_delimiter, concatenate_texts
from oats.nlp.small import remove_punctuation, remove_enclosing_brackets, add_prefix_safely

OUTPUT_DIR = "../reshaped/data"
mpl.rcParams["figure.dpi"] = 200
warnings.simplefilter('ignore')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


# In[2]:


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

# Creating and testing a lambda for finding gene model strings.
gene_model_pattern_1 = re.compile("grmzm.+")
gene_model_pattern_2 = re.compile("zm[0-9]+d[0-9]+")
is_gene_model = lambda s: bool(gene_model_pattern_1.match(s.lower() or gene_model_pattern_2.match(s.lower())))


# ### File with genes and phenotype descriptions (pheno_genes.txt)
# Note that fillna is being used here to replace missing values with an empty string. This is done so that the missing string will be quantified when checking for the number of occurences of unique values from different columns, see the analysis below. However this is not necessary as a preprocessing step because when the data is read in and appended to a dataset object later, any missing values or empty strings will be handled at that step.

# In[3]:


filename = "../databases/maizegdb/pheno_genes.txt"
usecols = ["phenotype_name", "phenotype_description", "locus_name", "alleles", "locus_synonyms", "v3_gene_model", "v4_gene_model", "uniprot_id", "ncbi_gene"]
df = pd.read_table(filename, usecols=usecols)
df.fillna("", inplace=True)
print(df[["phenotype_name","phenotype_description"]].head(10))
print(df.shape)


# In[4]:


df


# Text information about the phenotypes are contained in both the phenotype name and phenotype description for these data. The can be concatenated and retained together in a new description column that contains all this information, or just the phenotype description could be retained, depending on which data should be used downstream for making similarity comparisons. This is different than for most of the other sources of text used. The next cell looks at how many unique values there are in this data for each column.

# In[5]:


# Finding out how many unique values there are for each column.
unique_values = {col:len(pd.unique(df[col].values)) for col in df.columns}
for k,v in unique_values.items():
    print("{:24}{:8}".format(k,v))


# There are a fairly small number of distinct phenotype descriptions (379) compared to the number of lines that are in the complete dataset (3,616). This means that the same descriptions is occuring many times. Look at which descriptions are occuring most often.

# In[6]:


# Get a list sorted by number of occurences for each phenotype description.
description_counts = df["phenotype_description"].value_counts().to_dict()
sorted_tuples = sorted(description_counts.items(), key = lambda x: x[1], reverse=True)
for t in sorted_tuples[0:10]:
    print("{:6}    {:20}".format(t[1],t[0][:70]))


# The only description that occurs far more often than the next is an empty string, where this information is missing entirely. The next cell looks at how many phrases are included in the phenotype description values. Most have a single phrase, some have multiple. These look like they are mainly separated with semicolons.

# In[7]:


# Plotting distributions of number of phrases in each description.
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title("Phenotype Descriptions")
ax2.set_title("Phenotype Descriptions")
ax1.set_xlabel("Number of phrases")
ax2.set_xlabel("Number of words")
x1 = [len(sent_tokenize(x)) for x in df["phenotype_description"].values]
x2 = [len(word_tokenize(x)) for x in df["phenotype_description"].values]
ax1.hist(x1, bins=15, range=(0,15), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
ax2.hist(x2, bins=30, range=(0,150), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
fig.set_size_inches(15,4)
fig.tight_layout()
fig.show()
plt.close()


# In[8]:


# Restructuring the dataset to include all the expected column names.
combine_gene_columns = lambda row, columns: concatenate_with_delim("|", [row[column] for column in columns])
combine_text_columns = lambda row, columns: concatenate_with_delim("; ", [row[column] for column in columns])
df["text_unprocessed"] = df.apply(lambda x: combine_text_columns(x, ["phenotype_name", "phenotype_description"]), axis=1)
df["uniprot_id"] = df["uniprot_id"].apply(add_prefix_safely, prefix=UNIPROT_TAG)
df["ncbi_gene"] = df["ncbi_gene"].apply(add_prefix_safely, prefix=NCBI_TAG)
df["unique_gene_identifiers"] = df.apply(lambda x: combine_gene_columns(x, ["locus_name", "alleles", "v3_gene_model", "v4_gene_model", "uniprot_id", "ncbi_gene"]), axis=1)
df["other_gene_identifiers"] = df["locus_synonyms"]
df["gene_models"] = df.apply(lambda x: combine_gene_columns(x, ["v3_gene_model", "v4_gene_model"]), axis=1)
df["species_name"] = "maize"
df["species_code"] = "zma"
df["annotations"] = ""
df["reference_name"] = "MaizeGDB"
df["reference_link"] = "https://www.maizegdb.org/"
df["reference_file"] = "pheno_genes.txt"
df["other_gene_identifiers"] = df.apply(lambda row: subtract_string_lists("|", row["other_gene_identifiers"],row["unique_gene_identifiers"]), axis=1)
df = df[reshaped_columns]
df.head()


# In[9]:


df["text_unprocessed"].values[:1000]


# In[10]:


# Outputting the dataset of descriptions to a csv file.
path = os.path.join(OUTPUT_DIR,"maizegdb_phenotype_descriptions.csv")
df.to_csv(path, index=False)


# ### File with high confidence gene ontology annotations (maize_v3.gold.gaf)
# This file was generated as part of the [Maize GAMER](https://onlinelibrary.wiley.com/doi/full/10.1002/pld3.52)  publication (Wimalanathan et al., 2018). The annotations include all of the associations between maize genes and ontology terms from GO where the terms have been experimentally confirmed to represent correct functional annotations for those genes.

# In[11]:


filename = "../databases/maizegdb/maize_v3.gold.gaf"
df = pd.read_table(filename, skiprows=1)
df.fillna("", inplace=True)
df.head()


# In[13]:


# Restructuring the dataset to include all the expected column names.
df["text_unprocessed"] = ""
df["unique_gene_identifiers"] = df.apply(lambda x: combine_gene_columns(x, ["db_object_id", "db_object_symbol"]), axis=1)
df["other_gene_identifiers"] = df.apply(lambda x: combine_gene_columns(x, ["db_object_name", "db_object_synonym"]), axis=1)
df["gene_models"] =  df["unique_gene_identifiers"].map(lambda x: "".join([s for s in x.split("|") if is_gene_model(s)]))
df["species_name"] = "maize"
df["species_code"] = "zma"
df["annotations"] = df["term_accession"]
df["reference_name"] = "MaizeGDB"
df["reference_link"] = "https://www.maizegdb.org/"
df["reference_file"] = "pheno_genes.txt"
df = df[reshaped_columns]
df.head()


# In[14]:


# Outputting the dataset of annotations to a csv file.
path = os.path.join(OUTPUT_DIR,"maizegdb_curated_go_annotations.csv")
df.to_csv(path, index=False)

