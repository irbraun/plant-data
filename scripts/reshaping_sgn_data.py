#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available from the Sol Genomics Network (SGN) database
# The purpose of this notebook is to read in and do a preliminary survey of the data related to text descriptions obtained from the Sol Genomics Network. The data was provided in the form of the input file by request, rather than obtained through an already available file from the database. The data needs to be organized and restructured into a standard format that will allow it to be easily combined with datasets from other resources. This notebook takes the following input files that were obtained from SGN and produces a set of files that have standard columns that are listed and described below.
# 
# ### Files read
# ```
# plant-data/databases/sgn/sgn_tomato_phenotyped_loci.txt
# ```
# 
# 
# ### Files created
# ```
# plant-data/reshaped_data/sgn_phenotype_descriptions.csv
# ```
# 
# ### Columns in the created files
# * **species**: A string indicating what species the gene is in, currently uses the 3-letter codes from the KEGG database.
# * **unique_gene_identifiers**: Pipe delimited list of gene identifers, names, models, etc which must uniquely refer to this gene.
# * **other_gene_identifiers**: Pipe delimited list of other identifers, names, aliases, synonyms for the gene, which may but do not have to uniquely refer to it.
# * **gene_models**: Pipe delimited list of gene model names that map to this gene.
# * **descriptions**: A free text field for any descriptions of phenotyes associated with this gene.
# * **annotations**: Pipe delimited list of gene ontology term identifiers.
# * **sources**: Pipe delimited list of strings that indicate where this data comes from such as database names.

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
from nltk.tokenize import word_tokenize
from nltk.tokenize import sent_tokenize

sys.path.append("../utils")
from constants import ABBREVIATIONS_MAP

sys.path.append("../../oats")
from oats.nlp.preprocess import concatenate_with_delim, replace_delimiter
from oats.nlp.small import remove_punctuation, remove_enclosing_brackets

OUTPUT_DIR = "../reshaped_data"
mpl.rcParams["figure.dpi"] = 200
warnings.simplefilter('ignore')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


# In[2]:


# Columns that should be in the final reshaped files.
reshaped_columns = ["species", 
 "unique_gene_identifiers", 
 "other_gene_identifiers", 
 "gene_models", 
 "descriptions", 
 "annotations", 
 "sources"]


# ### File with genes and phenotype descriptions (sgn_tomato_phenotyped_loci.txt)
# Note that fillna is being used here to replace missing values with an empty string. This is done so that the missing string will be quantified when checking for the number of occurences of unique values from different columns, see the analysis below. However this is not necessary as a preprocessing step because when the data is read in and appended to a dataset object later, any missing values or empty strings will be handled at that step.
# 
# There are several columns that contain information about gene names and accessions. We need to know what type of information is in each in order to know which should be retained in the dataset we are preparing. We are interested in both gene names that should map to a specific accession (like cyp716A12 or Medtr3g021350) as well as gene names that are enzyme descriptions (like Ubiquitin-Specific Protease) that could map to more than one gene in a particular species. Each type of information is valuable, but needs to be differentiated so that when comparing whether two rows are specifying the same gene, this is not confused with specifying two different genes that have the same function. In the case of this dataset, we only want to considered the locus names in a single column, the rest of the columns are more ambiguous and as long as all the mapping can be done with the locus names the other names can be ignored for downstream analysis.
# 
# This section creates a set of columns that have standardized names and include data in a standardized format that other functions within the package expect. The species column contains strings which are KEGG abbreviations for particular species. The gene names column contains any strings we want to consider to be uniquely mapped to some particular gene.

# In[3]:


filename = "../databases/sgn/sgn_tomato_phenotyped_loci.txt"
df = pd.read_table(filename)
df.fillna("", inplace=True)
df.head(10)


# In[4]:


# Removing rows that have missing inforrmation in the columns we want to keep.
df = df[(df["locus"] != "") & (df["allele_phenotype"] != "")]


# In[5]:


# Finding out how many unique values there are for each column.
unique_values = {col:len(pd.unique(df[col].values)) for col in df.columns}
for k,v in unique_values.items():
    print("{:24}{:8}".format(k,v))


# In[6]:


# Plotting distributions of number of word and phrases in each description.
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title("Phenotype Descriptions")
ax2.set_title("Phenotype Descriptions")
ax1.set_xlabel("Number of phrases")
ax2.set_xlabel("Number of words")
x1 = [len(sent_tokenize(x)) for x in df["allele_phenotype"].values]
x2 = [len(word_tokenize(x)) for x in df["allele_phenotype"].values]
ax1.hist(x1, bins=10, range=(0,10), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
ax2.hist(x2, bins=25, range=(0,50), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
fig.set_size_inches(15,4)
fig.tight_layout()
fig.show()
plt.close()


# In[7]:


# Organizing the desired information into a standard set of column headers.
combine_columns = lambda row, columns: concatenate_with_delim("|", [row[column] for column in columns])
df["species"] = "sly"
df["descriptions"] = df["allele_phenotype"]
df["unique_gene_identifiers"] = df.apply(lambda x: combine_columns(x, ["locus", "locus_symbol"]), axis=1)
df["other_gene_identifiers"] = df.apply(lambda x: combine_columns(x, ["locus_name", "allele_symbol", "allele_name"]), axis=1)
df["gene_models"] = df["locus"]
df["annotations"] = ""
df["sources"] = "SGN"
df = df[reshaped_columns]
df.head(10)


# In[8]:


path = os.path.join(OUTPUT_DIR,"sgn_phenotype_descriptions.csv")
df.to_csv(path, index=False)

