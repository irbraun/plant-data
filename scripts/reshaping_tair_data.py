#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available through The Arabidopsis Information Resource (TAIR) database
# The purpose of this notebook is to read in and do a preliminary survey of the data related to text and ontology annotations of text that is available through TAIR. The datasets need to be organized and also restructured into a standard format that will allow it be combined with datasets from other resources. This notebook takes the following input files that were downloaded from the TAIR and produces a set of files that have standard columns, which are listed and described below.
# 
# ### Files read
# ```
# plant-data/databases/tair/Locus_Germplasm_Phenotype_20180702.txt
# plant-data/databases/tair/po_anatomy_gene_arabidopsis_tair.assoc
# plant-data/databases/tair/po_temporal_gene_arabidopsis_tair.assoc
# plant-data/databases/tair/ATH_GO_GOSLIM.txt
# ```
# 
# ### Files created
# ```
# plant-data/reshaped_data/tair_phenotype_descriptions.csv
# plant-data/reshaped_data/tair_all_go_annotations.csv
# plant-data/reshaped_data/tair_curated_go_annotations.csv
# plant-data/reshaped_data/tair_curated_po_annotations.csv
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
import re
import warnings
import pandas as pd
import numpy as np
import itertools
import re
import matplotlib.pyplot as plt
from nltk.tokenize import word_tokenize
from nltk.tokenize import sent_tokenize

sys.path.append("../../oats")
from oats.utils.constants import NCBI_TAG
from oats.utils.constants import EVIDENCE_CODES
from oats.utils.utils import to_abbreviation
from oats.nlp.preprocess import concatenate_with_bar_delim, remove_occurences_from_bar_delim_lists
from oats.nlp.preprocess import other_delim_to_bar_delim
from oats.nlp.preprocess import remove_punctuation
from oats.nlp.preprocess import remove_enclosing_brackets
from oats.nlp.preprocess import concatenate_descriptions
from oats.nlp.preprocess import add_prefix

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

# Creating and testing a lambda for finding gene model strings.
gene_model_pattern = re.compile("at[0-9]{1}g[0-9]+")
is_gene_model = lambda s: bool(gene_model_pattern.match(s.lower()))
assert is_gene_model("AT1G34534") == True
assert is_gene_model("at2g23452") == True
assert is_gene_model("ACAB1") == False


# ### File with genes and phenotype descriptions (Locus_Germplasm_Phenotype_20180702.txt)
# Reading in the dataset of phenotypic descriptions. There is only one value specified as the gene name (locus name) in the original dataset so this column does not need to be parsed further. The descriptions commmonly use semi-colons to separate phrases. The next cell gets the distribution of the number of phrases in each description field for the dataset of text descriptions, as determined by a sentence parser. The majority of the descriptions are a single sentence or phrase, but some contain more.

# In[3]:


filename = "../databases/tair/Locus_Germplasm_Phenotype_20180702.txt"
usecols = ["LOCUS_NAME", "PHENOTYPE"]
usenames = ["unique_gene_identifiers", "descriptions"]
renamed = {k:v for k,v in zip(usecols,usenames)}
df = pd.read_table(filename, usecols=usecols)
df.rename(columns=renamed, inplace=True)
df.sample(20)


# In[4]:


# Plotting distributions of number of phrases in each description.
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title("Phenotype Descriptions")
ax2.set_title("Phenotype Descriptions")
ax1.set_xlabel("Number of phrases")
ax2.set_xlabel("Number of words")
x1 = [len(sent_tokenize(x)) for x in df["descriptions"].values]
x2 = [len(word_tokenize(x)) for x in df["descriptions"].values]
ax1.hist(x1, bins=15, range=(0,15), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
ax2.hist(x2, bins=30, range=(0,200), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
fig.set_size_inches(15,4)
fig.tight_layout()
fig.show()
plt.close()


# In[5]:


# Restructuring the dataset to include all expected column names.
df["species"] = "ath"
df["other_gene_identifiers"] = ""
df["gene_models"] = df["unique_gene_identifiers"].map(lambda x: "|".join([s for s in x.split("|") if is_gene_model(s)]))
df["annotations"] = ""
df["sources"] = "TAIR"
df = df[reshaped_columns]

# Outputting the dataset of phenotype descriptions to csv file.
path = os.path.join(OUTPUT_DIR,"tair_phenotype_descriptions.csv")
df.to_csv(path, index=False)
df.head(20)


# ### File with gene ontology annotations (ATH_GO_GOSLIM.txt)
# Read in the file containing names of loci and corresponding information relating to gene ontology term annotation. Not all of the columns are used here, only a subset of them are read in. The relationship column refers to the relationships between the gene for that loci and the term mentioned on that given line. Evidence refer to the method of acquiring and the confidence in the annotation itself. This is retained so that we can subset that dataset based on whether the annotations are experimentally confirmed or simply predicted annotations. This section also looks at how many unique values are present for each field.
# 
# Each term annotation in this dataset is also associated with an evidence code specifying the method by which this annotation was made, which is related to the confidence that we can have in this annotation, and the tasks that the annotation should be used for. About half of the term annotations were made computationally, but there are also a high number of annotations available from high confidence annotations such as experimentally validated, curator statements, and author statements.

# In[6]:


filename = "../databases/tair/ATH_GO_GOSLIM.txt"
df_go = pd.read_table(filename, header=None, usecols=[0,2,3,4,5,9,12])
df_go.columns = ["locus","object","relationship","term_label","term_id","evidence_code","reference"]
unique_values = {col:len(pd.unique(df_go[col].values)) for col in df_go.columns}
print(df_go[["locus","object","term_id","evidence_code","reference"]].head(10))
print(df_go.shape)
for k,v in unique_values.items():
    print("{:18}{:8}".format(k,v))


# In[7]:


code_quantities = {c:len([x for x in df_go["evidence_code"] if EVIDENCE_CODES[x] in c]) 
             for c in list(set(EVIDENCE_CODES.values()))}
for k,v in code_quantities.items():
    print("{:25}{:8}".format(k,v))


# In[8]:


# Restructuring the dataset to include all the expected column names.
df_go["species"] = "ath"
df_go["unique_gene_identifiers"] = df_go["locus"]
df_go["other_gene_identifiers"] = ""
df_go["gene_models"] = df_go["unique_gene_identifiers"].map(lambda x: "".join([s for s in x.split("|") if is_gene_model(s)]))
df_go["descriptions"] = ""
df_go["annotations"] = df_go["term_id"]
df_go["sources"] = "TAIR"
high_confidence_categories = ["experimental","author_statement","curator_statement"]
df_go["high_confidence"] = df_go["evidence_code"].apply(lambda x: EVIDENCE_CODES[x] in high_confidence_categories)

# Subset to only include the high-quality GO annotations and ouptut the dataset to a csv file.
df_go_high_confidence = df_go[df_go["high_confidence"]==True]
path = os.path.join(OUTPUT_DIR,"tair_curated_go_annotations.csv")
df_go_high_confidence = df_go_high_confidence[reshaped_columns]
df_go_high_confidence.to_csv(path, index=False)

# Outputting the dataset of annotations to a csv file.
df_go = df_go[reshaped_columns]
path = os.path.join(OUTPUT_DIR,"tair_all_go_annotations.csv")
df_go.to_csv(path, index=False)
df_go.head(20)


# ### File with plant ontology term annotations (po_[termporal/anatomy]_gene_arabidopsis_tair.assoc)
# There are two separate files available that include annotations of PO terms. The files do not have headers so column names are added based on how the columns are described in the accompanying available readme files. One of the files contains annotations for PO terms that are spatial, or describe a specific part of plant anatomy or plant molecular structures. The other file contains annotations for PO terms that are temporal, or refer to a specific process or stage of development. These files are each read in separately, and the next cells look at the quantity of unique values in the columns of each dataset. There are more spatial annotations than temporal annotations, and a greater number of terms used to describe the spatial annotations.
# 
# The next field combines the two datasets of PO annotations and looks at the number of unique values for each column in the resulting dataset. Because there is no overlap in the terms between the two, the datasets are simply appended to one another and the total unique terms are a sum of the individual datasets.
# 
# Each term annotation in this dataset is also associated with an evidence code specifying the method by which this annotation was made, which is related to the confidence that we can have in this annotation, and the tasks that the annotation should be used for. Almost all of the PO term annotations are high confidence, they are experimentally validated, and only a few of them are derived from author statements.
# 
# The strings which are described in the synonyms column are included as references to each gene, and are combined with the gene name mentioned in the symbol column into a single bar delimited list.

# In[9]:


# Reading in the dataset of spatial PO term annotations.
filename = "../databases/tair/po_anatomy_gene_arabidopsis_tair.assoc"
df_po_spatial = pd.read_table(filename, header=None, skiprows=1, usecols=[2,4,5,6,9,10,11])
df_po_spatial.columns = ["symbol","term_id","references","evidence_code","name","synonyms","type"]
unique_values = {col:len(pd.unique(df_po_spatial[col].values)) for col in df_po_spatial.columns}
print(df_po_spatial.shape)
for k,v in unique_values.items():
    print("{:18}{:8}".format(k,v))


# In[10]:


# Reading in the dataset of temporal PO term annotations.
filename = "../databases/tair/po_temporal_gene_arabidopsis_tair.assoc"
df_po_temporal = pd.read_table(filename, header=None, skiprows=1, usecols=[2,4,5,6,9,10,11])
df_po_temporal.columns = ["symbol","term_id","references","evidence_code","name","synonyms","type"]
unique_values = {col:len(pd.unique(df_po_temporal[col].values)) for col in df_po_temporal.columns}
print(df_po_temporal.shape)
for k,v in unique_values.items():
    print("{:18}{:8}".format(k,v))


# In[11]:


# Looking at how many unique values each column has.
df_po = df_po_spatial.append(df_po_temporal, ignore_index=True)
unique_values = {col:len(pd.unique(df_po[col].values)) for col in df_po.columns}
print(df_po[["symbol","synonyms","evidence_code"]].head(10))
print(df_po.shape)
for k,v in unique_values.items():
    print("{:18}{:8}".format(k,v))


# In[12]:


# Quantifying the number of annotations of each type.
code_quantities = {c:len([x for x in df_po["evidence_code"] if EVIDENCE_CODES[x] in c]) 
             for c in list(set(EVIDENCE_CODES.values()))}
for k,v in code_quantities.items():
    print("{:25}{:8}".format(k,v))


# In[13]:


df_po['name'] = df_po['name'].astype(str)
df_po.dtypes


# In[14]:


# Restructuring the dataset to include all the expected column names.
df_po["species"] = "ath"
df_po["unique_gene_identifiers"] = np.vectorize(concatenate_with_bar_delim)(df_po["symbol"], df_po["name"])
df_po["other_gene_identifiers"] = df_po["synonyms"]
df_po["gene_model_strings_1"] =  df_po["unique_gene_identifiers"].map(lambda x: "|".join([s for s in x.split("|") if is_gene_model(s)]))
df_po["gene_model_strings_2"] =  df_po["unique_gene_identifiers"].map(lambda x: "|".join([s for s in x.split("|") if is_gene_model(s)]))
df_po["gene_models"] = np.vectorize(concatenate_with_bar_delim)(df_po["gene_model_strings_1"], df_po["gene_model_strings_2"])
df_po["descriptions"] = ""
df_po["annotations"] = df_po["term_id"]
df_po["sources"] = "TAIR"
df_po["other_gene_identifiers"] = df_po.apply(lambda row: remove_occurences_from_bar_delim_lists(row["other_gene_identifiers"],row["unique_gene_identifiers"]), axis=1)


df_po = df_po[reshaped_columns]

# Outputting the dataset of annotations to a csv file.
path = os.path.join(OUTPUT_DIR,"tair_curated_po_annotations.csv")
df_po.to_csv(path, index=False)
df_po.head(30)

