#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available from the Oellrich, Walls et al. (2015) paper
# The purpose of this notebook is to read in and do a preliminary analysis of the data that is present in the supplementary file of this paper. The dataset also needs to be converted to a standard set of columns containing information in a standard format. This notebook takes the following input files that were obtained from that study and produces a set of files that have standard columns that are listed and described below.
# 
# 
# 
# ### Files read
# ```
# plant-data/papers/oellrich_walls_et_al_2015/versions_cleaned_by_me/13007_2015_53_MOESM1_ESM.csv
# ```
# 
# ### Files created
# ```
# plant-data/reshaped-data/oellrich_walls_phenotypes_descriptions.csv
# plant-data/reshaped-data/oellrich_walls_phene_descriptions.csv
# plant-data/reshaped-data/oellrich_walls_annotations.csv
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

# Creating a list of lambdas for finding gene model strings.
gene_model_patterns = []
gene_model_patterns.append(re.compile("grmzm.+"))
gene_model_patterns.append(re.compile("zm[0-9]+d[0-9]+"))
gene_model_patterns.append(re.compile("at[0-9]{1}g[0-9]+"))
gene_model_patterns.append(re.compile("medtr[0-9]{1}g[0-9]+"))
gene_model_patterns.append(re.compile("os[0-9]+g[0-9]+"))
is_gene_model = lambda s: any([bool(pattern.match(s.lower())) for pattern in gene_model_patterns])


# ### Phenotypic Text Data (oellrich_walls_dataset_irb_cleaned.txt)
# This data contains the phenotype descriptions for dominant mutants of genes across six different plant species. The data is read in from a cleaned version that removed some small delimiter errors from the original dataset that is available as a supplemental file from that publication. The data itself is unchanged.
# 
# There are several columns that contain information about gene names and accessions. We need to know what type of information is in each in order to know which should be retained in the dataset we are preparing. We are interested in both gene names that should map to a specific accession (like cyp716A12 or Medtr3g021350) as well as gene names that are enzyme descriptions (like Ubiquitin-Specific Protease) that could map to more than one gene in a particular species. Each type of information is valuable, but needs to be differentiated so that when comparing whether two rows are specifying the same gene, this is not confused with specifying two different genes that have the same function. In the case of this dataset, the gene symbol and gene identifier columns contain strings that we want to consider to be unique to a particular gene for a particular species, meaning that we can use those strings to look for these gene objects in other resources such as databases of pathway membership. The strings in the gene name column could be unique (narrow sheath1), but they can also be generic descriptions of enzymes (Ubiquitin-Specific Protease). For this reason, this column is not used in downstream analysis.
# 
# This dataset includes both full phenotype descriptions in one field, and atomized statements (which are phene descriptions) in another field. Either or both of these can be used as a source of text annotations on which to calculate similarity between phenotypes, phenes, or assess a hypothesized connected between genes in a network. We will look at quantity and properties of each of these categerogies of descriptions available and save the restructured datasets separately for each type.
# 
# This section creates a set of columns that have standardized names and include data in a standardized format that other functions within the package expect. The species column contains strings which are KEGG abbreviations for particular species. The gene names column contains any strings we want to consider to be uniquely mapped to some particular gene.
# 
# When saving the dataset using the phenotype descriptions as the text description column, there will be duplicates with respect to the combination of that column and the gene names column. This is because for each phenotype description there can be one or more atomized statement that it is comprised of. However, merging these rows requires also merging the ontology term annotations that each was annotated with, and this requires logic that is applied later. At this step we're only concerned with getting the right information in the right columns, and any datset with that correct can be merged later.

# In[3]:


filename = "../papers/oellrich_walls_et_al_2015/versions_cleaned_by_me/13007_2015_53_MOESM1_ESM.csv"
usecols = ["Species", "gene symbol", "Gene Identifier", "allele (optional)", 
           "gene name", "phenotype name", "phenotype description", 'atomized statement', 
           'primary entity1 ID', 'primary entity1 text', 'relation_to (optional)', 
           'primary entity2 ID (optional)', 'primary entity2 text (optional)', 
           'quality ID', 'quality text', 'PATO Qualifier ID (optional)', 
           'PATO Qualifier text (optional)', 'secondary_entity1 ID (optional)', 
           'secondary_entity1 text (optional)', 'relation_to (optional)', 
           'secondary entity2 ID (optional)','secondary_entity2 text (opyional)',
           'developmental stage ID (optional)', 'developmental stage text (optional)', 
           'condition ID (optional)', 'condition text (optional)', 'Pubmed ID (optional)', 
           'Dominant, recessive, codominant, semi-dominant (optional)', 
           'Loss or gain of function (optional)', 'Comment on mode of inheritance (optional)']
df = pd.read_csv(filename, usecols=usecols)
df.fillna("", inplace=True)
print(df.shape)
print(df[["gene symbol","Gene Identifier","allele (optional)","gene name"]].head(15))


# In[4]:


# Plotting distributions of number of words in each class of description.
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title("Phenotype Descriptions")
ax2.set_title("Phene Descriptions")
ax1.set_xlabel("Number of words")
ax2.set_xlabel("Number of words")
x1 = [len(word_tokenize(x)) for x in df["phenotype description"].values]
x2 = [len(word_tokenize(x)) for x in df["atomized statement"].values]
ax1.hist(x1, bins=30, range=(0,150), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
ax2.hist(x2, bins=30, range=(0,150), density=False, alpha=0.8, histtype='stepfilled', color="black", edgecolor='none')
fig.set_size_inches(15,4)
fig.tight_layout()
fig.show()
plt.close()


# In[5]:


# Finding the number of unique descriptions in each class of text description.
print(len(pd.unique(df["phenotype description"])))
print(len(pd.unique(df["atomized statement"])))


# ### Ontology Term Annotations (oellrich_walls_dataset_irb_cleaned.txt)
# There are several columns in the original dataset which refer to ontology terms, and specify a particular aspect of the EQ statement structure that that particular term refers to. For this dataset we are constructing, we will treat ontology term annotations as a 'bag of terms', and ignore the context of multi-term structured annotations such as EQ statements. Therefore these columns can be combined and any mentioned terms can be combined into a new column (as a bar delimited list). Contex of these terms in their respective ontologies are ignored (more than just leaf terms are retained), because this is handled later when comparing term sets.

# In[6]:


# Combining the different components of the EQ statement into a single column.
combine_columns = lambda row, columns: concatenate_with_delim("|", [row[column] for column in columns])
df["annotations"] = df.apply(lambda x: combine_columns(x, [
    "primary entity1 ID",
    "primary entity2 ID (optional)",
    "quality ID",
    "PATO Qualifier ID (optional)",
    "secondary_entity1 ID (optional)",
    "secondary entity2 ID (optional)",
    "developmental stage ID (optional)","condition ID (optional)",
    ]), axis=1)
    
df[["annotations"]].head(15)


# In[7]:


# Organizing the desired information into a standard set of column headers.
df["species"] = df["Species"].map(ABBREVIATIONS_MAP)
df["unique_gene_identifiers"] = df.apply(lambda x: combine_columns(x,["gene symbol", "gene name", "Gene Identifier"]), axis=1)
df["other_gene_identifiers"] = df["allele (optional)"]
df["gene_models"] = df["unique_gene_identifiers"].map(lambda x: "".join([s for s in x.split("|") if is_gene_model(s)]))
df["sources"] = "Plant PhenomeNET"
df[["species","unique_gene_identifiers","other_gene_identifiers","gene_models"]].head(20)


# In[8]:


# Saving a version that uses the full phenotype descriptions.
df["descriptions"] = df["phenotype description"]
df_subset = df[reshaped_columns]
df_subset["annotations"] = ""
path = os.path.join(OUTPUT_DIR,"oellrich_walls_phenotype_descriptions.csv")
df_subset.to_csv(path, index=False)

# Saving a version that uses the individual phene descriptions.
df["descriptions"] = df["atomized statement"]
df_subset = df[reshaped_columns]
df_subset["annotations"] = ""
path = os.path.join(OUTPUT_DIR,"oellrich_walls_phene_descriptions.csv")
df_subset.to_csv(path, index=False)

# Saving a version that includes only the ontology term annotations.
df["descriptions"] = ""
df_subset = df[reshaped_columns]
path = os.path.join(OUTPUT_DIR,"oellrich_walls_annotations.csv")
df_subset.to_csv(path, index=False)

