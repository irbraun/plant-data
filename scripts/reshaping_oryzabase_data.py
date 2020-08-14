#!/usr/bin/env python
# coding: utf-8

# ## Reading in data available through Oryzabase (Integrated Rice Science Database)
# The purpose of this notebook is to read in and do a preliminary analysis of the data related to text descriptions and ontology term annotations that are available through Oryzabase. The data needs to be organized and also restructured into a standard format that will allow it to be easily combined with datasets from other resources. This notebook takes the following input files that were obtained from OryzaBase and produces a set of files that have standard columns that are listed and described below.
# 
# ### Files read
# ```
# plant-data/databases/oryzabase/OryzabaseGeneListEn_20190826010113.txt
# ```
# 
# ### Files created
# ```
# plant-data/reshaped_data/oryzabase_phenotype_descriptions_and_annotations.csv
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

# In[17]:


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
from constants import NCBI_TAG, EVIDENCE_CODES

sys.path.append("../../oats")
from oats.nlp.small import add_prefix_safely, get_ontology_ids, remove_punctuation, remove_enclosing_brackets
from oats.nlp.preprocess import concatenate_texts, concatenate_with_delim, replace_delimiter

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
gene_model_pattern_1 = re.compile("grmzm.+")
gene_model_pattern_2 = re.compile("zm[0-9]+d[0-9]+")
is_gene_model = lambda s: bool(gene_model_pattern_1.match(s.lower() or gene_model_pattern_2.match(s.lower())))


# ### Part 1: Phenotypic Text Data
# There are several columns that contain information about gene names and accessions. We need to know what type of information is in each in order to know which should be retained or parsed to form the desired cleaned dataset. We are interested in both gene names that should map to a specific accession (like cms-54257) as well as gene names that are enzyme descriptions (like Ubiquitin-Specific Protease) that could map to more than one gene in a particular species. Each type of information is valuable, but needs to be differentiated so that when comparing whether two rows are specifying the same gene, this is not confused with specifying two different genes that have the same function.
# 
# The gene symbols in this dataset are typically surrounded by square brackets but not always. If a second symbol for that same gene is mentioned in the same column, the second symbol might be enclosed in parentheses. The synonynms for the gene symbols are similarly sometimes enclosed in square brackets, and are typically separated by commas in cases where more than one are mentioned in this column. Also note that an underscore is being used to represent missing data, so this has to handled so that that character is not treated as a gene name that appears many times.

# In[3]:


filename = "../databases/oryzabase/OryzabaseGeneListEn_20190826010113.txt"
usecols = ["CGSNL Gene Symbol", "Gene symbol synonym(s)", "CGSNL Gene Name", "Gene name synonym(s)", "Protein Name", 
           "Allele", "Explanation", "Trait Class", "RAP ID", "MUS ID", "Gramene ID", "Gene Ontology", 
           "Trait Ontology", "Plant Ontology"]
df = pd.read_table(filename, usecols=usecols, sep="\t")
df.fillna("", inplace=True)
unique_values = {col:len(pd.unique(df[col].values)) for col in df.columns}
print(df.shape)
for k,v in unique_values.items():
    print("{:25}{:8}".format(k,v))


# In[4]:


print(df[["CGSNL Gene Symbol", "Gene symbol synonym(s)"]].head(5))
print(df[["CGSNL Gene Symbol", "Gene symbol synonym(s)"]].sample(5))


# The values in the gene name synonym(s) column can be comma delimited lists if more than one synonym for the gene name was known. Sometimes quotes are used. Empty strings and possibly underscores can be used to denote missing information.

# In[5]:


print(df[["CGSNL Gene Name"]].head(5))
print(df[["CGSNL Gene Name"]].sample(5))


# The gene name column has strings representing the full name of each gene rather than just the symbol. Note that an underscore is also being used to denote missing values in this column as well.

# In[6]:


print(df[["Gene name synonym(s)"]].head(5))
print(df[["Gene name synonym(s)"]].sample(5))


# The values in the gene name synonym(s) column can be comma delimited lists if more than one synonym for the gene name was known. Sometimes quotes are used. Empty strings and possibly underscores can be used to denote missing information.

# In[7]:


print(df[["Protein Name","Allele"]].sample(5))


# Both the protein name and allele columns are sparse within the dataset. Either can be a single value or a comma delimited list of values. These may not need to be retained for finding reference to genes in other resources because we already have more standardized representations of the genes in other columns.

# In[8]:


print(df[["RAP ID","MUS ID"]].sample(5))


# Both the RAP ID and the MUS ID can columns can contain multiple values for a given gene which are included as members of a comma delimited list. These values can also be missing using the same scheme for missing information as in the rest of the dataset.

# The following functions were created based on the needs following how the gene symbols, names, synonyms, and accessions are previously described in this dataset. These are a not guaranteed to be a perfectly accurate method of parsing in the information in this dataset but they are meant to approximate what is required based on going through the dataset by hand. The methods are meant to be applied only to specific columns within the dataset, and to make the code that later cleans the columns more readable by compressing multiple cleaning steps into a single function. Some of these rely on other very specific functions that are within the text preprocessing module and not shown here.

# In[9]:


def handle_synonym_in_parentheses(text, min_length):
    # Looks at a string that is suspected to be in a format like "name (othername)". If
    # that is the case then a list of strings is returned that looks like [name, othername].
    # This is useful when a column is specifying something like a gene name but a synonym
    # might be mentioned in the same column in parentheses, so the whole string in that 
    # column is not useful for searching against as whole. Does not consider text in
    # parentheses shorter than min_length to be a real synonym, but rather part of the 
    # name, such as gene_name(t) for example. 
    names = []
    pattern = r"\(.*?\)"
    results = re.findall(pattern, text)
    for result in results:
        enclosed_string = result[1:-1]
        if len(enclosed_string)>=min_length:
            text = text.replace(result, "")
            names.append(enclosed_string)
    names.append(text)
    names = [name.strip() for name in names]
    return(names)


def clean_oryzabase_symbol(string):
    # Should be applied to the gene symbol column in the dataset.
    # Returns a single string representing a bar delimited list of gene symbols.
    string = string.replace("*","")
    names = handle_synonym_in_parentheses(string, min_length=4)
    names = [remove_enclosing_brackets(name) for name in names]
    names = [name for name in names if len(name)>=2] # Retain only names that are atleast two characters.
    names_string = concatenate_with_delim("|", names)
    return(names_string)

def clean_oryzabase_symbol_synonyms(string):
    # Should be applied to the gene symbol synonym(s) column in the dataset.
    # Returns a single string representing a bar delimited list of gene symbols.
    string = string.replace("*","")
    names = string.split(",")
    names = [name.strip() for name in names]
    names = [remove_enclosing_brackets(name) for name in names]
    names_string = concatenate_with_delim("|", names)
    return(names_string)


# ### Part 2: Ontology Term Annotations
# Multiple columns within the dataset specify ontology term annotations that have been applied to the geen mentioned on that particular line. Ontology term annotations are separated into different columns based on which ontology the terms belong to, and both the term ID of each annotation and the accompanying label for that term and explicitly given. Columns for terms from the Gene Ontolgoy (GO), Plant Ontology (PO), and Plant Trait Ontology (TO) are all included. There is no information about what evidence codes these annotations are associated with in this dataset.

# In[10]:


print(df[["Gene Ontology"]].head(5))
print(df[["Plant Ontology"]].head(5))
print(df[["Trait Ontology"]].head(5))


# Both the term IDs and labels for each annotation are given. Multiple annotations from the same ontology for a given line are separated by commas. We want to parse out just the gene ontology IDs for the cleaned dataset so that they can be referenced later, all the other information is not needed. There is a function to the return a list of gene IDs present in a longer string of text that is in the preprocessing module.

# ### Part 3: Handling other text description and keyword information in this data
# This dataset does not contain any columns that consistly contain a natural language description of a phenotype associated with a given gene. But some text-based information is still present. The trait class column contains a value from a limited set of keyword descriptors for the trait a particular gene is associated with. The size of the vocabulary used is obtained here. Also see the specific description of this keyword vocabulary here (https://shigen.nig.ac.jp/rice/oryzabase/traitclass/). The explaination column also occasionally contains text information about a corresponding phenotype.

# In[11]:


# Get a list sorted by number of occurences for each trait class.
description_counts = df["Trait Class"].value_counts().to_dict()
sorted_tuples = sorted(description_counts.items(), key = lambda x: x[1], reverse=True)
for t in sorted_tuples[0:10]:
    print("{:6}  {:20}".format(t[1],t[0][:70]))


# The most common value in the trait class column is whitespace or an empty string indicating missing data. Another very common value though is 'Other' which has 2,566 occurences out of the 17,674 total instances. This needs to be handled if using this information as text descriptions because this contains no semantics relevant to the phenotype (two phenotypes with trait classes of 'Other' should not be considered similar).

# In[12]:


print(df[["Explanation"]].sample(30))


# The explanation column holds information potentially about the phenotype, but also sometimes contains redundant information about the gene names or identifiers and sometimes the ontology term annotations as well. Sometimes methods are mentioned as well. Some of this could be handled with parsing to remove the redundant information that already appears somewhere else in a particular column for this line, but this should be considered irregular text annotations or descriptions for the purposes of downstream analyses. The following cell contains a preliminary attempt at a function that cleans values in this column by removing some redundant information from other columns.

# In[13]:


def clean_oryzabase_explainations(string):
    # Should be applied to the explaniation column in the dataset.
    # Returns a version of the value in that column without some of the redundant information.
    ontology_ids = get_ontology_ids(string)
    for ontology_id in ontology_ids:
        string = string.replace(ontology_id,"")
        string = remove_punctuation(string)
    return(string)


# In[14]:


# Restructuring and combining columns that have gene name information.
combine_columns = lambda row, columns: concatenate_with_delim("|", [row[column] for column in columns])
df["CGSNL Gene Symbol"] = df["CGSNL Gene Symbol"].apply(clean_oryzabase_symbol)
df["Gene symbol synonym(s)"] = df["Gene symbol synonym(s)"].apply(clean_oryzabase_symbol_synonyms)
df["CGSNL Gene Name"] = df["CGSNL Gene Name"].apply(lambda x: x.replace("_","").strip())
df["Gene name synonym(s)"] = df["Gene name synonym(s)"].apply(lambda x: replace_delimiter(text=x, old_delim=",", new_delim="|"))
df["gene_names"] = df.apply(lambda x: combine_columns(x, ["RAP ID","MUS ID","CGSNL Gene Symbol", "Gene symbol synonym(s)", "CGSNL Gene Name", "Gene name synonym(s)"]), axis=1)

# Restructuring and combining columns that have ontology term annotations.
df["Gene Ontology"] = df["Gene Ontology"].apply(lambda x: concatenate_with_delim("|", get_ontology_ids(x)))
df["Trait Ontology"] = df["Trait Ontology"].apply(lambda x: concatenate_with_delim("|", get_ontology_ids(x))) 
df["Plant Ontology"] = df["Plant Ontology"].apply(lambda x: concatenate_with_delim("|", get_ontology_ids(x))) 
df["term_ids"] = df.apply(lambda x: combine_columns(x, ["Gene Ontology","Trait Ontology","Plant Ontology"]), axis=1)

# Adding other expected columns and subsetting the dataset.
df["species"] = "osa"
df["description"] = df["Explanation"].apply(clean_oryzabase_explainations)
df["pmid"] = ""
df = df[["species", "gene_names", "description", "term_ids"]]
print(df[["species","gene_names"]].head(10))
print(df[["species","gene_names"]].sample(10))


# In[15]:


# Outputting the dataset of descriptions to a csv file.
path = os.path.join(OUTPUT_DIR,"oryzabase_phenotype_descriptions_and_annotations.csv")
df.to_csv(path, index=False)

