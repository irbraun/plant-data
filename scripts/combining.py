#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import pandas as pd
import numpy as np
import glob
import os
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

sys.path.append("../../oats")
from oats.nlp.preprocess import concatenate_with_delim
from oats.annotation.ontology import Ontology
from oats.annotation.annotation import annotate_using_noble_coder


# In[2]:


# Which files are the reshaped ones that should be combined.
paths = [
    "oellrich_walls_phene_descriptions.csv", 
    "oellrich_walls_phenotype_descriptions.csv", 
    "oellrich_walls_annotations.csv",
    "sgn_phenotype_descriptions.csv", 
    "maizegdb_phenotype_descriptions.csv", 
    "maizegdb_curated_go_annotations.csv",
    "tair_phenotype_descriptions.csv",
    "tair_curated_go_annotations.csv", 
    "tair_curated_po_annotations.csv",
    "planteome_curated_annotations.csv"
]


# In[3]:


# Check to make sure that these are the columns present in each read in file.
expected_columns = [
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

# Stack all of those dataframes to create one new dataframe.
dfs_to_be_stacked = []
for path in paths:
    df = pd.read_csv(os.path.join("..","reshaped_data",path))
    assert set(expected_columns) == set(df.columns)
    dfs_to_be_stacked.append(df)
df = pd.concat(dfs_to_be_stacked, ignore_index=True)
dfs_to_be_stacked = None
df.head()


# In[4]:


# Need to use all the information in the gene identifier columns to create internal gene IDs
# Add a column that acts as an old index with one value for each existing row.
df = df.reset_index()
df.rename({"index":"old_id"},axis="columns",inplace=True)

# Create edges for the graph from the gene names.
# These edges should go from the old row IDs to each of the unique gene identifier strings.
def generate_edges(row,case_sensitive):
    if case_sensitive:
        names = row["unique_gene_identifiers"].split("|")
    else:
        names = row["unique_gene_identifiers"].lower().split("|")
    edges = [(str(row["old_id"]),"{}[SEP]{}".format(row["species_code"],name)) for name in names]
    return(edges)

# Create the network using those edges find the connected components and add those.
g = nx.Graph()
case_sensitive=False
edges = df.apply(generate_edges, case_sensitive=case_sensitive, axis=1)
edges = list(chain.from_iterable(edges.values))
g.add_edges_from(edges)

# Find the mapping between old IDs and connected components.
node_to_component = {}
component_index = 0
for node_set in nx.connected_components(g):
    for node in node_set:
        node_to_component[node] = component_index
    component_index = component_index+1
    
# The connected components number now serves as the gene ID.
df["old_id"] = df["old_id"].map(lambda x: str(x))
df["_gene_id"] = df["old_id"].map(node_to_component)
df.head()


# In[5]:


df._gene_id.value_counts()


# In[6]:


# Now the dataset looks clean, but the gene columns don't yet reflect all the information
# that was used in the network creation step. For example, we want all the gene identifiers
# found on any row to be present everywhere in the dataset for that given line.


# In[7]:


agg_df = df.groupby("_gene_id").agg({
    "unique_gene_identifiers": lambda x: concatenate_with_delim("|",x),
    "other_gene_identifiers": lambda x: concatenate_with_delim("|",x),
    "gene_models": lambda x: concatenate_with_delim("|",x)
})


# In[8]:


# This is only called by collapse_by_all_gene_names().
# A method necessary for cleaning up lists of gene identifiers after merging.
# This removes things from the other gene identifiers if they are already listed as a unique gene identifier.
# This could happen after merging if some string was unsure about being a unique identifier, but some other entry confirms that is is.
def remove_duplicate_names(row):
    gene_names = row["unique_gene_identifiers"].split("|")
    gene_synonyms = row["other_gene_identifiers"].split("|")
    updated_gene_synonyms = [x for x in gene_synonyms if x not in gene_names]
    gene_synonyms_str = concatenate_with_delim("|", updated_gene_synonyms)
    return(gene_synonyms_str)


# This is only called by collapse_by_all_gene_names().
# Another method necessary for cleaning up lists of gene identifiers after merging.
# This retains the order except for it puts anything that is also in the gene models column last.
def reorder_unique_gene_identifers(row):
    unique_identifiers = row["unique_gene_identifiers"].split("|")
    gene_models = row["gene_models"].split("|")
    reordered_unique_identifiers = [x for x in unique_identifiers if x not in gene_models]
    reordered_unique_identifiers.extend(gene_models)
    reordered_unique_identifiers_str = concatenate_with_delim("|", reordered_unique_identifiers)
    return(reordered_unique_identifiers_str)

  
agg_df["other_gene_identifiers"] = agg_df.apply(lambda x: remove_duplicate_names(x), axis=1)
agg_df["unique_gene_identifiers"] = agg_df.apply(lambda x: reorder_unique_gene_identifers(x), axis=1)


# In[9]:


cols_to_retain_from_old_df = ["_gene_id",
                              "species_name",
                              "species_code",
                              "text_unprocessed",
                              "annotations",
                              "reference_name",
                              "reference_link",
                              "reference_file"]
new_df = df[cols_to_retain_from_old_df].merge(right=agg_df, on="_gene_id", how="left")
new_df.head(10)


# In[10]:


new_df.shape


# In[11]:


df.shape


# In[ ]:





# In[12]:


new_df["text_unprocessed"].values


# In[13]:


new_df["text_unprocessed"].sample(5).values


# In[14]:


new_df["text_tokenized_sents"] = new_df["text_unprocessed"].map(lambda x: x.replace(";","."), na_action="ignore")
new_df["text_tokenized_sents"] = new_df["text_tokenized_sents"].map(sent_tokenize, na_action="ignore")
f = lambda sents: " ".join(["[SENT] {}".format(s) for s in sents])
new_df["text_tokenized_sents"] = new_df["text_tokenized_sents"].map(f, na_action="ignore")


# In[15]:


new_df["test"] = new_df["text_tokenized_sents"].map(lambda x: x.split("[SENT]"), na_action="ignore")
new_df["test"].sample(5).values


# 
# 
# 
# def preprocess(text):
#     sents = text.split("[SENT]")
#     sents = [" ".join(word_tokenize(s)) for s in sents]
# 
#     a = " [SENT] ".join(sents)
#     return(a)
#     
# 
# new_df["test"] = new_df["text_sent_tokenized"].map(preprocess, na_action="ignore")
# new_df["test"].sample(5).values                    
#                       
#                       

# for b in new_df["test"].sample(5).values:
#     print(b)
#     print()

# In[16]:


# The input should be one string with sentences separated by something
# This stems each word, removes puncutation and also all of the stopwords and lowercases.
def preprocess_sentences_full(text, sentence_delimiter):
    sentences = text.split(sentence_delimiter)
    sentences = [" ".join(preprocess_string(s)) for s in sentences]
    reformatted_text = " {} ".format(sentence_delimiter).join(sentences)
    reformatted_text = reformatted_text.strip()
    return(reformatted_text)


# The input should be one string with sentences separated by something
# This splits the strings into tokens but leaves the content of them alone.
def preprocess_sentences_partial(text, sentence_delimiter):
    sentences = text.split(sentence_delimiter)
    sentences = [" ".join(word_tokenize(s)) for s in sentences]
    reformatted_text = " {} ".format(sentence_delimiter).join(sentences)
    reformatted_text = reformatted_text.strip()
    return(reformatted_text)


    
SENT_DELIMITER = "[SENT]"
new_df["text_tokenized_stems"] = new_df["text_tokenized_sents"].map(lambda x: preprocess_sentences_full(x, SENT_DELIMITER), na_action="ignore")
new_df["text_tokenized_words"] = new_df["text_tokenized_sents"].map(lambda x: preprocess_sentences_partial(x, SENT_DELIMITER), na_action="ignore")

new_df[["text_tokenized_stems","text_tokenized_words"]].sample(5).values


# In[17]:


final_column_order = [
    "_gene_id",
    "species_name",
    "species_code",
    "unique_gene_identifiers",
    "other_gene_identifiers",
    "gene_models",
    "annotations",
    "text_unprocessed",
    "text_tokenized_sents",
    "text_tokenized_words",
    "text_tokenized_stems",
    "reference_name",
    "reference_link",
    "reference_file"
]


# In[18]:


df = new_df[final_column_order]
df.sort_values(by="_gene_id", ascending=True, inplace=True, ignore_index=True)
df.drop_duplicates(keep="first", inplace=True, ignore_index=True)
df.head(100)


# In[19]:


# Saving the full versions of the combined datasets.
#csv_path = "../final_data/genes_texts_annotsasdfa.csv"
#tsv_path = "../final_data/genes_texts_annotsasdfasdf.tsv"
#df.to_csv(tsv_path, sep="\t", index=False)
#df.to_csv(csv_path, index=False)


# ### Running NOBLE Coder on the text columns

# In[20]:


noblecoder_jarfile_path = "../lib/NobleCoder-1.0.jar"  


# In[21]:


# Small test for paths.
#texts = {1: "small plants had elongated leaves", 2:"short plants had dwarfism leaves"}
#annots = annotate_using_noble_coder(texts, noblecoder_jarfile_path, "pato", precise=1)
#annots


# In[22]:


# Create a mapping between the lines that include text and the row indices.
index_to_text = dict(zip(df[df["text_unprocessed"].notnull()].index, df[df["text_unprocessed"].notnull()]["text_unprocessed"]))


# In[23]:


# Create the set of precise NOBLE Coder annotations.
pato_annotations = annotate_using_noble_coder(index_to_text, noblecoder_jarfile_path, "pato", precise=1)
po_annotations = annotate_using_noble_coder(index_to_text, noblecoder_jarfile_path, "po", precise=1)
go_annotations = annotate_using_noble_coder(index_to_text, noblecoder_jarfile_path, "go", precise=1)
print("done running noble coder with precise parameter")


# In[24]:


# Combine those annotations and add them as a column in the dataset.
indices = []
list_of_annotation_lists = []
for index in index_to_text.keys():
    indices.append(index)
    annotations = []
    annotations.extend(pato_annotations[index])
    annotations.extend(po_annotations[index])
    annotations.extend(go_annotations[index])
    annotations_str = concatenate_with_delim("|", annotations)
    list_of_annotation_lists.append(annotations_str)
    
df["annotations_nc"] = np.nan    
df.loc[indices,"annotations_nc"] = list_of_annotation_lists
annotations_nc_col = df.pop("annotations_nc")
df.insert(7, "annotations_nc", annotations_nc_col)
df.head()


# In[25]:


df.sample(10)


# In[26]:


print(df.shape)


# In[27]:


# We don't need the genes that have ontology annotations but no text.
gene_ids_with_text = df[df["text_unprocessed"].notnull()]["_gene_id"].values
df = df[df["_gene_id"].isin(gene_ids_with_text)]
print(df.shape)


# ### Saving the combined datasets to new files

# In[28]:


# Saving the full versions of the combined datasets.
csv_path = "../final_data/genes_texts_annotations.csv"
tsv_path = "../final_data/genes_texts_annotations.tsv"
df.to_csv(tsv_path, sep="\t", index=False)
df.to_csv(csv_path, index=False)


# In[31]:


# Saving sample versions that should be viewable in the browser on GitHub.
# Prepare a sample file from that whole dataset that makes it easy to understand what the context is.

# Function to truncate strings for more readable sample files.
def truncate_string(text, char_limit):
    truncated_text = text[:char_limit]
    if len(text)>char_limit:
        truncated_text = "{}...".format(truncated_text)
    return(truncated_text)

def truncate_fields(sample_df):
    sample_df["unique_gene_identifiers"] = sample_df["unique_gene_identifiers"].map(lambda x: truncate_string(x, 30), na_action="ignore")
    sample_df["other_gene_identifiers"] = sample_df["other_gene_identifiers"].map(lambda x: truncate_string(x, 20), na_action="ignore")
    sample_df["gene_models"] = sample_df["gene_models"].map(lambda x: truncate_string(x, 30), na_action="ignore")
    sample_df["text_unprocessed"] = sample_df["text_unprocessed"].map(lambda x: truncate_string(x, 100), na_action="ignore")
    sample_df["text_tokenized_sents"] = sample_df["text_tokenized_sents"].map(lambda x: truncate_string(x, 100), na_action="ignore")
    sample_df["text_tokenized_words"] = sample_df["text_tokenized_words"].map(lambda x: truncate_string(x, 100), na_action="ignore")
    sample_df["text_tokenized_stems"] = sample_df["text_tokenized_stems"].map(lambda x: truncate_string(x, 100), na_action="ignore")
    sample_df["annotations"] = sample_df["annotations"].map(lambda x: truncate_string(x, 60), na_action="ignore")
    sample_df["annotations_nc"] = sample_df["annotations_nc"].map(lambda x: truncate_string(x, 60), na_action="ignore")
    return(sample_df)


csv_sample_path = "../final_samples/genes_texts_annotations.csv"
tsv_sample_path = "../final_samples/genes_texts_annotations.tsv"

# Taking only the first few rows and truncating values in some columns.
sample_df = df.head(100)
sample_df = truncate_fields(sample_df)
sample_df.to_csv(tsv_sample_path, sep="\t", index=False)
sample_df.to_csv(csv_sample_path, index=False)
print("done")


# ### Saving smaller subsets of the dataset to new files

# In[37]:


# Saving files for only the text fields and corresponding annotations.
df_subset = df[df["text_unprocessed"].notnull()]
subset_csv_path = "../final_data/genes_texts.csv"
subset_tsv_path = "../final_data/genes_texts.tsv"
df_subset.to_csv(subset_tsv_path, sep="\t", index=False)
df_subset.to_csv(subset_csv_path, index=False)

# Saving the sample versions.
sample_df_subset = df_subset.head(100)
sample_df_subset = truncate_fields(sample_df_subset)
cols_to_drop = ["annotations"]
sample_df_subset.drop(cols_to_drop, axis="columns", inplace=True)
subset_csv_path = "../final_samples/genes_texts.csv"
subset_tsv_path = "../final_samples/genes_texts.tsv"
sample_df_subset.to_csv(subset_tsv_path, sep="\t", index=False)
sample_df_subset.to_csv(subset_csv_path, index=False)


# In[38]:


# Saving files for only the annotation fields not the text ones.
df_subset = df[df["annotations"].notnull()]
subset_csv_path = "../final_data/genes_texts.csv"
subset_tsv_path = "../final_data/genes_texts.tsv"
df_subset.to_csv(subset_tsv_path, sep="\t", index=False)
df_subset.to_csv(subset_csv_path, index=False)

# Saving the sample versions.
sample_df_subset = df_subset.head(100)
sample_df_subset = truncate_fields(sample_df_subset)
cols_to_drop = ["annotations_nc","text_unprocessed", "text_tokenized_sents","text_tokenized_words","text_tokenized_stems"]
sample_df_subset.drop(cols_to_drop, axis="columns", inplace=True)
subset_csv_path = "../final_samples/genes_annotations.csv"
subset_tsv_path = "../final_samples/genes_annotations.tsv"
sample_df_subset.to_csv(subset_tsv_path, sep="\t", index=False)
sample_df_subset.to_csv(subset_csv_path, index=False)


# In[39]:


print("done with combining all files")

