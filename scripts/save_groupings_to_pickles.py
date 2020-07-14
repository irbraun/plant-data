import sys
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter('ignore')


sys.path.append("../../oats")
from oats.biology.groupings import Groupings
from oats.utils.utils import save_to_pickle
from oats.nlp.preprocess import other_delim_to_bar_delim, concatenate_with_bar_delim




# Mapping between species codes and relevant files, used in searching pathway databases.
species_dict = {
    "ath":"../databases/plantcyc/aracyc_pathways.20180702", 
    "zma":"../databases/plantcyc/corncyc_pathways.20180702", 
    "mtr":"../databases/plantcyc/mtruncatulacyc_pathways.20180702", 
    "osa":"../databases/plantcyc/oryzacyc_pathways.20180702", 
    "gmx":"../databases/plantcyc/soycyc_pathways.20180702",
    "sly":"../databases/plantcyc/tomatocyc_pathways.20180702"}

# Create and save the pathways object using PMN. This one uses the filenames in the dictionary.
pathways = Groupings(species_dict, source="pmn")
pathways.describe()
save_to_pickle(obj=pathways, path="../pickles/groupings_from_pmn_pathways.pickle")









# Mapping between species codes and relevant files, used in searching pathway databases.
# Also adding homo sapiens as a key here. 
# Note that the values in this dictionary are not actually used when the source used is KEGG as is the case here.
species_dict = {
    "ath":"", 
    "zma":"", 
    "mtr":"", 
    "osa":"", 
    "gmx":"",
    "sly":"",
    "hsa":""}

# Create and save a pathways object using KEGG. This one doesn't use the filenames in the dictionary.
# All the pathways asocations are looked up when the groupings object is created through the KEGG API.
pathways = Groupings(species_dict, source="kegg")
pathways.describe()
save_to_pickle(obj=pathways, path="../pickles/groupings_from_kegg_pathways.pickle")









# Some preprocessing necessary for the Lloyd three-level (group, class, subset) functional hierarchy dataset.
df = pd.read_csv("../papers/lloyd_meinke_2012/versions_cleaned_by_me/lloyd_gene_to_function_irb_cleaned.csv")
df.fillna("", inplace=True)
df["Alias Symbols"] = df["Alias Symbols"].apply(lambda x: other_delim_to_bar_delim(string=x, delim=";"))
df["gene_names"] = np.vectorize(concatenate_with_bar_delim)(df["Locus"], df["Gene Symbol"], df["Alias Symbols"], df["Full Gene Name"])
# Specific to classes (more general).
df_class = df[["Phenotype Classb", "gene_names"]]
df_class.columns = ["group_id", "gene_names"]
df_class.to_csv("../papers/lloyd_meinke_2012/versions_cleaned_by_me/arabidopsis_classes_irb_cleaned.csv", index=False)
# Specific to subsets (more specific).
df_subset = df[["Phenotype Subsetsb", "gene_names"]]
df_subset.columns = ["group_id", "gene_names"]
df_subset["group_id"] = df_subset["group_id"].apply(lambda x: x.replace("W:", "").replace("S:","").replace("(",",").replace(")",",").replace(";",","))
df_subset["group_id"] = df_subset["group_id"].apply(lambda x: other_delim_to_bar_delim(string=x, delim=","))
df_subset.to_csv("../papers/lloyd_meinke_2012/versions_cleaned_by_me/arabidopsis_subsets_irb_cleaned.csv", index=False)
# Provide a mapping from subset or class IDs to the longer names that define them.
df = pd.read_csv("../papers/lloyd_meinke_2012/versions_cleaned_by_me/lloyd_function_hierarchy_irb_cleaned.csv")
subset_id_to_name_dict = {row[5]:row[7] for row in df.itertuples()}
class_id_to_name_dict = {row[3]:row[4] for row in df.itertuples()}




# Create the save the grouping object using Lloyd function hierarchy dataset of subsets.
species_dict = {"ath":"../papers/lloyd_meinke_2012/versions_cleaned_by_me/arabidopsis_subsets_irb_cleaned.csv"}
subsets = Groupings(species_dict, source="csv", name_mapping=subset_id_to_name_dict)
subsets.describe()
save_to_pickle(obj=subsets, path="../pickles/groupings_from_lloyd_subsets.pickle")




# Create the save the grouping object using Lloyd function hierarchy dataset of classes.
species_dict = {"ath":"../papers/lloyd_meinke_2012/versions_cleaned_by_me/arabidopsis_classes_irb_cleaned.csv"}
classes = Groupings(species_dict, source="csv", name_mapping=class_id_to_name_dict)
classes.describe()
save_to_pickle(obj=classes, path="../pickles/groupings_from_lloyd_classes.pickle")






