import sys
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter('ignore')


sys.path.append("../../oats")
from oats.biology.groupings import Groupings
from oats.utils.utils import save_to_pickle
from oats.nlp.preprocess import other_delim_to_bar_delim, concatenate_with_bar_delim




# Paths to input files.
lloyd_meinke_cleaned_supplemental_table_path = "../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S2_Final_Revised.csv"



# Paths to the csv files that are created for each type of data for grouping genes.
lloyd_meinke_subsets_output_path = "../reshaped_data/lloyd_meinke_subsets.csv"
lloyd_meinke_classes_output_path = "../reshaped_data/lloyd_meinke_classes.csv"
kegg_pathways_output_path = "../reshaped_data/kegg_pathways.csv"
plantcyc_pathways_output_path = "../reshaped_data/plantcyc_pathways.csv"








# Mapping between species codes and files downloaded from PlantCyc.
plantcyc_paths_dictionary = {
    "ath":"../databases/plantcyc/aracyc_pathways.20180702", 
    "zma":"../databases/plantcyc/corncyc_pathways.20180702", 
    "mtr":"../databases/plantcyc/mtruncatulacyc_pathways.20180702", 
    "osa":"../databases/plantcyc/oryzacyc_pathways.20180702", 
    "gmx":"../databases/plantcyc/soycyc_pathways.20180702",
    "sly":"../databases/plantcyc/tomatocyc_pathways.20180702"}

# Create and save the pathways object using PlantCyc.
plantcyc_df = Groupings.get_dataframe_for_plantcyc(paths=plantcyc_paths_dictionary)
plantcyc_df.to_csv(plantcyc_pathways_output_path, index=False)







# Mapping between species codes and files saved using another script that uses the KEGG REST API.
kegg_paths_dictionary = {
    "ath":"/Users/irbraun/plant-data/databases/kegg/ath_pathway_files_from_api",
    "zma":"/Users/irbraun/plant-data/databases/kegg/zma_pathway_files_from_api",
    "osa":"/Users/irbraun/plant-data/databases/kegg/osa_pathway_files_from_api",
    "mtr":"/Users/irbraun/plant-data/databases/kegg/mtr_pathway_files_from_api",
    "gmx":"/Users/irbraun/plant-data/databases/kegg/gmx_pathway_files_from_api",
    "sly":"/Users/irbraun/plant-data/databases/kegg/sly_pathway_files_from_api",
    "hsa":"/Users/irbraun/plant-data/databases/kegg/hsa_pathway_files_from_api",
}
kegg_df = Groupings.get_dataframe_for_kegg(paths=kegg_paths_dictionary)
kegg_df.to_csv(kegg_pathways_output_path, index=False)




# Some preprocessing on the supplemental file from Lloyd and Meinke, 2012 paper to extrac the columns used.
df = pd.read_csv(lloyd_meinke_cleaned_supplemental_table_path)
df.fillna("", inplace=True)
df["Alias Symbols"] = df["Alias Symbols"].apply(lambda x: other_delim_to_bar_delim(string=x, delim=";"))
df["gene_identifiers"] = np.vectorize(concatenate_with_bar_delim)(df["Locus"], df["Gene Symbol"], df["Alias Symbols"], df["Full Gene Name"])

# Specific to classes (more general).
df_class = df[["Phenotype Classb", "gene_identifiers"]]
df_class["species"] = "ath"
df_class.columns = ["group_ids", "gene_identifiers","species"]
df_class = df_class[["species", "group_ids", "gene_identifiers"]]
df_class.to_csv(lloyd_meinke_classes_output_path, index=False)

# Specific to subsets (more specific).
df_subset = df[["Phenotype Subsetsb", "gene_identifiers"]]
df_subset["species"] = "ath"
df_subset.columns = ["group_ids", "gene_identifiers","species"]
df_subset = df_subset[["species", "group_ids", "gene_identifiers"]]
df_subset["group_ids"] = df_subset["group_ids"].apply(lambda x: x.replace("W:", "").replace("S:","").replace("(",",").replace(")",",").replace(";",","))
df_subset["group_ids"] = df_subset["group_ids"].apply(lambda x: other_delim_to_bar_delim(string=x, delim=","))
df_subset.to_csv(lloyd_meinke_subsets_output_path, index=False)



# Provide a mapping from subset or class IDs to the longer names that define them.
df = pd.read_csv("../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S1_Final.csv")
subset_id_to_name_dict = {row[5]:row[7] for row in df.itertuples()}
class_id_to_name_dict = {row[3]:row[4] for row in df.itertuples()}













# Create actual oats.Grouping objects using those CSV files that were created previously, and quick check of the contents.
print(Groupings(lloyd_meinke_subsets_output_path).describe())
print(Groupings(lloyd_meinke_classes_output_path).describe())
print(Groupings(kegg_pathways_output_path).describe())
print(Groupings(plantcyc_pathways_output_path).describe())











