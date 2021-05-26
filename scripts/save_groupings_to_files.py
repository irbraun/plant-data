import sys
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter('ignore')


sys.path.append("../../oats")
from oats.biology.groupings import Groupings
from oats.utils.utils import save_to_pickle
from oats.nlp.preprocess import replace_delimiter, concatenate_with_delim




# Paths to input files from databases or papers. Some others are specified below in the dictionaries as well.
lloyd_meinke_cleaned_supplemental_table_path_hierarchy = "../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S1_Final.csv"
lloyd_meinke_cleaned_supplemental_table_path_mappings = "../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S2_Final_Revised.csv"



# Paths to the csv files that are created for each type of data for grouping genes.
lloyd_meinke_subsets_output_path = "../reshaped_data/lloyd_meinke_subsets.csv"
lloyd_meinke_classes_output_path = "../reshaped_data/lloyd_meinke_classes.csv"
kegg_pathways_output_path = "../reshaped_data/kegg_pathways.csv"
plantcyc_pathways_output_path = "../reshaped_data/plantcyc_pathways.csv"



# Paths to the csv files that are created to specify the mappings between IDs and names for each group.
lloyd_meinke_subsets_name_mapping_path = "../reshaped_data/lloyd_meinke_subsets_name_map.csv"
lloyd_meinke_classes_name_mapping_path = "../reshaped_data/lloyd_meinke_classes_name_map.csv"
kegg_pathways_name_mapping_path = "../reshaped_data/kegg_pathways_name_map.csv"
plantcyc_pathways_name_mapping_path = "../reshaped_data/plantcyc_pathways_name_map.csv"







# PlantCyc


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
plantcyc_name_mapping = {row.pathway_id:row.pathway_name for row in plantcyc_df.itertuples()}
pd.DataFrame(plantcyc_name_mapping.items(), columns=["group_id","group_name"]).to_csv(plantcyc_pathways_name_mapping_path, index=False)











# KEGG


# Mapping between species codes and files saved using another script that uses the KEGG REST API.
kegg_paths_dictionary = {
    "ath":"../databases/kegg/ath_pathway_files_from_api",
    "zma":"../databases/kegg/zma_pathway_files_from_api",
    "osa":"../databases/kegg/osa_pathway_files_from_api",
    "mtr":"../databases/kegg/mtr_pathway_files_from_api",
    "gmx":"../databases/kegg/gmx_pathway_files_from_api",
    "sly":"../databases/kegg/sly_pathway_files_from_api",
    "hsa":"../kegg/hsa_pathway_files_from_api",
}
kegg_df = Groupings.get_dataframe_for_kegg(paths=kegg_paths_dictionary)
kegg_df.to_csv(kegg_pathways_output_path, index=False)
kegg_name_mapping = {row.pathway_id:row.pathway_name for row in kegg_df.itertuples()}
pd.DataFrame(kegg_name_mapping.items(), columns=["group_id","group_name"]).to_csv(kegg_pathways_name_mapping_path, index=False)










# Lloyd and Meinke et al., 2012


# Some preprocessing on the supplemental file from Lloyd and Meinke, 2012 paper to extrac the columns used.
df = pd.read_csv(lloyd_meinke_cleaned_supplemental_table_path_mappings)
df.fillna("", inplace=True)
combine_columns = lambda row, columns: concatenate_with_delim("|", [row[column] for column in columns])
df["Alias Symbols"] = df["Alias Symbols"].apply(lambda x: replace_delimiter(text=x, old_delim=";", new_delim="|"))
df["gene_identifiers"] = df.apply(lambda x: combine_columns(x, ["Locus", "Gene Symbol", "Alias Symbols", "Full Gene Name"]), axis=1)

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
df_subset["group_ids"] = df_subset["group_ids"].apply(lambda x: replace_delimiter(text=x, old_delim=",", new_delim="|"))
df_subset.to_csv(lloyd_meinke_subsets_output_path, index=False)



# Provide a mapping from subset or class IDs to the longer names that define them.
df = pd.read_csv(lloyd_meinke_cleaned_supplemental_table_path_hierarchy)
subset_id_to_name_dict = {row[5]:row[7] for row in df.itertuples()}
class_id_to_name_dict = {row[3]:row[4] for row in df.itertuples()}
pd.DataFrame(subset_id_to_name_dict.items(), columns=["group_id","group_name"]).to_csv(lloyd_meinke_subsets_name_mapping_path, index=False)
pd.DataFrame(class_id_to_name_dict.items(), columns=["group_id","group_name"]).to_csv(lloyd_meinke_classes_name_mapping_path, index=False)












# Briefly checking whether groupings object can be successfully built from the created files.
# Create actual oats.Grouping objects using those CSV files that were created previously, and quick check of the contents.
print(Groupings(path=lloyd_meinke_subsets_output_path, name_mapping=subset_id_to_name_dict).describe())
print(Groupings(path=lloyd_meinke_classes_output_path, name_mapping=class_id_to_name_dict).describe())
print(Groupings(path=kegg_pathways_output_path, name_mapping=kegg_name_mapping).describe())
print(Groupings(path=plantcyc_pathways_output_path, name_mapping=plantcyc_name_mapping).describe())











