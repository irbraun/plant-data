import sys
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter('ignore')


sys.path.append("../")
from utils.constants import ABBREVIATIONS_MAP


sys.path.append("../../oats")
from oats.biology.groupings import Groupings
from oats.utils.utils import save_to_pickle
from oats.nlp.preprocess import replace_delimiter, concatenate_with_delim






# Specifying input and output paths




# Paths to input files from databases or papers. Some others are specified below in the dictionaries as well.
lloyd_meinke_cleaned_supplemental_table_path_hierarchy = "../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S1_Final.csv"
lloyd_meinke_cleaned_supplemental_table_path_mappings = "../papers/lloyd_meinke_2012/versions_cleaned_by_me/192393Table_S2_Final_Revised.csv"


# Paths to the csv files that are created for each type of data for grouping genes.
lloyd_meinke_subsets_output_path = "../reshaped/data/lloyd_meinke_subsets.csv"
lloyd_meinke_classes_output_path = "../reshaped/data/lloyd_meinke_classes.csv"
kegg_pathways_output_path = "../reshaped/data/kegg_pathways.csv"
plantcyc_pathways_output_path = "../reshaped/data/plantcyc_pathways.csv"

# Paths to the csv files that are created to specify the mappings between IDs and names for each group.
lloyd_meinke_subsets_name_mapping_path = "../reshaped/data/lloyd_meinke_subsets_name_map.csv"
lloyd_meinke_classes_name_mapping_path = "../reshaped/data/lloyd_meinke_classes_name_map.csv"
kegg_pathways_name_mapping_path = "../reshaped/data/kegg_pathways_name_map.csv"
plantcyc_pathways_name_mapping_path = "../reshaped/data/plantcyc_pathways_name_map.csv"


# Paths to the final csv file that contains all the mapping information here.
final_groupings_path = "../final/data/groupings.csv"
final_groupings_sample_path = "../final/samples/groupings.csv"















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





# Putting the plantcyc data into a format that is standardized across the datasets.
plantcyc_df = plantcyc_df[["species", "pathway_id", "pathway_name", "gene_identifiers"]]
plantcyc_df["gene_identifiers"] = plantcyc_df["gene_identifiers"].map(lambda x: x.split("|"))
plantcyc_df = plantcyc_df.explode(column="gene_identifiers").reset_index(drop=True)
plantcyc_df = plantcyc_df.drop_duplicates(inplace=False, ignore_index=True)
plantcyc_df.columns = ["species_code", "group_id", "group_name", "gene_identifier"]
plantcyc_df["species_name"] = plantcyc_df["species_code"].map({v:k for k,v in ABBREVIATIONS_MAP.items()})
plantcyc_df["reference_name"] = "PlantCyc"
plantcyc_df["reference_file"] = plantcyc_df["species_code"].map(lambda x: plantcyc_paths_dictionary[x].split("/")[-1])
plantcyc_df["reference_link"] = "https://plantcyc.org/downloads"
plantcyc_df["group_type"] = "biochemical pathway"
plantcyc_df = plantcyc_df[["species_name",
         "species_code",
         "group_id",
         "group_name",
         "group_type",
         "gene_identifier",
         "reference_name",
         "reference_file",
         "reference_link"
         ]]


















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
#kegg_df = pd.read_csv(kegg_pathways_output_path)



# Putting the KEGG mappings dataframe into a shape that is standard across the datasets.
kegg_df = kegg_df[["species", "pathway_id", "pathway_name", "gene_identifiers"]]
kegg_df["gene_identifiers"] = kegg_df["gene_identifiers"].map(lambda x: x.split("|"))
kegg_df = kegg_df.explode(column="gene_identifiers").reset_index(drop=True)
kegg_df = kegg_df.drop_duplicates(inplace=False, ignore_index=True)
kegg_df["pathway_name"] = kegg_df["pathway_name"].map(lambda x: "".join(x.split("-")[0:-1]))
kegg_df.columns = ["species_code", "group_id", "group_name", "gene_identifier"]
kegg_df["species_name"] = kegg_df["species_code"].map({v:k for k,v in ABBREVIATIONS_MAP.items()})
kegg_df["reference_name"] = "KEGG"
kegg_df["reference_file"] = kegg_df["species_code"].map(lambda x: "{}/path_{}[].txt".format(kegg_paths_dictionary[x].split("/")[-1], x))
kegg_df["reference_link"] = "https://www.genome.jp/kegg/"
kegg_df["group_type"] = "biochemical pathway"
kegg_df = kegg_df[["species_name",
         "species_code",
         "group_id",
         "group_name",
         "group_type",
         "gene_identifier",
         "reference_name",
         "reference_file",
         "reference_link"
         ]]



















# Lloyd and Meinke, 2012




# Provide a mapping from subset or class IDs to the longer names that define them.
df = pd.read_csv(lloyd_meinke_cleaned_supplemental_table_path_hierarchy)
subset_id_to_name_dict = {row[5]:row[7] for row in df.itertuples()}
class_id_to_name_dict = {row[3]:row[4] for row in df.itertuples()}
pd.DataFrame(subset_id_to_name_dict.items(), columns=["group_id","group_name"]).to_csv(lloyd_meinke_subsets_name_mapping_path, index=False)
pd.DataFrame(class_id_to_name_dict.items(), columns=["group_id","group_name"]).to_csv(lloyd_meinke_classes_name_mapping_path, index=False)








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



# Putting the dataframe for the phenotype classes into a format that is standard across the data.
df_class = df_class[["species", "group_ids", "gene_identifiers"]]
df_class["gene_identifiers"] = df_class["gene_identifiers"].map(lambda x: x.split("|"))
df_class = df_class.explode(column="gene_identifiers").reset_index(drop=True)
df_class["group_ids"] = df_class["group_ids"].map(lambda x: x.split("|"))
df_class = df_class.explode(column="group_ids").reset_index(drop=True)
df_class = df_class.drop_duplicates(inplace=False, ignore_index=True)
df_class = df_class[["species", "group_ids", "gene_identifiers"]]
df_class.columns = ["species_code", "group_id", "gene_identifier"]
df_class["species_name"] = df_class["species_code"].map({v:k for k,v in ABBREVIATIONS_MAP.items()})
df_class["group_name"] = df_class["group_id"].map(class_id_to_name_dict)
df_class["reference_name"] = "Lloyd and Meinke, 2012"
df_class["reference_file"] = "192393Table_S2_Final_Revised.xls"
df_class["reference_link"] = "http://www.plantphysiol.org/highwire/filestream/122682/field_highwire_adjunct_files/2/192393Table_S2_Final_Revised.xls"
df_class["group_type"] = "phenotype class"
df_class = df_class[["species_name",
         "species_code",
         "group_id",
         "group_name",
         "group_type",
         "gene_identifier",
         "reference_name",
         "reference_file",
         "reference_link"
         ]]




# Specific to subsets (more specific).
df_subset = df[["Phenotype Subsetsb", "gene_identifiers"]]
df_subset["species"] = "ath"
df_subset.columns = ["group_ids", "gene_identifiers","species"]
df_subset = df_subset[["species", "group_ids", "gene_identifiers"]]
df_subset["group_ids"] = df_subset["group_ids"].apply(lambda x: x.replace("W:", "").replace("S:","").replace("(",",").replace(")",",").replace(";",","))
df_subset["group_ids"] = df_subset["group_ids"].apply(lambda x: replace_delimiter(text=x, old_delim=",", new_delim="|"))
df_subset.to_csv(lloyd_meinke_subsets_output_path, index=False)


# Putting the dataframe for the phenotype subsets into a standardized format across the data.
df_subset = df_subset[["species", "group_ids", "gene_identifiers"]]
df_subset["gene_identifiers"] = df_subset["gene_identifiers"].map(lambda x: x.split("|"))
df_subset = df_subset.explode(column="gene_identifiers").reset_index(drop=True)
df_subset["group_ids"] = df_subset["group_ids"].map(lambda x: x.split("|"))
df_subset = df_subset.explode(column="group_ids").reset_index(drop=True)
df_subset = df_subset.drop_duplicates(inplace=False, ignore_index=True)
df_subset = df_subset[["species", "group_ids", "gene_identifiers"]]
df_subset.columns = ["species_code", "group_id", "gene_identifier"]
df_subset["species_name"] = df_subset["species_code"].map({v:k for k,v in ABBREVIATIONS_MAP.items()})
df_subset["group_name"] = df_subset["group_id"].map(subset_id_to_name_dict)
df_subset["reference_name"] = "Lloyd and Meinke, 2012"
df_subset["reference_file"] = "192393Table_S2_Final_Revised.xls"
df_subset["reference_link"] = "http://www.plantphysiol.org/highwire/filestream/122682/field_highwire_adjunct_files/2/192393Table_S2_Final_Revised.xls"
df_subset["group_type"] = "phenotype subset"
df_subset = df_subset[["species_name",
         "species_code",
         "group_id",
         "group_name",
         "group_type",
         "gene_identifier",
         "reference_name",
         "reference_file",
         "reference_link"
         ]]














# Combine all four of the dataframes created above and put them into one file and save.
final_df = pd.concat([df_subset, df_class, kegg_df, plantcyc_df])
final_df.to_csv(final_groupings_path, index=False)
final_df.sample(100).sort_values(by="species_code").to_csv(final_groupings_sample_path, index=False)













# Briefly checking whether groupings object can be successfully built from the created files.
# Create actual oats.Grouping objects using those CSV files that were created previously, and quick check of the contents.
print(Groupings(path=lloyd_meinke_subsets_output_path, name_mapping=subset_id_to_name_dict).describe())
print(Groupings(path=lloyd_meinke_classes_output_path, name_mapping=class_id_to_name_dict).describe())
print(Groupings(path=kegg_pathways_output_path, name_mapping=kegg_name_mapping).describe())
print(Groupings(path=plantcyc_pathways_output_path, name_mapping=plantcyc_name_mapping).describe())











