# Convert the preprocessing notebooks to python scripts.
jupyter nbconvert --to script preprocessing/reshaping_maizegdb_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_sgn_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_tair_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_oryzabase_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_oellrich_walls_data.ipynb --output-dir scripts


# Run all the preprocessing scripts to generate files in the reshaped data directory.
cd scripts
python reshaping_maizegdb_data.py
python reshaping_oryzabase_data.py
python reshaping_tair_data.py
python reshaping_sgn_data.py
python reshaping_oellrich_walls_data.py


# Combine all those reshaped files into the combined dataset files for the genes, annotations, and phenotypes.
python save_datasets_to_files.py


# Create the files that map gene objects from the above resources to particular groupings.
# Commenting out running the script that first extracts the KEGG pathway files using the REST API, just retain those files.
# This way this pipeline can be stable and run without any of the underlying files changing, this way it is reproducible.
# python save_kegg_pathways_to_files.py
python save_groupings_to_files.py