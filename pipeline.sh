# Convert the preprocessing notebooks to python scripts.
jupyter nbconvert --to script preprocessing/reshaping_maizegdb_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_sgn_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_tair_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_oellrich_walls_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/reshaping_planteome_data.ipynb --output-dir scripts
jupyter nbconvert --to script preprocessing/combining.ipynb --output-dir scripts



# Run all the preprocessing scripts to generate files in the reshaped data directory.
cd scripts
python reshaping_maizegdb_data.py
python reshaping_tair_data.py
python reshaping_sgn_data.py
python reshaping_oellrich_walls_data.py
python reshaping_planteome_data.py



# Combine all those reshaped files into the combined dataset files for the genes, annotations, and phenotypes.
python combining.py #(used to be save_datasaets_to_files.py)

# Uncomment these lines again when ready to use or when replacing them with other steps.
#python save_datasets_to_json.py
#python check_dataset.py


# Create the files that map gene objects from the above resources to particular groupings.
# Commenting out running the script that first extracts the KEGG pathway files using the REST API, just retain those files.
# This way this pipeline can be stable and run without any of the underlying files changing, this way it is reproducible.
# python save_kegg_pathways_to_files.py


#uncomment this again when added
#python save_groupings_to_files.py


# Create samples of the first few lines of each of the created files into a samples directory for looking at the shape.
cd ..
for path in ./reshaped_data/*.csv; do
	filename=$(basename $path)
	head -100 $path > ./reshaped_samples/$filename
done