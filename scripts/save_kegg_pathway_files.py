sys.path.append("../../oats")
from oats.biology.groupings import Groupings


species_to_directories = {
	"ath":"/Users/irbraun/plant-data/databases/kegg/ath_pathway_files_from_api",
	"zma":"/Users/irbraun/plant-data/databases/kegg/zma_pathway_files_from_api",
	"osa":"/Users/irbraun/plant-data/databases/kegg/osa_pathway_files_from_api",
	"mtr":"/Users/irbraun/plant-data/databases/kegg/mtr_pathway_files_from_api",
	"gmx":"/Users/irbraun/plant-data/databases/kegg/gmx_pathway_files_from_api",
	"sly":"/Users/irbraun/plant-data/databases/kegg/sly_pathway_files_from_api",
	"hsa":"/Users/irbraun/plant-data/databases/kegg/hsa_pathway_files_from_api",
}

Groupings._save_all_kegg_pathway_files(species_to_directories)
