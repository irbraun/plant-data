## plant-data

### Project goal and description
The purpose of this repository is to document the collection of a dataset that organized information about plant genes, phenotypes, and annotations from a variety of data sources, and contains notebooks and scripts for preprocessing and combining this information. The intended use of the dataset is for looking into how well natural language processing (NLP) approaches work for recovering or predicting biological information from text descriptions of phenotypes, compared to using ontology-based annotations to make equivalent predictions.

### What's here?
1. All the files that were used for this dataset are listed in `file_descriptions.tsv`, which includes links to the original data source where applicable. These come from databases, papers, and other bioinformatics resources. Note that the files listed are not actually present in this repository, because some are only available through database subscriptions or requests. 

2. The `pipeline.sh` scripts runs each script (see `preprocessing` and `scripts` directories) for preprocessing and merging the information present in each of those files. This pipeline generates all the files in the `reshaped_data` directory, and the first few lines of each of those intermediate data files are available in the `reshaped_samples` directory. This is done to take the information from a variety of sources, and represent it with a standard set of columns so that it can be merged into a single dataset. Files that map genes to groups of any kind (e.g., pathways) use columns `species`, `group_ids`, and `gene_identifiers`. Files that specify mappings between those groups and full group names use columns `group_id` and `group_name`. Files that map genes to phenotype descriptions or annotations are use columns `species`, `unique_gene_identifiers`, `other_gene_identifiers`, `gene_models`, `descriptions`, `annotations`, and `sources`.

3. The primary dataset of interest that combines all this information is `genes_text_annots`, as csv, tsv, and json. 

4. A subset of this primary dataset with fewer genes and truncated fields (small enough to view on GitHub) is `genes_texts_annots_sample`, as tsv or json.

### Variable Names
Variable names share meaning across files and are listed and described here.
* `_gene_id`: A unique identifier used internally to reference each gene.
* `species_name`: A string referring to the name of the species for a given gene.
* `species_code`: A [three-letter organism code](https://www.genome.jp/kegg/catalog/org_list.html) identifier for each species.
* `unique_gene_identifiers`: A list of strings uniquely mapped to one particular gene.
* `other_gene_identifiers`: A list of strings that do not necessary (but might) uniquely refer to one gene.
* `gene_models`: A list of strings referencing gene models, a subset of `unique_gene_identifiers`.
* `annotations`: A list of ontology term IDs associated to a gene.
* `annotations_nc`: A list of ontology term IDs found using [NOBLE Coder](http://ties.dbmi.pitt.edu/noble-coder/) to process text descriptions.
* `text_unprocessed`: A string containing phenotype descriptions or similar text associated to a gene.
* `text_tokenized_sents`: The text string with sentence delimiters added but otherwise unprocessed.
* `text_tokenized_words`: The text string additionally split into individual word tokens.
* `text_tokenized_stems`: The text string additionally split into individual stemmed word tokens.
* `reference_name`: A string referring to the resource in which this information or mapping was found.
* `reference_link`: A link to that resource or file, where applicable.
* `reference_file`: A name for the file where the information was found, if applicable.
* `group_id`: A string name or code uniquely referring to a group of some type. 
* `group_name`: A string referring to a longer more specific name for a given group.
* `group_type`: A string stating the type of grouping (e.g., biochemical pathway).



### What does the data look like?

The primary dataset is a list of genes and accompanying fields that include information about the gene, descriptions of phenotypes associated with the gene, and ontology term anotations mapped to the gene.

A sample of the json dataset with truncated fields is shown here.
```
[
    ...
    {
        "_gene_id": 2,
        "species_code": "ath",
        "species_name": "Arabidopsis thaliana",
        "unique_gene_identifiers": ["UBP26", "At3g49600", "AT3G49600"],
        "other_gene_identifiers": ["Ubiquitin-Specific Protease", "SUP32", "ATUBP26", "ubiquitin-specific protease 26"],
        "gene_models": ["At3g49600", "AT3G49600"],
        "text_unprocessed": "Reduced fertility",
        "text_tokenized_sents": "[SENT] Reduced fertility",
        "text_tokenized_words": "[SENT] Reduced fertility",
        "text_tokenized_stems": "[SENT] reduc fertil",
        "annotations": [],
        "annotations_nc": ["PATO:0001997", "PATO:0000274", "PATO:0000587"],
        "reference_name": "Oellrich, Walls et al., 2015",
        "reference_file": "13007_2015_53_MOESM1_ESM.csv",
        "reference_link": "https://plantmethods.biomedcentral.com/articles/10.1186/s13007-015-0053-y"
    },
    ...
]
```

A sample of the tab-separated dataset with truncated fields is shown here. Lists are pipe-delimited.
```
id  species unique_gene_identifiers             other_gene_identifiers      gene_models                         descriptions                                                                                                annotations                                                         sources
...
40  ath     ZIG|Zigzag Stem|SGR4|SHOOT GRA...   VTI11|ATVTI1A|ATVTI1...     At5g39510|AT5G39510                 Abnormal hypocotyl gravitropism. Abnormal inflorescence gravitropism. Abnormal hypocotyl and inflore...     GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009959|P...     Plant PhenomeNET|TAIR
41  ath     SGR2|Shoot Gravitropism|SGR3|S...   T8E3.1|T8E3_1|ATVAM3...     At1g31480|At5g46860|At2g01940|...   Abnormal hypocotyl gravitropism. Abnormal inflorescence stem gravitropism. Late flowering. Male game...     GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009959|P...     Plant PhenomeNET|TAIR
42  ath     TIR1|Transport Inhibitor Respo...   DOC1|UMB1|ASA1|LPR1|...     At3g62980|At1g70560|At3g02260|...   Abnormal hypocotyl formation. Abnormal lateral root formation. Insensitive to auxin. Few lateral roo...     GO:0032502|PO:0020100|PATO:0001236|PATO:0000460|GO:0010311|P...     Plant PhenomeNET|TAIR
43  ath     NEK6|NIMA (Never In Mitosis, G...   ATNEK6|IBO1|NIMA-REL...     At3g44200|AT3G44200                 Abnormal hypocotyl protrusions. Abnormal petiole protrusions. Abnormal trichome branching. Disordere...     PO:0020100|PATO:0000051|PATO:0000460|PO:0020038|PO:0000282|P...     Plant PhenomeNET|TAIR
44  ath     ARG1|Altered Response to Gravi...   ALTERED RESPONSE TO ...     At1g68370|AT1G68370                 Abnormal hypocotyl gravitropism. Abnormal root gravitropism. Abnormal hypocotyl and root gravitropis...     GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009958|P...     Plant PhenomeNET|TAIR
45  ath     PIN5|Pin Formed|At5g16530|AT5G...   PIN-FORMED 5|MQK4.28...     At5g16530|AT5G16530                 Abnormal hypocotyl growth. Abnormal root growth. Insensitive to IAA. Abnormal hypocotyl and root gro...     GO:0040007|PO:0020100|PATO:0001236|PATO:0000460|GO:0048364|P...     Plant PhenomeNET|TAIR
46  ath     PIL1|Phytochrome Interacting F...   PIF1|phytochrome int...     At2g46970|At3g59060|At2g20180|...   Abnormal cotyledon size under far-red light. Abnormal cotyledon size under red light. Abnormal hypoc...     PO:0020030|PATO:0000117|PATO:0000460|EO:0007203|EO:0007207|P...     Plant PhenomeNET|TAIR
47  ath     GPAT5|Glycerol-3-Phosphate Acy...   ATGPAT1|sn-2-GPAT1|G...     At3g11430|At2g38110|At1g06520|...   Abnormal insoluble lipid polyester biosynthesis. Complete loss of cuticle on petals. Complete loss o...     GO:0009058|CHEBI:60160|PATO:0001236|PATO:0000460|PO:0009032|...     Plant PhenomeNET|TAIR
48  ath     INO|Inner No Outer|At1g23420|A...   INNER NO OUTER|F26F2...     At1g23420|AT1G23420                 Abnormal integuments. Strong allele. Ovule defective. Outer integument fails to initiate on the abax...     PO:0020021|PATO:0000051|PATO:0000460|GO:0048481|GO:0005515|G...     Plant PhenomeNET|TAIR
49  ath     BEL1|Bell-Shaped Ovules|At5g41...   BELL 1|MYC6.12|MYC6_...     At5g41410|AT5G41410                 Abnormal integuments. Female sterile. Malformed ovules. Only a single integument-like structure form...     PO:0020021|PATO:0000051|PATO:0000460|GO:0005515|GO:0009506|G...     Plant PhenomeNET|TAIR
50  ath     HLL|Huellenlos|At1g17560|AT1G1...   HUELLENLOS|F1L3.27|F...     At1g17560|AT1G17560                 Abnormal integuments. Female sterile. Female-sterile, integuments mostly absent, ovules senesce prio...     PO:0020021|PATO:0000051|PATO:0000460|PO:0009062|PATO:0000277...     Plant PhenomeNET|TAIR
...
```


### Publication
A publication that uses this dataset to look at how computational methods can be used to work with phenotype descriptions is in progress.



### Have any feedback?
Send any feedback, comments, questions, or suggestions to irbraun at iastate dot edu.
