## plant-data

### Description
This repository organizes information about plant genes, phenotypes, and annotations from a variety of data sources, and contains notebooks and scripts for preprocessing and combining this information. The objective was to obtain a dataset for looking into how well natural language processing (NLP) approaches work for recovering or predicting biological information from phenotype descriptions.



### What's here?
1. All the files that were used for this dataset are listed in `file_descriptions.csv`, which includes links to the original data source where applicable. Not all of the listed files are included in this repository because some are only available through database subscriptions or requests. 
2. The `pipeline.sh` scripts runs each script (see `preprocessing` and `scripts` directories) for preprocessing and mergining the information present in those files. This pipeline generates all the files in the `reshaped_data` directory, and the first few lines of each of those intermediate data files are available in the `reshaped_samples` directory.
3. The primary dataset of interest that is produced is available here as `genes_text_annots`, and is available as csv, tsv, and json.
4. The primary dataset including just a small subset of the genes in the full dataset, and truncated fields, is available as `genes_texts_annots_sample`, as csv, tsv, or json, and is viewable on GitHub.



### What does the data look like?

The primary dataset is a list of genes and accompanying fields that include information about the gene, descriptions of phenotypes associated with the gene, and ontology term anotations mapped to the gene.

A sample of the json dataset with truncated fields is shown here.
```
[
    {
        "id": 95,
        "species": "ath", "unique_gene_identifiers": ["CSN5B", "COP9-Signalosome", "At1g71230", "AT1G71230"],
        "other_gene_identifiers": ["CSN5", "AJH2", "COP9-signalosome 5B", "F3I17.12", "..."],
        "gene_models": ["At1g71230", "AT1G71230"],
        "descriptions": "Abnormal root hair development. Altered response to auxin. Altered response to light. Altered respon...",
        "annotations": ["GO:0080147", "PATO:0001236", "PATO:0000460", "PO:0000003", "..."],
        "sources": ["Plant PhenomeNET", "TAIR"]
    },
    {
        "id": 96,
        "species": "ath", "unique_gene_identifiers": ["OXI1", "Oxidative Signal-Inducible", "AGC2-1", "At3g25250", "..."],
        "other_gene_identifiers": ["AGC2", "AtOXI1", "AGC2 kinase 1", "oxidative signal-inducible1", "..."],
        "gene_models": ["At3g25250", "AT3G25250"],
        "descriptions": "Abnormal root hair development. Susceptible to avirulent bacteria. Susceptible to virulent bacteria....",
        "annotations": ["GO:0080147", "PATO:0001236", "PATO:0000460", "GO:0009617", "..."],
        "sources": ["Plant PhenomeNET", "TAIR"]
    },
    ...
]
```

A sample of the tab-separated dataset with truncated fields is shown here.
```
id	species	unique_gene_identifiers	other_gene_identifiers	gene_models	descriptions	annotations	sources
40	ath	ZIG|Zigzag Stem|SGR4|SHOOT GRA...	VTI11|ATVTI1A|ATVTI1...	At5g39510|AT5G39510	Abnormal hypocotyl gravitropism. Abnormal inflorescence gravitropism. Abnormal hypocotyl and inflore...	GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009959|P...	Plant PhenomeNET|TAIR
41	ath	SGR2|Shoot Gravitropism|SGR3|S...	T8E3.1|T8E3_1|ATVAM3...	At1g31480|At5g46860|At2g01940|...	Abnormal hypocotyl gravitropism. Abnormal inflorescence stem gravitropism. Late flowering. Male game...	GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009959|P...	Plant PhenomeNET|TAIR
42	ath	TIR1|Transport Inhibitor Respo...	DOC1|UMB1|ASA1|LPR1|...	At3g62980|At1g70560|At3g02260|...	Abnormal hypocotyl formation. Abnormal lateral root formation. Insensitive to auxin. Few lateral roo...	GO:0032502|PO:0020100|PATO:0001236|PATO:0000460|GO:0010311|P...	Plant PhenomeNET|TAIR
43	ath	NEK6|NIMA (Never In Mitosis, G...	ATNEK6|IBO1|NIMA-REL...	At3g44200|AT3G44200	Abnormal hypocotyl protrusions. Abnormal petiole protrusions. Abnormal trichome branching. Disordere...	PO:0020100|PATO:0000051|PATO:0000460|PO:0020038|PO:0000282|P...	Plant PhenomeNET|TAIR
44	ath	ARG1|Altered Response to Gravi...	ALTERED RESPONSE TO ...	At1g68370|AT1G68370	Abnormal hypocotyl gravitropism. Abnormal root gravitropism. Abnormal hypocotyl and root gravitropis...	GO:0009630|PO:0020100|PATO:0001236|PATO:0000460|GO:0009958|P...	Plant PhenomeNET|TAIR
45	ath	PIN5|Pin Formed|At5g16530|AT5G...	PIN-FORMED 5|MQK4.28...	At5g16530|AT5G16530	Abnormal hypocotyl growth. Abnormal root growth. Insensitive to IAA. Abnormal hypocotyl and root gro...	GO:0040007|PO:0020100|PATO:0001236|PATO:0000460|GO:0048364|P...	Plant PhenomeNET|TAIR
46	ath	PIL1|Phytochrome Interacting F...	PIF1|phytochrome int...	At2g46970|At3g59060|At2g20180|...	Abnormal cotyledon size under far-red light. Abnormal cotyledon size under red light. Abnormal hypoc...	PO:0020030|PATO:0000117|PATO:0000460|EO:0007203|EO:0007207|P...	Plant PhenomeNET|TAIR
47	ath	GPAT5|Glycerol-3-Phosphate Acy...	ATGPAT1|sn-2-GPAT1|G...	At3g11430|At2g38110|At1g06520|...	Abnormal insoluble lipid polyester biosynthesis. Complete loss of cuticle on petals. Complete loss o...	GO:0009058|CHEBI:60160|PATO:0001236|PATO:0000460|PO:0009032|...	Plant PhenomeNET|TAIR
48	ath	INO|Inner No Outer|At1g23420|A...	INNER NO OUTER|F26F2...	At1g23420|AT1G23420	Abnormal integuments. Strong allele. Ovule defective. Outer integument fails to initiate on the abax...	PO:0020021|PATO:0000051|PATO:0000460|GO:0048481|GO:0005515|G...	Plant PhenomeNET|TAIR
49	ath	BEL1|Bell-Shaped Ovules|At5g41...	BELL 1|MYC6.12|MYC6_...	At5g41410|AT5G41410	Abnormal integuments. Female sterile. Malformed ovules. Only a single integument-like structure form...	PO:0020021|PATO:0000051|PATO:0000460|GO:0005515|GO:0009506|G...	Plant PhenomeNET|TAIR
50	ath	HLL|Huellenlos|At1g17560|AT1G1...	HUELLENLOS|F1L3.27|F...	At1g17560|AT1G17560	Abnormal integuments. Female sterile. Female-sterile, integuments mostly absent, ovules senesce prio...	PO:0020021|PATO:0000051|PATO:0000460|PO:0009062|PATO:0000277...	Plant PhenomeNET|TAIR
...
```


### Publication
A publication that uses this dataset to look at how computational methods can handle phenotype descriptions is in progress.



### Feedback
Send any feedback, questions, or suggestions to irbraun at iastate dot edu.
