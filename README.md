# SplitVep
The INFO/CSQ field in the Ensembl VEP annotated VCF file are separated with pipe "|", which make the annotations hard to extract by bcftools. Also, for multi-allelic variants, VEP stores the annotations in the same INFO/CSQ line with comma "," as seperator， which exacerbate the situation. This software is to convert the VEP annotated VCF file into bgzip compressed a TSV formatted file.

# Dependency
SplitVep requires htslib installed in the system.

# Input
A compressed VCF file with "##INFO=<ID=CSQ" in the header.

# Run the software
The software supports **Linux** system only currently.
```console
./SplitVep <input.vcf.gz> <output_prefix> <tmp_folder>
```

# Output
CSQ information is outputed in duplication mode. The different CSQ annotations for one variant will be outputed into different rows.

| CHROM  | POS | REF | ALT | Allele | Consequence | IMPACT | SYMBOL | Gene | Feature_type | Feature | BIOTYPE | ... | 
| ------ | --- | --- | --- | ------ | ----------- | ------ | ------ | ---- | ------------ | ------- | ------- | --- |
| chr1 | 10146 | AC | A | - | upstream_gene_variant | MODIFIER | DDX11L1 | ENSG00000223972 | Transcript | ENST00000450305 | transcribed_unprocessed_pseudogene | ... |
| chr1 | 10146 | AC | A | - | upstream_gene_variant | MODIFIER | DDX11L16 | ENSG00000290825 | Transcript | ENST00000456328 | lncRNA | ... |
| chr1 | 10146 | AC | A | - | downstream_gene_variant | MODIFIER | WASH7P | ENSG00000223972 | Transcript | ENST00000488147 | transcribed_unprocessed_pseudogene | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
