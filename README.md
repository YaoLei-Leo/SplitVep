# SplitVep
The INFO/CSQ field annotated by VEP is separated with | and hard to read. This software is to convert the VEP annotated VCF file into CSV formatted file.

# Input
A compressed VCF file with "##INFO=<ID=CSQ" in the header.

# Run the software
The software currently support Linux system only.
```console
./SplitVep <input.vcf.gz> <output_prefix>
```

# Output
CSQ information is outputed in duplication mode. The different CSQ annotations for one variant will be outputed into different rows.

| CHROM  | POS | REF | ALT | Allele | Consequence | IMPACT | SYMBOL | Gene | Feature_type | Feature | BIOTYPE | ... | 
| ------ | --- | --- | --- | ------ | ----------- | ------ | ------ | ---- | ------------ | ------- | ------- | --- |
| chr1 | 10146 | AC | A | - | upstream_gene_variant | MODIFIER | DDX11L1 | ENSG00000223972 | Transcript | ENST00000450305 | transcribed_unprocessed_pseudogene | ... |
| chr1 | 10146 | AC | A | - | upstream_gene_variant | MODIFIER | DDX11L16 | ENSG00000290825 | Transcript | ENST00000456328 | lncRNA | ... |
| chr1 | 10146 | AC | A | - | downstream_gene_variant | MODIFIER | WASH7P | ENSG00000223972 | Transcript | ENST00000488147 | transcribed_unprocessed_pseudogene | ... |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
