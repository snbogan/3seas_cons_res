# Title of Dataset: Data processing and visualization files and scripts for Zostera marina WGS.

=====================================

The text below summarizes the code and files used to process Zostera marina WGS data associated with the manuscript, "Genomic responses to parallel selection in the eelgrass Zostera marina in adjacent bays" by Lauren M. Schiebelhut, Rachael A. Bay, Richard K. Grosberg, & John J. Stachowicz.

## Sharing/access Information

Raw reads are deposited in the NCBI SRA under BioProject PRJNA887384.

The reads were mapped to the publically available Zostera marina 3.1 genome (PRJNA701932)

The following files are available on Dryad under https://doi.org/10.6071/M3DD4F

	The main filtered vcf file generated for further pruning and/or analyses is:
	Zm_TomBod_MAF01MM85INDM30_AllChr.recode.vcf.gz

	The accompanying metadata file with relevant site information and temperature data is:
	Zm_TomBod_MAF01MM85INDM30_AllChr.metadata.txt
	
	Temperature data from loggers:
	temperature_data_shared_dates.txt
	
	Code and associated files that are described in the remainder of this readme file:
	code_and_associated_files.tar.gz
	

=====================================

## Description of the Data and file structure. The headings match the directory names.

=====================================

### 1_alignment_and_snp_calling

	# 1
	bbduk from the BBTools suite v38.73 for adapter trimming and quality and length filtering

	# 2
	aligned whole genome sequences to the Zostera marina 3.1 genome (NCBI accession number PRJNA701932) with bwa-mem in bwa v0.7.13

	# 3
	Called SNPs using GATK v. 4.1.0.0, converted sam files to bam format and sorted using samtools v1.9. We marked duplicates, called haplotypes, combined g.vcf files and genotyped individuals across batches of 50 using GATK. 

	# 4
	After retaining only SNPs, we applied additional hard filters for mapping quality, strand bias, variant confidence, and variants with excessive depth following the GATK best practices guidelines. We further filtered sites to have a minimum depth of 10 and a minimum genotype quality score of 30. 

	# 5_identify_clones
	We also used the SNP data to identify and remove clones. To filter clones from the data set, we used Rclone v1.0.2 in R v4.0.3 with a reduced set of SNPs without any missing data (as required by the program) for each geographic location separately. 

	# 6_remove_clones_final_filter
	After removing clones we used vcftools v0.1.16 to filter the final data set to include only bi-allelic SNPs with a minor allele frequency of at least 0.01 and a genotype call rate of at least 85% of individuals.

=====================================

### 2_population_genetic_analyses

	# 1_prune_pca
	We thinned the SNPs based on linkage disequilibrium (LD) using SNPRelate with an LD threshold of 0.5 based on SNPRelate’s “composite” measure of LD. This thinned set was used for clustering analysis, PCA, and FST. We used SNPRelate to conduct principal components analysis.

	# 2_fst_clustering_analysis
	We used SNPRelate to calculate pairwise FST between all pairs of sites. We conducted clustering analysis and estimated ancestry proportions using the R package tess3r. For each K value (1-6) we ran five replicate runs with the lambda parameter set to 0 to ignore priors based on the spatial distribution of samples. 


=====================================

### 3_environmental_association_analyses

We used (1) latent factor mixed models (LFMM) and (2) FST outliers using OutFLANK v0.2 (grouping warmer vs. cooler sites) to search for SNPs associated with environmental gradients in Tomales Bay and Bodega Harbor separately correlating individual genotypes at each SNP with mean temperature at the sampling location. Each analysis was done after first filtering for SNPs with a minor allele frequency greater than 5% for samples within a given subset of sampling locations.

	# 1_imputed_missing_data_lfmm_outflank
	Using lfmm2 in LEA v4.0.0 genotypes were first imputed in LEA using ancestry coefficients estimated by LEA (K=3). Filtered for SNPs with a minor allele frequency greater than 5%. LFMM. The temperature gradient in Tomales Bay reaches much higher temperatures than that in Bodega Harbor, so to capture a similar range of temperatures for comparisons between bays, we also repeated the analyses with a subset of four sites in Tomales Bay (Lawson’s Landing, Pita Beach, Nick’s Cove, and Sacramento Landing; figure 1) that more closely matched the temperature range of Bodega Harbor. Using FST outliers using OutFLANK v0.2, grouping warmer vs. cooler sites. We also used OutFLANK to identify SNPs associated with higher versus lower intertidal habitat in Bodega. 

	# 2_identifying_parallel_associations
	To characterize parallel associations with temperature between the two bays:

		## 1_snp
		At the SNP level, we directly compared the genomic positions of outlier SNPs from LFMM and OutFLANK analyses. 

		## 2_gene_ldannot
		At the gene level, we used LD-Annot v0.4 with r2 = 0.9 to first identify the genes in linkage disequilibrium with candidate SNPs and then compare lists to detect overlapping genes. 

		## 3_function_topgo
		Finally, to determine whether there was overlap at the functional level, we tested whether the outlier-associated gene sets were enriched for particular gene ontology (GO) terms using TopGO v2.40.0 in R. We used a possible gene universe of all genes within linkage disequilibrium (r2 = 0.9) of the full set of SNPs, rather than the full set of annotated genes. To identify enriched GO terms we used a significance cutoff of p<0.05 and required that more than two genes be significant for a particular GO term.

	# 3_pgs
	To test predictability of genotype-environment associations within and across bays, we created polygenic scores using the R-package randomForest v4.6-14. We used the mean temperature for each site as the response variable and candidate SNPs derived from LFMM or OutFLANK analyses as predictors. We first thinned the candidate SNP sets based on linkage disequilibrium in SNPRelate (ld.threshold=0.5). Because randomForest cannot accept missing data, we used genotypes that were imputed with LEA. We ran separate random forest models with the sets of candidate SNPs from LFMM and OutFLANK for each bay. For each SNP set, we conducted runs using either all samples to train the model, or only samples from the bay from which the candidate set was derived (e.g. when candidate SNPs were identified using analysis of Bodega Harbor samples, we trained using either all samples or Bodega Harbor samples only). For validation, we used either all samples, Bodega Harbor only, or Tomales Bay only. To limit bias in the training set, we conducted cross-validation at the level of the sampling site, training the random forest on the dataset with a single sampling site removed, then predicting temperature for that site. For each random forest run, this procedure gave us a random-forest predicted temperature for each individual sample based on candidate SNPs. We then compared predicted and observed temperatures using Spearman’s correlation coefficient. For each run, we report the percent variance explained in the training model and the correlation coefficient when comparing observed and predicted temperatures. We included in our random forest predictors the sample loadings on the first two PC axes and we ran both training and prediction with 100 random sets of SNPs (in addition to the first two PC axes). We used these randomizations to estimate a 90% confidence interval for a null distribution, showing how well our predictions performed when based on neutral genetic variation alone.

=====================================


