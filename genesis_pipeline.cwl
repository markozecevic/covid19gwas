{
    "class": "Workflow",
    "cwlVersion": "v1.2",
    "label": "GENESIS pipeline",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "inputs": [
        {
            "id": "rename_chrs",
            "sbg:fileTypes": "TXT",
            "type": "File?",
            "label": "Rename chromosomes",
            "doc": "Rename chromosomes according to the map in file, with OLD_NAME NEW_NAME pairs separated by whitespaces, each on a separate line.",
            "sbg:x": -1215.6956787109375,
            "sbg:y": -349.7339172363281
        },
        {
            "id": "phenotype_file",
            "sbg:fileTypes": "RDATA",
            "type": "File",
            "label": "Phenotype file",
            "doc": "RData file with an AnnotatedDataFrame of phenotypes and covariates. Sample identifiers must be in column named “sample.id”.",
            "sbg:x": -634.0674438476562,
            "sbg:y": -397
        },
        {
            "id": "outcome",
            "type": "string",
            "label": "Outcome",
            "doc": "Name of column in Phenotype file containing outcome variable.",
            "sbg:exposed": true
        },
        {
            "id": "covars",
            "type": "string[]?",
            "label": "Covariates",
            "doc": "Names of columns in Phenotype file containing covariates.",
            "sbg:exposed": true
        },
        {
            "id": "out_prefix",
            "type": "string",
            "label": "Output prefix",
            "doc": "Prefix that will be included in all output files.",
            "sbg:x": -156.25,
            "sbg:y": -534.5
        },
        {
            "id": "n_pcs",
            "type": "int?",
            "label": "Number of PCs to include as covariates",
            "doc": "Number of PCs from PCA file to include as covariates.",
            "sbg:exposed": true
        },
        {
            "id": "in_database_compressed",
            "sbg:fileTypes": "TAR.GZ",
            "type": "File?",
            "label": "Database bundle",
            "doc": "Compressed database directory. Available as LZ_Database.tar.gz in the public file gallery. Required if Database Directory is not provided.",
            "sbg:x": 1863.8502197265625,
            "sbg:y": -572.1030883789062
        },
        {
            "id": "samples_file",
            "sbg:fileTypes": "TXT",
            "type": "File?",
            "label": "Samples file",
            "doc": "File of samples to include (or exclude with \"^\" prefix).",
            "sbg:x": -1388.005859375,
            "sbg:y": -218.0272674560547
        },
        {
            "id": "in_variants",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
            "type": "File",
            "label": "Input variants file",
            "doc": "Input variants file.",
            "sbg:x": -1383.5166015625,
            "sbg:y": 18.786296844482422
        },
        {
            "id": "field",
            "type": "string",
            "label": "Case/control indicator variable",
            "doc": "A variable in the phenotype which indicates samples belonging to either cases or control group.",
            "sbg:x": -613.9481811523438,
            "sbg:y": -9.818663597106934
        },
        {
            "id": "autosome_only",
            "type": "boolean?",
            "label": "Autosomes only",
            "doc": "Only include variants on the autosomes.",
            "sbg:exposed": true
        },
        {
            "id": "controls_threshold",
            "type": "float?",
            "label": "HWE test p.value threshold for controls",
            "doc": "HWE test p.value threshold for controls",
            "sbg:exposed": true
        },
        {
            "id": "cases_threshold",
            "type": "float?",
            "label": "HWE threshold for cases",
            "doc": "HWE test p.value threshold for cases",
            "sbg:exposed": true
        },
        {
            "id": "chain",
            "type": "File",
            "label": "Chain file",
            "doc": "Chain file for remapping genomic coordinates",
            "sbg:x": 1835.43408203125,
            "sbg:y": 277.6538391113281
        },
        {
            "id": "conditional_variant_file",
            "sbg:fileTypes": "RDATA",
            "type": "File?",
            "label": "Conditional variant file",
            "doc": "RData file with a data.frame of identifiers for variants to be included as covariates for conditional analysis. Columns should include “chromosome” and “variant.id” that match the variant.id in the GDS files. The alternate allele dosage of these variants will be included as covariates in the analysis.",
            "sbg:x": 473.1524658203125,
            "sbg:y": 92.6304931640625
        },
        {
            "id": "custom_input_1",
            "type": "boolean",
            "label": "Remap to hg38",
            "sbg:x": 1638.8797607421875,
            "sbg:y": 79.59199523925781
        },
        {
            "id": "outname",
            "type": "string",
            "label": "Output name",
            "doc": "Name to be given to the merged CSV file.",
            "sbg:exposed": true
        },
        {
            "id": "threshold",
            "type": "float?",
            "label": "Threshold for results filtering",
            "doc": "Only the results with p values below this threshold will be kept in the final output file.",
            "sbg:exposed": true
        },
        {
            "id": "lz_threshold",
            "type": "float?",
            "sbg:exposed": true
        },
        {
            "id": "out_prefix_1",
            "type": "string?",
            "label": "Output prefix",
            "doc": "Prefix for files created by this script.",
            "sbg:exposed": true
        },
        {
            "id": "track_threshold",
            "type": "float?",
            "label": "Track threshold",
            "doc": "P-value threshold for selecting regions to display.",
            "sbg:exposed": true
        },
        {
            "id": "maf_threshold",
            "type": "float?",
            "label": "MAF threshold",
            "doc": "Minimum minor allele frequency for variants to include in test. Only used if MAC threshold is NA.",
            "sbg:exposed": true
        },
        {
            "id": "mac_threshold",
            "type": "float?",
            "label": "MAC threshold",
            "doc": "Minimum minor allele count for variants to include in test. Recommend to use a higher threshold when outcome is binary or count data. To disable it set it to NA.",
            "sbg:exposed": true
        },
        {
            "id": "signif_line_fixed",
            "type": "float?",
            "label": "Significance line",
            "doc": "P-value for the significance line. Only used if `signif_type = fixed`.",
            "sbg:exposed": true
        }
    ],
    "outputs": [
        {
            "id": "html_reports",
            "outputSource": [
                "null_model/html_reports"
            ],
            "sbg:fileTypes": "html",
            "type": "File[]?",
            "label": "Null model reports",
            "doc": "HTML Reports generated by the tool.",
            "sbg:x": 1355.7774658203125,
            "sbg:y": 116.25656127929688
        },
        {
            "id": "assoc_plots",
            "outputSource": [
                "single_variant_association_testing/assoc_plots"
            ],
            "sbg:fileTypes": "PNG",
            "type": "File[]?",
            "label": "Single variant association test plots",
            "doc": "QQ and Manhattan Plots of p-values in association test results.",
            "sbg:x": 1740.6304931640625,
            "sbg:y": -735.7452392578125
        },
        {
            "id": "king_robust_plots",
            "outputSource": [
                "king_robust/king_robust_plots"
            ],
            "sbg:fileTypes": "PDF",
            "type": "File[]?",
            "label": "KING kinship plots",
            "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
            "sbg:x": 527.9237670898438,
            "sbg:y": -624.5982055664062
        },
        {
            "id": "pcair_plots",
            "outputSource": [
                "pc_air/pcair_plots"
            ],
            "type": "File[]?",
            "label": "PC-Air plots",
            "doc": "PC plots",
            "sbg:x": 849.4087524414062,
            "sbg:y": 139.45037841796875
        },
        {
            "id": "pcrelate_plots",
            "outputSource": [
                "pc_relate/pcrelate_plots"
            ],
            "sbg:fileTypes": "PDF",
            "type": "File[]?",
            "label": "PC-Relate kinship plots",
            "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
            "sbg:x": 1119.14697265625,
            "sbg:y": -613.5075073242188
        },
        {
            "id": "out_pdf_reports",
            "outputSource": [
                "genesis_locuszoom/out_pdf_reports"
            ],
            "sbg:fileTypes": "PDF",
            "type": "File[]?",
            "label": "Locuszoom plots",
            "doc": "One LZ plot per locus.",
            "sbg:x": 2583.903076171875,
            "sbg:y": -682.5721435546875
        },
        {
            "id": "pcair_output",
            "outputSource": [
                "pc_air/pcair_output"
            ],
            "sbg:fileTypes": "RDATA",
            "type": "File?",
            "label": "RData file with PC-AiR PCs for all samples",
            "sbg:x": 893.217041015625,
            "sbg:y": -164.47802734375
        },
        {
            "id": "pc_correlation_plots",
            "outputSource": [
                "pc_air/pc_correlation_plots"
            ],
            "sbg:fileTypes": "PNG",
            "type": "File[]?",
            "label": "PC-variant correlation plots",
            "doc": "PC-variant correlation plots",
            "sbg:x": 1019.3255615234375,
            "sbg:y": 67.15247344970703
        },
        {
            "id": "hg38_out",
            "outputSource": [
                "remap/hg38_out"
            ],
            "sbg:fileTypes": "CSV",
            "type": "File?",
            "label": "Single variant association test results - hg38",
            "sbg:x": 2647.8916015625,
            "sbg:y": -154.48696899414062
        },
        {
            "id": "hg19_out",
            "outputSource": [
                "remap/hg19_out"
            ],
            "sbg:fileTypes": "CSV",
            "type": "File?",
            "label": "Single variant association test results - hg19",
            "sbg:x": 2651.98828125,
            "sbg:y": 194.73007202148438
        },
        {
            "id": "merged",
            "outputSource": [
                "sbg_merge_and_filter_genesis_results_cwl2/merged"
            ],
            "sbg:fileTypes": "CSV",
            "type": "File?",
            "label": "Merged CSV",
            "sbg:x": 2243.531005859375,
            "sbg:y": -235.0034942626953
        }
    ],
    "steps": [
        {
            "id": "null_model",
            "in": [
                {
                    "id": "relatedness_matrix_file",
                    "source": "pc_relate/pcrelate_matrix"
                },
                {
                    "id": "phenotype_file",
                    "source": "phenotype_file"
                },
                {
                    "id": "pca_file",
                    "source": "pc_air/pcair_output"
                },
                {
                    "id": "gds_files",
                    "source": [
                        "vcf_to_gds_1/unique_variant_id_gds_per_chr"
                    ]
                },
                {
                    "id": "conditional_variant_file",
                    "source": "conditional_variant_file"
                },
                {
                    "id": "outcome",
                    "source": "outcome"
                },
                {
                    "id": "n_pcs",
                    "source": "n_pcs"
                },
                {
                    "id": "cpu",
                    "default": 8
                },
                {
                    "id": "covars",
                    "source": [
                        "covars"
                    ]
                },
                {
                    "id": "family",
                    "default": "binomial"
                }
            ],
            "out": [
                {
                    "id": "null_model_phenotypes"
                },
                {
                    "id": "rmd_files"
                },
                {
                    "id": "html_reports"
                },
                {
                    "id": "null_model_file"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/null-model/0",
                "doc": "**Null Model** workflow fits the regression or mixed effects model under the null hypothesis of no genotype effects. i.e., The outcome variable is regressed on the specified fixed effect covariates and random effects. The output of this null model is then used in the association tests.\n\nQuantitative and binary outcomes are both supported. Set parameter **family** to gaussian, binomial or poisson depending on the outcome type. Fixed effect covariates from the **Phenotype file** are specified using the **Covariates** parameter, and ancestry principal components can be included as fixed effects using the **PCA Files** and **Number of PCs to include as covariates** parameters. A kinship matrix (KM) or genetic relationship matrix (GRM) can be provided using the **Relatedness matrix file** parameter to account for genetic similarity among samples as a random effect. \n\nWhen no **Relatedness matrix file** is provided, standard linear regression is used if the parameter **family** is set to gaussian, logistic regression is used if the parameter **family** is set to binomial and poisson regression is used in case when **family** is set to poisson. When **Relatedness matrix file** is provided, a linear mixed model (LMM) is fit if **family** is set to gaussian, or a generalized linear mixed model (GLMM) is fit using the [GMMAT method](https://doi.org/10.1016/j.ajhg.2016.02.012) if **family** is set to binomial or poisson. For either the LMM or GLMM, the [AI-REML algorithm](https://doi.org/10.1111/jbg.12398) is used to estimate the variance components and fixed effects parameters.\n \nWhen samples come from multiple groups (e.g., study or ancestry group), it is common to observe different variances by group for quantitative traits. It is recommended to allow the null model to fit heterogeneous residual variances by group using the parameter group_var. This often provides better control of false positives in subsequent association testing. Note that this only applies when **family** is set to gaussian.\n\nRank-based inverse Normal transformation is supported for quantitative outcomes via the inverse_normal parameter. This parameter is TRUE by default. When **inverse normal** parameter is set to TRUE, (1) the null model is fit using the original outcome values, (2) the marginal residuals are rank-based inverse Normal transformed, and (3) the null model is fit again using the transformed residuals as the outcome; fixed effects and random effects are included both times the null model is fit. It has been shown that this fully adjusted two-stage procedure provides better false positive control and power in association tests than simply inverse Normalizing the outcome variable prior to analysis [(**reference**)](https://doi.org/10.1002/gepi.22188).\n\nThis workflow utilizes the *fitNullModel* function from the [GENESIS](doi.org/10.1093/bioinformatics/btz567) software.\n\nWorkflow consists of two steps. First step fits the null model, and the second one generates reports based on data. Reports are available both in RMD and HTML format. If **inverse normal** is TRUE, reports are generated for the model both before and after the transformation.\nReports contain the following information: Config info, phenotype distributions, covariate effect size estimates, marginal residuals, adjusted phenotype values and session information.\n\n### Common use cases:\nThis workflow is the first step in association testing. This workflow fits the null model and produces several files which are used in the association testing workflows:\n* Null model file which contains adjusted outcome values, the model matrix, the estimated covariance structure, and other parameters required for downstream association testing. This file will be input in association testing workflows.\n* Phenotype file which is a subset of the provided phenotype file, containing only the outcome and covariates used in fitting null model.\n* *Reportonly* null model file which is used to generate the report for the association test\n\n\nThis workflow can be used for trait heritability estimation.\n\nIndividual genetic variants or groups of genetic variants can be directly tested for association with this workflow by including them as fixed effect covariates in the model (via the **Conditional Variant File** parameter). This would be extremely inefficient genome-wide, but is useful for follow-up analyses testing variants of interest.\n\n\n### Common issues and important notes:\n* If **PCA File** is not provided, the **Number of PCs to include as covariates** parameter **must** be set to 0.\n\n* **PCA File** must be an RData object output from the *pcair* function in the GENESIS package.\n\n* The null model job can be very computationally demanding in large samples (e.g. > 20K). GENESIS supports using sparse representations of matrices in the **Relatedness matrix file** via the R Matrix package, and this can substantially reduce memory usage and CPU time.\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running time and cost on spot instances. \n      \n\n| Sample Count | Relatedness matrix   | Duration   | Cost - spot ($)  |  Instance (AWS)  |\n|--------------------|----------------------------|-------------|---------------------|------------------------|\n| 2.5k samples  |                 | 4 min                 | $0.01   |   1x c4.xlarge |\n| 2.5k samples  | sparse     | 5 min                 | $0.01   |   1x c4.xlarge |\n| 2.5k samples  | dense      | 5 min                 | $0.01   |   1x c4.xlarge |\n| 10k samples   |                 | 6 min                 | $0.06   |   1x r4.8xlarge |\n| 10k samples   | sparse     | 7 min                 | $0.07   |   1x r4.8xlarge |\n| 10k samples   | dense      | 16 min               | $0.13    |   1x r4.8xlarge |\n| 36k samples  |                 | 7 min                 | $0.06    |   1x r4.8xlarge |\n| 36k samples  | sparse     | 24 min               | $0.27    |   1x r4.8xlarge |\n| 36k samples  | dense      | 52 min               | $0.56    |   1x r4.8xlarge |\n| 54k samples  |                 | 7 min                 | $0.07     |   1x r4.8xlarge |\n| 54k samples  | sparse     | 32 min               | $0.36    |   1x r4.8xlarge |\n| 54k samples  | dense      | 2 h                     | $1.5       | 1x r4.8xlarge |\n\n\n*Cost shown here were obtained with **spot instances** enabled. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*",
                "label": "GENESIS Null Model",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "sample_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample include file",
                        "doc": "RData file with a vector of sample.id to include. If not provided, all samples in the Phenotype file will be included in the analysis.",
                        "sbg:x": -608.4528198242188,
                        "sbg:y": -792.2285766601562
                    },
                    {
                        "id": "relatedness_matrix_file",
                        "sbg:fileTypes": "GDS, RDATA",
                        "type": "File?",
                        "label": "Relatedness matrix file",
                        "doc": "RData or GDS file with a kinship matrix or genetic relatedness matrix (GRM). For RData files, R object type may be “matrix” or “Matrix”. For very large sample sizes, a block diagonal sparse Matrix object from the “Matrix” package is recommended.",
                        "sbg:x": -611,
                        "sbg:y": -661
                    },
                    {
                        "id": "phenotype_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File",
                        "label": "Phenotype file",
                        "doc": "RData file with an AnnotatedDataFrame of phenotypes and covariates. Sample identifiers must be in column named “sample.id”.",
                        "sbg:x": -779.132080078125,
                        "sbg:y": -731.9622802734375
                    },
                    {
                        "id": "pca_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "PCA File",
                        "doc": "RData file containing principal components for ancestry adjustment. R object type may be “pcair”, data.frame, or matrix. Row names must contain sample identifiers.",
                        "sbg:x": -877.1132202148438,
                        "sbg:y": -620.0188598632812
                    },
                    {
                        "id": "output_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Base for all output file names. By default it is null_model.",
                        "sbg:x": -693.75,
                        "sbg:y": -475.75
                    },
                    {
                        "id": "gds_files",
                        "sbg:fileTypes": "GDS",
                        "type": "File[]?",
                        "label": "GDS Files",
                        "doc": "List of gds files with genotype data for variants to be included as covariates for conditional analysis. Only required if Conditional Variant file is specified.",
                        "sbg:x": -726.8867797851562,
                        "sbg:y": -351.88677978515625
                    },
                    {
                        "id": "inverse_normal",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "TRUE",
                                    "FALSE"
                                ],
                                "name": "inverse_normal"
                            }
                        ],
                        "label": "Two stage model",
                        "doc": "TRUE if a two-stage model should be implemented. Stage 1: a null model is fit using the original outcome variable. Stage 2: a second null model is fit using the inverse-normal transformed residuals from Stage 1 as the outcome variable. When FALSE, only the Stage 1 model is fit.  Only applies when Family is “gaussian”.",
                        "sbg:toolDefaultValue": "TRUE",
                        "sbg:x": -565.0943603515625,
                        "sbg:y": -251.3773651123047
                    },
                    {
                        "id": "conditional_variant_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Conditional variant file",
                        "doc": "RData file with a data.frame of identifiers for variants to be included as covariates for conditional analysis. Columns should include “chromosome” and “variant.id” that match the variant.id in the GDS files. The alternate allele dosage of these variants will be included as covariates in the analysis.",
                        "sbg:x": -365.25,
                        "sbg:y": -231.75
                    },
                    {
                        "id": "rescale_variance",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "marginal",
                                    "varcomp",
                                    "none"
                                ],
                                "name": "rescale_variance"
                            }
                        ],
                        "label": "Rescale residuals",
                        "doc": "Applies only if Two stage model is TRUE. Controls whether to rescale the inverse-normal transformed residuals before fitting the Stage 2 null model, restoring the values to their original scale before the transform. “Marginal” rescales by the standard deviation of the marginal residuals from the Stage 1 model. “Varcomp” rescales by an estimate of the standard deviation based on the Stage 1 model variance component estimates; this can only be used if Norm by group is TRUE. “None” does not rescale.",
                        "sbg:toolDefaultValue": "Marginal",
                        "sbg:x": -576,
                        "sbg:y": -43.5
                    },
                    {
                        "id": "outcome",
                        "type": "string",
                        "label": "Outcome",
                        "doc": "Name of column in Phenotype file containing outcome variable.",
                        "sbg:x": -671.5,
                        "sbg:y": 1
                    },
                    {
                        "id": "norm_bygroup",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "TRUE",
                                    "FALSE"
                                ],
                                "name": "norm_bygroup"
                            }
                        ],
                        "label": "Norm by group",
                        "doc": "Applies only if Two stage model is TRUE and Group variate is provided. If TRUE,the inverse-normal transformation (and rescaling) is done on each group separately. If FALSE, this is done on all samples jointly.",
                        "sbg:toolDefaultValue": "FALSE",
                        "sbg:x": -562.75,
                        "sbg:y": 133
                    },
                    {
                        "id": "n_pcs",
                        "type": "int?",
                        "label": "Number of PCs to include as covariates",
                        "doc": "Number of PCs from PCA file to include as covariates.",
                        "sbg:toolDefaultValue": "0",
                        "sbg:x": -651.81103515625,
                        "sbg:y": 180.75
                    },
                    {
                        "id": "group_var",
                        "type": "string?",
                        "label": "Group variate",
                        "doc": "Name of column in Phenotype file providing groupings for heterogeneous residual error variances in the model. Only applies when Family is “gaussian”.",
                        "sbg:x": -449.409423828125,
                        "sbg:y": 285.0613098144531
                    },
                    {
                        "id": "cpu",
                        "type": "int?",
                        "label": "CPU",
                        "doc": "Number of CPUs to use per job.",
                        "sbg:toolDefaultValue": "1",
                        "sbg:x": -538,
                        "sbg:y": 329.25
                    },
                    {
                        "id": "covars",
                        "type": "string[]?",
                        "label": "Covariates",
                        "doc": "Names of columns in Phenotype file containing covariates.",
                        "sbg:x": -619.25,
                        "sbg:y": 377.5
                    },
                    {
                        "id": "family",
                        "type": {
                            "type": "enum",
                            "symbols": [
                                "gaussian",
                                "poisson",
                                "binomial"
                            ],
                            "name": "family"
                        },
                        "label": "Family",
                        "doc": "The distribution used to fit the model. Select “gaussian” for continuous outcomes, “binomial” for binary or case/control outcomes, or “poisson” for count outcomes.",
                        "sbg:x": -463.0585021972656,
                        "sbg:y": 69.44185638427734
                    },
                    {
                        "id": "n_categories_boxplot",
                        "type": "int?",
                        "label": "Number of categories in boxplot",
                        "doc": "If a covariate has fewer than the specified value, boxplots will be used instead of scatter plots for that covariate in the null model report.",
                        "sbg:x": 293.9433898925781,
                        "sbg:y": 314.99371337890625
                    }
                ],
                "outputs": [
                    {
                        "id": "null_model_phenotypes",
                        "outputSource": [
                            "null_model_r/null_model_phenotypes"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Null model Phenotypes file",
                        "doc": "Phenotype file containing all covariates used in the model. This file should be used as the “Phenotype file” input for the GENESIS association testing workflows.",
                        "sbg:x": 529.75,
                        "sbg:y": 81
                    },
                    {
                        "id": "rmd_files",
                        "outputSource": [
                            "null_model_report/rmd_files"
                        ],
                        "sbg:fileTypes": "Rmd",
                        "type": "File[]?",
                        "label": "Rmd files",
                        "doc": "R markdown files used to generate the HTML reports.",
                        "sbg:x": 834.340576171875,
                        "sbg:y": -492.55523681640625
                    },
                    {
                        "id": "html_reports",
                        "outputSource": [
                            "null_model_report/html_reports"
                        ],
                        "sbg:fileTypes": "html",
                        "type": "File[]?",
                        "label": "HTML Reports",
                        "doc": "HTML Reports generated by the tool.",
                        "sbg:x": 822.75,
                        "sbg:y": 24.75
                    },
                    {
                        "id": "null_model_file",
                        "outputSource": [
                            "null_model_r/null_model_output"
                        ],
                        "type": "File[]?",
                        "label": "Null model file",
                        "sbg:x": 501.5128173828125,
                        "sbg:y": 291.2054138183594
                    }
                ],
                "steps": [
                    {
                        "id": "null_model_r",
                        "in": [
                            {
                                "id": "outcome",
                                "source": "outcome"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "gds_files",
                                "source": [
                                    "gds_files"
                                ]
                            },
                            {
                                "id": "pca_file",
                                "source": "pca_file"
                            },
                            {
                                "id": "relatedness_matrix_file",
                                "source": "relatedness_matrix_file"
                            },
                            {
                                "id": "family",
                                "source": "family"
                            },
                            {
                                "id": "conditional_variant_file",
                                "source": "conditional_variant_file"
                            },
                            {
                                "id": "covars",
                                "source": [
                                    "covars"
                                ]
                            },
                            {
                                "id": "group_var",
                                "source": "group_var"
                            },
                            {
                                "id": "inverse_normal",
                                "source": "inverse_normal"
                            },
                            {
                                "id": "n_pcs",
                                "source": "n_pcs"
                            },
                            {
                                "id": "rescale_variance",
                                "source": "rescale_variance"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "cpu",
                                "source": "cpu"
                            },
                            {
                                "id": "output_prefix",
                                "source": "output_prefix"
                            },
                            {
                                "id": "norm_bygroup",
                                "source": "norm_bygroup"
                            }
                        ],
                        "out": [
                            {
                                "id": "configs"
                            },
                            {
                                "id": "null_model_phenotypes"
                            },
                            {
                                "id": "null_model_files"
                            },
                            {
                                "id": "null_model_params"
                            },
                            {
                                "id": "null_model_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/null-model-r/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Configs",
                                    "id": "outcome",
                                    "type": "string",
                                    "label": "Outcome",
                                    "doc": "Name of column in Phenotype File containing outcome variable."
                                },
                                {
                                    "sbg:category": "Categories",
                                    "id": "phenotype_file",
                                    "type": "File",
                                    "label": "Phenotype file",
                                    "doc": "RData file with AnnotatedDataFrame of phenotypes.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:altPrefix": "GDS file.",
                                    "sbg:category": "Configs",
                                    "id": "gds_files",
                                    "type": "File[]?",
                                    "label": "GDS Files",
                                    "doc": "List of gds files. Required if conditional_variant_file is specified.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "pca_file",
                                    "type": "File?",
                                    "label": "PCA File",
                                    "doc": "RData file with PCA results created by PC-AiR.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Categories",
                                    "id": "relatedness_matrix_file",
                                    "type": "File?",
                                    "label": "Relatedness matrix file",
                                    "doc": "RData or GDS file with a kinship matrix or GRM.",
                                    "sbg:fileTypes": "GDS, RDATA, RData"
                                },
                                {
                                    "sbg:category": "Configs",
                                    "sbg:toolDefaultValue": "gaussian",
                                    "id": "family",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "gaussian",
                                            "poisson",
                                            "binomial"
                                        ],
                                        "name": "family"
                                    },
                                    "label": "Family",
                                    "doc": "Depending on the output type (quantitative or qualitative) one of possible values should be chosen: Gaussian, Binomial, Poisson."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "conditional_variant_file",
                                    "type": "File?",
                                    "label": "Conditional Variant File",
                                    "doc": "RData file with data frame of of conditional variants. Columns should include chromosome and variant.id. The alternate allele dosage of these variants will be included as covariates in the analysis.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "covars",
                                    "type": "string[]?",
                                    "label": "Covariates",
                                    "doc": "Names of columns phenotype_file containing covariates."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "group_var",
                                    "type": "string?",
                                    "label": "Group variate",
                                    "doc": "Name of covariate to provide groupings for heterogeneous residual error variances in the mixed model."
                                },
                                {
                                    "sbg:toolDefaultValue": "TRUE",
                                    "sbg:category": "Configs",
                                    "id": "inverse_normal",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "TRUE",
                                                "FALSE"
                                            ],
                                            "name": "inverse_normal"
                                        }
                                    ],
                                    "label": "Inverse normal",
                                    "doc": "TRUE if an inverse-normal transform should be applied to the outcome variable. If Group variate is provided, the transform is done on each group separately."
                                },
                                {
                                    "sbg:toolDefaultValue": "3",
                                    "sbg:category": "Configs",
                                    "id": "n_pcs",
                                    "type": "int?",
                                    "label": "Number of PCs to include as covariates",
                                    "doc": "Number of PCs to include as covariates."
                                },
                                {
                                    "sbg:toolDefaultValue": "marginal",
                                    "sbg:category": "Configs",
                                    "id": "rescale_variance",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "marginal",
                                                "varcomp",
                                                "none"
                                            ],
                                            "name": "rescale_variance"
                                        }
                                    ],
                                    "label": "Rescale variance",
                                    "doc": "Applies only if Inverse normal is TRUE and Group variate is provided. Controls whether to rescale the variance for each group after inverse-normal transform, restoring it to the original variance before the transform. Options are marginal, varcomp, or none."
                                },
                                {
                                    "sbg:toolDefaultValue": "TRUE",
                                    "sbg:category": "Configs",
                                    "id": "resid_covars",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "TRUE",
                                                "FALSE"
                                            ],
                                            "name": "resid_covars"
                                        }
                                    ],
                                    "label": "Residual covariates",
                                    "doc": "Applies only if Inverse normal is TRUE. Logical for whether covariates should be included in the second null model using the residuals as the outcome variable."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample include file",
                                    "doc": "RData file with vector of sample.id to include.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:toolDefaultValue": "1",
                                    "sbg:category": "Input Options",
                                    "id": "cpu",
                                    "type": "int?",
                                    "label": "CPU",
                                    "doc": "Number of CPUs for each tool job. Default value: 1."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "sbg:toolDefaultValue": "null_model",
                                    "id": "output_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Base for all output file names. By default it is null_model."
                                },
                                {
                                    "sbg:category": "General",
                                    "sbg:toolDefaultValue": "FALSE",
                                    "id": "norm_bygroup",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "TRUE",
                                                "FALSE"
                                            ],
                                            "name": "norm_bygroup"
                                        }
                                    ],
                                    "label": "Norm by group",
                                    "doc": "If TRUE and group_var is provided, the inverse normal transform is done on each group separately."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "configs",
                                    "doc": "Config files.",
                                    "label": "Config files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.config*"
                                    }
                                },
                                {
                                    "id": "null_model_phenotypes",
                                    "doc": "Phenotypes file",
                                    "label": "Null model Phenotypes file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${\n    if(inputs.null_model_file)\n    {\n        return inputs.phenotype_file.basename\n    }\n    else\n    {\n        return \"*phenotypes.RData\"\n    }\n}",
                                        "outputEval": "$(inheritMetadata(self, inputs.phenotype_file))"
                                    },
                                    "sbg:fileTypes": "RData"
                                },
                                {
                                    "id": "null_model_files",
                                    "doc": "Null model file.",
                                    "label": "Null model file",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "${\n    if(inputs.null_model_file)\n    {\n        return inputs.null_model_file.basename\n    }\n    else\n    {\n        if(inputs.output_prefix)\n        {\n            return inputs.output_prefix + '_null_model*RData'\n        }\n        return \"*null_model*RData\"\n    }\n}"
                                    },
                                    "sbg:fileTypes": "RData"
                                },
                                {
                                    "id": "null_model_params",
                                    "doc": "Parameter file",
                                    "label": "Parameter file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.params"
                                    },
                                    "sbg:fileTypes": "params"
                                },
                                {
                                    "id": "null_model_output",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "${\n    if(inputs.null_model_file)\n    {\n        return inputs.null_model_file.basename\n    }\n    else\n    {\n        if(inputs.output_prefix)\n        {\n            return inputs.output_prefix + '_null_model*RData'\n        }\n        return \"*null_model*RData\"\n    }\n}",
                                        "outputEval": "${\n    var result = []\n    var len = self.length\n    var i;\n\n    for(i=0; i<len; i++){\n        if(!self[i].path.split('/')[self[0].path.split('/').length-1].includes('reportonly')){\n            result.push(self[i])\n        }\n    }\n    return result\n}"
                                    }
                                }
                            ],
                            "label": "null_model.R",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "${\n        return \"Rscript /usr/local/analysis_pipeline/R/null_model.R null_model.config\"\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    if (inputs.cpu)\n        return 'export NSLOTS=' + inputs.cpu + ' &&'\n    else\n        return ''\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "coresMin": "${\n    if(inputs.cpu){\n        return inputs.cpu\n    }\n    else{\n        return 1\n    }\n}"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "null_model.config",
                                            "entry": "${  \n\n    var arg = [];\n    if(inputs.output_prefix){\n        var filename = inputs.output_prefix + \"_null_model\";\n        arg.push('out_prefix \\\"' + filename + '\\\"');\n        var phenotype_filename = inputs.output_prefix + \"_phenotypes.RData\";\n        arg.push('out_phenotype_file \\\"' + phenotype_filename + '\\\"');\n    }\n    else{\n        arg.push('out_prefix \"null_model\"');\n        arg.push('out_phenotype_file \"phenotypes.RData\"');\n    }\n    arg.push('outcome ' + inputs.outcome);\n    arg.push('phenotype_file \"' + inputs.phenotype_file.basename + '\"');\n    if(inputs.gds_files){\n        \n       function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n       }    \n        \n       var gds = inputs.gds_files[0].path.split('/').pop();    \n       var right = gds.split('chr')[1];\n       var chr = [];\n       \n       if(isNumeric(parseInt(right.charAt(1)))) chr.push(right.substr(0,2))\n       else chr.push(right.substr(0,1))\n       \n        arg.push('gds_file \"' + inputs.gds_files[0].basename.split(\"chr\")[0] + \"chr \"+gds.split(\"chr\"+chr)[1] +'\"')\n        \n        \n    }\n    if(inputs.pca_file){\n        arg.push('pca_file \"' + inputs.pca_file.basename + '\"')\n    }\n    if(inputs.relatedness_matrix_file){\n        arg.push('relatedness_matrix_file \"' + inputs.relatedness_matrix_file.basename + '\"')\n    }\n    if(inputs.family){\n        arg.push('family ' + inputs.family)\n    }\n    if(inputs.conditional_variant_file){\n        arg.push('conditional_variant_file \"' + inputs.conditional_variant_file.basename + '\"')\n    }\n    if(inputs.covars){\n        var temp = [];\n        for(var i=0; i<inputs.covars.length; i++){\n            temp.push(inputs.covars[i])\n        }\n        arg.push('covars \"' + temp.join(' ') + '\"')\n    }\n    if(inputs.group_var){\n        arg.push('group_var \"' + inputs.group_var + '\"')\n    }\n    if(inputs.inverse_normal){\n        arg.push('inverse_normal ' + inputs.inverse_normal)\n    }\n    if(inputs.n_pcs){\n        if(inputs.n_pcs > 0)\n            arg.push('n_pcs ' + inputs.n_pcs)\n    }\n    if(inputs.rescale_variance){\n        arg.push('rescale_variance \"' + inputs.rescale_variance + '\"')\n    }\n    if(inputs.resid_covars){\n        arg.push('resid_covars ' + inputs.resid_covars)\n    }\n    if(inputs.sample_include_file){\n        arg.push('sample_include_file \"' + inputs.sample_include_file.basename + '\"')\n    }\n    if(inputs.norm_bygroup){\n        arg.push('norm_bygroup ' + inputs.norm_bygroup)\n    }\n    return arg.join('\\n')\n}",
                                            "writable": false
                                        },
                                        {
                                            "entry": "${\n    return inputs.phenotype_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.gds_files\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.pca_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.relatedness_matrix_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.conditional_variant_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.sample_include_file\n}",
                                            "writable": true
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement",
                                    "expressionLib": [
                                        "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                                    ]
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:AWSInstanceType",
                                    "value": "r4.8xlarge;ebs-gp2;500"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105816,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/null-model-r/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105816,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105816,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a29ebdb2dd77dd3a1b21a60950dc736723d2cd50f33c9702ef7ea9b21de79e471"
                        },
                        "label": "Fit Null Model",
                        "sbg:x": 112.75,
                        "sbg:y": 112.25
                    },
                    {
                        "id": "null_model_report",
                        "in": [
                            {
                                "id": "family",
                                "source": "family"
                            },
                            {
                                "id": "inverse_normal",
                                "source": "inverse_normal"
                            },
                            {
                                "id": "null_model_params",
                                "source": "null_model_r/null_model_params"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "pca_file",
                                "source": "pca_file"
                            },
                            {
                                "id": "relatedness_matrix_file",
                                "source": "relatedness_matrix_file"
                            },
                            {
                                "id": "null_model_files",
                                "source": [
                                    "null_model_r/null_model_files"
                                ]
                            },
                            {
                                "id": "output_prefix",
                                "source": "output_prefix"
                            },
                            {
                                "id": "conditional_variant_file",
                                "source": "conditional_variant_file"
                            },
                            {
                                "id": "gds_files",
                                "source": [
                                    "gds_files"
                                ]
                            },
                            {
                                "id": "n_categories_boxplot",
                                "source": "n_categories_boxplot"
                            }
                        ],
                        "out": [
                            {
                                "id": "html_reports"
                            },
                            {
                                "id": "rmd_files"
                            },
                            {
                                "id": "null_model_report_config"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/null-model-report/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "General",
                                    "id": "family",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "gaussian",
                                            "binomial",
                                            "poisson"
                                        ],
                                        "name": "family"
                                    },
                                    "label": "Family",
                                    "doc": "Possible values: Gaussian, Binomial, Poisson"
                                },
                                {
                                    "sbg:category": "General",
                                    "sbg:toolDefaultValue": "TRUE",
                                    "id": "inverse_normal",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "TRUE",
                                                "FALSE"
                                            ],
                                            "name": "inverse_normal"
                                        }
                                    ],
                                    "label": "Inverse normal",
                                    "doc": "TRUE if an inverse-normal transform should be applied to the outcome variable. If Group variate is provided, the transform is done on each group separately. Default value is TRUE."
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "null_model_params",
                                    "type": "File",
                                    "label": "Null Model Params",
                                    "doc": ".params file generated by null model.",
                                    "sbg:fileTypes": "params"
                                },
                                {
                                    "sbg:category": "Null model",
                                    "id": "phenotype_file",
                                    "type": "File?",
                                    "label": "Phenotype file",
                                    "doc": "RData file with AnnotatedDataFrame of phenotypes.",
                                    "sbg:fileTypes": "RData"
                                },
                                {
                                    "sbg:category": "Null model",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample include file",
                                    "doc": "RData file with vector of sample.id to include.",
                                    "sbg:fileTypes": "RData"
                                },
                                {
                                    "sbg:category": "Null model",
                                    "id": "pca_file",
                                    "type": "File?",
                                    "label": "PCA File",
                                    "doc": "RData file with PCA results created by PC-AiR.",
                                    "sbg:fileTypes": "RData"
                                },
                                {
                                    "sbg:category": "Null model",
                                    "id": "relatedness_matrix_file",
                                    "type": "File?",
                                    "label": "Relatedness Matrix File",
                                    "doc": "RData or GDS file with a kinship matrix or GRM.",
                                    "sbg:fileTypes": "GDS, RData, gds, RDATA"
                                },
                                {
                                    "id": "null_model_files",
                                    "type": "File[]?",
                                    "label": "Null model files",
                                    "doc": "Null model files."
                                },
                                {
                                    "id": "output_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Output prefix."
                                },
                                {
                                    "sbg:category": "General",
                                    "id": "conditional_variant_file",
                                    "type": "File?",
                                    "label": "Conditional variant file",
                                    "doc": "Conditional variant file",
                                    "sbg:fileTypes": "RData, RDATA"
                                },
                                {
                                    "sbg:category": "General",
                                    "id": "gds_files",
                                    "type": "File[]?",
                                    "label": "GDS Files",
                                    "doc": "GDS files",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:toolDefaultValue": "10",
                                    "id": "n_categories_boxplot",
                                    "type": "int?",
                                    "label": "Number of categories in boxplot",
                                    "doc": "Number of categories in boxplot",
                                    "default": 10
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "html_reports",
                                    "doc": "HTML Reports generated by the tool.",
                                    "label": "HTML Reports",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.html"
                                    },
                                    "sbg:fileTypes": "html"
                                },
                                {
                                    "id": "rmd_files",
                                    "doc": "R markdown files used to generate the HTML reports.",
                                    "label": "Rmd files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.Rmd"
                                    },
                                    "sbg:fileTypes": "Rmd"
                                },
                                {
                                    "id": "null_model_report_config",
                                    "doc": "Null model report config",
                                    "label": "Null model report config",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*config"
                                    }
                                }
                            ],
                            "label": "null_model_report",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "${\n    return \" Rscript /usr/local/analysis_pipeline/R/null_model_report.R null_model_report.config\"\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "null_model_report.config",
                                            "entry": "${\n    var config = \"\";\n    if(inputs.family)\n    {\n        config += \"family \" + inputs.family + \"\\n\";\n    }\n    \n    if(inputs.inverse_normal)\n    {\n        config += \"inverse_normal \" + inputs.inverse_normal + \"\\n\";\n    }\n\n    if(inputs.output_prefix)\n    {\n        config += \"out_prefix \\\"\" + inputs.output_prefix + \"\\\"\\n\";\n    }\n    else\n    {\n        config += \"out_prefix \\\"null_model\\\"\" + \"\\n\";\n    }\n    \n    if(inputs.n_categories_boxplot){\n        \n        config += \"n_categories_boxplot \" + inputs.n_categories_boxplot + \"\\n\";\n    }\n    return config\n}",
                                            "writable": false
                                        },
                                        {
                                            "entry": "${\n    return inputs.null_model_params\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.phenotype_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.sample_include_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.pca_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.relatedness_matrix_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.null_model_files\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.conditional_variant_file\n}",
                                            "writable": true
                                        },
                                        {
                                            "entry": "${\n    return inputs.gds_files\n}",
                                            "writable": true
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105817,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/null-model-report/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105817,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105817,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a50b396a1895f88f2ab57b95903db3710b39e88c9789447c03220746cff82936f"
                        },
                        "label": "Null Model Report",
                        "sbg:x": 537,
                        "sbg:y": -212.5
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:AWSInstanceType",
                        "value": "c5.2xlarge;ebs-gp2;512"
                    },
                    {
                        "class": "sbg:AzureInstanceType",
                        "value": "Standard_D8s_v4;PremiumSSD;512"
                    }
                ],
                "requirements": [
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105817,
                        "sbg:revisionNotes": "Workflow decomposed"
                    }
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/null-model/0.png",
                "sbg:toolAuthor": "TOPMed DCC",
                "sbg:license": "MIT",
                "sbg:links": [
                    {
                        "id": "https://github.com/UW-GAC/analysis_pipeline",
                        "label": "Source Code, Download"
                    },
                    {
                        "id": "doi.org/10.1093/bioinformatics/btz567",
                        "label": "Publication"
                    },
                    {
                        "id": "https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test.html",
                        "label": "Home Page"
                    },
                    {
                        "id": "https://bioconductor.org/packages/devel/bioc/manuals/GENESIS/man/GENESIS.pdf",
                        "label": "Documentation"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "CWL1.0",
                    "Genomics"
                ],
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2",
                    "v1.1"
                ],
                "sbg:id": "markoz/hgi/null-model/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": "Workflow decomposed",
                "sbg:modifiedOn": 1637105817,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105817,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a0a5f6fe9fcbebf8f52196d8342a3d37de43b93df6e5141868aab18eab4509460"
            },
            "label": "GENESIS Null Model",
            "sbg:x": 1054.6588134765625,
            "sbg:y": -323.4139404296875
        },
        {
            "id": "single_variant_association_testing",
            "in": [
                {
                    "id": "genome_build",
                    "default": "hg19"
                },
                {
                    "id": "test_type",
                    "default": "score.spa"
                },
                {
                    "id": "phenotype_file",
                    "source": "null_model/null_model_phenotypes"
                },
                {
                    "id": "null_model_file",
                    "source": "null_model/null_model_file",
                    "pickValue": "first_non_null"
                },
                {
                    "id": "maf_threshold",
                    "source": "maf_threshold"
                },
                {
                    "id": "mac_threshold",
                    "source": "mac_threshold"
                },
                {
                    "id": "out_prefix",
                    "default": "single_variant"
                },
                {
                    "id": "input_gds_files",
                    "linkMerge": "merge_flattened",
                    "source": [
                        "vcf_to_gds_1/unique_variant_id_gds_per_chr"
                    ]
                },
                {
                    "id": "variant_include_files",
                    "source": [
                        "hwe_filter_cwl1/variants_kept"
                    ]
                },
                {
                    "id": "signif_line_fixed",
                    "source": "signif_line_fixed"
                }
            ],
            "out": [
                {
                    "id": "assoc_combined"
                },
                {
                    "id": "assoc_plots"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/single-variant-association-testing/4",
                "doc": "**Single Variant workflow** runs single-variant association tests. It consists of several steps. Define Segments divides genome into segments, either by a number of segments, or by a segment length. Note that number of segments refers to whole genome, not a number of segments per chromosome. Association test is then performed for each segment in parallel, before combining results on chromosome level. Final step produces QQ and Manhattan plots.\n\nThis workflow uses the output from a model fit using the null model workflow to perform score tests for all variants individually. The reported effect estimate is for the alternate allele, and multiple alternate alleles for a single variant are tested separately.\n\nWhen testing a binary outcome, the saddlepoint approximation (SPA) for p-values [1][2] can be used by specifying **Test type** = ‘score.spa’; this is generally recommended. SPA will provide better calibrated p-values, particularly for rarer variants in samples with case-control imbalance. \n\nWhen testing a binary outcome, the BinomiRare test is also available[3]. This is a “carriers only” exact test that compares the observed number of variant carriers who are cases to the expected number based on the probability of being a case under the null hypothesis of no association between outcome and variant. This test may be useful when testing association of very rare variants with rare outcomes.\n\nWhen using the test_type ‘score‘ or ‘score.spa‘, the fast approximation to the score standard error (SE) implemented by Zhou et al. (2018) in their SAIGE software [5] is available by using a null_model_file prepared with the Update Null Model for Fast Score Test workflow. This approximation may be much faster than computing the true score SE in large samples, as it replaces the full covariance matrix in the calculation with the product of a diagonal matrix and a scalar correction factor (se.correction) in the updated null model output.\n\nIf your genotype data has sporadic missing values, they are mean imputed using the allele frequency observed in the sample.\n\nOn the X chromosome, males have genotype values coded as 0/2 (females as 0/1/2).\n\nThis workflow utilizes the *assocTestSingle* function from the GENESIS software [4].\n\n### Common Use Cases\n\nSingle Variant Association Testing workflow is designed to run single-variant association tests using GENESIS software. Set of variants on which to run association testing can be reduced by providing **Variant Include Files** - One file per chromosome containing variant IDs for variants on which association testing will be performed.\n\n\n### Common issues and important notes\n* Association Testing - Single job can be very memory demanding, depending on number of samples and null model used. We suggest running with at least 5GB of memory allocated for small studies, and to use approximation of 0.5GB per thousand samples for larger studies (with more than 10k samples), but this again depends on complexity of null model. If a run fails with *error 137*, and with message killed, most likely cause is lack of memory. Memory can be allocated using the **memory GB** parameter.\n\n* This workflow expects **GDS** files split by chromosome, and will not work otherwise. If provided, **variant include** files must also be split in the same way. Also GDS and Variant include files should be properly named. It is expected that chromosome is included in the filename in following format: chr## , where ## is the name of the chromosome (1-24 or X, Y). Chromosome can be included at any part of the filename. Examples for GDS: data_subset_chr1.gds, data_chr1_subset.gds. Examples for Variant include files: variant_include_chr1.RData, chr1_variant_include.RData.\n\n* Some input arguments are mutually exclusive, for more information, please visit workflow [github page](https://github.com/UW-GAC/analysis_pipeline/tree/v2.5.0)\n\n### Changes introduced by Seven Bridges\nThere are no changes introduced by Seven Bridges.\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running time and cost. \n        \n\n| Samples &nbsp; &nbsp; |    | Rel. Matrix &nbsp; &nbsp;|Parallel instances &nbsp; &nbsp; | Instance type  &nbsp; &nbsp; &nbsp; &nbsp;| Spot/On Dem. &nbsp; &nbsp; |CPU &nbsp; &nbsp; | RAM &nbsp; &nbsp; | Time  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Cost |\n|--------------------|---------------|----------------|------------------------|--------------------|--------------------|--------|--------|---------|-------|\n| 2.5k   |                 |   w/o          | 8                           |  r4.8xlarge | Spot     |1  | 2   | 1 h, 8 min   | $5  |\n| 2.5k   |               |   Dense     | 8                           |  r4.8xlarge | Spot     |1  | 2   | 1 h, 8 min   | $5  |\n| 10k   |                 |   w/o           | 8                           |  c5.18xlarge | On Demand     |1  | 2   | 50 min   | $10  |\n| 10k   |                |   Sparse     | 8                           |  c5.18xlarge | On Demand     |1  | 2   | 58 min   | $11  |\n| 10k   |                |   Sparse     | 8                           |  r4.8xlarge | On Demand     |1  | 2   | 1 h, 30 min   | $11  |\n| 10k   |                 |   Dense      | 8                           |  r5.4xlarge | On Demand     |1  | 8   | 3 h   | $24  |\n| 36k  |                  |   w/o           | 8                           |  r5.4xlarge | On Demand     |1  | 5   | 3 h, 20 min   | $27  |\n| 36k  |                  |   Sparse     | 8                           |  r5.4xlarge | On Demand     |1  | 5   | 4 h   | $32  |\n| 36k   |                  |   Sparse     | 8                           |  r5.12xlarge | On Demand     |1  | 5   | 1 h, 20 min   | $32  |\n| 36k   |                  |   Dense      | 8                           |  r5.12xlarge | On Demand     |1  | 50   | 1 d, 15 h   | $930  |\n| 36k   |                 |   Dense      | 8                           |  r5.24xlarge | On Demand     |1  | 50   | 17 h   | $800  |\n| 50k   |                  |   w/o           | 8                           |  r5.12xlarge | On Demand     |1  | 8   | 2 h   | $44  |\n| 50k   |                  |   Sparse     | 8                           |  r5.12xlarge | On Demand     |1  | 8   | 2 h   | $48 |\n| 50k   |                  |   Dense      | 8                           |  r5.24xlarge | On Demand     |48  | 100   | 11 d   | $13500  |\n| 2.5k   |                  |   w/o          | 8                           |  n1-standard-64 | Preemptible    |1  | 2   | 1 h   | $7  |\n| 2.5k   |                  |   Dense     | 8                           |  n1-standard-64 | Preemptible    |1  | 2   | 1 h   | $7  |\n| 10k   |                  |   w/o           | 8                           |  n1-standard-4 | On Demand     |1  | 2   | 1 h, 12 min  | $13  |\n| 10k   |                  |   Sparse     | 8                           |  n1-standard-4 | On Demand     |1  | 2   | 1 h, 13  min   | $14 |\n| 10k  |                  |   Dense      | 8                           |  n1-highmem-32 | On Demand     |1  | 8   | 2 h, 20  min   | $30  |\n| 36k   |                  |   w/o           | 8                           |  n1-standard-64 | On Demand     |1  | 5   | 1 h, 30  min   | $35  |\n| 36k   |                 |   Sparse     | 8                           |  n1-highmem-16 | On Demand     |1  | 5   | 4 h, 30  min   | $35  |\n| 36k   |                  |   Sparse     | 8                           |  n1-standard-64 | On Demand     |1  | 5   | 1 h, 30  min   | $35  |\n| 36k   |                  |   Dense      | 8                           |  n1-highmem-96 | On Demand     |1  | 50   | 1 d, 6  h   | $1300  |\n| 50k   |                  |   w/o           | 8                           |  n1-standard-96 | On Demand     |1  | 8    | 2  h   | $73  |\n| 50k   |                  |   Sparse     | 8                           |  n1-standard-96 | On Demand     |1  | 8    | 2  h   | $73  |\n| 50k   |                  |   Dense      | 8                           |  n1-highmem-96 | On Demand     |16  | 100    | 6  d   | $6600  |\n\nIn tests performed we used 1000G (tasks with 2.5k participants) and TOPMed freeze5 datasets (tasks with 10k or more participants). \nAll these tests are done with applied **MAF >= 1%** filter. The number of variants that have been tested is **14 mio in 1000G** and **12 mio in TOPMed freeze 5** dataset. \n\nAlso, a common filter in these analysis is **MAC>=5**. In that case the number of variants would be **32 mio for 1000G** and **92 mio for TOPMed freeze5** data. Since for single variant testing, the compute time grows linearly with the number of variants tested the execution time and price can be easily estimated from the results above.\n\n*For more details on **spot/preemptible instances** please visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances).*   \n\n\n### API Python Implementation\n\nThe app's draft task can also be submitted via the **API**. In order to learn how to get your **Authentication token** and **API endpoint** for the corresponding Platform visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).\n\n```python\nfrom sevenbridges import Api\n\nauthentication_token, api_endpoint = \"enter_your_token\", \"enter_api_endpoint\"\napi = Api(token=authentication_token, url=api_endpoint)\n# Get project_id/app_id from your address bar. Example: https://f4c.sbgenomics.com/u/your_username/project/app\nproject_id, app_id = \"your_username/project\", \"your_username/project/app\"\n# Get file names from files in your project. Example: Names are taken from Data/Public Reference Files.\ninputs = {\n    \"input_gds_files\": api.files.query(project=project_id, names=[\"basename_chr1.gds\", \"basename_chr2.gds\", ..]),\n    \"phenotype_file\": api.files.query(project=project_id, names=[\"name_of_phenotype_file\"])[0],\n    \"null_model_file\": api.files.query(project=project_id, names=[\"name_of_null_model_file\"])[0]\n}\ntask = api.tasks.create(name='Single Variant Association Testing - API Run', project=project_id, app=app_id, inputs=inputs, run=False)\n```\nInstructions for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation). For more information about using the API Python client, consult [sevenbridges-python documentation](http://sevenbridges-python.readthedocs.io/en/latest/). **More examples** are available [here](https://github.com/sbg/okAPI).\n\nAdditionally, [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java) clients are available. To learn more about using these API clients please refer to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).\n\n\n### References\n - [1] [SaddlePoint Approximation (SPA)](https://doi.org/10.1016/j.ajhg.2017.05.014)  \n - [2] [SPA - additional reference](https://doi.org/10.1038/s41588-018-0184-y)  \n - [3] [BinomiRare](https://pubmed.ncbi.nlm.nih.gov/28393384/)  \n - [4] [GENESIS toolkit](doi.org/10.1093/bioinformatics/btz567/) \n - [5] [Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies](https://www.nature.com/articles/s41588-018-0184-y)",
                "label": "GENESIS Single Variant Association Testing",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "segment_length",
                        "type": "int?",
                        "label": "Segment length",
                        "doc": "Segment length in kb, used for parallelization.",
                        "sbg:toolDefaultValue": "10000kb",
                        "sbg:x": -361,
                        "sbg:y": -204
                    },
                    {
                        "id": "n_segments",
                        "type": "int?",
                        "label": "Number of segments",
                        "doc": "Number of segments, used for parallelization (overrides Segment length). Note that this parameter defines the number of segments for the entire genome, so using this argument with selected chromosomes may result in fewer segments than you expect (and the minimum is one segment per chromosome).",
                        "sbg:x": -484,
                        "sbg:y": -88
                    },
                    {
                        "id": "genome_build",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "hg19",
                                    "hg38"
                                ],
                                "name": "genome_build"
                            }
                        ],
                        "label": "Genome build",
                        "doc": "Genome build for the genotypes in the GDS file (hg19 or hg38). Used to divide the genome into segments for parallel processing.",
                        "sbg:toolDefaultValue": "hg38",
                        "sbg:x": -363,
                        "sbg:y": 6
                    },
                    {
                        "id": "variant_block_size",
                        "type": "int?",
                        "label": "Variant block size",
                        "doc": "Number of variants to read from the GDS file in a single block. For smaller sample sizes, increasing this value will reduce the number of iterations in the code. For larger sample sizes, values that are too large will result in excessive memory use.",
                        "sbg:toolDefaultValue": "1024",
                        "sbg:x": 58.5833740234375,
                        "sbg:y": 42.84904098510742
                    },
                    {
                        "id": "test_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "score",
                                    "score.spa",
                                    "BinomiRare"
                                ],
                                "name": "test_type"
                            }
                        ],
                        "label": "Test type",
                        "doc": "Type of association test to perform. “Score” performs a score test and can be used with any null model. “Score.spa” uses the saddle point approximation (SPA) to provide more accurate p-values, especially for rare variants, in samples with unbalanced case:control ratios; “score.spa” can only be used if the null model family is “binomial”. “BinomiRare” is a carriers only exact test that may perform better when testing very rare variants with rare outcomes; “BinomiRare” can only be used if the null model family is “binomial”.",
                        "sbg:toolDefaultValue": "score",
                        "sbg:x": -50,
                        "sbg:y": 93
                    },
                    {
                        "id": "phenotype_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File",
                        "label": "Phenotype file",
                        "doc": "RData file with an AnnotatedDataFrame of phenotypes and covariates. Sample identifiers must be in column named “sample.id”. It is recommended to use the phenotype file output by the GENESIS Null Model app.",
                        "sbg:x": 57,
                        "sbg:y": 150
                    },
                    {
                        "id": "pass_only",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "TRUE",
                                    "FALSE"
                                ],
                                "name": "pass_only"
                            }
                        ],
                        "label": "Pass only",
                        "doc": "TRUE to select only variants with FILTER=PASS. If FALSE, variants that failed the quality filter will be included in the test.",
                        "sbg:toolDefaultValue": "TRUE",
                        "sbg:x": -49,
                        "sbg:y": 202
                    },
                    {
                        "id": "null_model_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File",
                        "label": "Null model file",
                        "doc": "RData file containing a null model object. Run the GENESIS Null Model app to create this file. A null model object created with the GENESIS Update Null Model for Fast Score Test app can also be used.",
                        "sbg:x": 60,
                        "sbg:y": 268
                    },
                    {
                        "id": "memory_gb",
                        "type": "float?",
                        "label": "memory GB",
                        "doc": "Memory in GB per job.",
                        "sbg:toolDefaultValue": "8",
                        "sbg:x": -44,
                        "sbg:y": 319
                    },
                    {
                        "id": "maf_threshold",
                        "type": "float?",
                        "label": "MAF threshold",
                        "doc": "Minimum minor allele frequency for variants to include in test. Only used if MAC threshold is NA.",
                        "sbg:toolDefaultValue": "0.001",
                        "sbg:x": 59,
                        "sbg:y": 383
                    },
                    {
                        "id": "mac_threshold",
                        "type": "float?",
                        "label": "MAC threshold",
                        "doc": "Minimum minor allele count for variants to include in test. Recommend to use a higher threshold when outcome is binary or count data. To disable it set it to NA.",
                        "sbg:toolDefaultValue": "5",
                        "sbg:x": -42,
                        "sbg:y": 432
                    },
                    {
                        "id": "cpu",
                        "type": "int?",
                        "label": "CPU",
                        "doc": "Number of CPUs for each job.",
                        "sbg:toolDefaultValue": "1",
                        "sbg:x": -46.285701751708984,
                        "sbg:y": 548.0167846679688
                    },
                    {
                        "id": "disable_thin",
                        "type": "boolean?",
                        "label": "Disable Thin",
                        "doc": "Logical for whether to thin points in the QQ and Manhattan plots. By default, points are thinned in dense regions to reduce plotting time. If this parameter is set to TRUE, all variant p-values will be included in the plots, and the plotting will be very long and memory intensive.",
                        "sbg:toolDefaultValue": "TRUE",
                        "sbg:x": 952.9915771484375,
                        "sbg:y": 432.899169921875
                    },
                    {
                        "id": "known_hits_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Known hits file",
                        "doc": "RData file with data.frame containing columns chr and pos. If provided, 1 Mb regions surrounding each variant listed will be omitted from the QQ and manhattan plots.",
                        "sbg:x": 1070,
                        "sbg:y": 407
                    },
                    {
                        "id": "thin_nbins",
                        "type": "int?",
                        "label": "Thin N bins",
                        "doc": "Number of bins to use for thinning.",
                        "sbg:toolDefaultValue": "10",
                        "sbg:x": 1052.4117431640625,
                        "sbg:y": 253.57142639160156
                    },
                    {
                        "id": "thin_npoints",
                        "type": "int?",
                        "label": "Thin N points",
                        "doc": "Number of points in each bin after thinning.",
                        "sbg:toolDefaultValue": "10000",
                        "sbg:x": 934.3193359375,
                        "sbg:y": 234
                    },
                    {
                        "id": "out_prefix",
                        "type": "string",
                        "label": "Output prefix",
                        "doc": "Prefix that will be included in all output files.",
                        "sbg:x": 74.62183380126953,
                        "sbg:y": 502.3025207519531
                    },
                    {
                        "id": "input_gds_files",
                        "sbg:fileTypes": "GDS",
                        "type": "File[]",
                        "label": "GDS files",
                        "doc": "GDS files with genotype data for variants to be tested for association. If multiple files are selected, they will be run in parallel. Files separated by chromosome are expected to have ‘chr##’ strings indicating chromosome number, where ‘##’ can be (1-24, X, Y). Output files for each chromosome will include the corresponding chromosome number.",
                        "sbg:x": -358.317626953125,
                        "sbg:y": -345.1597595214844
                    },
                    {
                        "id": "truncate_pval_threshold",
                        "type": "float?",
                        "label": "Truncate pval threshold",
                        "doc": "Maximum p-value to display in truncated QQ and manhattan plots.",
                        "sbg:x": 1158.9296875,
                        "sbg:y": 598.3650512695312
                    },
                    {
                        "id": "plot_mac_threshold",
                        "type": "int?",
                        "label": "Plot MAC threshold",
                        "doc": "Minimum minor allele count for variants or aggregate units to include in plots (if different from MAC threshold).",
                        "sbg:x": 1044.307861328125,
                        "sbg:y": 560.8524169921875
                    },
                    {
                        "id": "variant_include_files",
                        "sbg:fileTypes": "RData",
                        "type": "File[]?",
                        "label": "Variant Include Files",
                        "doc": "RData file containing ids of variants to be included.",
                        "sbg:x": 13.8739652633667,
                        "sbg:y": -437.260498046875
                    },
                    {
                        "id": "signif_line_fixed",
                        "type": "float?",
                        "label": "Significance line",
                        "doc": "P-value for the significance line. Only used if `signif_type = fixed`.",
                        "sbg:exposed": true
                    }
                ],
                "outputs": [
                    {
                        "id": "assoc_combined",
                        "outputSource": [
                            "assoc_combine_r/assoc_combined"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File[]?",
                        "label": "Association test results",
                        "doc": "RData file with data.frame of association test results (test statistic, p-value, etc.) See the documentation of the GENESIS R package for detailed description of output.",
                        "sbg:x": 1370.5833740234375,
                        "sbg:y": -3.150960683822632
                    },
                    {
                        "id": "assoc_plots",
                        "outputSource": [
                            "assoc_plots_r/assoc_plots"
                        ],
                        "sbg:fileTypes": "PNG",
                        "type": "File[]?",
                        "label": "Association test plots",
                        "doc": "QQ and Manhattan Plots of p-values in association test results.",
                        "sbg:x": 1577.5833740234375,
                        "sbg:y": 331.8490295410156
                    }
                ],
                "steps": [
                    {
                        "id": "define_segments_r",
                        "in": [
                            {
                                "id": "segment_length",
                                "source": "segment_length"
                            },
                            {
                                "id": "n_segments",
                                "source": "n_segments"
                            },
                            {
                                "id": "genome_build",
                                "source": "genome_build"
                            }
                        ],
                        "out": [
                            {
                                "id": "config"
                            },
                            {
                                "id": "define_segments_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/define-segments-r/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:altPrefix": "-s",
                                    "sbg:toolDefaultValue": "10000",
                                    "sbg:category": "Optional parameters",
                                    "id": "segment_length",
                                    "type": "int?",
                                    "inputBinding": {
                                        "prefix": "--segment_length",
                                        "shellQuote": false,
                                        "position": 1
                                    },
                                    "label": "Segment length",
                                    "doc": "Segment length in kb, used for paralelization."
                                },
                                {
                                    "sbg:altPrefix": "-n",
                                    "sbg:category": "Optional parameters",
                                    "id": "n_segments",
                                    "type": "int?",
                                    "inputBinding": {
                                        "prefix": "--n_segments",
                                        "shellQuote": false,
                                        "position": 2
                                    },
                                    "label": "Number of segments",
                                    "doc": "Number of segments, used for paralelization (overrides Segment length). Note that this parameter defines the number of segments for the entire genome, so using this argument with selected chromosomes may result in fewer segments than you expect (and the minimum is one segment per chromosome)."
                                },
                                {
                                    "sbg:toolDefaultValue": "hg38",
                                    "sbg:category": "Configs",
                                    "id": "genome_build",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "hg19",
                                                "hg38"
                                            ],
                                            "name": "genome_build"
                                        }
                                    ],
                                    "label": "Genome build",
                                    "doc": "Genome build for the genotypes in the GDS file (hg19 or hg38). Used to divide the genome into segments for parallel processing.",
                                    "default": "hg38"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "config",
                                    "doc": "Config file.",
                                    "label": "Config file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.config"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                },
                                {
                                    "id": "define_segments_output",
                                    "doc": "Segments txt file.",
                                    "label": "Segments file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "segments.txt"
                                    },
                                    "sbg:fileTypes": "TXT"
                                }
                            ],
                            "label": "define_segments.R",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "separate": false,
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "define_segments.config"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "Rscript /usr/local/analysis_pipeline/R/define_segments.R"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "define_segments.config",
                                            "entry": "${\n    var argument = [];\n    argument.push('out_file \"segments.txt\"')\n    if(inputs.genome_build){\n         argument.push('genome_build \"' + inputs.genome_build + '\"')\n    }\n    return argument.join('\\n')\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105820,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/define-segments-r/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105820,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105820,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "abcb7884a4e9f96eab06afefcfd6ac9a971605d3a26b810578009f05e0f63455d"
                        },
                        "label": "define_segments.R",
                        "sbg:x": -199.3984375,
                        "sbg:y": -88
                    },
                    {
                        "id": "assoc_single_r",
                        "in": [
                            {
                                "id": "gds_file",
                                "linkMerge": "merge_flattened",
                                "source": [
                                    "sbg_prepare_segments_1/gds_output"
                                ],
                                "valueFrom": "$(self ? [].concat(self)[0] : self)"
                            },
                            {
                                "id": "null_model_file",
                                "source": "null_model_file"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "mac_threshold",
                                "source": "mac_threshold"
                            },
                            {
                                "id": "maf_threshold",
                                "source": "maf_threshold"
                            },
                            {
                                "id": "pass_only",
                                "source": "pass_only"
                            },
                            {
                                "id": "segment_file",
                                "linkMerge": "merge_flattened",
                                "source": [
                                    "define_segments_r/define_segments_output"
                                ],
                                "valueFrom": "$(self ? [].concat(self)[0] : self)"
                            },
                            {
                                "id": "test_type",
                                "source": "test_type"
                            },
                            {
                                "id": "variant_include_file",
                                "linkMerge": "merge_flattened",
                                "source": [
                                    "sbg_prepare_segments_1/variant_include_output"
                                ],
                                "valueFrom": "$(self ? [].concat(self)[0] : self)"
                            },
                            {
                                "id": "segment",
                                "linkMerge": "merge_flattened",
                                "source": [
                                    "sbg_prepare_segments_1/segments"
                                ],
                                "valueFrom": "$(self ? [].concat(self)[0] : self)"
                            },
                            {
                                "id": "memory_gb",
                                "default": 8,
                                "source": "memory_gb"
                            },
                            {
                                "id": "cpu",
                                "source": "cpu"
                            },
                            {
                                "id": "variant_block_size",
                                "source": "variant_block_size"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "genome_build",
                                "source": "genome_build"
                            }
                        ],
                        "out": [
                            {
                                "id": "configs"
                            },
                            {
                                "id": "assoc_single"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/assoc-single-r/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Configs",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS file",
                                    "doc": "GDS file.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "null_model_file",
                                    "type": "File",
                                    "label": "Null model file",
                                    "doc": "Null model file.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "phenotype_file",
                                    "type": "File",
                                    "label": "Phenotype file",
                                    "doc": "RData file with AnnotatedDataFrame of phenotypes.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:toolDefaultValue": "5",
                                    "sbg:category": "Configs",
                                    "id": "mac_threshold",
                                    "type": "float?",
                                    "label": "MAC threshold",
                                    "doc": "Minimum minor allele count for variants to include in test. Use a higher threshold when outcome is binary. To disable it set it to NA. Tool default: 5."
                                },
                                {
                                    "sbg:toolDefaultValue": "0.001",
                                    "sbg:category": "Configs",
                                    "id": "maf_threshold",
                                    "type": "float?",
                                    "label": "MAF threshold",
                                    "doc": "Minimum minor allele frequency for variants to include in test. Only used if MAC threshold is NA. Tool default: 0.001."
                                },
                                {
                                    "sbg:toolDefaultValue": "TRUE",
                                    "sbg:category": "Configs",
                                    "id": "pass_only",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "TRUE",
                                                "FALSE"
                                            ],
                                            "name": "pass_only"
                                        }
                                    ],
                                    "label": "Pass only",
                                    "doc": "TRUE to select only variants with FILTER=PASS."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "segment_file",
                                    "type": "File?",
                                    "label": "Segment file",
                                    "doc": "Segment file.",
                                    "sbg:fileTypes": "TXT"
                                },
                                {
                                    "sbg:toolDefaultValue": "score",
                                    "sbg:category": "Configs",
                                    "id": "test_type",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "score",
                                                "score.spa",
                                                "BinomiRare"
                                            ],
                                            "name": "test_type"
                                        }
                                    ],
                                    "label": "Test type",
                                    "doc": "Type of test to perform. If samples are related (mixed model), options are score if binary is FALSE, score and score.spa if binary is TRUE. Default value: score."
                                },
                                {
                                    "sbg:category": "Configs",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant include file",
                                    "doc": "RData file with vector of variant.id to include.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:altPrefix": "-c",
                                    "sbg:category": "Optional inputs",
                                    "id": "chromosome",
                                    "type": "string?",
                                    "inputBinding": {
                                        "prefix": "--chromosome",
                                        "shellQuote": false,
                                        "position": 1
                                    },
                                    "label": "Chromosome",
                                    "doc": "Chromosome (1-24 or X,Y)."
                                },
                                {
                                    "sbg:category": "Optional parameters",
                                    "id": "segment",
                                    "type": "int?",
                                    "inputBinding": {
                                        "prefix": "--segment",
                                        "shellQuote": false,
                                        "position": 2
                                    },
                                    "label": "Segment number",
                                    "doc": "Segment number."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "8",
                                    "id": "memory_gb",
                                    "type": "float?",
                                    "label": "memory GB",
                                    "doc": "Memory in GB per job. Default value: 8."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "cpu",
                                    "type": "int?",
                                    "label": "CPU",
                                    "doc": "Number of CPUs for each tool job. Default value: 1."
                                },
                                {
                                    "sbg:category": "General",
                                    "sbg:toolDefaultValue": "1024",
                                    "id": "variant_block_size",
                                    "type": "int?",
                                    "label": "Variant block size",
                                    "doc": "Number of variants to read in a single block. Default: 1024"
                                },
                                {
                                    "sbg:toolDefaultValue": "assoc_single",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Output prefix"
                                },
                                {
                                    "id": "genome_build",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "hg19",
                                                "hg38"
                                            ],
                                            "name": "genome_build"
                                        }
                                    ],
                                    "label": "Genome build",
                                    "doc": "Genome build for the genotypes in the GDS file (hg19 or hg38). Used to divide the genome into segments for parallel processing.",
                                    "default": "hg38"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "configs",
                                    "doc": "Config files.",
                                    "label": "Config files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*config*"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                },
                                {
                                    "id": "assoc_single",
                                    "doc": "Assoc single output.",
                                    "label": "Assoc single output",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.RData",
                                        "outputEval": "$(inheritMetadata(self, inputs.gds_file))"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "assoc_single.R",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "separate": false,
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    if(inputs.is_unrel)\n    {\n        return \"assoc_single_unrel.config\"\n    }\n    else\n    {\n        return \"assoc_single.config\"\n    }\n    \n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "${\n    if(inputs.is_unrel)\n    {\n        return \"Rscript /usr/local/analysis_pipeline/R/assoc_single_unrel.R\"\n    }\n    else\n    {\n        return \"Rscript /usr/local/analysis_pipeline/R/assoc_single.R\"\n    }\n    \n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    if (inputs.cpu)\n        return 'export NSLOTS=' + inputs.cpu + ' &&'\n    else\n        return ''\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": "${\n    if(inputs.memory_gb)\n        return parseFloat(inputs.memory_gb * 1024)\n    else\n        return 8*1024\n}",
                                    "coresMin": "${ if(inputs.cpu)\n        return inputs.cpu \n    else \n        return 1\n}"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "assoc_single.config",
                                            "entry": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n        \n    var chr = find_chromosome(inputs.gds_file.path)\n    \n    var argument = [];\n    if(!inputs.is_unrel)\n    {   \n        if(inputs.out_prefix){\n            argument.push(\"out_prefix \\\"\" + inputs.out_prefix + \"_chr\"+chr + \"\\\"\");\n        }\n        if(!inputs.out_prefix){\n        var data_prefix = inputs.gds_file.basename.split('chr');\n        var data_prefix2 = inputs.gds_file.basename.split('.chr');\n        if(data_prefix.length == data_prefix2.length)\n            argument.push('out_prefix \"' + data_prefix2[0] + '_single_chr' + chr + inputs.gds_file.basename.split('chr'+chr)[1].split('.gds')[0] +'\"');\n        else\n            argument.push('out_prefix \"' + data_prefix[0] + 'single_chr' + chr +inputs.gds_file.basename.split('chr'+chr)[1].split('.gds')[0]+'\"');}\n        argument.push('gds_file \"' + inputs.gds_file.path +'\"');\n        argument.push('null_model_file \"' + inputs.null_model_file.path + '\"');\n        argument.push('phenotype_file \"' + inputs.phenotype_file.path + '\"');\n        if(inputs.mac_threshold){\n            argument.push('mac_threshold ' + inputs.mac_threshold);\n        }\n        if(inputs.maf_threshold){\n            argument.push('maf_threshold ' + inputs.maf_threshold);\n        }\n        if(inputs.pass_only){\n            argument.push('pass_only ' + inputs.pass_only);\n        }\n        if(inputs.segment_file){\n            argument.push('segment_file \"' + inputs.segment_file.path + '\"');\n        }\n        if(inputs.test_type){\n            argument.push('test_type \"' + inputs.test_type + '\"') ;\n        }\n        if(inputs.variant_include_file){\n            argument.push('variant_include_file \"' + inputs.variant_include_file.path + '\"');\n        }\n        if(inputs.variant_block_size){\n            argument.push('variant_block_size ' + inputs.variant_block_size);\n        }\n        if(inputs.genome_build){\n            argument.push('genome_build ' + inputs.genome_build);\n        }\n        \n        argument.push('');\n        return argument.join('\\n');\n    }\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement",
                                    "expressionLib": [
                                        "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                                    ]
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:AWSInstanceType",
                                    "value": "r4.8xlarge;ebs-gp2;500"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105821,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/assoc-single-r/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105821,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105821,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a94935ca3d553186c6d4ba57c35b4ad86ea00b2a2aeda3c3d3186763f473f4aa8"
                        },
                        "label": "Association Testing Single",
                        "scatter": [
                            "gds_file",
                            "variant_include_file",
                            "segment"
                        ],
                        "scatterMethod": "dotproduct",
                        "sbg:x": 418.52001953125,
                        "sbg:y": 79.32788848876953
                    },
                    {
                        "id": "assoc_combine_r",
                        "in": [
                            {
                                "id": "chromosome",
                                "source": [
                                    "sbg_group_segments_1/chromosome"
                                ],
                                "valueFrom": "$(self ? [].concat(self) : self)"
                            },
                            {
                                "id": "assoc_type",
                                "default": "single"
                            },
                            {
                                "id": "assoc_files",
                                "source": [
                                    "sbg_group_segments_1/grouped_assoc_files"
                                ],
                                "valueFrom": "$(self ? [].concat(self) : self)"
                            },
                            {
                                "id": "memory_gb",
                                "default": 8
                            },
                            {
                                "id": "cpu",
                                "default": 2
                            }
                        ],
                        "out": [
                            {
                                "id": "assoc_combined"
                            },
                            {
                                "id": "configs"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/assoc-combine-r/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:altPrefix": "-c",
                                    "sbg:category": "Optional inputs",
                                    "id": "chromosome",
                                    "type": "string[]?",
                                    "inputBinding": {
                                        "prefix": "--chromosome",
                                        "shellQuote": false,
                                        "position": 10
                                    },
                                    "label": "Chromosome",
                                    "doc": "Chromosome (1-24 or X,Y)."
                                },
                                {
                                    "id": "assoc_type",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "single",
                                            "aggregate",
                                            "window"
                                        ],
                                        "name": "assoc_type"
                                    },
                                    "label": "Association Type",
                                    "doc": "Type of association test: single, window or aggregate."
                                },
                                {
                                    "id": "assoc_files",
                                    "type": "File[]",
                                    "label": "Association files",
                                    "doc": "Association files to be combined.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Out Prefix",
                                    "doc": "Output prefix."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "4",
                                    "id": "memory_gb",
                                    "type": "float?",
                                    "label": "memory GB",
                                    "doc": "Memory in GB per one job. Default value: 4GB."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "cpu",
                                    "type": "int?",
                                    "label": "CPU",
                                    "doc": "Number of CPUs for each tool job. Default value: 1."
                                },
                                {
                                    "sbg:category": "General",
                                    "id": "conditional_variant_file",
                                    "type": "File?",
                                    "label": "Conditional variant file",
                                    "doc": "RData file with data frame of of conditional variants. Columns should include chromosome (or chr) and variant.id. The alternate allele dosage of these variants will be included as covariates in the analysis.",
                                    "sbg:fileTypes": "RData, RDATA"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "assoc_combined",
                                    "doc": "Assoc combined.",
                                    "label": "Assoc combined",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${\n    \n    //var input_files = [].concat(inputs.assoc_files);\n    //var first_filename = input_files[0].basename;\n    \n    //var chr = first_filename.split('_chr')[1].split('_')[0].split('.RData')[0];\n    \n    //return first_filename.split('chr')[0]+'chr'+chr+'.RData';\n    \n    return '*.RData'\n}",
                                        "outputEval": "$(inheritMetadata(self, inputs.assoc_files))"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "configs",
                                    "doc": "Config files.",
                                    "label": "Config files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*config*"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                }
                            ],
                            "label": "assoc_combine.R",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "assoc_combine.config"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 5,
                                    "valueFrom": "Rscript /usr/local/analysis_pipeline/R/assoc_combine.R"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "${\n    var command = '';\n    var i;\n    for(i=0; i<inputs.assoc_files.length; i++)\n        command += \"ln -s \" + inputs.assoc_files[i].path + \" \" + inputs.assoc_files[i].path.split(\"/\").pop() + \" && \"\n    \n    return command\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": "${\n    if(inputs.memory_gb)\n        return parseInt(inputs.memory_gb * 1024)\n    else\n        return 4*1024\n}",
                                    "coresMin": "${ if(inputs.cpu)\n        return inputs.cpu \n    else \n        return 1\n}"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "assoc_combine.config",
                                            "entry": "${\n    var argument = [];\n    argument.push('assoc_type \"'+ inputs.assoc_type + '\"');\n    var data_prefix = inputs.assoc_files[0].basename.split('_chr')[0];\n    if (inputs.out_prefix)\n    {\n        argument.push('out_prefix \"' + inputs.out_prefix+ '\"');\n    }\n    else\n    {\n        argument.push('out_prefix \"' + data_prefix+ '\"');\n    }\n    \n    if(inputs.conditional_variant_file){\n        argument.push('conditional_variant_file \"' + inputs.conditional_variant_file.path + '\"');\n    }\n    //if(inputs.assoc_files)\n    //{\n    //    arguments.push('assoc_files \"' + inputs.assoc_files[0].path + '\"')\n    //}\n    return argument.join('\\n') + '\\n'\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement",
                                    "expressionLib": [
                                        "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                                    ]
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105822,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/assoc-combine-r/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105822,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105822,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a9441836c8bc986fc185a4d0cacafb79eee2380bb33c56e1b49de6a4cabdbf4b8"
                        },
                        "label": "Association Combine",
                        "scatter": [
                            "chromosome",
                            "assoc_files"
                        ],
                        "scatterMethod": "dotproduct",
                        "sbg:x": 1087,
                        "sbg:y": 113
                    },
                    {
                        "id": "assoc_plots_r",
                        "in": [
                            {
                                "id": "assoc_files",
                                "linkMerge": "merge_flattened",
                                "source": [
                                    "assoc_combine_r/assoc_combined"
                                ],
                                "valueFrom": "$(self ? [].concat(self) : self)"
                            },
                            {
                                "id": "assoc_type",
                                "default": "single"
                            },
                            {
                                "id": "plots_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "disable_thin",
                                "source": "disable_thin"
                            },
                            {
                                "id": "known_hits_file",
                                "source": "known_hits_file"
                            },
                            {
                                "id": "thin_npoints",
                                "source": "thin_npoints"
                            },
                            {
                                "id": "thin_nbins",
                                "source": "thin_nbins"
                            },
                            {
                                "id": "plot_mac_threshold",
                                "source": "plot_mac_threshold"
                            },
                            {
                                "id": "truncate_pval_threshold",
                                "source": "truncate_pval_threshold"
                            },
                            {
                                "id": "signif_type",
                                "default": "fixed"
                            },
                            {
                                "id": "signif_line_fixed",
                                "source": "signif_line_fixed"
                            }
                        ],
                        "out": [
                            {
                                "id": "assoc_plots"
                            },
                            {
                                "id": "configs"
                            },
                            {
                                "id": "Lambdas"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/assoc-plots-r/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "assoc_files",
                                    "type": "File[]",
                                    "label": "Results from association testing",
                                    "doc": "Rdata files. Results from association testing workflow.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input options",
                                    "id": "assoc_type",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "single",
                                            "window",
                                            "aggregate"
                                        ],
                                        "name": "assoc_type"
                                    },
                                    "label": "Association Type",
                                    "doc": "Type of association test: single, window or aggregate"
                                },
                                {
                                    "sbg:toolDefaultValue": "1-23",
                                    "sbg:category": "Input options",
                                    "id": "chromosomes",
                                    "type": "string?",
                                    "label": "Chromosomes",
                                    "doc": "List of chromosomes. If not provided, in case of multiple files, it will be automatically generated with assumtion that files are in format *chr*.RData\nExample: 1 2 3"
                                },
                                {
                                    "sbg:toolDefaultValue": "plots",
                                    "sbg:category": "Input Options",
                                    "id": "plots_prefix",
                                    "type": "string?",
                                    "label": "Plots prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "disable_thin",
                                    "type": "boolean?",
                                    "label": "Disable Thin",
                                    "doc": "Logical for whether to thin points in the QQ and Manhattan plots. By default, points are thinned in dense regions to reduce plotting time. If this parameter is set to TRUE, all variant p-values will be included in the plots, and the plotting will be very long and memory intensive."
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "known_hits_file",
                                    "type": "File?",
                                    "label": "Known hits file",
                                    "doc": "RData file with data.frame containing columns chr and pos. If provided, 1 Mb regions surrounding each variant listed will be omitted from the QQ and manhattan plots.",
                                    "sbg:fileTypes": "RData, RDATA"
                                },
                                {
                                    "sbg:category": "General",
                                    "sbg:toolDefaultValue": "10000",
                                    "id": "thin_npoints",
                                    "type": "int?",
                                    "label": "Number of points in each bin after thinning",
                                    "doc": "Number of points in each bin after thinning."
                                },
                                {
                                    "sbg:toolDefaultValue": "10",
                                    "sbg:category": "General",
                                    "id": "thin_nbins",
                                    "type": "int?",
                                    "label": "Thin N binsNumber of bins to use for thinning",
                                    "doc": "Number of bins to use for thinning."
                                },
                                {
                                    "id": "plot_mac_threshold",
                                    "type": "int?",
                                    "label": "Plot MAC threshold",
                                    "doc": "Minimum minor allele count for variants or Minimum cumulative minor allele count for aggregate units to include in plots (if different from threshold used to run tests; see `mac_threshold`)."
                                },
                                {
                                    "id": "truncate_pval_threshold",
                                    "type": "float?",
                                    "label": "Truncate pval threshold",
                                    "doc": "Truncate pval threshold."
                                },
                                {
                                    "sbg:toolDefaultValue": "FALSE",
                                    "id": "plot_qq_by_chrom",
                                    "type": "boolean?",
                                    "label": "Plot qq by chromosome",
                                    "doc": "Logical indicator for whether to generate QQ plots faceted by chromosome."
                                },
                                {
                                    "id": "plot_include_file",
                                    "type": "File?",
                                    "label": "Plot include file",
                                    "doc": "RData file with vector of ids to include. See `TopmedPipeline::assocFilterByFile` for format requirements.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "signif_type",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "fixed",
                                                "bonferroni",
                                                "none"
                                            ],
                                            "name": "signif_type"
                                        }
                                    ],
                                    "label": "Significance type",
                                    "doc": "`fixed`, `bonferroni`, or `none`; character string for how to calculate the significance threshold. Default is `fixed` for single variant analysis and `bonferroni` for other analysis types."
                                },
                                {
                                    "sbg:toolDefaultValue": "5e-9",
                                    "id": "signif_line_fixed",
                                    "type": "float?",
                                    "label": "Significance line",
                                    "doc": "P-value for the significance line. Only used if `signif_type = fixed`."
                                },
                                {
                                    "id": "qq_mac_bins",
                                    "type": "string?",
                                    "label": "QQ MAC bins",
                                    "doc": "Space separated string of integers (e.g., `\"5 20 50\"`). If set, generate a QQ plot binned by the specified MAC thresholds. 0 and Infinity will automatically be added."
                                },
                                {
                                    "id": "qq_maf_bins",
                                    "type": "string?",
                                    "label": "QQ MAF bins",
                                    "doc": "Space separated string of minor allele frequencies (e.g., \"0.01 0.05 0.1\"). If set, generate a QQ plot binned by the specified minor allele frequencies. 0 and Infinity will automatically be added. Single variant tests only."
                                },
                                {
                                    "id": "lambda_quantiles",
                                    "type": "string?",
                                    "label": "Lambda quantiles",
                                    "doc": "Space separated string of quantiles at which to calculate genomic inflation lambda (e.g., “0.25 0.5 0.75”). If set, create a text file with lambda calculated at the specified quantiles stored in `out_file_lambdas`."
                                },
                                {
                                    "sbg:toolDefaultValue": "lambda.txt",
                                    "id": "out_file_lambdas",
                                    "type": "string?",
                                    "label": "Lambda outfile name",
                                    "doc": "File name of file to store lambda calculated at different quantiles. The default is `lambda.txt`."
                                },
                                {
                                    "sbg:toolDefaultValue": "1",
                                    "id": "plot_max_p",
                                    "type": "float?",
                                    "label": "Plot max p",
                                    "doc": "Maximum p-value to plot in QQ and Manhattan plots. Expected QQ values are still calculated using the full set of p-values."
                                },
                                {
                                    "id": "plot_maf_threshold",
                                    "type": "float?",
                                    "label": "Plot MAF threshold",
                                    "doc": "Minimum minor allele frequency for variants to include in plots. Ignored if `plot_mac_threshold` is specified. Single variant association tests only."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "assoc_plots",
                                    "doc": "QQ and Manhattan Plots generated by assoc_plots.R script.",
                                    "label": "Assoc plots",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.png"
                                    },
                                    "sbg:fileTypes": "PNG"
                                },
                                {
                                    "id": "configs",
                                    "doc": "Config files.",
                                    "label": "Config files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*config*"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                },
                                {
                                    "id": "Lambdas",
                                    "doc": "File to store lambda calculated at different quantiles.",
                                    "label": "File to store lambda calculated at different quantiles",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.txt"
                                    },
                                    "sbg:fileTypes": "TXT"
                                }
                            ],
                            "doc": "### Description\n\nThe UW-GAC GENESIS Association Result Plotting standalone app creates Manhattan and QQ plots from GENESIS association test results with additional filtering and stratification options available. This app is run automatically with default options set by the GENESIS Association Testing Workflows. Users can fine-tune the Manhattan and QQ plots by running this app separately, after one of the association testing workflows. The available options are:\n - Create QQ plots by chromosome.\n - Include a user-specified subset of the results in the plots.\n - Filter results to only those with MAC or MAF greater than a specified threshold.\n - Calculate genomic inflation lambda at various quantiles.\n - Specify the significance type and level.\n - Create QQ plots stratified by MAC or MAF.\n - Specify a maximum p-value to display on the plots.\n\n### Common use cases\n\nThe UW-GAC GENESIS Association Result Plotting standalone app creates Manhattan and QQ plots from GENESIS association test results with additional filtering and stratification options available.\n\n### Changes introduced by Seven Bridges\n\nNo changes introduced by Seven Bridges.",
                            "label": "GENESIS Association results plotting",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 5,
                                    "valueFrom": "assoc_file.config"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "Rscript /usr/local/analysis_pipeline/R/assoc_plots.R"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "${\n    var command = '';\n    var i;\n    for(i=0; i<inputs.assoc_files.length; i++)\n        command += \"ln -s \" + inputs.assoc_files[i].path + \" \" + inputs.assoc_files[i].path.split(\"/\").pop() + \" && \"\n    \n    return command\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": 64000,
                                    "coresMin": 1
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "assoc_file.config",
                                            "entry": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    var argument = [];\n    argument.push('out_prefix \"assoc_single\"');\n    var a_file = [].concat(inputs.assoc_files)[0];\n    var chr = find_chromosome(a_file.basename);\n    var path = a_file.path.split('chr'+chr);\n    var extension = path[1].split('.')[1];\n    \n     \n    if(inputs.plots_prefix){\n        argument.push('plots_prefix ' + inputs.plots_prefix);\n        argument.push('out_file_manh ' + inputs.plots_prefix + '_manh.png');\n        argument.push('out_file_qq ' + inputs.plots_prefix + '_qq.png');\n    }\n    else{\n        var data_prefix = path[0].split('/').pop();\n        argument.push('out_file_manh ' + data_prefix + 'manh.png');\n        argument.push('out_file_qq ' + data_prefix + 'qq.png');\n        argument.push('plots_prefix \"plots\"')\n    }\n    if(inputs.assoc_type){\n        argument.push('assoc_type ' + inputs.assoc_type)\n    }\n    \n    argument.push('assoc_file ' + '\"' + path[0].split('/').pop() + 'chr ' +path[1] + '\"')\n\n    if(inputs.chromosomes){\n        argument.push('chromosomes \"' + inputs.chromosomes + '\"')\n    }\n    else {\n        var chr_array = [];\n        var chrom_num;\n        var i;\n        for (var i = 0; i < inputs.assoc_files.length; i++) \n        {\n            chrom_num = inputs.assoc_files[i].path.split(\"/\").pop()\n            chrom_num = find_chromosome(chrom_num)\n            \n            chr_array.push(chrom_num)\n        }\n        \n        chr_array = chr_array.sort(function(a, b) { a.localeCompare(b, 'en', {numeric: true, ignorePunctuation: true})})\n        \n        var chrs = \"\";\n        for (var i = 0; i < chr_array.length; i++) \n        {\n            chrs += chr_array[i] + \" \"\n        }\n        argument.push('chromosomes \"' + chrs + '\"')\n    }\n    if(inputs.disable_thin){\n        argument.push('thin FALSE')\n    }\n    if(inputs.thin_npoints)\n        argument.push('thin_npoints ' + inputs.thin_npoints)\n    if(inputs.thin_npoints)\n        argument.push('thin_nbins ' + inputs.thin_nbins)\n    if(inputs.known_hits_file)\n        argument.push('known_hits_file \"' + inputs.known_hits_file.path + '\"')\n    if(inputs.plot_mac_threshold)\n        argument.push('plot_mac_threshold ' + inputs.plot_mac_threshold)  \n    if(inputs.truncate_pval_threshold)\n        argument.push('truncate_pval_threshold ' + inputs.truncate_pval_threshold)    \n    if(inputs.plot_qq_by_chrom){\n        argument.push('plot_qq_by_chrom ' + inputs.plot_qq_by_chrom)\n    }\n    if(inputs.plot_include_file){\n        argument.push('plot_include_file ' + '\"'+ inputs.plot_include_file.path + '\"')\n    }\n    if(inputs.signif_type){\n        argument.push('signif_type ' + inputs.signif_type)\n    }    \n    if(inputs.signif_line_fixed){\n        argument.push('signif_line_fixed ' + inputs.signif_line_fixed)\n    } \n    if(inputs.qq_mac_bins){\n        argument.push('qq_mac_bins ' + inputs.qq_mac_bins)\n    }\n    if(inputs.qq_maf_bins){\n        argument.push('qq_maf_bins ' + inputs.qq_maf_bins)\n    }    \n    if(inputs.lambda_quantiles){\n        argument.push('lambda_quantiles ' + inputs.lambda_quantiles)\n    }    \n    if(inputs.out_file_lambdas){\n        argument.push('out_file_lambdas ' + inputs.out_file_lambdas)\n    } \n    if(inputs.plot_max_p){\n        argument.push('plot_max_p ' + inputs.plot_max_p)\n    } \n    if(inputs.plot_maf_threshold){\n        argument.push('plot_maf_threshold ' + inputs.plot_maf_threshold)\n    }\n        \n        \n    argument.push('\\n')\n    return argument.join('\\n')\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:AzureInstanceType",
                                    "value": "Standard_D8s_v4;PremiumSSD;512"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105823,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/assoc-plots-r/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105823,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105823,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "ac69853666a83b66464719ab99d0f3422d6986b55831a7556e56b6d8311c6d61b"
                        },
                        "label": "Association Plots",
                        "sbg:x": 1367,
                        "sbg:y": 306
                    },
                    {
                        "id": "sbg_gds_renamer",
                        "in": [
                            {
                                "id": "in_variants",
                                "source": "input_gds_files"
                            }
                        ],
                        "out": [
                            {
                                "id": "renamed_variants"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/sbg-gds-renamer/0",
                            "baseCommand": [
                                "cp"
                            ],
                            "inputs": [
                                {
                                    "id": "in_variants",
                                    "type": "File",
                                    "label": "GDS input",
                                    "doc": "This tool removes suffix after 'chr##' in GDS filename. ## stands for chromosome name and can be (1-22,X,Y).",
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "renamed_variants",
                                    "doc": "Renamed GDS file.",
                                    "label": "Renamed GDS",
                                    "type": "File",
                                    "outputBinding": {
                                        "glob": "${\n    return '*'+inputs.in_variants.nameext\n}"
                                    },
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "doc": "This tool renames GDS file in GENESIS pipelines if they contain suffixes after chromosome (chr##) in the filename.\nFor example: If GDS file has name data_chr1_subset.gds the tool will rename GDS file to data_chr1.gds.",
                            "label": "SBG GDS renamer",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    if(inputs.in_variants){\n    return inputs.in_variants.path}\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n     function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n      function find_chromosome(file){\n         var chr_array = [];\n         var chrom_num = file.split(\"chr\")[1];\n        \n         if(isNumeric(chrom_num.charAt(1)))\n         {\n            chr_array.push(chrom_num.substr(0,2))\n         }\n         else\n         {\n            chr_array.push(chrom_num.substr(0,1))\n         }\n         return chr_array.toString()\n    }\n    \n    var chr = find_chromosome(inputs.in_variants.nameroot)\n    var base = inputs.in_variants.nameroot.split('chr'+chr)[0]\n    \n    return base+'chr' + chr + inputs.in_variants.nameext\n  \n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log && chmod -R 777 .' \n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105824,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/sbg-gds-renamer/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105824,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105824,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "aab167dd7a4bb3a5384595723197895ae486fdb06edfd1f4c39eec2299eced6c5"
                        },
                        "label": "SBG GDS renamer",
                        "scatter": [
                            "in_variants"
                        ],
                        "sbg:x": -138.8903045654297,
                        "sbg:y": -234.21176147460938
                    },
                    {
                        "id": "sbg_flatten_lists",
                        "in": [
                            {
                                "id": "input_list",
                                "source": [
                                    "assoc_single_r/assoc_single"
                                ],
                                "valueFrom": "${     var out = [];     for (var i = 0; i<self.length; i++){         if (self[i])    out.push(self[i])     }     return out }"
                            }
                        ],
                        "out": [
                            {
                                "id": "output_list"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/sbg-flatten-lists/0",
                            "baseCommand": [
                                "echo"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "File inputs",
                                    "id": "input_list",
                                    "type": "File[]?",
                                    "label": "Input list of files and lists",
                                    "doc": "List of inputs, can be any combination of lists of files and single files, it will be combined into a single list of files at the output."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "output_list",
                                    "doc": "Single list of files that combines all files from all inputs.",
                                    "label": "Output list of files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "outputEval": "${\n    function flatten(files) {\n        var a = [];\n        for (var i = 0; i < files.length; i++) {\n            if (files[i]) {\n                if (files[i].constructor == Array) a = a.concat(flatten(files[i]))\n                else a = a.concat(files[i])\n            }\n        }\n        return a\n    }\n\n    {\n        if (inputs.input_list) {\n            var arr = [].concat(inputs.input_list);\n            var return_array = [];\n            return_array = flatten(arr)\n            return return_array\n        }\n    }\n}"
                                    }
                                }
                            ],
                            "doc": "###**Overview** \n\nSBG FlattenLists is used to merge any combination of single file and list of file inputs into a single list of files. This is important because most tools and the CWL specification doesn't allow array of array types, and combinations of single file and array need to be converted into a single list for tools that can process a list of files.\n\n###**Input** \n\nAny combination of input nodes that are of types File or array of File, and any tool outputs that produce types File or array of File.\n\n###**Output** \n\nSingle array of File list containing all Files from all inputs combined, provided there are no duplicate files in those lists.\n\n###**Usage example** \n\nExample of usage is combining the outputs of two tools, one which produces a single file, and the other that produces an array of files, so that the next tool, which takes in an array of files, can process them together.",
                            "label": "SBG FlattenLists",
                            "arguments": [
                                {
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "\"Output"
                                },
                                {
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "is"
                                },
                                {
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "now"
                                },
                                {
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "a"
                                },
                                {
                                    "shellQuote": false,
                                    "position": 4,
                                    "valueFrom": "single"
                                },
                                {
                                    "shellQuote": false,
                                    "position": 5,
                                    "valueFrom": "list\""
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": 1000,
                                    "coresMin": 1
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.8.1"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        "$(inputs.input_list)"
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement",
                                    "expressionLib": [
                                        "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};"
                                    ]
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:toolAuthor": "Seven Bridges",
                            "sbg:cmdPreview": "echo \"Output is now a single list\"",
                            "sbg:image_url": null,
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105825,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:license": "Apache License 2.0",
                            "sbg:categories": [
                                "Other"
                            ],
                            "sbg:toolkit": "SBGTools",
                            "sbg:toolkitVersion": "1.0",
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/sbg-flatten-lists/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105825,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105825,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a8ab04a2a11a3f02f5cb29025dbeebbe3bb71cc8f1eb7caafb6e2140373cc62f3"
                        },
                        "label": "SBG FlattenLists",
                        "sbg:x": 684.666015625,
                        "sbg:y": 128.0019073486328
                    },
                    {
                        "id": "sbg_prepare_segments_1",
                        "in": [
                            {
                                "id": "input_gds_files",
                                "source": [
                                    "sbg_gds_renamer/renamed_variants"
                                ]
                            },
                            {
                                "id": "segments_file",
                                "source": "define_segments_r/define_segments_output"
                            },
                            {
                                "id": "variant_include_files",
                                "source": [
                                    "variant_include_files"
                                ]
                            }
                        ],
                        "out": [
                            {
                                "id": "gds_output"
                            },
                            {
                                "id": "segments"
                            },
                            {
                                "id": "aggregate_output"
                            },
                            {
                                "id": "variant_include_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/sbg-prepare-segments-1/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Inputs",
                                    "id": "input_gds_files",
                                    "type": "File[]",
                                    "label": "GDS files",
                                    "doc": "GDS files.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "segments_file",
                                    "type": "File",
                                    "label": "Segments file",
                                    "doc": "Segments file.",
                                    "sbg:fileTypes": "TXT"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "aggregate_files",
                                    "type": "File[]?",
                                    "label": "Aggregate files",
                                    "doc": "Aggregate files.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "variant_include_files",
                                    "type": "File[]?",
                                    "label": "Variant Include Files",
                                    "doc": "RData file containing ids of variants to be included.",
                                    "sbg:fileTypes": "RData"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "gds_output",
                                    "doc": "GDS files.",
                                    "label": "GDS files",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "loadContents": true,
                                        "glob": "*.txt",
                                        "outputEval": "${\n     function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    \n    \n    function pair_chromosome_gds(file_array){\n        var gdss = {};\n        for(var i=0; i<file_array.length; i++){\n            gdss[find_chromosome(file_array[i].path)] = file_array[i]\n        }\n        return gdss\n    }\n\n    var input_gdss = pair_chromosome_gds(inputs.input_gds_files)\n    var output_gdss = [];\n    var segments = self[0].contents.split('\\n');\n    var chr;\n    \n    segments = segments.slice(1)\n    for(var i=0;i<segments.length;i++){\n        chr = segments[i].split('\\t')[0]\n        if(chr in input_gdss){\n            output_gdss.push(input_gdss[chr])\n        }\n    }\n    return output_gdss\n}"
                                    },
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "id": "segments",
                                    "doc": "Segments.",
                                    "label": "Segments",
                                    "type": "int[]?",
                                    "outputBinding": {
                                        "loadContents": true,
                                        "glob": "*.txt",
                                        "outputEval": "${\n     function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    function pair_chromosome_gds(file_array){\n        var gdss = {};\n        for(var i=0; i<file_array.length; i++){\n            gdss[find_chromosome(file_array[i].path)] = file_array[i]\n        }\n        return gdss\n    }\n    \n    var input_gdss = pair_chromosome_gds(inputs.input_gds_files)\n    var output_segments = []\n    var segments = self[0].contents.split('\\n');\n    segments = segments.slice(1)\n    var chr;\n    \n    for(var i=0;i<segments.length;i++){\n        chr = segments[i].split('\\t')[0]\n        if(chr in input_gdss){\n            output_segments.push(i+1)\n        }\n    }\n    return output_segments\n    \n}"
                                    }
                                },
                                {
                                    "id": "aggregate_output",
                                    "doc": "Aggregate output.",
                                    "label": "Aggregate output",
                                    "type": [
                                        "null",
                                        {
                                            "type": "array",
                                            "items": [
                                                "null",
                                                "File"
                                            ]
                                        }
                                    ],
                                    "outputBinding": {
                                        "loadContents": true,
                                        "glob": "*.txt",
                                        "outputEval": "${\n     function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    function pair_chromosome_gds(file_array){\n        var gdss = {};\n        for(var i=0; i<file_array.length; i++){\n            gdss[find_chromosome(file_array[i].path)] = file_array[i]\n        }\n        return gdss\n    }\n    function pair_chromosome_gds_special(file_array, agg_file){\n        var gdss = {};\n        for(var i=0; i<file_array.length; i++){\n            gdss[find_chromosome(file_array[i].path)] = agg_file\n        }\n        return gdss\n    }\n    var input_gdss = pair_chromosome_gds(inputs.input_gds_files)\n    var segments = self[0].contents.split('\\n');\n    segments = segments.slice(1)\n    var chr;\n    \n    if(inputs.aggregate_files){\n        if (inputs.aggregate_files[0] != null){\n            if (inputs.aggregate_files[0].basename.includes('chr'))\n                var input_aggregate_files = pair_chromosome_gds(inputs.aggregate_files);\n            else\n                var input_aggregate_files = pair_chromosome_gds_special(inputs.input_gds_files, inputs.aggregate_files[0].path);\n            var output_aggregate_files = []\n            for(var i=0;i<segments.length;i++){\n                chr = segments[i].split('\\t')[0]\n                if(chr in input_aggregate_files){\n                    output_aggregate_files.push(input_aggregate_files[chr])\n                }\n                else if(chr in input_gdss){\n                    output_aggregate_files.push(null)\n                }\n            }\n            return output_aggregate_files\n        }\n    }\n    else{\n        var null_outputs = []\n        for(var i=0; i<segments.length; i++){\n            chr = segments[i].split('\\t')[0]\n            if(chr in input_gdss){\n                null_outputs.push(null)\n            }\n        }\n        return null_outputs\n    }\n}"
                                    }
                                },
                                {
                                    "id": "variant_include_output",
                                    "doc": "Variant Include Output",
                                    "label": "Variant Include Output",
                                    "type": [
                                        "null",
                                        {
                                            "type": "array",
                                            "items": [
                                                "null",
                                                "File"
                                            ]
                                        }
                                    ],
                                    "outputBinding": {
                                        "loadContents": true,
                                        "glob": "*.txt",
                                        "outputEval": "${\n     function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    function pair_chromosome_gds(file_array){\n        var gdss = {};\n        for(var i=0; i<file_array.length; i++){\n            gdss[find_chromosome(file_array[i].path)] = file_array[i]\n        }\n        return gdss\n    }\n    var input_gdss = pair_chromosome_gds(inputs.input_gds_files)\n    var segments = self[0].contents.split('\\n');\n    segments = segments.slice(1)\n    var chr;\n    \n    if(inputs.variant_include_files){\n        if (inputs.variant_include_files[0] != null){\n            var input_variant_files = pair_chromosome_gds(inputs.variant_include_files)\n            var output_variant_files = []\n            for(var i=0;i<segments.length;i++){\n                chr = segments[i].split('\\t')[0]\n                if(chr in input_variant_files){\n                    output_variant_files.push(input_variant_files[chr])\n                }\n                else if(chr in input_gdss){\n                    output_variant_files.push(null)\n                }\n            }\n            return output_variant_files\n        }\n    }\n    else{\n        var null_outputs = [];\n        for(var i=0; i<segments.length; i++){\n            chr = segments[i].split('\\t')[0]\n            if(chr in input_gdss){\n                null_outputs.push(null)\n            }\n        }\n        return null_outputs\n    }\n}"
                                    }
                                }
                            ],
                            "label": "SBG Prepare Segments",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    return \"cp \" + inputs.segments_file.path + \" .\"\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.8.1"
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105826,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/sbg-prepare-segments-1/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105826,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105826,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "af5431cfdc789d53445974b82b534a1ba1c6df2ac79d7b39af88dce65def8cb34"
                        },
                        "label": "SBG Prepare Segments",
                        "sbg:x": 76.38661193847656,
                        "sbg:y": -183.02523803710938
                    },
                    {
                        "id": "sbg_group_segments_1",
                        "in": [
                            {
                                "id": "assoc_files",
                                "source": [
                                    "sbg_flatten_lists/output_list"
                                ]
                            }
                        ],
                        "out": [
                            {
                                "id": "grouped_assoc_files"
                            },
                            {
                                "id": "chromosome"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.1",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/sbg-group-segments-1/0",
                            "baseCommand": [
                                "echo",
                                "\"Grouping\""
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Inputs",
                                    "id": "assoc_files",
                                    "type": "File[]",
                                    "label": "Assoc files",
                                    "doc": "Assoc files.",
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "grouped_assoc_files",
                                    "type": [
                                        "null",
                                        {
                                            "type": "array",
                                            "items": [
                                                {
                                                    "type": "array",
                                                    "items": [
                                                        "File",
                                                        "null"
                                                    ]
                                                },
                                                "null"
                                            ]
                                        }
                                    ],
                                    "outputBinding": {
                                        "outputEval": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"/\").pop();\n        chrom_num = chrom_num.substr(0,chrom_num.lastIndexOf(\".\")).split('_').slice(0,-1).join('_')\n        if(isNumeric(chrom_num.charAt(chrom_num.length-2)))\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 1))\n        }\n        return chr_array.toString()\n    }\n    \n    var assoc_files_dict = {};\n    var grouped_assoc_files = [];\n    var chr;\n    for(var i=0; i<inputs.assoc_files.length; i++){\n        chr = find_chromosome(inputs.assoc_files[i].path)\n        if(chr in assoc_files_dict){\n            assoc_files_dict[chr].push(inputs.assoc_files[i])\n        }\n        else{\n            assoc_files_dict[chr] = [inputs.assoc_files[i]]\n        }\n    }\n    for(var key in assoc_files_dict){\n        grouped_assoc_files.push(assoc_files_dict[key])\n    }\n    return grouped_assoc_files\n    \n}"
                                    }
                                },
                                {
                                    "id": "chromosome",
                                    "doc": "Chromosomes.",
                                    "label": "Chromosomes",
                                    "type": "string[]?",
                                    "outputBinding": {
                                        "outputEval": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"/\").pop();\n        chrom_num = chrom_num.substr(0,chrom_num.lastIndexOf(\".\")).split('_').slice(0,-1).join('_')\n        if(isNumeric(chrom_num.charAt(chrom_num.length-2)))\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 1))\n        }\n        return chr_array.toString()\n    }\n    \n    var assoc_files_dict = {};\n    var output_chromosomes = [];\n    var chr;\n    for(var i=0; i<inputs.assoc_files.length; i++){\n        chr = find_chromosome(inputs.assoc_files[i].path)\n        if(chr in assoc_files_dict){\n            assoc_files_dict[chr].push(inputs.assoc_files[i])\n        }\n        else{\n            assoc_files_dict[chr] = [inputs.assoc_files[i]]\n        }\n    }\n    for(var key in assoc_files_dict){\n        output_chromosomes.push(key)\n    }\n    return output_chromosomes\n    \n}"
                                    }
                                }
                            ],
                            "label": "SBG Group Segments",
                            "requirements": [
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.8.1"
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105827,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:appVersion": [
                                "v1.1"
                            ],
                            "sbg:id": "markoz/hgi/sbg-group-segments-1/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105827,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105827,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a515be0f5124c62e65c743e3ca9940a2d4d90f71217b08949ce69537195ad562c"
                        },
                        "label": "SBG Group Segments",
                        "sbg:x": 855.9915771484375,
                        "sbg:y": 119.47896575927734
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:maxNumberOfParallelInstances",
                        "value": "8"
                    },
                    {
                        "class": "sbg:AWSInstanceType",
                        "value": "c5.2xlarge;ebs-gp2;1024"
                    },
                    {
                        "class": "sbg:AzureInstanceType",
                        "value": "Standard_D8s_v4;PremiumSSD;1024"
                    }
                ],
                "requirements": [
                    {
                        "class": "ScatterFeatureRequirement"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105828,
                        "sbg:revisionNotes": "Workflow decomposed"
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648835486,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648835625,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1653923070,
                        "sbg:revisionNotes": "revert (qq plot)"
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1654733934,
                        "sbg:revisionNotes": "back to 0"
                    }
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/single-variant-association-testing/4.png",
                "sbg:toolAuthor": "TOPMed DCC",
                "sbg:license": "MIT",
                "sbg:categories": [
                    "GWAS",
                    "CWL1.0",
                    "Genomics"
                ],
                "sbg:links": [
                    {
                        "id": "https://github.com/UW-GAC/analysis_pipeline",
                        "label": "Source Code, Download"
                    },
                    {
                        "id": "https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz567/5536872?redirectedFrom=fulltext",
                        "label": "Publication"
                    },
                    {
                        "id": "https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test.html",
                        "label": "Home Page"
                    },
                    {
                        "id": "https://bioconductor.org/packages/devel/bioc/manuals/GENESIS/man/GENESIS.pdf",
                        "label": "Documentation"
                    }
                ],
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2",
                    "v1.1"
                ],
                "sbg:id": "markoz/hgi/single-variant-association-testing/4",
                "sbg:revision": 4,
                "sbg:revisionNotes": "back to 0",
                "sbg:modifiedOn": 1654733934,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105828,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 4,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a1633c5ae14c2277956607cd0bdaf5440537deca76393f02e05b064fae21d838c",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "GENESIS Single Variant Association Testing",
            "hints": [
                {
                    "class": "sbg:AWSInstanceType",
                    "value": "r4.8xlarge;ebs-gp2;1024"
                }
            ],
            "sbg:x": 1538.51025390625,
            "sbg:y": -461.8770751953125
        },
        {
            "id": "ld_pruning",
            "in": [
                {
                    "id": "gds_file",
                    "linkMerge": "merge_flattened",
                    "source": [
                        "vcf_to_gds_1/unique_variant_id_gds_per_chr"
                    ]
                },
                {
                    "id": "variant_include_file",
                    "source": "merge_rds_arrays/output"
                },
                {
                    "id": "genome_build",
                    "default": "hg19"
                },
                {
                    "id": "autosome_only",
                    "source": "autosome_only"
                }
            ],
            "out": [
                {
                    "id": "pruned_gds_output"
                },
                {
                    "id": "ld_pruning_output"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/ld-pruning/1",
                "doc": "This workflow LD prunes variants and creates a new GDS file containing only the pruned variants. Linkage disequilibrium (LD) is a measure of correlation of genotypes between a pair of variants. LD-pruning is the process filtering variants so that those that remain have LD measures below a given threshold. This procedure is typically used to identify a (nearly) independent subset of variants. This is often the first step in evaluating relatedness and population structure to avoid having results driven by clusters of variants in high LD regions of the genome.",
                "label": "LD Pruning",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "gds_file",
                        "sbg:fileTypes": "GDS",
                        "type": "File[]",
                        "label": "GDS files",
                        "doc": "Input GDS files, one per chromosome with string \"chr*\" in the file name.",
                        "sbg:x": -440.613037109375,
                        "sbg:y": 21.550386428833008
                    },
                    {
                        "id": "sample_include_file_pruning",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample Include file for LD pruning",
                        "doc": "RData file with vector of sample.id to use for LD pruning (unrelated samples are recommended). If not provided, all samples in the GDS files are included.",
                        "sbg:x": -482.9205627441406,
                        "sbg:y": -174.98597717285156
                    },
                    {
                        "id": "variant_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Variant Include file for LD pruning",
                        "doc": "RData file with vector of variant.id to consider for LD pruning. If not provided, all variants in the GDS files are included.",
                        "sbg:x": -426.8399658203125,
                        "sbg:y": -303.2926025390625
                    },
                    {
                        "id": "sample_include_file_gds",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample include file for output GDS",
                        "doc": "RData file with vector of sample.id to include in the output GDS. If not provided, all samples in the GDS files are included.",
                        "sbg:x": -258.4799499511719,
                        "sbg:y": 59.013526916503906
                    },
                    {
                        "id": "out_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Prefix for output files.",
                        "sbg:x": -535.4487915039062,
                        "sbg:y": -42.224388122558594
                    },
                    {
                        "id": "ld_r_threshold",
                        "type": "float?",
                        "label": "LD |r| threshold",
                        "doc": "|r| threshold for LD pruning.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "0.32 (r^2 = 0.1)"
                    },
                    {
                        "id": "ld_win_size",
                        "type": "float?",
                        "label": "LD window size",
                        "doc": "Sliding window size in Mb for LD pruning.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "10"
                    },
                    {
                        "id": "maf_threshold",
                        "type": "float?",
                        "label": "MAF threshold",
                        "doc": "Minimum MAF for variants used in LD pruning. Variants below this threshold are removed.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "0.01"
                    },
                    {
                        "id": "missing_threshold",
                        "type": "float?",
                        "label": "Missing call rate threshold",
                        "doc": "Maximum missing call rate for variants used in LD pruning. Variants above this threshold are removed.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "0.01"
                    },
                    {
                        "id": "exclude_pca_corr",
                        "type": "boolean?",
                        "label": "Exclude PCA corr",
                        "doc": "Exclude variants in genomic regions known to result in high PC-variant correlations when included (HLA, LCT, inversions).",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "TRUE"
                    },
                    {
                        "id": "genome_build",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "hg18",
                                    "hg19",
                                    "hg38"
                                ],
                                "name": "genome_build"
                            }
                        ],
                        "label": "Genome build",
                        "doc": "Genome build, used to define genomic regions to filter for PC-variant correlation.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "hg38"
                    },
                    {
                        "id": "autosome_only",
                        "type": "boolean?",
                        "label": "Autosomes only",
                        "doc": "Only include variants on the autosomes.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "TRUE"
                    }
                ],
                "outputs": [
                    {
                        "id": "pruned_gds_output",
                        "outputSource": [
                            "merge_gds/merged_gds_output"
                        ],
                        "sbg:fileTypes": "GDS",
                        "type": "File?",
                        "label": "Pruned GDS output file",
                        "doc": "GDS output file containing sample genotypes at pruned variants from all chromosomes.",
                        "sbg:x": 510.43572998046875,
                        "sbg:y": -266.7860412597656
                    },
                    {
                        "id": "ld_pruning_output",
                        "outputSource": [
                            "ld_pruning/ld_pruning_output"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Pruned output file",
                        "doc": "RData file with variant.id of pruned variants.",
                        "sbg:x": -61.60548782348633,
                        "sbg:y": -287.0057373046875
                    }
                ],
                "steps": [
                    {
                        "id": "ld_pruning",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "ld_r_threshold",
                                "source": "ld_r_threshold"
                            },
                            {
                                "id": "ld_win_size",
                                "source": "ld_win_size"
                            },
                            {
                                "id": "maf_threshold",
                                "source": "maf_threshold"
                            },
                            {
                                "id": "missing_threshold",
                                "source": "missing_threshold"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file_pruning"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "variant_include_file"
                            },
                            {
                                "id": "exclude_pca_corr",
                                "source": "exclude_pca_corr"
                            },
                            {
                                "id": "genome_build",
                                "source": "genome_build"
                            },
                            {
                                "id": "autosome_only",
                                "source": "autosome_only"
                            }
                        ],
                        "out": [
                            {
                                "id": "ld_pruning_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/ld-pruning/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input files",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS file",
                                    "doc": "Input GDS file.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "0.32 (r^2 = 0.1)",
                                    "id": "ld_r_threshold",
                                    "type": "float?",
                                    "label": "LD |r| threshold",
                                    "doc": "|r| threshold for LD pruning."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "10",
                                    "id": "ld_win_size",
                                    "type": "float?",
                                    "label": "LD window size",
                                    "doc": "Sliding window size in Mb for LD pruning."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "0.01",
                                    "id": "maf_threshold",
                                    "type": "float?",
                                    "label": "MAF threshold",
                                    "doc": "Minimum MAF for variants used in LD pruning. Variants below this threshold are removed."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "0.01",
                                    "id": "missing_threshold",
                                    "type": "float?",
                                    "label": "Missing call rate threshold",
                                    "doc": "Maximum missing call rate for variants used in LD pruning. Variants above this threshold are removed."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant Include file",
                                    "doc": "RData file with vector of variant.id to consider for LD pruning. If not provided, all variants in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "true",
                                    "id": "exclude_pca_corr",
                                    "type": "boolean?",
                                    "label": "Exclude PCA corr",
                                    "doc": "Exclude variants in genomic regions known to result in high PC-variant correlations when included (HLA, LCT, inversions).",
                                    "default": true
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "hg38",
                                    "id": "genome_build",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "hg18",
                                                "hg19",
                                                "hg38"
                                            ],
                                            "name": "genome_build"
                                        }
                                    ],
                                    "label": "Genome build",
                                    "doc": "Genome build, used to define genomic regions to filter for PC-variant correlation.",
                                    "default": "hg38"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "chromosome",
                                    "type": "string?",
                                    "inputBinding": {
                                        "prefix": "--chromosome",
                                        "shellQuote": false,
                                        "position": 20
                                    },
                                    "label": "Chromosome",
                                    "doc": "Chromosome range of gds file (1-24 or X,Y)."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "true",
                                    "id": "autosome_only",
                                    "type": "boolean?",
                                    "label": "Autosomes only",
                                    "doc": "Only include variants on the autosomes.",
                                    "default": true
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "ld_pruning_output",
                                    "doc": "RData file with variant.id of pruned variants.",
                                    "label": "Pruned output file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.RData"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "ld_pruning",
                            "arguments": [
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "ld_pruning.config"
                                },
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/ld_pruning.R"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "ld_pruning.config",
                                            "entry": "${\n\n    var cmd_line = \"\"\n    \n    if(inputs.gds_file)\n        cmd_line += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n    if(!inputs.exclude_pca_corr)\n        cmd_line += \"exclude_pca_corr FALSE\" + \"\\n\"\n    if(!inputs.autosome_only)\n        cmd_line += \"autosome_only FALSE\" + \"\\n\"\n    if(inputs.genome_build)\n        cmd_line += \"genome_build \\\"\" + inputs.genome_build + \"\\\"\\n\"\n    if(inputs.ld_r_threshold)\n        cmd_line += \"ld_r_threshold \" + inputs.ld_r_threshold + \"\\n\"\n    if(inputs.ld_win_size)\n        cmd_line += \"ld_win_size \" + inputs.ld_win_size + \"\\n\"\n    if(inputs.maf_threshold)\n        cmd_line += \"maf_threshold \" + inputs.maf_threshold + \"\\n\"\n    if(inputs.missing_threshold)\n        cmd_line += \"missing_threshold \" + inputs.missing_threshold + \"\\n\"\n    \n    if(inputs.sample_include_file)\n        cmd_line += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\"\n        \n        \n        \n    if(inputs.variant_include_file){\n        cmd_line += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n        //var chr = inputs.variant_include_file.path.split('/').pop()\n        //cmd_line += 'out_file \"' + outfile_temp + '_' + chr + '\"\\n'\n    }\n    if (inputs.gds_file.nameroot.includes('chr'))\n    {\n        var parts = inputs.gds_file.nameroot.split('chr')\n        var outfile_temp = 'pruned_variants_chr' + parts[1] + '.RData'\n    } else {\n        var outfile_temp = 'pruned_variants.RData'\n    }\n    if(inputs.out_prefix){\n        outfile_temp = inputs.out_prefix + '_' + outfile_temp\n    }\n    \n    cmd_line += 'out_file \"' + outfile_temp + '\"\\n'\n        \n    return cmd_line\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "ld_pruning.config"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105830,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/ld-pruning/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105830,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105830,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a6c9034d53f3d8edea0d82b0594ec0b34ce1f64878600cb2a29f73e4e1ee58aa4"
                        },
                        "label": "ld_pruning",
                        "scatter": [
                            "gds_file"
                        ],
                        "sbg:x": -259.3580627441406,
                        "sbg:y": -159.14634704589844
                    },
                    {
                        "id": "subset_gds",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file_gds"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "ld_pruning/ld_pruning_output"
                            }
                        ],
                        "out": [
                            {
                                "id": "output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/subset-gds/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Inputs",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS file",
                                    "doc": "Input GDS file.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample include file",
                                    "doc": "RData file with vector of sample.id to include. All samples in the GDS file are included when not provided.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "variant_include_file",
                                    "type": "File",
                                    "label": "Variant include file",
                                    "doc": "RData file with vector of variant.id to include. All variants in the GDS files are used when not provided.",
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "output",
                                    "doc": "GDS file with subset of variants from original file",
                                    "label": "Subset file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.gds"
                                    },
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "label": "subset_gds",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/subset_gds.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "subset_gds.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "subset_gds.config",
                                            "entry": "${\n    var config = \"\";\n  \tif (inputs.gds_file) \n  \t{\n  \t    config += 'gds_file \"' + inputs.gds_file.path  + \"\\\"\\n\";\n  \t}\n  \tif (inputs.variant_include_file)\n  \t{\n  \t    config += 'variant_include_file \"' + inputs.variant_include_file.path + \"\\\"\\n\"\n  \t}\n    \n    if (inputs.sample_include_file) \n      config += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\";\n      \n    if (inputs.out_prefix)\n    {\n        var chromosome = inputs.gds_file.nameroot.split('chr')[1].split('.')[0];\n        config += 'subset_gds_file \"' + inputs.out_prefix + '_chr' + chromosome + '.gds\"\\n';\n    }\n    else\n    {\n        var chromosome = inputs.gds_file.nameroot.split('chr')[1].split('.')[0]\n        var basename = inputs.gds_file.nameroot.split('chr')[0]\n        config += 'subset_gds_file \"' + basename + 'subset_chr' + chromosome + '.gds\"\\n'\n    }\n\n    return config\n}\n",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "subset_gds.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105831,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/subset-gds/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105831,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105831,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "aa863df6a1ecf0a4b8b39a397e6383efbfdd0e3ca86652cf09e115e001e1569bb"
                        },
                        "label": "subset_gds",
                        "scatter": [
                            "gds_file",
                            "variant_include_file"
                        ],
                        "scatterMethod": "dotproduct",
                        "sbg:x": -71.82777404785156,
                        "sbg:y": -48.61606216430664
                    },
                    {
                        "id": "merge_gds",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": [
                                    "subset_gds/output"
                                ]
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            }
                        ],
                        "out": [
                            {
                                "id": "merged_gds_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/merge-gds/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "id": "gds_file",
                                    "type": "File[]",
                                    "label": "GDS files",
                                    "doc": "Input GDS files.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "merged_gds_output",
                                    "doc": "GDS output file",
                                    "label": "Merged GDS output file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.gds"
                                    },
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "label": "merge_gds",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 10,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/merge_gds.R"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    var cmd_line = \"\"\n    for (var i=0; i<inputs.gds_file.length; i++)\n        cmd_line += \"ln -s \" + inputs.gds_file[i].path + \" \" + inputs.gds_file[i].basename + \" && \"\n    \n    return cmd_line\n}"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "merge_gds.config"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "R -q --vanilla"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "merge_gds.config",
                                            "entry": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    var a_file = inputs.gds_file[0]\n    var chr = find_chromosome(a_file.basename);\n    var path = a_file.path.split('chr'+chr);\n\n    var chr_array = [];\n        var chrom_num;\n        for (var i = 0; i < inputs.gds_file.length; i++) \n        {\n            chrom_num = find_chromosome(inputs.gds_file[i].nameroot)\n            \n            chr_array.push(chrom_num)\n        }\n        \n        chr_array = chr_array.sort((a, b) => a.localeCompare(b, 'en', {numeric: true, ignorePunctuation: true}))\n        \n        var chrs = chr_array.join(' ')\n    \n    if(inputs.out_prefix)\n    {\n        var merged_gds_file_name = inputs.out_prefix + \".gds\"\n    }\n    else\n    {\n        var merged_gds_file_name = \"merged.gds\"\n    }\n    var arguments = []\n    arguments.push('gds_file ' + '\"' + path[0].split('/').pop() + 'chr ' + path[1] + '\"')\n    arguments.push('merged_gds_file \"' + merged_gds_file_name + '\"')\n    arguments.push('chromosomes \"' + chrs + '\"')\n    \n    return arguments.join('\\n')\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "merge_gds.config"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105832,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/merge-gds/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105832,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105832,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "acee209113ca200e23c2e40bdc8ec609216f01279f609080a27d425e4c69ca48a"
                        },
                        "label": "merge_gds",
                        "sbg:x": 170.64193725585938,
                        "sbg:y": -207.07080078125
                    },
                    {
                        "id": "check_merged_gds",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "subset_gds/output"
                            },
                            {
                                "id": "merged_gds_file",
                                "source": "merge_gds/merged_gds_output"
                            }
                        ],
                        "out": [],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/check-merged-gds/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Inputs",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS File",
                                    "doc": "Base reference file for comparison.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "merged_gds_file",
                                    "type": "File",
                                    "label": "Merged GDS file",
                                    "doc": "Output of merge_gds script. This file is being checked against starting gds file.",
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "outputs": [],
                            "label": "check_merged_gds",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/check_merged_gds.R"
                                },
                                {
                                    "prefix": "--chromosome",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "${\n    return inputs.gds_file.path.split('chr')[1].split('.')[0]\n}"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "check_merged_gds.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "check_merged_gds.config",
                                            "entry": "${\n    var config = \"\";\n    \n    if(inputs.gds_file)\n        var gds_name = \"\";\n        var gds_file = [].concat(inputs.gds_file)\n        var gds = gds_file[0].path;\n        var gds_first_part = gds.split('chr')[0];\n        var gds_second_part = gds.split('chr')[1].split('.');\n        gds_second_part.shift();\n        gds_name = gds_first_part + 'chr .' + gds_second_part.join('.');\n        config += \"gds_file \\\"\" + gds_name + \"\\\"\\n\"\n\n    if (inputs.merged_gds_file)\n        config += \"merged_gds_file \\\"\" + inputs.merged_gds_file.path + \"\\\"\\n\"\n\n    return config\n}\n",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "check_merged_gds.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105833,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/check-merged-gds/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105833,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105833,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a36acbf3d55cee804b31eaaa59f8049a1127973b7f0f90212b180041e1f4c44a6"
                        },
                        "label": "check_merged_gds",
                        "scatter": [
                            "gds_file"
                        ],
                        "sbg:x": 419.7127380371094,
                        "sbg:y": -23.686861038208008
                    }
                ],
                "requirements": [
                    {
                        "class": "ScatterFeatureRequirement"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "Ancestry and Relatedness"
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/ld-pruning/1.png",
                "sbg:original_source": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/smgogarten/genesis-relatedness/ld-pruning-pipeline/32/raw/",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105830,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105834,
                        "sbg:revisionNotes": "Workflow decomposed"
                    }
                ],
                "sbg:toolkit": "UW-GAC Ancestry and Relatedness",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/ld-pruning/1",
                "sbg:revision": 1,
                "sbg:revisionNotes": "Workflow decomposed",
                "sbg:modifiedOn": 1637105834,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105830,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 1,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "ad7a840aeadc00c2f0e0ab044621726a7acd6f57dd8e13e978a76e494cf920924"
            },
            "label": "LD Pruning",
            "hints": [
                {
                    "class": "sbg:AWSInstanceType",
                    "value": "m5.xlarge;ebs-gp2;1024"
                }
            ],
            "sbg:x": 59.864078521728516,
            "sbg:y": 40.21495819091797
        },
        {
            "id": "king_robust",
            "in": [
                {
                    "id": "gds_file",
                    "source": "ld_pruning/pruned_gds_output"
                },
                {
                    "id": "variant_include_file",
                    "source": "merge_rds_arrays/output"
                },
                {
                    "id": "out_prefix",
                    "source": "out_prefix"
                },
                {
                    "id": "phenotype_file",
                    "source": "phenotype_file"
                }
            ],
            "out": [
                {
                    "id": "king_robust_output"
                },
                {
                    "id": "king_robust_plots"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/king-robust/1",
                "doc": "This workflow uses the KING-robust method to estimate kinship coefficients, and returns results for all pairs of samples. Due to the negative bias of these kinship estimates for samples of different ancestry, they can be used as a measure of ancestry divergence in PC-AiR.",
                "label": "KING robust",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "gds_file",
                        "sbg:fileTypes": "GDS",
                        "type": "File",
                        "label": "GDS file",
                        "doc": "Input GDS file. It is recommended to use an LD pruned file with all chromosomes.",
                        "sbg:x": -289,
                        "sbg:y": 70
                    },
                    {
                        "id": "sample_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample Include file",
                        "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                        "sbg:x": -353.2032470703125,
                        "sbg:y": -120
                    },
                    {
                        "id": "variant_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Variant Include file",
                        "doc": "RData file with vector of variant.id to use for kinship estimation. If not provided, all variants in the GDS file are included.",
                        "sbg:x": -294.2032470703125,
                        "sbg:y": -242.32879638671875
                    },
                    {
                        "id": "out_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Prefix for output files.",
                        "sbg:x": -397.2032470703125,
                        "sbg:y": 0.6711986660957336
                    },
                    {
                        "id": "phenotype_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Phenotype File",
                        "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by group.",
                        "sbg:x": -117,
                        "sbg:y": 117
                    },
                    {
                        "id": "kinship_plot_threshold",
                        "type": "float?",
                        "label": "Kinship plotting threshold",
                        "doc": "Minimum kinship for a pair to be included in the plot.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)"
                    },
                    {
                        "id": "group",
                        "type": "string?",
                        "label": "Group column name",
                        "doc": "Name of column in phenotype_file containing group variable (e.g., study) for plotting.",
                        "sbg:exposed": true
                    }
                ],
                "outputs": [
                    {
                        "id": "king_robust_output",
                        "outputSource": [
                            "king_robust/king_robust_output"
                        ],
                        "sbg:fileTypes": "GDS",
                        "type": "File?",
                        "label": "KING robust output",
                        "doc": "GDS file with matrix of pairwise kinship estimates.",
                        "sbg:x": 98,
                        "sbg:y": -235
                    },
                    {
                        "id": "king_robust_plots",
                        "outputSource": [
                            "kinship_plots/kinship_plots"
                        ],
                        "sbg:fileTypes": "PDF",
                        "type": "File[]?",
                        "label": "Kinship plots",
                        "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
                        "sbg:x": 312,
                        "sbg:y": -108
                    }
                ],
                "steps": [
                    {
                        "id": "king_robust",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "variant_include_file"
                            }
                        ],
                        "out": [
                            {
                                "id": "king_robust_output"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/king-robust/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input files",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS file",
                                    "doc": "Input GDS file.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant Include file",
                                    "doc": "RData file with vector of variant.id to use for kinship estimation. If not provided, all variants in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "king_robust_output",
                                    "doc": "GDS file with matrix of pairwise kinship estimates.",
                                    "label": "KING robust output",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.gds"
                                    },
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "label": "king_robust",
                            "arguments": [
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "ibd_king.config"
                                },
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/ibd_king.R"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "coresMin": 8
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "ibd_king.config",
                                            "entry": "${\n\n    var cmd_line = \"\"\n    \n    if(inputs.gds_file)\n        cmd_line += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n        \n    if(inputs.sample_include_file)\n        cmd_line += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\"\n        \n    if(inputs.variant_include_file){\n        cmd_line += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n    }\n    \n    if(inputs.out_prefix){\n        cmd_line += 'out_file \"' + inputs.out_prefix + '_king_robust.gds\"\\n'\n    } else {\n        cmd_line += 'out_file \\\"king_robust.gds\\\"\\n'\n    }\n    \n        \n    return cmd_line\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "EnvVarRequirement",
                                    "envDef": [
                                        {
                                            "envName": "NSLOTS",
                                            "envValue": "${ return runtime.cores }"
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "ibd_king.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105836,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/king-robust/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105836,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105836,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a8bdc1361577029af7e2a464064115e6ba24f3fedd1b9d850c1033415a1237ac2"
                        },
                        "label": "king_robust",
                        "sbg:x": -168,
                        "sbg:y": -76
                    },
                    {
                        "id": "kinship_plots",
                        "in": [
                            {
                                "id": "kinship_file",
                                "source": "king_robust/king_robust_output"
                            },
                            {
                                "id": "kinship_method",
                                "default": "king_robust"
                            },
                            {
                                "id": "kinship_plot_threshold",
                                "source": "kinship_plot_threshold"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "group",
                                "source": "group"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix",
                                "valueFrom": "${ return inputs.out_prefix + \"_king_robust\" }"
                            }
                        ],
                        "out": [
                            {
                                "id": "kinship_plots"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/kinship-plots/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input",
                                    "id": "kinship_file",
                                    "type": "File",
                                    "label": "Kinship File",
                                    "doc": "Kinship file",
                                    "sbg:fileTypes": "RDATA, SEG, KIN, GDS"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "kinship_method",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "king_ibdseg",
                                            "pcrelate",
                                            "king_robust"
                                        ],
                                        "name": "kinship_method"
                                    },
                                    "label": "Kinship method",
                                    "doc": "Method used to generate kinship estimates."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)",
                                    "id": "kinship_plot_threshold",
                                    "type": "float?",
                                    "label": "Kinship plotting threshold",
                                    "doc": "Minimum kinship for a pair to be included in the plot."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "phenotype_file",
                                    "type": "File?",
                                    "label": "Phenotype File",
                                    "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by group.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "NA",
                                    "id": "group",
                                    "type": "string?",
                                    "label": "Group column name",
                                    "doc": "Name of column in phenotype_file containing group variable (e.g., study) for plotting."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include File",
                                    "doc": "RData file with vector of sample.id to include.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "kinship",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "kinship_plots",
                                    "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
                                    "label": "Kinship plots",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.pdf"
                                    },
                                    "sbg:fileTypes": "PDF"
                                }
                            ],
                            "label": "kinship_plots",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/kinship_plots.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "kinship_plots.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "kinship_plots.config",
                                            "entry": "${\nvar cmd_line = \"\"\n\nif(inputs.kinship_file)\n    cmd_line += \"kinship_file \\\"\" + inputs.kinship_file.path + \"\\\"\\n\"\n\nif(inputs.kinship_method){\n    if(inputs.kinship_method == \"king_robust\"){\n        cmd_line +='kinship_method \"king\"\\n'\n    } else {\n        cmd_line +='kinship_method \"' + inputs.kinship_method + '\"\\n'\n    }\n}\n\nif(inputs.kinship_plot_threshold)\n    cmd_line +='kinship_threshold \"' + inputs.kinship_plot_threshold + '\\n'\n\nif(inputs.out_prefix) {\n    cmd_line += 'out_file_all \"' + inputs.out_prefix + '_all.pdf\"\\n'\n    cmd_line += 'out_file_cross \"' + inputs.out_prefix + '_cross_group.pdf\"\\n'\n    cmd_line += 'out_file_study \"' + inputs.out_prefix + '_within_group.pdf\"\\n'\n}\n\nif(inputs.phenotype_file)\n    cmd_line += 'phenotype_file \"' + inputs.phenotype_file.path + '\"\\n'\n\nif(inputs.group)\n    cmd_line += 'study \"' + inputs.group + '\"\\n'\n\nif(inputs.sample_include_file)\n    cmd_line += 'sample_include_file \"' + inputs.sample_include_file.path + '\"\\n'\n\nreturn cmd_line\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "kinship_plots.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105837,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/kinship-plots/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105837,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105837,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "ac07d7f46f1402f82f9271881eca181d1a4e44613e62de62f9242aff85413f411"
                        },
                        "label": "kinship_plots",
                        "sbg:x": 30,
                        "sbg:y": -9
                    }
                ],
                "requirements": [
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "Ancestry and Relatedness"
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/king-robust/1.png",
                "sbg:original_source": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/smgogarten/genesis-relatedness/king-robust-1/15/raw/",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105836,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105838,
                        "sbg:revisionNotes": "Workflow decomposed"
                    }
                ],
                "sbg:toolkit": "UW-GAC Ancestry and Relatedness",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/king-robust/1",
                "sbg:revision": 1,
                "sbg:revisionNotes": "Workflow decomposed",
                "sbg:modifiedOn": 1637105838,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105836,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 1,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a2841c7c47a9dfaac5c2ce2e6cb6f1fcb6ebb6936bacb2e72c2ade728e1536b3d"
            },
            "label": "KING robust",
            "sbg:x": 261,
            "sbg:y": -408.25
        },
        {
            "id": "pc_air",
            "in": [
                {
                    "id": "kinship_file",
                    "source": "king_robust/king_robust_output"
                },
                {
                    "id": "gds_file",
                    "source": "ld_pruning/pruned_gds_output"
                },
                {
                    "id": "variant_include_file",
                    "source": "merge_rds_arrays/output"
                },
                {
                    "id": "phenotype_file",
                    "source": "phenotype_file"
                },
                {
                    "id": "gds_file_full",
                    "linkMerge": "merge_flattened",
                    "source": [
                        "vcf_to_gds_1/unique_variant_id_gds_per_chr"
                    ]
                },
                {
                    "id": "run_correlation",
                    "default": true
                },
                {
                    "id": "pruned_variant_file",
                    "source": [
                        "ld_pruning/ld_pruning_output"
                    ]
                }
            ],
            "out": [
                {
                    "id": "out_unrelated_file"
                },
                {
                    "id": "out_related_file"
                },
                {
                    "id": "pcair_output"
                },
                {
                    "id": "pcair_plots"
                },
                {
                    "id": "pc_correlation_plots"
                },
                {
                    "id": "pca_corr_gds"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/pc-air/1",
                "doc": "This workflow uses the PC-AiR algorithm to compute ancestry principal components (PCs) while accounting for kinship.\n\nStep 1 uses pairwise kinship estimates to assign samples to an unrelated set that is representative of all ancestries in the sample. Step 2 performs Principal Component Analysis (PCA) on the unrelated set, then projects relatives onto the resulting set of PCs. Step 3 plots the PCs, optionally color-coding by a grouping variable. Step 4 (optional) calculates the correlation between each PC and variants in the dataset, then plots this correlation to allow screening for PCs that are driven by particular genomic regions.",
                "label": "PC-AiR",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "kinship_file",
                        "sbg:fileTypes": "RDATA, GDS",
                        "type": "File",
                        "label": "Kinship File",
                        "doc": "Pairwise kinship matrix used to identify unrelated and related sets of samples in Step 1. It is recommended to use KING-IBDseg or PC-Relate estimates.",
                        "sbg:x": -566.5120849609375,
                        "sbg:y": 182.8510284423828
                    },
                    {
                        "id": "out_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Prefix for output files.",
                        "sbg:x": -665.81005859375,
                        "sbg:y": 56.52513885498047
                    },
                    {
                        "id": "gds_file",
                        "sbg:fileTypes": "GDS",
                        "type": "File",
                        "label": "Pruned GDS File",
                        "doc": "Input GDS file for PCA. It is recommended to use an LD pruned file with all chromosomes.",
                        "sbg:x": -381.15753173828125,
                        "sbg:y": 127.77397155761719
                    },
                    {
                        "id": "divergence_file",
                        "sbg:fileTypes": "RDATA, GDS",
                        "type": "File?",
                        "label": "Divergence File",
                        "doc": "Pairwise matrix used to identify ancestrally divergent pairs of samples in Step 1. It is recommended to use KING-robust estimates.",
                        "sbg:x": -590.9608764648438,
                        "sbg:y": 319.02606201171875
                    },
                    {
                        "id": "sample_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample Include file",
                        "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                        "sbg:x": -628.7839965820312,
                        "sbg:y": -77.0744857788086
                    },
                    {
                        "id": "variant_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Variant include file",
                        "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included.",
                        "sbg:x": -537.6592407226562,
                        "sbg:y": -221.6741180419922
                    },
                    {
                        "id": "phenotype_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Phenotype file",
                        "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for color-coding PCA plots by group.",
                        "sbg:x": -14.366241455078125,
                        "sbg:y": -274.2578430175781
                    },
                    {
                        "id": "kinship_threshold",
                        "type": "float?",
                        "label": "Kinship threshold",
                        "doc": "Minimum kinship estimate to use for identifying relatives.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)"
                    },
                    {
                        "id": "divergence_threshold",
                        "type": "float?",
                        "label": "Divergence threshold",
                        "doc": "Maximum divergence estimate to use for identifying ancestrally divergent pairs of samples.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "-2^(-9/2)"
                    },
                    {
                        "id": "n_pairs",
                        "type": "int?",
                        "label": "Number of PCs",
                        "doc": "Number of PCs to include in the pairs plot.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "6"
                    },
                    {
                        "id": "group",
                        "type": "string?",
                        "label": "Group",
                        "doc": "Name of column in phenotype_file containing group variable for color-coding plots.",
                        "sbg:exposed": true
                    },
                    {
                        "id": "n_pcs",
                        "type": "int?",
                        "label": "Number of PCs",
                        "doc": "Number of PCs (Principal Components) to return.",
                        "sbg:toolDefaultValue": "32",
                        "sbg:x": -427.1446533203125,
                        "sbg:y": -333.2216796875
                    },
                    {
                        "id": "gds_file_full",
                        "sbg:fileTypes": "GDS",
                        "type": "File[]",
                        "label": "Full GDS Files",
                        "doc": "GDS files (one per chromosome) used to calculate PC-variant correlations.",
                        "sbg:x": -283.20037841796875,
                        "sbg:y": 376.7168273925781
                    },
                    {
                        "id": "run_correlation",
                        "type": "boolean",
                        "label": "Run PC-variant correlation",
                        "doc": "For pruned variants as well as a random sample of additional variants, compute correlation between the variants and PCs, and generate plots. This step can be computationally intensive, but is useful for verifying that PCs are not driven by small regions of the genome.",
                        "sbg:x": -405.8658447265625,
                        "sbg:y": 275.1279296875
                    },
                    {
                        "id": "pruned_variant_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File[]",
                        "label": "Pruned variant files",
                        "doc": "RData files (one per chromosome) with vector of variant.id produced by the LD pruning workflow. These variants will be added to the set of randomly selected variants.",
                        "sbg:x": -241.30958557128906,
                        "sbg:y": 195.14608764648438
                    },
                    {
                        "id": "n_corr_vars",
                        "type": "int?",
                        "label": "Number of variants to select",
                        "doc": "Randomly select this number of variants distributed across the entire genome to use for PC-variant correlation. If running on a single chromosome, the variants returned will be scaled by the proportion of that chromosome in the genome.",
                        "sbg:exposed": true
                    },
                    {
                        "id": "n_pcs_plot",
                        "type": "int?",
                        "label": "Number of PCs to plot",
                        "doc": "Number of PCs to plot.",
                        "sbg:exposed": true
                    },
                    {
                        "id": "n_perpage",
                        "type": "int?",
                        "label": "Number of plots per page",
                        "doc": "Number of PC-variant correlation plots to stack in a single page. The number of png files generated will be ceiling(n_pcs_plot/n_perpage).",
                        "sbg:exposed": true
                    }
                ],
                "outputs": [
                    {
                        "id": "out_unrelated_file",
                        "outputSource": [
                            "find_unrelated/out_unrelated_file"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Unrelated file",
                        "doc": "RData file with vector of sample.id of unrelated samples identified in Step 1",
                        "sbg:x": -223.16201782226562,
                        "sbg:y": -296.2215881347656
                    },
                    {
                        "id": "out_related_file",
                        "outputSource": [
                            "find_unrelated/out_related_file"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Related file",
                        "doc": "RData file with vector of sample.id of samples related to the set of unrelated samples identified in Step 1",
                        "sbg:x": -200.21041870117188,
                        "sbg:y": -189.47299194335938
                    },
                    {
                        "id": "pcair_output",
                        "outputSource": [
                            "pca_byrel/pcair_output"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "RData file with PC-AiR PCs for all samples",
                        "sbg:x": 372.57098388671875,
                        "sbg:y": -50.2413215637207
                    },
                    {
                        "id": "pcair_plots",
                        "outputSource": [
                            "pca_plots/pca_plots"
                        ],
                        "type": "File[]?",
                        "label": "PC plots",
                        "doc": "PC plots",
                        "sbg:x": 471.4736022949219,
                        "sbg:y": -205.1788787841797
                    },
                    {
                        "id": "pc_correlation_plots",
                        "outputSource": [
                            "pc_variant_correlation/pc_correlation_plots"
                        ],
                        "sbg:fileTypes": "PNG",
                        "type": "File[]?",
                        "label": "PC-variant correlation plots",
                        "doc": "PC-variant correlation plots",
                        "sbg:x": 332.8373107910156,
                        "sbg:y": 208.7924346923828
                    },
                    {
                        "id": "pca_corr_gds",
                        "outputSource": [
                            "pc_variant_correlation/pca_corr_gds"
                        ],
                        "sbg:fileTypes": "GDS",
                        "type": "File[]?",
                        "label": "PC-variant correlation",
                        "doc": "GDS file with PC-variant correlation results",
                        "sbg:x": 285.73321533203125,
                        "sbg:y": 85.907958984375
                    }
                ],
                "steps": [
                    {
                        "id": "find_unrelated",
                        "in": [
                            {
                                "id": "kinship_file",
                                "source": "kinship_file"
                            },
                            {
                                "id": "divergence_file",
                                "source": "divergence_file"
                            },
                            {
                                "id": "kinship_threshold",
                                "source": "kinship_threshold"
                            },
                            {
                                "id": "divergence_threshold",
                                "source": "divergence_threshold"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            }
                        ],
                        "out": [
                            {
                                "id": "out_related_file"
                            },
                            {
                                "id": "out_unrelated_file"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "smgogarten/genesis-relatedness/find-unrelated/21",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "kinship_file",
                                    "type": "File",
                                    "label": "Kinship File",
                                    "doc": "Pairwise kinship matrix used to identify unrelated and related sets of samples. It is recommended to use KING-IBDseg or PC-Relate estimates.",
                                    "sbg:fileTypes": "RDATA, GDS"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "divergence_file",
                                    "type": "File?",
                                    "label": "Divergence File",
                                    "doc": "Pairwise matrix used to identify ancestrally divergent pairs of samples. It is recommended to use KING-robust estimates.",
                                    "sbg:fileTypes": "RDATA, GDS"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)",
                                    "id": "kinship_threshold",
                                    "type": "float?",
                                    "label": "Kinship threshold",
                                    "doc": "Minimum kinship estimate to use for identifying relatives."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "-2^(-9/2)",
                                    "id": "divergence_threshold",
                                    "type": "float?",
                                    "label": "Divergence threshold",
                                    "doc": "Maximum divergence estimate to use for identifying ancestrally divergent pairs of samples."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the kinship file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "out_related_file",
                                    "doc": "RData file with vector of sample.id of samples related to the set of unrelated samples",
                                    "label": "Related file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${   \n    if (!inputs.out_prefix) {     \n        var comm = \"related.RData\"   \n    } else {     \t\n        var comm = inputs.out_prefix + \"_related.RData\"   \n    }   \n    return comm \n}"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "out_unrelated_file",
                                    "doc": "RData file with vector of sample.id of unrelated samples",
                                    "label": "Unrelated file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${       \n    if (!inputs.out_prefix) {              \n        var comm = \"unrelated.RData\"        \n    } else {     \t         \n        var comm = inputs.out_prefix + \"_unrelated.RData\"       \n    }        \n    return comm  \n}"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "find_unrelated",
                            "arguments": [
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "find_unrelated.config"
                                },
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/find_unrelated.R"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "find_unrelated.config",
                                            "entry": "${\n    var config = \"\"\n    \n    if(inputs.kinship_file)\n        config += \"kinship_file \\\"\" + inputs.kinship_file.path + \"\\\"\\n\"\n\n    if(inputs.divergence_file)\n        config += \"divergence_file \\\"\" + inputs.divergence_file.path + \"\\\"\\n\"\n\n  \tif (inputs.kinship_threshold) \n      config += \"kinship_threshold \" + inputs.kinship_threshold  + \"\\n\"\n    \n  \tif (inputs.divergence_threshold) \n      config += \"divergence_threshold \" + inputs.kinship_threshold  + \"\\n\"\n    \n  \tif (inputs.out_prefix) {\n      config += \"out_related_file \\\"\" + inputs.out_prefix  + \"_related.RData\\\"\\n\"\n    \n      config += \"out_unrelated_file \\\"\" + inputs.out_prefix  + \"_unrelated.RData\\\"\\n\"\n  \t}\n  \n    if (inputs.sample_include_file) \n      config += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\"\n\n    return config\n}\n",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "find_unrelated.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:content_hash": "a299286715c93f8dc19ab52236662f813bcf12f4593a415132c3e4f070eb6c824",
                            "sbg:contributors": [
                                "smgogarten"
                            ],
                            "sbg:createdBy": "smgogarten",
                            "sbg:createdOn": 1581464601,
                            "sbg:id": "smgogarten/genesis-relatedness/find-unrelated/21",
                            "sbg:image_url": null,
                            "sbg:latestRevision": 21,
                            "sbg:modifiedBy": "smgogarten",
                            "sbg:modifiedOn": 1623274853,
                            "sbg:project": "smgogarten/genesis-relatedness",
                            "sbg:projectName": "GENESIS relatedness",
                            "sbg:publisher": "sbg",
                            "sbg:revision": 21,
                            "sbg:revisionNotes": "update docker image",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1581464601,
                                    "sbg:revision": 0,
                                    "sbg:revisionNotes": "Copy of boris_majic/topmed-optimization/find-unrelated/2"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1583954906,
                                    "sbg:revision": 1,
                                    "sbg:revisionNotes": "import changes from RC"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604987303,
                                    "sbg:revision": 2,
                                    "sbg:revisionNotes": "use out_prefix instead of separate file names"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604987972,
                                    "sbg:revision": 3,
                                    "sbg:revisionNotes": "fix config file"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604988821,
                                    "sbg:revision": 4,
                                    "sbg:revisionNotes": "use master docker image"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604990537,
                                    "sbg:revision": 5,
                                    "sbg:revisionNotes": "must keep separate output ports for use in workflow"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604990841,
                                    "sbg:revision": 6,
                                    "sbg:revisionNotes": "allow for missing output prefix"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604990983,
                                    "sbg:revision": 7,
                                    "sbg:revisionNotes": "using prefix in output ports"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604991264,
                                    "sbg:revision": 8,
                                    "sbg:revisionNotes": "update category"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604991568,
                                    "sbg:revision": 9,
                                    "sbg:revisionNotes": "update description names"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605548620,
                                    "sbg:revision": 10,
                                    "sbg:revisionNotes": "try stdout redirect"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605570523,
                                    "sbg:revision": 11,
                                    "sbg:revisionNotes": "separate arguments, move config to logs"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605571246,
                                    "sbg:revision": 12,
                                    "sbg:revisionNotes": "try redirecting stdout and stderr to same file"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605571681,
                                    "sbg:revision": 13,
                                    "sbg:revisionNotes": "redirecting stdout and stderr to same file just loses stderr"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605572061,
                                    "sbg:revision": 14,
                                    "sbg:revisionNotes": "see what happens when we pipe to R instead of using Rscript"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605572744,
                                    "sbg:revision": 15,
                                    "sbg:revisionNotes": "rearrange base command and args"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1606355852,
                                    "sbg:revision": 16,
                                    "sbg:revisionNotes": "update input and output labels"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1608621762,
                                    "sbg:revision": 17,
                                    "sbg:revisionNotes": "default label"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609370869,
                                    "sbg:revision": 18,
                                    "sbg:revisionNotes": "update descriptions"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609447255,
                                    "sbg:revision": 19,
                                    "sbg:revisionNotes": "update descriptions"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1615936345,
                                    "sbg:revision": 20,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1623274853,
                                    "sbg:revision": 21,
                                    "sbg:revisionNotes": "update docker image"
                                }
                            ],
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": []
                        },
                        "label": "find_unrelated",
                        "sbg:x": -343.6015625,
                        "sbg:y": -80.5
                    },
                    {
                        "id": "pca_byrel",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "related_file",
                                "source": "find_unrelated/out_related_file"
                            },
                            {
                                "id": "unrelated_file",
                                "source": "find_unrelated/out_unrelated_file"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "variant_include_file"
                            },
                            {
                                "id": "n_pcs",
                                "source": "n_pcs"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            }
                        ],
                        "out": [
                            {
                                "id": "pcair_output"
                            },
                            {
                                "id": "pcair_output_unrelated"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "smgogarten/genesis-relatedness/pca-byrel/21",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS File",
                                    "doc": "Input GDS file used for PCA. It is recommended to use an LD pruned file with all chromosomes.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "related_file",
                                    "type": "File",
                                    "label": "Related file",
                                    "doc": "RData file with related subjects.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "unrelated_file",
                                    "type": "File",
                                    "label": "Unrelated file",
                                    "doc": "RData file with unrelated subjects.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant include file",
                                    "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "32",
                                    "id": "n_pcs",
                                    "type": "int?",
                                    "label": "Number of PCs",
                                    "doc": "Number of PCs (Principal Components) to return.",
                                    "default": 32
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Output prefix."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "4",
                                    "id": "cpu",
                                    "type": "int?",
                                    "label": "cpu",
                                    "doc": "Number of CPUs used for each job.",
                                    "default": 4
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "pcair_output",
                                    "label": "RData file with PC-AiR PCs for all samples",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${\n  if (!inputs.out_prefix) {\n    var comm = \"pca.RData\"\n  } else {\n    var comm = inputs.out_prefix + \"_pca.RData\"\n  }\n  return comm\n}"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "pcair_output_unrelated",
                                    "doc": "RData file with PC-AiR PCs for unrelated samples",
                                    "label": "PCA byrel unrelated",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "${\n  if (!inputs.out_prefix) {\n    var comm = \"pca_unrel.RData\"\n  } else {\n    var comm = inputs.out_prefix + \"_pca_unrel.RData\"\n  }\n  return comm\n}"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "pca_byrel",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/pca_byrel.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "pca_byrel.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "coresMin": "${ if(inputs.cpu)\n        return inputs.cpu \n    else \n        return 4\n}"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "pca_byrel.config",
                                            "entry": "${\n\n    var cmd_line = ''\n    var arguments = []\n    \n    if(inputs.gds_file)\n     arguments.push(\"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\")\n    if(inputs.related_file)\n        arguments.push(\"related_file \\\"\" + inputs.related_file.path + \"\\\"\")\n    \n    if(inputs.unrelated_file)\n        arguments.push(\"unrelated_file \\\"\" + inputs.unrelated_file.path + \"\\\"\")\n        \n    if(inputs.variant_include_file)\n        arguments.push(\"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\")\n        \n    if(inputs.n_pcs)\n        arguments.push('n_pcs ' + inputs.n_pcs)\n    \n    if(inputs.out_prefix) {\n        arguments.push(cmd_line +='out_file \"' + inputs.out_prefix + '_pca.RData\"')\n    \n        arguments.push('out_file_unrel \"' + inputs.out_prefix + '_pca_unrel.RData\"')\n    }\n    \n    if(inputs.sample_include_file)\n        arguments.push('sample_include_file \"' + inputs.sample_include_file.path + '\"')\n    \n    return arguments.join(\"\\n\")\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "EnvVarRequirement",
                                    "envDef": [
                                        {
                                            "envName": "NSLOTS",
                                            "envValue": "${ return inputs.cpu }"
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "pca_byrel.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:content_hash": "aaf17382bd8e01b67f608549fea0af9c64b2feaab8ce0868b3727c4f84341b36e",
                            "sbg:contributors": [
                                "smgogarten"
                            ],
                            "sbg:createdBy": "smgogarten",
                            "sbg:createdOn": 1604988650,
                            "sbg:id": "smgogarten/genesis-relatedness/pca-byrel/21",
                            "sbg:image_url": null,
                            "sbg:latestRevision": 21,
                            "sbg:modifiedBy": "smgogarten",
                            "sbg:modifiedOn": 1623274896,
                            "sbg:project": "smgogarten/genesis-relatedness",
                            "sbg:projectName": "GENESIS relatedness",
                            "sbg:publisher": "sbg",
                            "sbg:revision": 21,
                            "sbg:revisionNotes": "update docker image",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604988650,
                                    "sbg:revision": 0,
                                    "sbg:revisionNotes": "Copy of boris_majic/topmed-optimization/pca-byrel/0"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604990205,
                                    "sbg:revision": 1,
                                    "sbg:revisionNotes": "base command ->  arguments"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604991094,
                                    "sbg:revision": 2,
                                    "sbg:revisionNotes": "fix output ports"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604991605,
                                    "sbg:revision": 3,
                                    "sbg:revisionNotes": "use out_prefix"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605570979,
                                    "sbg:revision": 4,
                                    "sbg:revisionNotes": "try redirecting stdout and stderr to same file"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605571421,
                                    "sbg:revision": 5,
                                    "sbg:revisionNotes": "variant_include should not be required"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605571692,
                                    "sbg:revision": 6,
                                    "sbg:revisionNotes": "redirecting stdout and stderr to same file just loses stderr"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605572819,
                                    "sbg:revision": 7,
                                    "sbg:revisionNotes": "use R as base command instead of Rscript"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605577005,
                                    "sbg:revision": 8,
                                    "sbg:revisionNotes": "add cpu as an argument"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605639326,
                                    "sbg:revision": 9,
                                    "sbg:revisionNotes": "set NSLOTS for parallelization"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605640696,
                                    "sbg:revision": 10,
                                    "sbg:revisionNotes": "set default n_pcs to 32"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605647083,
                                    "sbg:revision": 11,
                                    "sbg:revisionNotes": "try using EnvVarRequirement"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605648960,
                                    "sbg:revision": 12,
                                    "sbg:revisionNotes": "try EnvironmentDef"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605721746,
                                    "sbg:revision": 13,
                                    "sbg:revisionNotes": "fix EnvVarRequirement"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605726503,
                                    "sbg:revision": 14,
                                    "sbg:revisionNotes": "set NSLOTS automatically instead of as input"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605737241,
                                    "sbg:revision": 15,
                                    "sbg:revisionNotes": "add resource requirement"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1606356018,
                                    "sbg:revision": 16,
                                    "sbg:revisionNotes": "update input and output labels"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609448527,
                                    "sbg:revision": 17,
                                    "sbg:revisionNotes": "update descriptions"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609461048,
                                    "sbg:revision": 18,
                                    "sbg:revisionNotes": "fix missing quotation mark"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609461797,
                                    "sbg:revision": 19,
                                    "sbg:revisionNotes": "update output names and descriptions"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1615936346,
                                    "sbg:revision": 20,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1623274896,
                                    "sbg:revision": 21,
                                    "sbg:revisionNotes": "update docker image"
                                }
                            ],
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": []
                        },
                        "label": "pca_byrel",
                        "sbg:x": -94.60113525390625,
                        "sbg:y": -19.5
                    },
                    {
                        "id": "pca_plots",
                        "in": [
                            {
                                "id": "pca_file",
                                "source": "pca_byrel/pcair_output"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "n_pairs",
                                "source": "n_pairs"
                            },
                            {
                                "id": "group",
                                "source": "group"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            }
                        ],
                        "out": [
                            {
                                "id": "pca_plots"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "smgogarten/genesis-relatedness/pca-plots/7",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "pca_file",
                                    "type": "File",
                                    "label": "PCA File",
                                    "doc": "RData file containing pcair object (output by pca_byrel tool)",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "phenotype_file",
                                    "type": "File?",
                                    "label": "Phenotype file",
                                    "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for color-coding PCA plots by group.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "6",
                                    "id": "n_pairs",
                                    "type": "int?",
                                    "label": "Number of PCs",
                                    "doc": "Number of PCs to include in the pairs plot."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "NA",
                                    "id": "group",
                                    "type": "string?",
                                    "label": "Group",
                                    "doc": "Name of column in phenotype_file containing group variable for color-coding plots."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "pca_plots",
                                    "doc": "PC plots",
                                    "label": "PC plots",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.p*"
                                    }
                                }
                            ],
                            "label": "pca_plots",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/pca_plots.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "pca_plots.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "pca_plots.config",
                                            "entry": "${\n var cmd_line = \"\"\n\nif(inputs.pca_file)\n    cmd_line += \"pca_file \\\"\" + inputs.pca_file.path + \"\\\"\\n\"\n\nif(inputs.n_pairs)\n    cmd_line += \"n_pairs \\\"\" + inputs.n_pairs + \"\\\"\\n\"\n    \n\nif(inputs.out_prefix) {\n     cmd_line += \"out_file_scree \\\"\" + inputs.out_prefix + \"_pca_scree.pdf\\\"\\n\"\n    \n     cmd_line += \"out_file_pc12 \\\"\" + inputs.out_prefix + \"_pca_pc12.pdf\\\"\\n\"\n   \n     cmd_line += \"out_file_parcoord \\\"\" + inputs.out_prefix + \"_pca_parcoord.pdf\\\"\\n\"\n   \n    cmd_line += 'out_file_pairs \"' + inputs.out_prefix + '_pca_pairs.png\"\\n'\n}\n\nif(inputs.group)\n    cmd_line += 'group \"' + inputs.group + '\"\\n'\n\nif(inputs.phenotype_file)\n    cmd_line += 'phenotype_file \"' + inputs.phenotype_file.path + '\"\\n'\n\n\n return cmd_line \n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "pca_plots.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:content_hash": "ab4a2bebd6f68a0a465fc6a6487ce68dee2924ef24c5a48eb5df6be60dd3c8538",
                            "sbg:contributors": [
                                "smgogarten"
                            ],
                            "sbg:createdBy": "smgogarten",
                            "sbg:createdOn": 1604988656,
                            "sbg:id": "smgogarten/genesis-relatedness/pca-plots/7",
                            "sbg:image_url": null,
                            "sbg:latestRevision": 7,
                            "sbg:modifiedBy": "smgogarten",
                            "sbg:modifiedOn": 1623274941,
                            "sbg:project": "smgogarten/genesis-relatedness",
                            "sbg:projectName": "GENESIS relatedness",
                            "sbg:publisher": "sbg",
                            "sbg:revision": 7,
                            "sbg:revisionNotes": "update docker image",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1604988656,
                                    "sbg:revision": 0,
                                    "sbg:revisionNotes": "Copy of boris_majic/topmed-optimization/pca-plots/0"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605574521,
                                    "sbg:revision": 1,
                                    "sbg:revisionNotes": "update to new format"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1605576107,
                                    "sbg:revision": 2,
                                    "sbg:revisionNotes": "specify docker version"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1606356113,
                                    "sbg:revision": 3,
                                    "sbg:revisionNotes": "update input and output labels"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609370862,
                                    "sbg:revision": 4,
                                    "sbg:revisionNotes": "update descriptions"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1609449160,
                                    "sbg:revision": 5,
                                    "sbg:revisionNotes": "use devel docker image to allow data.frame input for phenotype file"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1615936347,
                                    "sbg:revision": 6,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1623274941,
                                    "sbg:revision": 7,
                                    "sbg:revisionNotes": "update docker image"
                                }
                            ],
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": []
                        },
                        "label": "pca_plots",
                        "sbg:x": 206.48231506347656,
                        "sbg:y": -174.49720764160156
                    },
                    {
                        "id": "pc_variant_correlation",
                        "in": [
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "gds_file_full",
                                "source": [
                                    "gds_file_full"
                                ]
                            },
                            {
                                "id": "variant_include_file",
                                "source": [
                                    "pruned_variant_file"
                                ]
                            },
                            {
                                "id": "pca_file",
                                "source": "pca_byrel/pcair_output_unrelated"
                            },
                            {
                                "id": "n_corr_vars",
                                "source": "n_corr_vars"
                            },
                            {
                                "id": "n_pcs",
                                "source": "n_pcs"
                            },
                            {
                                "id": "n_pcs_plot",
                                "source": "n_pcs_plot"
                            },
                            {
                                "id": "n_perpage",
                                "source": "n_perpage"
                            },
                            {
                                "id": "run_correlation",
                                "source": "run_correlation"
                            }
                        ],
                        "out": [
                            {
                                "id": "pc_correlation_plots"
                            },
                            {
                                "id": "pca_corr_gds"
                            }
                        ],
                        "run": {
                            "class": "Workflow",
                            "cwlVersion": "v1.2",
                            "id": "smgogarten/genesis-relatedness-pre-build/pc-variant-correlation/2",
                            "label": "PC-variant correlation",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "inputs": [
                                {
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files.",
                                    "sbg:x": -58,
                                    "sbg:y": 30
                                },
                                {
                                    "id": "gds_file_full",
                                    "sbg:fileTypes": "GDS",
                                    "type": "File[]",
                                    "label": "Full GDS Files",
                                    "doc": "GDS files (one per chromosome) used to calculate PC-variant correlations.",
                                    "sbg:x": -25.0925350189209,
                                    "sbg:y": 187.01194763183594
                                },
                                {
                                    "id": "variant_include_file",
                                    "sbg:fileTypes": "RDATA",
                                    "type": "File[]",
                                    "label": "Variant include file",
                                    "doc": "RData files (one per chromosome) with vector of variant.id to include. These variants will be added to the set of randomly selected variants. It is recommended to provide the set of pruned variants used for PCA.",
                                    "sbg:x": -106,
                                    "sbg:y": -85
                                },
                                {
                                    "id": "pca_file",
                                    "sbg:fileTypes": "RDATA",
                                    "type": "File",
                                    "label": "PCA file",
                                    "doc": "RData file with PCA results for unrelated samples",
                                    "sbg:x": 11.241790771484375,
                                    "sbg:y": -188.33731079101562
                                },
                                {
                                    "id": "n_corr_vars",
                                    "type": "int?",
                                    "label": "Number of variants to select",
                                    "doc": "Randomly select this number of variants distributed across the entire genome to use for PC-variant correlation. If running on a single chromosome, the variants returned will be scaled by the proportion of that chromosome in the genome.",
                                    "sbg:exposed": true
                                },
                                {
                                    "id": "n_pcs",
                                    "type": "int?",
                                    "label": "Number of PCs",
                                    "doc": "Number of PCs (Principal Components) to return.",
                                    "sbg:x": -67.57342529296875,
                                    "sbg:y": -311.5232849121094
                                },
                                {
                                    "id": "n_pcs_plot",
                                    "type": "int?",
                                    "label": "Number of PCs to plot",
                                    "doc": "Number of PCs to plot.",
                                    "sbg:exposed": true
                                },
                                {
                                    "id": "n_perpage",
                                    "type": "int?",
                                    "label": "Number of plots per page",
                                    "doc": "Number of PC-variant correlation plots to stack in a single page. The number of png files generated will be ceiling(n_pcs_plot/n_perpage).",
                                    "sbg:exposed": true
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "pc_correlation_plots",
                                    "outputSource": [
                                        "pca_corr_plots/pca_corr_plots"
                                    ],
                                    "sbg:fileTypes": "PNG",
                                    "type": "File[]?",
                                    "label": "PC-variant correlation plots",
                                    "doc": "PC-variant correlation plots",
                                    "sbg:x": 678.9522094726562,
                                    "sbg:y": -9.483582496643066
                                },
                                {
                                    "id": "pca_corr_gds",
                                    "outputSource": [
                                        "pca_corr/pca_corr_gds"
                                    ],
                                    "sbg:fileTypes": "GDS",
                                    "type": "File[]?",
                                    "label": "PC-variant correlation",
                                    "doc": "GDS file with PC-variant correlation results",
                                    "sbg:x": 581.6796875,
                                    "sbg:y": 159
                                }
                            ],
                            "steps": [
                                {
                                    "id": "pca_corr_vars",
                                    "in": [
                                        {
                                            "id": "gds_file",
                                            "source": "gds_file_full"
                                        },
                                        {
                                            "id": "variant_include_file",
                                            "source": "variant_include_file"
                                        },
                                        {
                                            "id": "out_prefix",
                                            "source": "out_prefix"
                                        },
                                        {
                                            "id": "n_corr_vars",
                                            "source": "n_corr_vars"
                                        },
                                        {
                                            "id": "chromosome",
                                            "valueFrom": "${\n    if (inputs.gds_file.nameroot.includes('chr')) {\n        var parts = inputs.gds_file.nameroot.split('chr')\n        var chrom = parts[1]\n    } else {\n        var chrom = \"NA\"\n    }\n    return chrom\n}"
                                        }
                                    ],
                                    "out": [
                                        {
                                            "id": "pca_corr_vars"
                                        }
                                    ],
                                    "run": {
                                        "class": "CommandLineTool",
                                        "cwlVersion": "v1.2",
                                        "$namespaces": {
                                            "sbg": "https://sevenbridges.com"
                                        },
                                        "id": "smgogarten/genesis-relatedness/pca-corr-vars/7",
                                        "baseCommand": [
                                            "R -q --vanilla"
                                        ],
                                        "inputs": [
                                            {
                                                "sbg:category": "Input Files",
                                                "id": "gds_file",
                                                "type": "File",
                                                "label": "GDS File",
                                                "doc": "Input GDS file",
                                                "sbg:fileTypes": "GDS"
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "variant_include_file",
                                                "type": "File?",
                                                "label": "Variant include file",
                                                "doc": "RData file with vector of variant.id to include. These variants will be added to the set of randomly selected variants. It is recommended to provide the set of pruned variants used for PCA.",
                                                "sbg:fileTypes": "RDATA"
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "out_prefix",
                                                "type": "string?",
                                                "label": "Output prefix",
                                                "doc": "Prefix for output files."
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "sbg:toolDefaultValue": "10e6",
                                                "id": "n_corr_vars",
                                                "type": "int?",
                                                "label": "Number of variants to select",
                                                "doc": "Randomly select this number of variants distributed across the entire genome to use for PC-variant correlation. If running on a single chromosome, the variants returned will be scaled by the proportion of that chromosome in the genome."
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "chromosome",
                                                "type": "int",
                                                "inputBinding": {
                                                    "prefix": "--chromosome",
                                                    "shellQuote": false,
                                                    "position": 3
                                                },
                                                "label": "Chromosome",
                                                "doc": "Run on this chromosome only. 23=X, 24=Y"
                                            }
                                        ],
                                        "outputs": [
                                            {
                                                "id": "pca_corr_vars",
                                                "doc": "RData file with a randomly selected set of variant.ids distributed across the genome, plus any variants from variant_include_file.",
                                                "label": "Variants to use for PC correlation",
                                                "type": "File?",
                                                "outputBinding": {
                                                    "glob": "*.RData"
                                                },
                                                "sbg:fileTypes": "RDATA"
                                            }
                                        ],
                                        "label": "pca_corr_vars",
                                        "arguments": [
                                            {
                                                "prefix": "--args",
                                                "shellQuote": false,
                                                "position": 2,
                                                "valueFrom": "pca_corr_vars.config"
                                            },
                                            {
                                                "prefix": "<",
                                                "shellQuote": false,
                                                "position": 4,
                                                "valueFrom": "/usr/local/analysis_pipeline/R/pca_corr_vars.R"
                                            }
                                        ],
                                        "requirements": [
                                            {
                                                "class": "ShellCommandRequirement"
                                            },
                                            {
                                                "class": "DockerRequirement",
                                                "dockerPull": "uwgac/topmed-master:2.8.1"
                                            },
                                            {
                                                "class": "InitialWorkDirRequirement",
                                                "listing": [
                                                    {
                                                        "entryname": "pca_corr_vars.config",
                                                        "entry": "${\n    var config = \"\"\n    //if (inputs.gds_file.nameroot.includes('chr'))\n    //{\n    //    var parts = inputs.gds_file.nameroot.split('chr')\n    //    var outfile_temp = 'pca_corr_vars_chr' + parts[1] + '.RData'\n    //} else {\n    //    var outfile_temp = 'pca_corr_vars.RData'\n    //}\n    var outfile_temp = 'pca_corr_vars_chr .RData'\n    if(inputs.out_prefix){\n        outfile_temp = inputs.out_prefix + '_' + outfile_temp\n    }\n    config += 'out_file \"' + outfile_temp + '\"\\n'\n\n    \n    if (inputs.variant_include_file)\n      config += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n      \n    if(inputs.gds_file)\n        config += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n\n    config += \"segment_file \\\"/usr/local/analysis_pipeline/segments_hg38.txt\\\"\\n\"\n\n    return config\n}",
                                                        "writable": false
                                                    }
                                                ]
                                            },
                                            {
                                                "class": "InlineJavascriptRequirement"
                                            }
                                        ],
                                        "hints": [
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "job.out.log"
                                            },
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "pca_corr_vars.config"
                                            }
                                        ],
                                        "stdout": "job.out.log",
                                        "sbg:appVersion": [
                                            "v1.2"
                                        ],
                                        "sbg:content_hash": "a365e9e704d2f30712910d5e8d7ad9ae8db50b827df7f666f738757c362e99602",
                                        "sbg:contributors": [
                                            "smgogarten"
                                        ],
                                        "sbg:createdBy": "smgogarten",
                                        "sbg:createdOn": 1608671560,
                                        "sbg:id": "smgogarten/genesis-relatedness/pca-corr-vars/7",
                                        "sbg:image_url": null,
                                        "sbg:latestRevision": 7,
                                        "sbg:modifiedBy": "smgogarten",
                                        "sbg:modifiedOn": 1623457859,
                                        "sbg:project": "smgogarten/genesis-relatedness",
                                        "sbg:projectName": "GENESIS relatedness",
                                        "sbg:publisher": "sbg",
                                        "sbg:revision": 7,
                                        "sbg:revisionNotes": "update to CWL 1.2",
                                        "sbg:revisionsInfo": [
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1608671560,
                                                "sbg:revision": 0,
                                                "sbg:revisionNotes": null
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1608672869,
                                                "sbg:revision": 1,
                                                "sbg:revisionNotes": "tool for selecting variant for pc correlation"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1608674522,
                                                "sbg:revision": 2,
                                                "sbg:revisionNotes": "add newline in config, chromosome arg is required"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1608681003,
                                                "sbg:revision": 3,
                                                "sbg:revisionNotes": "need blank space for chromosome number"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609370853,
                                                "sbg:revision": 4,
                                                "sbg:revisionNotes": "update descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609448724,
                                                "sbg:revision": 5,
                                                "sbg:revisionNotes": "update descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1615936348,
                                                "sbg:revision": 6,
                                                "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1623457859,
                                                "sbg:revision": 7,
                                                "sbg:revisionNotes": "update to CWL 1.2"
                                            }
                                        ],
                                        "sbg:sbgMaintained": false,
                                        "sbg:validationErrors": []
                                    },
                                    "label": "pca_corr_vars",
                                    "scatter": [
                                        "gds_file",
                                        "variant_include_file"
                                    ],
                                    "scatterMethod": "dotproduct",
                                    "sbg:x": 121,
                                    "sbg:y": -8
                                },
                                {
                                    "id": "pca_corr",
                                    "in": [
                                        {
                                            "id": "gds_file",
                                            "source": "gds_file_full"
                                        },
                                        {
                                            "id": "variant_include_file",
                                            "source": "pca_corr_vars/pca_corr_vars"
                                        },
                                        {
                                            "id": "pca_file",
                                            "source": "pca_file"
                                        },
                                        {
                                            "id": "n_pcs_corr",
                                            "source": "n_pcs"
                                        },
                                        {
                                            "id": "out_prefix",
                                            "source": "out_prefix"
                                        }
                                    ],
                                    "out": [
                                        {
                                            "id": "pca_corr_gds"
                                        }
                                    ],
                                    "run": {
                                        "class": "CommandLineTool",
                                        "cwlVersion": "v1.2",
                                        "$namespaces": {
                                            "sbg": "https://sevenbridges.com"
                                        },
                                        "id": "smgogarten/genesis-relatedness/pca-corr/15",
                                        "baseCommand": [
                                            "R -q --vanilla"
                                        ],
                                        "inputs": [
                                            {
                                                "sbg:category": "Input Files",
                                                "id": "gds_file",
                                                "type": "File",
                                                "label": "GDS File",
                                                "doc": "Input GDS file",
                                                "sbg:fileTypes": "GDS"
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "variant_include_file",
                                                "type": "File?",
                                                "label": "Variant include file",
                                                "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included.",
                                                "sbg:fileTypes": "RDATA"
                                            },
                                            {
                                                "id": "pca_file",
                                                "type": "File",
                                                "label": "PCA file",
                                                "doc": "RData file with PCA results for unrelated samples",
                                                "sbg:fileTypes": "RDATA"
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "sbg:toolDefaultValue": "32",
                                                "id": "n_pcs_corr",
                                                "type": "int?",
                                                "label": "Number of PCs",
                                                "doc": "Number of PCs (Principal Components) to use for PC-variant correlation",
                                                "default": 32
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "out_prefix",
                                                "type": "string?",
                                                "label": "Output prefix",
                                                "doc": "Prefix for output files."
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "id": "chromosome",
                                                "type": "int?",
                                                "inputBinding": {
                                                    "prefix": "--chromosome",
                                                    "shellQuote": false,
                                                    "position": 3
                                                },
                                                "label": "Chromosome",
                                                "doc": "Run on this chromosome only. 23=X, 24=Y"
                                            }
                                        ],
                                        "outputs": [
                                            {
                                                "id": "pca_corr_gds",
                                                "doc": "GDS file with PC-SNP correlation results",
                                                "label": "PC-SNP correlation",
                                                "type": "File?",
                                                "outputBinding": {
                                                    "glob": "*.gds"
                                                },
                                                "sbg:fileTypes": "GDS"
                                            }
                                        ],
                                        "label": "pca_corr",
                                        "arguments": [
                                            {
                                                "prefix": "--args",
                                                "shellQuote": false,
                                                "position": 2,
                                                "valueFrom": "pca_corr.config"
                                            },
                                            {
                                                "prefix": "<",
                                                "shellQuote": false,
                                                "position": 4,
                                                "valueFrom": "/usr/local/analysis_pipeline/R/pca_corr.R"
                                            }
                                        ],
                                        "requirements": [
                                            {
                                                "class": "ShellCommandRequirement"
                                            },
                                            {
                                                "class": "ResourceRequirement",
                                                "coresMin": 4
                                            },
                                            {
                                                "class": "DockerRequirement",
                                                "dockerPull": "uwgac/topmed-master:2.10.0"
                                            },
                                            {
                                                "class": "InitialWorkDirRequirement",
                                                "listing": [
                                                    {
                                                        "entryname": "pca_corr.config",
                                                        "entry": "${\n    var config = \"\"\n    if (inputs.gds_file.nameroot.includes('chr'))\n    {\n        var parts = inputs.gds_file.nameroot.split('chr')\n        var outfile_temp = 'pca_corr_chr' + parts[1] + '.gds'\n    } else {\n        var outfile_temp = 'pca_corr.gds'\n    }\n    if(inputs.out_prefix){\n        outfile_temp = inputs.out_prefix + '_' + outfile_temp\n    }\n    config += 'out_file \"' + outfile_temp + '\"\\n'\n\n    \n  \tif (inputs.n_pcs_corr) {\n        config += \"n_pcs \" + inputs.n_pcs_corr + \"\\n\"\n  \t}\n    \n    if (inputs.variant_include_file)\n      config += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n      \n    if(inputs.gds_file)\n        config += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n\n    if(inputs.pca_file)\n        config += \"pca_file \\\"\" + inputs.pca_file.path + \"\\\"\\n\"\n        \n    return config\n}",
                                                        "writable": false
                                                    }
                                                ]
                                            },
                                            {
                                                "class": "EnvVarRequirement",
                                                "envDef": [
                                                    {
                                                        "envName": "NSLOTS",
                                                        "envValue": "${ return runtime.cores }"
                                                    }
                                                ]
                                            },
                                            {
                                                "class": "InlineJavascriptRequirement"
                                            }
                                        ],
                                        "hints": [
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "pca_corr.config"
                                            },
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "job.out.log"
                                            }
                                        ],
                                        "stdout": "job.out.log",
                                        "sbg:appVersion": [
                                            "v1.2"
                                        ],
                                        "sbg:content_hash": "a9156cc90b7882de0c565ed53f8ebda666d0d8b17a6a9b98a1e12d12808da7d50",
                                        "sbg:contributors": [
                                            "smgogarten"
                                        ],
                                        "sbg:createdBy": "smgogarten",
                                        "sbg:createdOn": 1604988662,
                                        "sbg:id": "smgogarten/genesis-relatedness/pca-corr/15",
                                        "sbg:image_url": null,
                                        "sbg:latestRevision": 15,
                                        "sbg:modifiedBy": "smgogarten",
                                        "sbg:modifiedOn": 1623457838,
                                        "sbg:project": "smgogarten/genesis-relatedness",
                                        "sbg:projectName": "GENESIS relatedness",
                                        "sbg:publisher": "sbg",
                                        "sbg:revision": 15,
                                        "sbg:revisionNotes": "update docker image",
                                        "sbg:revisionsInfo": [
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1604988662,
                                                "sbg:revision": 0,
                                                "sbg:revisionNotes": "Copy of boris_majic/topmed-optimization/pca-corr/0"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605577131,
                                                "sbg:revision": 1,
                                                "sbg:revisionNotes": "save complete R output"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605639828,
                                                "sbg:revision": 2,
                                                "sbg:revisionNotes": "add NSLOTS for parallelization"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605640837,
                                                "sbg:revision": 3,
                                                "sbg:revisionNotes": "set default n_pcs to 32"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605726814,
                                                "sbg:revision": 4,
                                                "sbg:revisionNotes": "use runtime.cores to set NSLOTS"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605729508,
                                                "sbg:revision": 5,
                                                "sbg:revisionNotes": "try checking the multi-thread box"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605730462,
                                                "sbg:revision": 6,
                                                "sbg:revisionNotes": "try setting min CPUs"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606335780,
                                                "sbg:revision": 7,
                                                "sbg:revisionNotes": "put chromosome number in output file name"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606341005,
                                                "sbg:revision": 8,
                                                "sbg:revisionNotes": "fix output file name"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606355536,
                                                "sbg:revision": 9,
                                                "sbg:revisionNotes": "update input and output labels"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606356416,
                                                "sbg:revision": 10,
                                                "sbg:revisionNotes": "update labels"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1608672691,
                                                "sbg:revision": 11,
                                                "sbg:revisionNotes": "fix position of chromosome argument"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609370838,
                                                "sbg:revision": 12,
                                                "sbg:revisionNotes": "update descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609448619,
                                                "sbg:revision": 13,
                                                "sbg:revisionNotes": "update descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1615936349,
                                                "sbg:revision": 14,
                                                "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1623457838,
                                                "sbg:revision": 15,
                                                "sbg:revisionNotes": "update docker image"
                                            }
                                        ],
                                        "sbg:sbgMaintained": false,
                                        "sbg:validationErrors": []
                                    },
                                    "label": "pca_corr",
                                    "scatter": [
                                        "gds_file",
                                        "variant_include_file"
                                    ],
                                    "scatterMethod": "dotproduct",
                                    "sbg:x": 350,
                                    "sbg:y": 30
                                },
                                {
                                    "id": "pca_corr_plots",
                                    "in": [
                                        {
                                            "id": "n_pcs_plot",
                                            "source": "n_pcs_plot"
                                        },
                                        {
                                            "id": "corr_file",
                                            "source": [
                                                "pca_corr/pca_corr_gds"
                                            ]
                                        },
                                        {
                                            "id": "n_perpage",
                                            "source": "n_perpage"
                                        },
                                        {
                                            "id": "out_prefix",
                                            "source": "out_prefix"
                                        }
                                    ],
                                    "out": [
                                        {
                                            "id": "pca_corr_plots"
                                        }
                                    ],
                                    "run": {
                                        "class": "CommandLineTool",
                                        "cwlVersion": "v1.2",
                                        "$namespaces": {
                                            "sbg": "https://sevenbridges.com"
                                        },
                                        "id": "smgogarten/genesis-relatedness/pca-corr-plots/10",
                                        "baseCommand": [],
                                        "inputs": [
                                            {
                                                "sbg:category": "Input Options",
                                                "sbg:toolDefaultValue": "20",
                                                "id": "n_pcs_plot",
                                                "type": "int?",
                                                "label": "Number of PCs to plot",
                                                "doc": "Number of PCs to plot.",
                                                "default": 20
                                            },
                                            {
                                                "sbg:category": "Input File",
                                                "id": "corr_file",
                                                "type": "File[]",
                                                "label": "PC correlation file",
                                                "doc": "PC correlation file",
                                                "sbg:fileTypes": "GDS"
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "sbg:toolDefaultValue": "4",
                                                "id": "n_perpage",
                                                "type": "int?",
                                                "label": "Number of plots per page",
                                                "doc": "Number of PC-variant correlation plots to stack in a single page. The number of png files generated will be ceiling(n_pcs_plot/n_perpage).",
                                                "default": 4
                                            },
                                            {
                                                "sbg:category": "Input Options",
                                                "sbg:toolDefaultValue": "pca_corr",
                                                "id": "out_prefix",
                                                "type": "string?",
                                                "label": "Output prefix",
                                                "doc": "Prefix for output files."
                                            }
                                        ],
                                        "outputs": [
                                            {
                                                "id": "pca_corr_plots",
                                                "doc": "PC-variant correlation plots",
                                                "label": "PC-variant correlation plots",
                                                "type": "File[]?",
                                                "outputBinding": {
                                                    "glob": "*.png"
                                                },
                                                "sbg:fileTypes": "PNG"
                                            }
                                        ],
                                        "label": "pca_corr_plots",
                                        "arguments": [
                                            {
                                                "prefix": "<",
                                                "shellQuote": false,
                                                "position": 3,
                                                "valueFrom": "/usr/local/analysis_pipeline/R/pca_corr_plots.R"
                                            },
                                            {
                                                "prefix": "",
                                                "shellQuote": false,
                                                "position": 0,
                                                "valueFrom": "${\n    var cmd_line = \"\"\n    for (var i=0; i<inputs.corr_file.length; i++)\n        cmd_line += \"ln -s \" + inputs.corr_file[i].path + \" \" + inputs.corr_file[i].basename + \" && \"\n    \n    return cmd_line\n}"
                                            },
                                            {
                                                "prefix": "--args",
                                                "shellQuote": false,
                                                "position": 2,
                                                "valueFrom": "pca_corr_plots.config"
                                            },
                                            {
                                                "prefix": "",
                                                "shellQuote": false,
                                                "position": 1,
                                                "valueFrom": "R -q --vanilla"
                                            }
                                        ],
                                        "requirements": [
                                            {
                                                "class": "ShellCommandRequirement"
                                            },
                                            {
                                                "class": "DockerRequirement",
                                                "dockerPull": "uwgac/topmed-master:2.10.0"
                                            },
                                            {
                                                "class": "InitialWorkDirRequirement",
                                                "listing": [
                                                    {
                                                        "entryname": "pca_corr_plots.config",
                                                        "entry": "${ \n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s))\n    }\n\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    var a_file = inputs.corr_file[0]\n    var chr = find_chromosome(a_file.basename);\n    var path = a_file.path.split('chr'+chr);\n\n    var chr_array = [];\n    var chrom_num;\n    for (var i = 0; i < inputs.corr_file.length; i++) \n    {\n        chrom_num = find_chromosome(inputs.corr_file[i].nameroot)\n            \n        chr_array.push(chrom_num)\n    }\n        \n    chr_array = chr_array.sort((a, b) => a.localeCompare(b, 'en', {numeric: true, ignorePunctuation: true}))\n        \n    var chrs = chr_array.join(' ')\n    \n    var arguments = []\n    arguments.push('corr_file ' + '\"' + path[0].split('/').pop() + 'chr ' + path[1] + '\"')\n    arguments.push('chromosomes \"' + chrs + '\"')\n    \n    arguments.push('n_pcs ' + inputs.n_pcs_plot)\n    arguments.push('n_perpage ' + inputs.n_perpage)\n    \n    if(inputs.out_prefix) {\n        arguments.push('out_prefix \"' + inputs.out_prefix + '\"')\n    }\n    \n    arguments.push('\\n')\n    return arguments.join('\\n')\n}",
                                                        "writable": false
                                                    }
                                                ]
                                            },
                                            {
                                                "class": "InlineJavascriptRequirement"
                                            }
                                        ],
                                        "hints": [
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "job.out.log"
                                            },
                                            {
                                                "class": "sbg:SaveLogs",
                                                "value": "pca_corr_plots.config"
                                            }
                                        ],
                                        "stdout": "job.out.log",
                                        "sbg:appVersion": [
                                            "v1.2"
                                        ],
                                        "sbg:content_hash": "a8f8ac8c3e89b92e8988df33f4627f59b9161386b97798abaeb41f81d667a8eb1",
                                        "sbg:contributors": [
                                            "smgogarten"
                                        ],
                                        "sbg:createdBy": "smgogarten",
                                        "sbg:createdOn": 1604988665,
                                        "sbg:id": "smgogarten/genesis-relatedness/pca-corr-plots/10",
                                        "sbg:image_url": null,
                                        "sbg:latestRevision": 10,
                                        "sbg:modifiedBy": "smgogarten",
                                        "sbg:modifiedOn": 1623457909,
                                        "sbg:project": "smgogarten/genesis-relatedness",
                                        "sbg:projectName": "GENESIS relatedness",
                                        "sbg:publisher": "sbg",
                                        "sbg:revision": 10,
                                        "sbg:revisionNotes": "update docker image",
                                        "sbg:revisionsInfo": [
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1604988665,
                                                "sbg:revision": 0,
                                                "sbg:revisionNotes": "Copy of boris_majic/topmed-optimization/pca-corr-plots/0"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605730277,
                                                "sbg:revision": 1,
                                                "sbg:revisionNotes": "initial update"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1605765307,
                                                "sbg:revision": 2,
                                                "sbg:revisionNotes": "copy config javascript from assoc_plots tool"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606334463,
                                                "sbg:revision": 3,
                                                "sbg:revisionNotes": "update javascript for config file"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606340439,
                                                "sbg:revision": 4,
                                                "sbg:revisionNotes": "fix filenames in config"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606355355,
                                                "sbg:revision": 5,
                                                "sbg:revisionNotes": "update input and output labels"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1606780666,
                                                "sbg:revision": 6,
                                                "sbg:revisionNotes": "use basename instead of path.split"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609370845,
                                                "sbg:revision": 7,
                                                "sbg:revisionNotes": "update descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1609448682,
                                                "sbg:revision": 8,
                                                "sbg:revisionNotes": "upate descriptions"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1615936369,
                                                "sbg:revision": 9,
                                                "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: (uncommitted file)"
                                            },
                                            {
                                                "sbg:modifiedBy": "smgogarten",
                                                "sbg:modifiedOn": 1623457909,
                                                "sbg:revision": 10,
                                                "sbg:revisionNotes": "update docker image"
                                            }
                                        ],
                                        "sbg:sbgMaintained": false,
                                        "sbg:validationErrors": []
                                    },
                                    "label": "pca_corr_plots",
                                    "sbg:x": 489.2178955078125,
                                    "sbg:y": -117.57015228271484
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ScatterFeatureRequirement"
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                },
                                {
                                    "class": "StepInputExpressionRequirement"
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:content_hash": "a07e18736ce639c6ec7fd7e582b1398500b4313e0ddd4d2584eaa6fa7516ce056",
                            "sbg:contributors": [
                                "smgogarten"
                            ],
                            "sbg:createdBy": "smgogarten",
                            "sbg:createdOn": 1623705933,
                            "sbg:id": "smgogarten/genesis-relatedness-pre-build/pc-variant-correlation/2",
                            "sbg:image_url": "https://platform.sb.biodatacatalyst.nhlbi.nih.gov/ns/brood/images/smgogarten/genesis-relatedness-pre-build/pc-variant-correlation/2.png",
                            "sbg:latestRevision": 2,
                            "sbg:modifiedBy": "smgogarten",
                            "sbg:modifiedOn": 1637104100,
                            "sbg:original_source": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/smgogarten/genesis-relatedness-pre-build/pc-variant-correlation/2/raw/",
                            "sbg:project": "smgogarten/genesis-relatedness-pre-build",
                            "sbg:projectName": "GENESIS relatedness - Pre-build",
                            "sbg:publisher": "sbg",
                            "sbg:revision": 2,
                            "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: c2eb59b",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1623705933,
                                    "sbg:revision": 0,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: 07a9d69"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1623710906,
                                    "sbg:revision": 1,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: 2097da0"
                                },
                                {
                                    "sbg:modifiedBy": "smgogarten",
                                    "sbg:modifiedOn": 1637104100,
                                    "sbg:revision": 2,
                                    "sbg:revisionNotes": "Uploaded using sbpack v2020.10.05. \nSource: \nrepo: git@github.com:UW-GAC/analysis_pipeline_cwl.git\nfile: \ncommit: c2eb59b"
                                }
                            ],
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": []
                        },
                        "label": "PC-variant correlation",
                        "when": "$(inputs.run_correlation)",
                        "sbg:x": 85.90345764160156,
                        "sbg:y": 199.70834350585938
                    }
                ],
                "requirements": [
                    {
                        "class": "SubworkflowFeatureRequirement"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "Ancestry and Relatedness"
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/pc-air/1.png",
                "sbg:original_source": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/smgogarten/genesis-relatedness/pc-air/32/raw/",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105848,
                        "sbg:revisionNotes": "Workflow decomposed"
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1637106003,
                        "sbg:revisionNotes": "updated PC-AiR (nov 17th)"
                    }
                ],
                "sbg:toolkit": "UW-GAC Ancestry and Relatedness",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/pc-air/1",
                "sbg:revision": 1,
                "sbg:revisionNotes": "updated PC-AiR (nov 17th)",
                "sbg:modifiedOn": 1637106003,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105848,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 1,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "af68ea9db03682b5a14e3ff34c56a09dfa778542170f5f83a72137779f6230ed5"
            },
            "label": "PC-AiR",
            "hints": [
                {
                    "class": "sbg:AWSInstanceType",
                    "value": "m5.4xlarge;ebs-gp2;1024"
                }
            ],
            "sbg:x": 517.511474609375,
            "sbg:y": -222.83387756347656
        },
        {
            "id": "pc_relate",
            "in": [
                {
                    "id": "pca_file",
                    "source": "pc_air/pcair_output"
                },
                {
                    "id": "gds_file",
                    "source": "ld_pruning/pruned_gds_output"
                },
                {
                    "id": "n_pcs",
                    "default": 8
                },
                {
                    "id": "variant_include_file",
                    "source": "merge_rds_arrays/output"
                },
                {
                    "id": "phenotype_file",
                    "source": "phenotype_file"
                },
                {
                    "id": "ibd_probs",
                    "default": true
                }
            ],
            "out": [
                {
                    "id": "pcrelate_plots"
                },
                {
                    "id": "pcrelate_output"
                },
                {
                    "id": "pcrelate_matrix"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/pc-relate/0",
                "doc": "This workflow estimates kinship and IBD sharing probabilities between all pairs of samples using the PC-Relate method, which accounts for population structure by conditioning on ancestry PCs.",
                "label": "PC-Relate",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "pca_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File",
                        "label": "PCA file",
                        "doc": "RData file with PCA results from PC-AiR workflow; used to adjust for population structure.",
                        "sbg:x": -189.6508026123047,
                        "sbg:y": 53.30672836303711
                    },
                    {
                        "id": "gds_file",
                        "sbg:fileTypes": "GDS",
                        "type": "File",
                        "label": "GDS File",
                        "doc": "Input GDS file. It is recommended to use an LD pruned file with all chromosomes.",
                        "sbg:x": -205.0176544189453,
                        "sbg:y": 181.32745361328125
                    },
                    {
                        "id": "sample_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Sample Include file",
                        "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                        "sbg:x": -175.52853393554688,
                        "sbg:y": -72.8362808227539
                    },
                    {
                        "id": "n_pcs",
                        "type": "int?",
                        "label": "Number of PCs",
                        "doc": "Number of PCs (Principal Components) in `pca_file` to use in adjusting for ancestry.",
                        "sbg:toolDefaultValue": "3",
                        "sbg:x": -61.75441360473633,
                        "sbg:y": 284.0424499511719
                    },
                    {
                        "id": "variant_include_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Variant include file",
                        "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included.",
                        "sbg:x": -169.77308654785156,
                        "sbg:y": -226.9329071044922
                    },
                    {
                        "id": "variant_block_size",
                        "type": "int?",
                        "label": "Variant block size",
                        "doc": "Number of variants to read in a single block.",
                        "sbg:toolDefaultValue": "1024",
                        "sbg:x": -60.2248649597168,
                        "sbg:y": -132.0818634033203
                    },
                    {
                        "id": "out_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Prefix for output files.",
                        "sbg:x": -65.81759643554688,
                        "sbg:y": 105.9605941772461
                    },
                    {
                        "id": "n_sample_blocks",
                        "type": "int?",
                        "label": "Number of sample blocks",
                        "doc": "Number of blocks to divide samples into for parallel computation. Adjust depending on computer memory and number of samples in the analysis.",
                        "sbg:toolDefaultValue": "1",
                        "sbg:x": 49.20663070678711,
                        "sbg:y": -341.2738037109375
                    },
                    {
                        "id": "phenotype_file",
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Phenotype File",
                        "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by group.",
                        "sbg:x": 645.658203125,
                        "sbg:y": 169.88719177246094
                    },
                    {
                        "id": "kinship_plot_threshold",
                        "type": "float?",
                        "label": "Kinship plotting threshold",
                        "doc": "Minimum kinship for a pair to be included in the plot.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)"
                    },
                    {
                        "id": "group",
                        "type": "string?",
                        "label": "Group column name",
                        "doc": "Name of column in phenotype_file containing group variable (e.g., study) for plotting.",
                        "sbg:exposed": true
                    },
                    {
                        "id": "sparse_threshold",
                        "type": "float?",
                        "label": "Sparse threshold",
                        "doc": "Threshold for making the output kinship matrix sparse. A block diagonal matrix will be created such that any pair of samples with a kinship estimate greater than the threshold is in the same block; all pairwise estimates within a block are kept, and pairwise estimates between blocks are set to 0.",
                        "sbg:exposed": true,
                        "sbg:toolDefaultValue": "2^(-11/2) (~0.022, 4th degree)"
                    },
                    {
                        "id": "ibd_probs",
                        "type": "boolean?",
                        "label": "Return IBD probabilities?",
                        "doc": "Set this to FALSE to skip computing pairwise IBD probabilities (k0, k1, k2). If FALSE, the plottng step is also skipped, as it requires values for k0.",
                        "default": true,
                        "sbg:x": 210.16744995117188,
                        "sbg:y": 258.3512268066406
                    }
                ],
                "outputs": [
                    {
                        "id": "pcrelate_plots",
                        "outputSource": [
                            "kinship_plots/kinship_plots"
                        ],
                        "sbg:fileTypes": "PDF",
                        "type": "File[]?",
                        "label": "Kinship plots",
                        "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
                        "sbg:x": 1051.3870849609375,
                        "sbg:y": 166.4375457763672
                    },
                    {
                        "id": "pcrelate_output",
                        "outputSource": [
                            "pcrelate_correct/pcrelate_output"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "PC-Relate output file",
                        "doc": "PC-Relate output file with all samples",
                        "sbg:x": 931.4324951171875,
                        "sbg:y": -263.4626159667969
                    },
                    {
                        "id": "pcrelate_matrix",
                        "outputSource": [
                            "pcrelate_correct/pcrelate_matrix"
                        ],
                        "sbg:fileTypes": "RDATA",
                        "type": "File?",
                        "label": "Kinship matrix",
                        "doc": "A block diagonal matrix of pairwise kinship estimates with sparsity set by sparse_threshold. Samples are clustered into blocks of relatives based on `sparse_threshold`; all kinship estimates within a block are kept, and kinship estimates between blocks are set to 0. When `sparse_threshold` is 0, this is a dense matrix with all pairwise kinship estimates.",
                        "sbg:x": 927.3430786132812,
                        "sbg:y": -108.06597900390625
                    }
                ],
                "steps": [
                    {
                        "id": "pcrelate_beta",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "pca_file",
                                "source": "pca_file"
                            },
                            {
                                "id": "n_pcs",
                                "source": "n_pcs"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "variant_include_file"
                            },
                            {
                                "id": "variant_block_size",
                                "source": "variant_block_size"
                            }
                        ],
                        "out": [
                            {
                                "id": "beta"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/pcrelate-beta/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS File",
                                    "doc": "Input GDS file. It is recommended to use an LD pruned file with all chromosomes.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "id": "pca_file",
                                    "type": "File",
                                    "label": "PCA file",
                                    "doc": "RData file with PCA results from PC-AiR workflow; used to adjust for population structure.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "3",
                                    "id": "n_pcs",
                                    "type": "int?",
                                    "label": "Number of PCs",
                                    "doc": "Number of PCs (Principal Components) to use in adjusting for ancestry.",
                                    "default": 3
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant include file",
                                    "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1024",
                                    "id": "variant_block_size",
                                    "type": "int?",
                                    "label": "Variant block size",
                                    "doc": "Number of variants to read in a single block.",
                                    "default": 1024
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "beta",
                                    "doc": "RData file with ISAF beta values",
                                    "label": "ISAF beta values",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.RData"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "pcrelate_beta",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/pcrelate_beta.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "pcrelate_beta.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "pcrelate_beta.config",
                                            "entry": "${\n    var config = \"\"\n\n    if (inputs.gds_file) \n      config += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n      \n    if (inputs.pca_file) \n      config += \"pca_file \\\"\" + inputs.pca_file.path + \"\\\"\\n\"\n      \n    if (inputs.variant_include_file) \n      config += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n\n    if (inputs.out_prefix)\n      config += \"out_prefix \\\"\" + inputs.out_prefix  + \"_pcrelate_beta\\\"\\n\"\n  \t\n  \tif (inputs.n_pcs)\n      config += \"n_pcs \" + inputs.n_pcs + \"\\n\"\n    \n    if (inputs.sample_include_file) \n      config += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\"\n      \n  \tif (inputs.variant_block_size) \n      config += \"variant_block_size \" + inputs.variant_block_size + \"\\n\"\n      \n    return config\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "pcrelate_beta.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105851,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/pcrelate-beta/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105851,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105851,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a19bbabe8d172104f6906390276eaadcbcb8fa090c1d9a9485b8c286cb6ed7d20"
                        },
                        "label": "pcrelate_beta",
                        "sbg:x": 111.58096313476562,
                        "sbg:y": 18
                    },
                    {
                        "id": "sample_blocks_to_segments",
                        "in": [
                            {
                                "id": "n_sample_blocks",
                                "source": "n_sample_blocks"
                            }
                        ],
                        "out": [
                            {
                                "id": "segments"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/sample-blocks-to-segments/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "n_sample_blocks",
                                    "type": "int?",
                                    "label": "Number of sample blocks",
                                    "doc": "Number of blocks to divide samples into for parallel computation. Adjust depending on computer memory and number of samples in the analysis.",
                                    "default": 1
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "segments",
                                    "type": "int[]?",
                                    "outputBinding": {
                                        "outputEval": "${ \n    var blocks = [];\n    var N = inputs.n_sample_blocks;\n    for (var i = 1; i <= N; i++) {\n        for (var j = i; j <= N; j++) {\n            blocks.push([i, j]);\n        }\n    }\n    \n    var segments = [];\n    for (var i = 1; i <= blocks.length; i++)\n        segments.push(i)\n        \n    return segments;\n}"
                                    }
                                }
                            ],
                            "label": "sample blocks to segments",
                            "requirements": [
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105852,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/sample-blocks-to-segments/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105852,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105852,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a3feda07c7f070233aef62be7905b945e1ca9a2cd1feb1c472117026d555f5902"
                        },
                        "label": "sample blocks to segments",
                        "sbg:x": 195.43441772460938,
                        "sbg:y": -185.92086791992188
                    },
                    {
                        "id": "pcrelate",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": "gds_file"
                            },
                            {
                                "id": "pca_file",
                                "source": "pca_file"
                            },
                            {
                                "id": "beta_file",
                                "source": "pcrelate_beta/beta"
                            },
                            {
                                "id": "n_pcs",
                                "source": "n_pcs"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix"
                            },
                            {
                                "id": "variant_include_file",
                                "source": "variant_include_file"
                            },
                            {
                                "id": "variant_block_size",
                                "source": "variant_block_size"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "n_sample_blocks",
                                "source": "n_sample_blocks"
                            },
                            {
                                "id": "segment",
                                "source": "sample_blocks_to_segments/segments"
                            },
                            {
                                "id": "ibd_probs",
                                "source": "ibd_probs"
                            }
                        ],
                        "out": [
                            {
                                "id": "pcrelate"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/pcrelate/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS File",
                                    "doc": "Input GDS file. It is recommended to use an LD pruned file with all chromosomes.",
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "pca_file",
                                    "type": "File",
                                    "label": "PCA file",
                                    "doc": "RData file with PCA results from PC-AiR workflow; used to adjust for population structure.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "beta_file",
                                    "type": "File",
                                    "label": "ISAF beta values",
                                    "doc": "RData file with output from pcrelate_beta tool.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "3",
                                    "id": "n_pcs",
                                    "type": "int?",
                                    "label": "Number of PCs",
                                    "doc": "Number of PCs to use in adjusting for ancestry."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "pcrelate",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "variant_include_file",
                                    "type": "File?",
                                    "label": "Variant include file",
                                    "doc": "RData file with vector of variant.id to include. If not provided, all variants in the GDS file are included."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1024",
                                    "id": "variant_block_size",
                                    "type": "int?",
                                    "label": "Variant block size",
                                    "doc": "Number of variants to read in a single block.",
                                    "default": 1024
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include file",
                                    "doc": "RData file with vector of sample.id to include. If not provided, all samples in the GDS file are included.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "n_sample_blocks",
                                    "type": "int?",
                                    "label": "Number of sample blocks",
                                    "doc": "Number of blocks to divide samples into for parallel computation. Adjust depending on computer memory and number of samples in the analysis.",
                                    "default": 1
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "segment",
                                    "type": "int?",
                                    "label": "Sample block combination",
                                    "doc": "If number of sample blocks is > 1, run on this combination of sample blocks. Allowed values are 1:N where N is the number of possible combinations of sample blocks [i, j].",
                                    "default": 1
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "TRUE",
                                    "id": "ibd_probs",
                                    "type": "boolean?",
                                    "label": "Return IBD probabilities?",
                                    "doc": "Set this to FALSE to skip computing pairwise IBD probabilities (k0, k1, k2). If FALSE, the plottng step is also skipped, as it requires values for k0."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "pcrelate",
                                    "doc": "RData files with PC-Relate results for each sample block.",
                                    "label": "PC-Relate results",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.RData"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "pcrelate",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/pcrelate.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "pcrelate.config"
                                },
                                {
                                    "prefix": "--segment",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "${ return inputs.segment }"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "pcrelate.config",
                                            "entry": "${\n    var config = \"\"\n\n    if (inputs.gds_file) \n      config += \"gds_file \\\"\" + inputs.gds_file.path + \"\\\"\\n\"\n      \n    if (inputs.pca_file) \n      config += \"pca_file \\\"\" + inputs.pca_file.path + \"\\\"\\n\"\n            \n    if (inputs.beta_file) \n      config += \"beta_file \\\"\" + inputs.beta_file.path + \"\\\"\\n\"\n      \n    if (inputs.variant_include_file) \n      config += \"variant_include_file \\\"\" + inputs.variant_include_file.path + \"\\\"\\n\"\n\n    if (inputs.out_prefix)\n      config += \"out_prefix \\\"\" + inputs.out_prefix  + \"\\\"\\n\"\n  \t\n  \tif (inputs.n_pcs)\n      config += \"n_pcs \" + inputs.n_pcs + \"\\n\"\n    \n    if (inputs.sample_include_file) \n      config += \"sample_include_file \\\"\" + inputs.sample_include_file.path + \"\\\"\\n\"\n      \n  \tif (inputs.n_sample_blocks) \n      config += \"n_sample_blocks \" + inputs.n_sample_blocks + \"\\n\"\n      \n    if (!inputs.ibd_probs)\n      config += \"ibd_probs FALSE\" + \"\\n\"\n      \n    return config\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "pcrelate.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105853,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/pcrelate/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105853,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105853,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "ab16fa7e832289b8a4186790302020ada10ac708454fd6ac410f9fd80e3f81e25"
                        },
                        "label": "pcrelate",
                        "scatter": [
                            "segment"
                        ],
                        "sbg:x": 406.9648742675781,
                        "sbg:y": 4.164716720581055
                    },
                    {
                        "id": "pcrelate_correct",
                        "in": [
                            {
                                "id": "n_sample_blocks",
                                "source": "n_sample_blocks"
                            },
                            {
                                "id": "pcrelate_block_files",
                                "source": [
                                    "pcrelate/pcrelate"
                                ]
                            },
                            {
                                "id": "sparse_threshold",
                                "source": "sparse_threshold"
                            }
                        ],
                        "out": [
                            {
                                "id": "pcrelate_output"
                            },
                            {
                                "id": "pcrelate_matrix"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/pcrelate-correct/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "n_sample_blocks",
                                    "type": "int",
                                    "label": "Number of sample blocks",
                                    "doc": "Number of blocks to divide samples into for parallel computation. Adjust depending on computer memory and number of samples in the analysis.",
                                    "default": 1
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "pcrelate_block_files",
                                    "type": "File[]",
                                    "label": "PCRelate files for all sample blocks",
                                    "doc": "PCRelate files for all sample blocks",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "2^(-11/2) (~0.022, 4th degree)",
                                    "id": "sparse_threshold",
                                    "type": "float?",
                                    "label": "Sparse threshold",
                                    "doc": "Threshold for making the output kinship matrix sparse. A block diagonal matrix will be created such that any pair of samples with a kinship estimate greater than the threshold is in the same block; all pairwise estimates within a block are kept, and pairwise estimates between blocks are set to 0.",
                                    "default": 0.02209709
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "pcrelate_output",
                                    "doc": "PC-Relate output file with all samples",
                                    "label": "PC-Relate output file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*_pcrelate.RData"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "id": "pcrelate_matrix",
                                    "doc": "A block diagonal matrix of pairwise kinship estimates with sparsity set by sparse_threshold. Samples are clustered into blocks of relatives based on `sparse_threshold`; all kinship estimates within a block are kept, and kinship estimates between blocks are set to 0. When `sparse_threshold` is 0, this is a dense matrix with all pairwise kinship estimates.",
                                    "label": "Kinship matrix",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*_pcrelate_Matrix.RData"
                                    },
                                    "sbg:fileTypes": "RDATA"
                                }
                            ],
                            "label": "pcrelate_correct",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 3,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/pcrelate_correct.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "pcrelate_correct.config"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n    var cmd_line = \"\"\n    \n    for (var i=0; i<inputs.pcrelate_block_files.length; i++)\n        cmd_line += \"ln -s \" + inputs.pcrelate_block_files[i].path + \" \" + inputs.pcrelate_block_files[i].basename + \" && \"\n    \n    return cmd_line\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "R -q --vanilla"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "pcrelate_correct.config",
                                            "entry": "${\n    var config = \"\"\n\n    if (inputs.pcrelate_block_files) {\n        var prefix = inputs.pcrelate_block_files[0].nameroot.split(\"_block_\")[0]\n        config += \"pcrelate_prefix \\\"\" + prefix + \"\\\"\\n\"\n    }\n      \n  \tif (inputs.n_sample_blocks) \n        config += \"n_sample_blocks \" + inputs.n_sample_blocks + \"\\n\"\n      \n    if (inputs.sparse_threshold) \n        config += \"sparse_threshold \" + inputs.sparse_threshold + \"\\n\"  \n    \n    return config\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "pcrelate_correct.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105854,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/pcrelate-correct/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105854,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105854,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a1086737b1e69c4e0f5d69ae897150a558b603820eff7c298b547bfebeb236b83"
                        },
                        "label": "pcrelate_correct",
                        "sbg:x": 610.6256713867188,
                        "sbg:y": -144.12290954589844
                    },
                    {
                        "id": "kinship_plots",
                        "in": [
                            {
                                "id": "kinship_file",
                                "source": "pcrelate_correct/pcrelate_output"
                            },
                            {
                                "id": "kinship_method",
                                "default": "pcrelate"
                            },
                            {
                                "id": "kinship_plot_threshold",
                                "source": "kinship_plot_threshold"
                            },
                            {
                                "id": "phenotype_file",
                                "source": "phenotype_file"
                            },
                            {
                                "id": "group",
                                "source": "group"
                            },
                            {
                                "id": "sample_include_file",
                                "source": "sample_include_file"
                            },
                            {
                                "id": "out_prefix",
                                "source": "out_prefix",
                                "valueFrom": "${ return inputs.out_prefix + \"_pcrelate\" }"
                            },
                            {
                                "id": "run_plots",
                                "source": "ibd_probs"
                            }
                        ],
                        "out": [
                            {
                                "id": "kinship_plots"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/kinship-plots/0",
                            "baseCommand": [
                                "R -q --vanilla"
                            ],
                            "inputs": [
                                {
                                    "sbg:category": "Input",
                                    "id": "kinship_file",
                                    "type": "File",
                                    "label": "Kinship File",
                                    "doc": "Kinship file",
                                    "sbg:fileTypes": "RDATA, SEG, KIN, GDS"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "kinship_method",
                                    "type": {
                                        "type": "enum",
                                        "symbols": [
                                            "king_ibdseg",
                                            "pcrelate",
                                            "king_robust"
                                        ],
                                        "name": "kinship_method"
                                    },
                                    "label": "Kinship method",
                                    "doc": "Method used to generate kinship estimates."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "2^(-9/2) (third-degree relatives and closer)",
                                    "id": "kinship_plot_threshold",
                                    "type": "float?",
                                    "label": "Kinship plotting threshold",
                                    "doc": "Minimum kinship for a pair to be included in the plot."
                                },
                                {
                                    "sbg:category": "Input Files",
                                    "id": "phenotype_file",
                                    "type": "File?",
                                    "label": "Phenotype File",
                                    "doc": "RData file with data.frame or AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by group.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "NA",
                                    "id": "group",
                                    "type": "string?",
                                    "label": "Group column name",
                                    "doc": "Name of column in phenotype_file containing group variable (e.g., study) for plotting."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "id": "sample_include_file",
                                    "type": "File?",
                                    "label": "Sample Include File",
                                    "doc": "RData file with vector of sample.id to include.",
                                    "sbg:fileTypes": "RDATA"
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "kinship",
                                    "id": "out_prefix",
                                    "type": "string?",
                                    "label": "Output prefix",
                                    "doc": "Prefix for output files."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "kinship_plots",
                                    "doc": "Hexbin plots of estimated kinship coefficients vs. IBS0. If \"group\" is provided, additional plots will be generated within each group and across groups.",
                                    "label": "Kinship plots",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*.pdf"
                                    },
                                    "sbg:fileTypes": "PDF"
                                }
                            ],
                            "label": "kinship_plots",
                            "arguments": [
                                {
                                    "prefix": "<",
                                    "shellQuote": false,
                                    "position": 2,
                                    "valueFrom": "/usr/local/analysis_pipeline/R/kinship_plots.R"
                                },
                                {
                                    "prefix": "--args",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "kinship_plots.config"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.10.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "kinship_plots.config",
                                            "entry": "${\nvar cmd_line = \"\"\n\nif(inputs.kinship_file)\n    cmd_line += \"kinship_file \\\"\" + inputs.kinship_file.path + \"\\\"\\n\"\n\nif(inputs.kinship_method){\n    if(inputs.kinship_method == \"king_robust\"){\n        cmd_line +='kinship_method \"king\"\\n'\n    } else {\n        cmd_line +='kinship_method \"' + inputs.kinship_method + '\"\\n'\n    }\n}\n\nif(inputs.kinship_plot_threshold)\n    cmd_line +='kinship_threshold \"' + inputs.kinship_plot_threshold + '\\n'\n\nif(inputs.out_prefix) {\n    cmd_line += 'out_file_all \"' + inputs.out_prefix + '_all.pdf\"\\n'\n    cmd_line += 'out_file_cross \"' + inputs.out_prefix + '_cross_group.pdf\"\\n'\n    cmd_line += 'out_file_study \"' + inputs.out_prefix + '_within_group.pdf\"\\n'\n}\n\nif(inputs.phenotype_file)\n    cmd_line += 'phenotype_file \"' + inputs.phenotype_file.path + '\"\\n'\n\nif(inputs.group)\n    cmd_line += 'study \"' + inputs.group + '\"\\n'\n\nif(inputs.sample_include_file)\n    cmd_line += 'sample_include_file \"' + inputs.sample_include_file.path + '\"\\n'\n\nreturn cmd_line\n}",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                },
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "kinship_plots.config"
                                }
                            ],
                            "stdout": "job.out.log",
                            "sbg:image_url": null,
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105837,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/kinship-plots/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105837,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105837,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "ac07d7f46f1402f82f9271881eca181d1a4e44613e62de62f9242aff85413f411"
                        },
                        "label": "kinship_plots",
                        "when": "$(inputs.run_plots)",
                        "sbg:x": 801.525146484375,
                        "sbg:y": 40.709495544433594
                    }
                ],
                "requirements": [
                    {
                        "class": "ScatterFeatureRequirement"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "Ancestry and Relatedness"
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/pc-relate/0.png",
                "sbg:original_source": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/smgogarten/genesis-relatedness/pc-relate/17/raw/",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105855,
                        "sbg:revisionNotes": "Workflow decomposed"
                    }
                ],
                "sbg:toolkit": "UW-GAC Ancestry and Relatedness",
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/pc-relate/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": "Workflow decomposed",
                "sbg:modifiedOn": 1637105855,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105855,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "abfc19052511dbf4f1b8f8a60d7ff2652a3562350fc2e7850a207d4bc8b26ba5e"
            },
            "label": "PC-Relate",
            "sbg:x": 818.5,
            "sbg:y": -504
        },
        {
            "id": "genesis_locuszoom",
            "in": [
                {
                    "id": "out_prefix",
                    "source": "out_prefix_1"
                },
                {
                    "id": "in_assoc_files",
                    "source": [
                        "single_variant_association_testing/assoc_combined"
                    ]
                },
                {
                    "id": "in_loci_file",
                    "source": "sbg_merge_and_filter_genesis_results_cwl2/out_lz"
                },
                {
                    "id": "locus_type",
                    "default": "variant"
                },
                {
                    "id": "in_gds_files",
                    "linkMerge": "merge_flattened"
                },
                {
                    "id": "genome_build",
                    "default": "hg19"
                },
                {
                    "id": "track_file_type",
                    "default": "window"
                },
                {
                    "id": "track_label",
                    "default": "SKATO"
                },
                {
                    "id": "track_threshold",
                    "default": 0.000001,
                    "source": "track_threshold"
                },
                {
                    "id": "database_directory",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "in_database_compressed",
                    "source": "in_database_compressed"
                }
            ],
            "out": [
                {
                    "id": "out_pdf_reports"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/genesis-locuszoom/20",
                "baseCommand": [],
                "inputs": [
                    {
                        "sbg:category": "Configuration",
                        "id": "out_prefix",
                        "type": "string?",
                        "label": "Output prefix",
                        "doc": "Prefix for files created by this script."
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "in_assoc_files",
                        "type": "File[]",
                        "label": "Association results files",
                        "doc": "Files with single-variant association test results. File has to follow same naming specification as GENESIS Association Test files (basename_chr##.RData).",
                        "sbg:fileTypes": "RDATA"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "in_loci_file",
                        "type": "File",
                        "label": "Locus file",
                        "doc": "Text file with columns chr, pop and either variantID (for locus_type=variant) or start, end (for locus_type=region)",
                        "sbg:fileTypes": "TXT, BED"
                    },
                    {
                        "sbg:toolDefaultValue": "variant",
                        "sbg:category": "Configuration",
                        "id": "locus_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "variant",
                                    "region"
                                ],
                                "name": "locus_type"
                            }
                        ],
                        "label": "Locus type",
                        "doc": "Type of region to plot (variant with flanking region, or region)"
                    },
                    {
                        "sbg:toolDefaultValue": "500",
                        "sbg:category": "Configuration",
                        "id": "flanking_region",
                        "type": "int?",
                        "label": "Flanking region",
                        "doc": "Flanking region in kb. Default is 500kb."
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "in_gds_files",
                        "type": "File[]?",
                        "label": "GDS files",
                        "doc": "GDS file to use for calculating LD",
                        "sbg:fileTypes": "GDS"
                    },
                    {
                        "sbg:toolDefaultValue": "hg38",
                        "sbg:category": "Configuration",
                        "id": "genome_build",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "hg19",
                                    "hg38"
                                ],
                                "name": "genome_build"
                            }
                        ],
                        "label": "Genome build",
                        "doc": "Genome build (hg19 or hg38). Default: hg38",
                        "default": "hg38"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "ld_sample_include",
                        "type": "File?",
                        "label": "LD sample include",
                        "doc": "RData file with vector of sample.id to include when calculating LD.",
                        "sbg:fileTypes": "RDATA"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "in_track_files",
                        "type": "File[]?",
                        "label": "Track files",
                        "doc": "File with aggregate or window association test results. Regions will be displayed in a track in the LocusZoom plot. Include a space to insert chromosome.",
                        "sbg:fileTypes": "RDATA"
                    },
                    {
                        "sbg:toolDefaultValue": "window",
                        "sbg:category": "Configuration",
                        "id": "track_file_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "window",
                                    "aggregate"
                                ],
                                "name": "track_file_type"
                            }
                        ],
                        "label": "Track file type",
                        "doc": "Type of association regions in track_file (window or aggregate)."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "track_label",
                        "type": "string?",
                        "label": "Track label",
                        "doc": "Label to display to the right of the track in the plot."
                    },
                    {
                        "sbg:toolDefaultValue": "5e-8",
                        "sbg:category": "Configuration",
                        "id": "track_threshold",
                        "type": "float?",
                        "label": "Track threshold",
                        "doc": "P-value threshold for selecting regions to display."
                    },
                    {
                        "loadListing": "deep_listing",
                        "sbg:category": "General",
                        "id": "database_directory",
                        "type": "Directory?",
                        "label": "Database directory",
                        "doc": "Directory containing databases used for LD calculation and annotation. It can be obtained by decompressing LD_Database.tar.gz file found in the public files gallery, or can be downloaded from LocusZoom Standalone github page."
                    },
                    {
                        "sbg:toolDefaultValue": "1",
                        "id": "segment",
                        "type": "int?",
                        "label": "Segment",
                        "doc": "If provided, only a single segment from Locus File will be used. Otherwise, one PDF file will be provided per line in the Locus File."
                    },
                    {
                        "sbg:category": "General",
                        "sbg:toolDefaultValue": "4",
                        "id": "threads",
                        "type": "int?",
                        "label": "Number of CPUs",
                        "doc": "Number of plots to create in parallel. Must not be larger than number of threads on the instance. Default: 4"
                    },
                    {
                        "sbg:toolDefaultValue": "5e-8",
                        "sbg:category": "Configuration",
                        "id": "signif_line",
                        "type": "float?",
                        "label": "Significance line",
                        "doc": "Where to draw significance line on plots. Default: 5e-8"
                    },
                    {
                        "sbg:category": "General",
                        "id": "in_database_compressed",
                        "type": "File?",
                        "label": "Database bundle",
                        "doc": "Compressed database directory. Available as LZ_Database.tar.gz in the public file gallery. Required if Database Directory is not provided.",
                        "sbg:fileTypes": "TAR.GZ"
                    }
                ],
                "outputs": [
                    {
                        "id": "out_pdf_reports",
                        "doc": "One LZ plot per locus.",
                        "label": "Locuszoom plots",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "*.pdf",
                            "outputEval": "$(inheritMetadata(self, inputs.assoc_file))"
                        },
                        "sbg:fileTypes": "PDF"
                    }
                ],
                "doc": "**LocusZoom for GENESIS** visualizes association testing results using the LocusZoom standalone software. This App is a wrapper around LocusZoom standalone software to enable it to work with outputs of GENESIS association pipelines [1]. Main goal of this App is to visualize results of **GENESIS Single Variant Association Test**, however regions from sliding window or aggregate tests with p-values below a certain threshold can be displayed in a separate track. A list of all inputs and parameters with corresponding descriptions can be found at the bottom of this page.\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common Use Cases \n**LocusZoom for GENESIS** is used to make annotated Manhattan plots on specific regions from association files generated by GENESIS workflows. **LocusZoom for GENESIS** can work with outputs from **GENESIS Single Variant Association Testing**, **GENESIS Aggregate Association Testing** and **GENESIS Sliding Window Association testing** workflows. **Association results file** originating from **GENESIS Single Variant Association Testing** is required even when results from Sliding Window or Aggregate tests are visualized. Sliding Window or Aggregate results are optional inputs provided via the **Track file** input. One PDF file is generated for each locus defined in the **Loci file**.\n\nLoci to plot are specified in the **Loci file**, with chromosome *chr* and either *variantID* (to specify the reference variant) or *start end* (to indicate a region to plot, in which case the variant with the smallest p-value will be the reference). Population (pop) is either TOPMED or one of the 1000 Genomes populations (hg19:AFR, AMR, ASN, EUR; hg38: AFR, AMR, EUR, EAS, SAS). If pop = TOPMED, LD is computed from the TOPMed data using the sample set in the **LD Sample Include** file. Example of region - defined **Loci file**:\n\n| chr | start     | end       | pop |\n|-----|-----------|-----------|-----|\n| 1   | 69423922  | 70423922  | AFR |\n| 1   | 94393630  | 95393630  | AMR |\n| 1   | 193326139 | 194326139 | EUR |\n| 2   | 2009400   | 2009430   | EUR |\n\nExample of variant - defined **Loci file**:\n\n| variant.id | chr | pop    |\n|------------|-----|--------|\n| 339        | 1   | AFR    |\n| 831        | 1   | topmed |\n\nFor visualization, LD calculation and additional information, **LocusZoom for GENESIS** uses several databases. Due to the size of these databases, they cannot be included in the Docker image, but have to be provided every time this App is run. There are two ways to provide these database to the **GENESIS LocusZoom** App:\n1. Via **Database bundle** input: \n    1. Copy *LZ_Database.tar.gz* file from the Public Files Gallery.\n    2. Use the copied file as **Database bundle** input.\n2. Via **Database directory** input:\n    1. Copy \"*LZ_Database.tar.gz*\" from the Public Files Gallery.\n    2. Copy the **SBG Decompressor CWL1.0** App from the Public Apps Gallery into the current project.\n    3. Run **SBG Decompressor CWL1.0** with *LZ_Database.tar.gz* as an input file, with **Flatten outputs** set to **False**.\n    4. Resulting directory should be named *data*. If it is not, rename it so.\n    5. Use the *data* directory as the **Database directory** input\n\nBoth approaches come with pros and cons. First approach is very simple and does not lead to additional storage cost as copies of files from the Public Apps Gallery are not billed. However, the resulting task will be longer and more expensive. If **LocusZoom for GENESIS** is run multiple times in the span of several days, it is recommended to use the second approach, and delete the *data* directory when it is no longer used. Time difference between the two approaches is around 10 minutes and around 0.10$ per task (directory approach is cheaper), while the cost of storing the uncompressed directory is around 0.06$/day. \n\n### Changes Introduced by Seven Bridges\n\n * Multiple loci are run in parallel if the **segment** argument is not provided. If the **segment** argument is provided, only the *n*-th segment from the **Loci file** will be plotted. Number of plots generated in parallel is determined using the **Number of CPUs** argument.\n\n### Common Issues and Important Notes\n\n * Either **Database bundle** or **Database directory** must be provided.\n * All **Association result files**, **GDS files** and **Track files** must follow the same naming standard as in other **GENESIS** workflows, meaning that these files must be provided as one file per chromosome, with the same name except for the number of chromosome, and part of the filename must include *_chr##.* substring, where *##* is replaced by the chromosome number (1-22, X, Y). Note that basename must be identical for all the files in these three input groups, but only within the group, meaning that **GDS files** can have a different basename from **Association result files** and **Track files**.\n * There must be an **Association result file** for each unique chromosome in the **Loci file**\n * **Loci file** must contain exclusively regions or variants, and the type of loci in the file must match the type defined by the **Locus type** app setting.\n * **Loci file** must follow its file specification:\n1. Values and header should be space-delimited.\n2. Header must contain [chr, start, end, population] in case of **Locus type**=*region*, or [variant.id, chr, population] in case **Locus type**=*variant*.\n3. Population must be one of [*AFR, AMR, ASN, EUR, TOPMED*] if **Genome build**=*hg19*, or [*AFR, AMR, EUR, EAS, SAS, TOPMED*] if **Genome build**=*hg38*. Population is not case sensitive.\n * If **Track files** are provided, **Track file type** setting must be set to the appropriate value\n * If **Track files** are provided, **Track threshold** value must be higher than the value of lowest p-value of any variant in all regions. If there is a region without any variants with p-value lower than **Track threshold** the task will fail! It is recommended to set a high value for **Track threshold**, for example *0.1* if it is not known which values to expect. However, some knowledge of variants and regions plotted and setting an appropriate threshold value is recommended.\n * **Number of CPUs** setting must not be larger than the number of threads on the instance. However, it is advised to leave at least one idle thread as additional jobs unrelated to plotting might use it. \n\n### Performance Benchmarking \nUploading and unpacking files and databases is a fixed cost for running **LocusZoom for GENESIS** and takes around 10 minutes when working with **Database directory**, or 20 minutes when working with **Database bundle**. The process of generating plots takes minutes and can impact final task execution time and cost only if a large number of plots is generated (more than 100). Default instance cost is 0.48$/h and can be lowered by using [spot instances](https://docs.sevenbridges.com/docs/about-spot-instances). \n\n### References\n[1] [LocusZoom Standalone](https://github.com/UW-GAC/locuszoom-standalone)",
                "label": "LocusZoom for GENESIS",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    if(inputs.database_directory)\n    {\n        return \"mkdir /usr/local/src/locuszoom-standalone/data && mkdir /usr/local/src/locuszoom-standalone/data/database && export PATH=/usr/local/locuszoom-standalone/bin:$PATH && cp -r \" + inputs.database_directory.path + \"/. /usr/local/src/locuszoom-standalone/data/ && \"\n    }\n    else\n    {\n        if(inputs.in_database_compressed)\n            return \"mkdir /usr/local/src/locuszoom-standalone/data && mkdir /usr/local/src/locuszoom-standalone/data/database && export PATH=/usr/local/locuszoom-standalone/bin:$PATH && tar -xzf \" + inputs.in_database_compressed.path + \" -C /usr/local/src/locuszoom-standalone/ && \"\n        else\n        {\n            var error = \"Database not provided! Please provide LZ database either as Database bundle or Database directory.\";\n            throw error;\n        }\n            \n    }\n    \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 10,
                        "valueFrom": "${\n    if (inputs.segment)\n        return \"Rscript locuszoom.R locuszoom2.config --segment \" + inputs.segment\n    else\n    {\n        var num_of_threads = 4;\n        if(inputs.threads)\n            num_of_threads = inputs.threads\n        return \"seq 1 $(expr $(wc -l \" + inputs.in_loci_file.path + \"| awk '{ print $1 }') - 1) | xargs -P\" + num_of_threads + \" -n1 -t -I % Rscript locuszoom.R locuszoom2.config --segment %\"\n    }\n        \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 5,
                        "valueFrom": "${\n    var command = '';\n    for(var i=0; i<inputs.in_assoc_files.length; i++)\n        command += \"ln -s \" + inputs.in_assoc_files[i].path + \" \" + inputs.in_assoc_files[i].path.split(\"/\").pop() + \" && \"\n    if(inputs.in_gds_files)\n        for(var i=0; i<inputs.in_gds_files.length; i++)\n            command += \"ln -s \" + inputs.in_gds_files[i].path + \" \" + inputs.in_gds_files[i].path.split(\"/\").pop() + \" && \"\n    if(inputs.in_track_files)\n        for(var i=0; i<inputs.in_track_files.length; i++)\n            command += \"ln -s \" + inputs.in_track_files[i].path + \" \" + inputs.in_track_files[i].path.split(\"/\").pop() + \" && \"\n    \n    return command\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 4000,
                        "coresMin": "${\n    if(inputs.threads)\n        return inputs.threads\n    else\n        return 1\n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "uwgac/topmed-master:2.8.1"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "locuszoom2.config",
                                "entry": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    \n    var argument = [];\n    var a_file = [].concat(inputs.in_assoc_files)[0];\n    var chr = find_chromosome(a_file.basename);\n    var path = a_file.path.split('chr' + chr);\n    var extension = path[1].split('.')[1];\n    \n    argument.push('assoc_file ' + '\"' + path[0].split('/').pop() + 'chr ' +path[1] + '\"');\n    \n    if(inputs.locus_type)\n        argument.push('locus_type \"' + inputs.locus_type + '\"')\n    if(inputs.flanking_region)\n        argument.push('flanking_region \"' + inputs.flanking_region + '\"')\n    if(inputs.genome_build)\n        argument.push('genome_build \"' + inputs.genome_build + '\"')\n    if(inputs.in_gds_files)\n    {\n        var g_file = [].concat(inputs.in_gds_files)[0];\n        var gchr = find_chromosome(g_file.basename);\n        var gpath = g_file.path.split('chr' + gchr);\n        argument.push('gds_file ' + '\"' + gpath[0].split('/').pop() + 'chr ' +gpath[1] + '\"');\n    }\n    if(inputs.ld_sample_include)\n        argument.push('ld_sample_include \"' + inputs.ld_sample_include.path + '\"')\n    if(inputs.in_loci_file)\n        argument.push('locus_file \"' + inputs.in_loci_file.path + '\"')\n        \n    if(inputs.in_track_files)\n    {\n        var t_file = [].concat(inputs.in_track_files)[0];\n        var tchr = find_chromosome(t_file.basename);\n        var tpath = t_file.path.split('chr' + tchr);\n        argument.push('track_file ' + '\"' + tpath[0].split('/').pop() + 'chr ' +tpath[1] + '\"');\n    }\n    if(inputs.track_file_type)\n        argument.push('track_file_type \"' + inputs.track_file_type + '\"')\n    if(inputs.track_label)\n        argument.push('track_label \"' + inputs.track_label + '\"')\n    if(inputs.track_threshold)\n        argument.push('track_threshold \"' + inputs.track_threshold + '\"')\n    if(inputs.signif_line)\n        argument.push('signif_line \"' + inputs.signif_line + '\"')\n    if(inputs.out_prefix)\n        argument.push('out_prefix \"' + inputs.out_prefix + '\"')\n    \n    argument.push('\\n');\n    return argument.join('\\n')\n    \n}\n",
                                "writable": false
                            },
                            {
                                "entryname": "locuszoom.R",
                                "entry": "library(argparser)\nlibrary(TopmedPipeline)\nlibrary(SeqVarTools)\nlibrary(dplyr)\nsessionInfo()\n\nargp <- arg_parser(\"LocusZoom plots\")\nargp <- add_argument(argp, \"config\", help=\"path to config file\")\nargp <- add_argument(argp, \"--segment\", help=\"row in locus_file to plot\", default=1, type=\"integer\")\nargp <- add_argument(argp, \"--version\", help=\"pipeline version number\")\nargv <- parse_args(argp)\ncat(\">>> TopmedPipeline version \", argv$version, \"\\n\")\nconfig <- readConfig(argv$config)\nsegment <- argv$segment\n\nrequired <- c(\"assoc_file\",\n              \"locus_file\")\noptional <- c(\"flanking_region\"=500,\n              \"gds_file\"=NA,\n              \"genome_build\"=\"hg38\",\n              \"ld_sample_include\"=NA,\n              \"locus_type\"=\"variant\",\n              \"out_prefix\"=\"locuszoom\",\n              \"signif_line\"=5e-8,\n              \"track_file\"=NA,\n              \"track_file_type\"=\"window\",\n              \"track_label\"=\"\",\n              \"track_threshold\"=5e-8)\nconfig <- setConfigDefaults(config, required, optional)\nprint(config)\n\nstopifnot(config[\"locus_type\"] %in% c(\"variant\", \"region\"))\n\n# read selected locus\nlocus <- read.table(config[\"locus_file\"], header=TRUE, as.is=TRUE)[segment,]\nstopifnot(all(c(\"chr\", \"pop\") %in% names(locus)))\nprint(locus)\n\n# population for LD\npop <- toupper(locus$pop)\nstopifnot(pop %in% c(\"TOPMED\", \"AFR\", \"AMR\", \"ASN\", \"EUR\", \"EAS\", \"SAS\"))\n\n## get association test results\nvar.chr <- locus$chr\nassocfile <- insertChromString(config[\"assoc_file\"], var.chr)\nassoc <- getobj(assocfile)\n\nif (config[\"locus_type\"] == \"variant\") {\n    stopifnot(\"variant.id\" %in% names(locus))\n    variant <- locus$variant.id\n    flank <- as.numeric(config[\"flanking_region\"]) * 1000\n    var.pos <- assoc$pos[assoc$variant.id == variant]\n    start <- var.pos - flank\n    if (start < 1) start <- 1\n    end <- var.pos + flank\n    \n    lz.name <- paste0(\"chr\", var.chr, \":\", var.pos)\n    ld.region <- paste0(\"--refsnp \\\"\", lz.name, \"\\\"\", \" --flank \", config[\"flanking_region\"], \"kb\")\n    prefix <- paste0(config[\"out_prefix\"], \"_var\", variant, \"_ld_\", pop)\n    freq <- assoc$freq[assoc$variant.id == variant]\n    maf <- min(freq, 1-freq)\n    mac <- assoc$MAC[assoc$variant.id == variant]\n    title <- paste(lz.name, \"- MAF:\", formatC(maf, digits=3), \"- MAC:\", mac)\n    \n} else if (config[\"locus_type\"] == \"region\") {\n    stopifnot(all(c(\"start\", \"end\") %in% names(locus)))\n    start <- locus$start\n    end <- locus$end\n\n    ld.region <- paste(\"--chr\", var.chr, \"--start\", start, \"--end\", end)\n    prefix <- paste0(config[\"out_prefix\"], \"_ld_\", pop)\n    title <- \"\"\n}\n\n## construct METAL-format file\nassoc <- assoc %>%\n    filter(chr == var.chr, pos > start, pos < end) %>%\n    select(variant.id, chr, pos, ends_with(\"pval\"))\nnames(assoc)[4] <- \"pval\"\n\n##### NOTE: i added the following two lines of code #####\n## remove duplicate chr:pos rows, removing the variant with the less significant (higher) pvalue\nassoc <- assoc[order(assoc$pval,decreasing=FALSE),]\nassoc <- assoc[!duplicated(assoc$pos),]\nassoc <- assoc[order(assoc$pos),]\n#####\n\nassoc.filename <- tempfile()\nwriteMETAL(assoc, file=assoc.filename)\n\n# LD\nif (pop != \"TOPMED\") {\n    ld.cmd <- paste(\"--pop\", pop, \"--source 1000G_Nov2014\")\n    ld.title <- paste(\"LD: 1000G\", pop)\n} else {\n    if (!is.na(config[\"ld_sample_include\"])) {\n        sample.id <- getobj(config[\"ld_sample_include\"])\n    } else {\n        sample.id <- NULL\n    }\n    if (config[\"locus_type\"] == \"variant\") {\n        ref.var <- variant\n    } else {\n        ref.var <- filter(assoc, pval == min(pval))$variant.id\n        if (length(ref.var) > 1) {\n            message(\"Multiple variants with minimum pval; selecting the first one as reference\")\n            ref.var <- ref.var[1]\n        }\n    }\n    gdsfile <- insertChromString(config[\"gds_file\"], var.chr)\n    ld <- calculateLD(gdsfile, variant.id=assoc$variant.id, ref.var=ref.var, sample.id=sample.id)\n    ld.filename <- tempfile()\n    writeLD(assoc, ld, ref.var, file=ld.filename)\n\n    ld.cmd <- paste(\"--ld\", ld.filename)\n    ld.title <- \"LD: TOPMed\"\n}\ntitle <- if (title == \"\") ld.title else paste(ld.title, title, sep=\" - \")\n\n## construct BED track file\nif (!is.na(config[\"track_file\"])) {\n    trackfile <- insertChromString(config[\"track_file\"], var.chr)\n    track <- getAssoc(trackfile, config[\"track_file_type\"]) %>%\n        filter(pval < as.numeric(config[\"track_threshold\"]))\n    track.filename <- tempfile()\n    writeBED(track, file=track.filename, track.label=config[\"track_label\"])\n    track.cmd <- paste(\"--bed-tracks\", track.filename)\n} else {\n    track.cmd <- \"\"\n}\n\nsignif <- as.numeric(config[\"signif_line\"])\n\ncommand <- paste(\"locuszoom\",\n                 \"theme=publication\",\n                 \"--cache None\",\n                 \"--no-date\",\n                 \"--plotonly\",\n                 \"--gene-table gencode\",\n                 \"--build\", config[\"genome_build\"],\n                 \"--chr\", var.chr,\n                 \"--metal\", assoc.filename,\n                 track.cmd,\n                 ld.cmd,\n                 ld.region,\n                 \"--prefix \", prefix,\n                 paste0(\"title=\\\"\", title, \"\\\"\"),\n                 \"width=10\", \n                 \"height=11\",\n                 \"signifLine=\\\"5,7.3\\\"\",\n                 \"signifLineColor=\\\"orange,orange\\\"\",\n                 \"signifLineWidth=\\\"2,3\\\"\",\n                 \"rfrows=9\", \n                 \"ldColors=\\\"gray50, steelblue4, steelblue1, palegreen2, chocolate1, brown1, darkorchid2\\\"\",\n                 #paste0(\"signifLine=\\\"\", -log10(signif), \"\\\" signifLineColor=\\\"gray\\\" signifLineWidth=\\\"2\\\"\"),\n                 \"ylab=\\\"-log10(p-value) from single variant test\\\"\")\n\ncat(paste(command, \"\\n\"))\nsystem(command)\n\nunlink(assoc.filename)\nif (exists(\"track.filename\")) unlink(track.filename)\nif (exists(\"ld.filename\")) unlink(ld.filename)",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:AWSInstanceType",
                        "value": "c5.2xlarge;ebs-gp2;512"
                    },
                    {
                        "class": "sbg:GoogleInstanceType",
                        "value": "n1-standard-8;pd-ssd;512"
                    }
                ],
                "sbg:image_url": null,
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105869,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647079075,
                        "sbg:revisionNotes": "MZ edit"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647085849,
                        "sbg:revisionNotes": "mz edit2"
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647087365,
                        "sbg:revisionNotes": "mz edit3"
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647088983,
                        "sbg:revisionNotes": "mz edit4"
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647090584,
                        "sbg:revisionNotes": "mz update5"
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647092406,
                        "sbg:revisionNotes": "revert to original script"
                    },
                    {
                        "sbg:revision": 7,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647108067,
                        "sbg:revisionNotes": "mz edit11"
                    },
                    {
                        "sbg:revision": 8,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647162181,
                        "sbg:revisionNotes": "color palette edit"
                    },
                    {
                        "sbg:revision": 9,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647167732,
                        "sbg:revisionNotes": "hex colors in uppercase"
                    },
                    {
                        "sbg:revision": 10,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647174658,
                        "sbg:revisionNotes": "mz edit (significance lines)"
                    },
                    {
                        "sbg:revision": 11,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647178912,
                        "sbg:revisionNotes": "mz update"
                    },
                    {
                        "sbg:revision": 12,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647284524,
                        "sbg:revisionNotes": "black and white added to LD colors"
                    },
                    {
                        "sbg:revision": 13,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647287055,
                        "sbg:revisionNotes": "colors updated"
                    },
                    {
                        "sbg:revision": 14,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647287177,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 15,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648482209,
                        "sbg:revisionNotes": "rfrows=6"
                    },
                    {
                        "sbg:revision": 16,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648492404,
                        "sbg:revisionNotes": "8 rows"
                    },
                    {
                        "sbg:revision": 17,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648761045,
                        "sbg:revisionNotes": "15 11 to 10 7"
                    },
                    {
                        "sbg:revision": 18,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648762908,
                        "sbg:revisionNotes": "10 14"
                    },
                    {
                        "sbg:revision": 19,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648765966,
                        "sbg:revisionNotes": "10 13"
                    },
                    {
                        "sbg:revision": 20,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648766004,
                        "sbg:revisionNotes": "10 11"
                    }
                ],
                "sbg:toolkit": "LocusZoom Standalone",
                "sbg:links": [
                    {
                        "id": "https://genome.sph.umich.edu/wiki/LocusZoom_Standalone",
                        "label": "Wiki"
                    },
                    {
                        "id": "https://github.com/statgen/locuszoom-standalone",
                        "label": "github"
                    },
                    {
                        "id": "https://genome.sph.umich.edu/wiki/LocusZoom_Standalone#License",
                        "label": "License"
                    }
                ],
                "sbg:categories": [
                    "Plotting and Rendering"
                ],
                "sbg:license": "GNU General Public License v3",
                "sbg:wrapperAuthor": "SBG",
                "sbg:toolAuthor": "TopMed DCC",
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/genesis-locuszoom/20",
                "sbg:revision": 20,
                "sbg:revisionNotes": "10 11",
                "sbg:modifiedOn": 1648766004,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105869,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 20,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a79f28b327ead573a5554158074caeaf4e60915e993c00434012af1c5525c40fe",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "LocusZoom for GENESIS",
            "sbg:x": 2170.1318359375,
            "sbg:y": -497.8125305175781
        },
        {
            "id": "bcftools_annotate_cwl1",
            "in": [
                {
                    "id": "input_file",
                    "source": [
                        "bcftools_view_1_10_1/out_variants"
                    ]
                },
                {
                    "id": "rename_chrs",
                    "source": "rename_chrs"
                }
            ],
            "out": [
                {
                    "id": "output_file"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/bcftools-annotate-cwl1/0",
                "baseCommand": [],
                "inputs": [
                    {
                        "sbg:category": "File Input",
                        "id": "input_file",
                        "type": "File[]",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 43,
                            "valueFrom": "${\n    var files_array = [].concat(inputs.input_file)\n    var fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() == 'gz') {\n        fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n        return fname + \".gz\"\n    } else {\n        return fname + \".gz\"\n    }\n}"
                        },
                        "label": "Input file",
                        "doc": "Input file.",
                        "sbg:fileTypes": "VCF, VCF.GZ, BED"
                    },
                    {
                        "sbg:altPrefix": "-i",
                        "sbg:category": "Configuration",
                        "id": "include_expression",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--include",
                            "shellQuote": false,
                            "position": 16
                        },
                        "label": "include expression",
                        "doc": "Include only sites for which the expression is true."
                    },
                    {
                        "sbg:altPrefix": "-e",
                        "sbg:category": "Configuration",
                        "id": "exclude_expression",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--exclude",
                            "shellQuote": false,
                            "position": 13
                        },
                        "label": "Exclude expression",
                        "doc": "Exclude sites for which the expression is true."
                    },
                    {
                        "sbg:altPrefix": "-m",
                        "sbg:category": "Configuration",
                        "id": "mark_sites",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--mark-sites",
                            "shellQuote": false,
                            "position": 17
                        },
                        "label": "Filter mode",
                        "doc": "Annotate sites which are present (\"+\") or absent (\"-\") in the -a file with a new INFO/TAG flag."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "output_name",
                        "type": "string?",
                        "label": "Output file name",
                        "doc": "Output file name."
                    },
                    {
                        "sbg:toolDefaultValue": "Uncompressed VCF",
                        "sbg:altPrefix": "-O",
                        "sbg:category": "Configuration",
                        "id": "output_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "CompressedBCF",
                                    "UncompressedBCF",
                                    "CompressedVCF",
                                    "UncompressedVCF"
                                ],
                                "name": "output_type"
                            }
                        ],
                        "label": "Output type",
                        "doc": "b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]."
                    },
                    {
                        "sbg:altPrefix": "-r",
                        "sbg:category": "Configuration",
                        "id": "regions",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--regions",
                            "shellQuote": false,
                            "position": 20
                        },
                        "label": "Regions for processing",
                        "doc": "Restrict to comma-separated list of regions."
                    },
                    {
                        "sbg:altPrefix": "-R",
                        "sbg:category": "Configuration",
                        "id": "regions_from_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--regions-file",
                            "shellQuote": false,
                            "position": 21
                        },
                        "label": "Processing regions from file",
                        "doc": "Restrict to regions listed in a file."
                    },
                    {
                        "sbg:altPrefix": "-s",
                        "sbg:category": "Configuration",
                        "id": "samples",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--samples",
                            "shellQuote": false,
                            "position": 23
                        },
                        "label": "Samples list",
                        "doc": "Subset of samples to annotate."
                    },
                    {
                        "sbg:altPrefix": "-S",
                        "sbg:category": "Configuration",
                        "id": "samples_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--samples-file",
                            "shellQuote": false,
                            "position": 25
                        },
                        "label": "Samples file",
                        "doc": "Subset of samples to annotate. If the samples are named differently in the target VCF and the -a, --annotations VCF, the name mapping can be given as \"src_name dst_name\\n\", separated by whitespaces, each pair on a separate line."
                    },
                    {
                        "sbg:toolDefaultValue": "0",
                        "sbg:category": "Execution",
                        "id": "threads",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--threads",
                            "shellQuote": false,
                            "position": 28
                        },
                        "label": "Threads",
                        "doc": "Number of output compression threads to use in addition to main thread. Only used when output type is CompressedBCF CompressedVCF."
                    },
                    {
                        "sbg:altPrefix": "-a",
                        "sbg:category": "Configuration",
                        "id": "annotations",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--annotations",
                            "shellQuote": false,
                            "position": 7
                        },
                        "label": "Annotation file",
                        "doc": "Bgzip-compressed and tabix-indexed file with annotations. The file can be VCF, BED, or a tab-delimited file with mandatory columns CHROM, POS (or, alternatively, FROM and TO), optional columns REF and ALT, and arbitrary number of annotation columns. BED files are expected to have the \".bed.gz\" suffix (case-insensitive), otherwise a tab-delimited file is assumed. Note that in case of tab-delimited file, the coordinates POS, FROM and TO are one-based and inclusive. When REF and ALT are present, only matching VCF records will be annotated. When multiple ALT alleles are present in the annotation file (given as comma-separated list of alleles), at least one must match one of the alleles in the corresponding VCF record. Similarly, at least one alternate allele from a multi-allelic VCF record must be present in the annotation file. Note that flag types, such as \"INFO/FLAG\", can be annotated by including a field with the value \"1\" to set the flag, \"0\" to remove it, or \".\" to keep existing flags.",
                        "sbg:fileTypes": "VCF.GZ, BED.GZ",
                        "secondaryFiles": [
                            {
                                "pattern": ".tbi",
                                "required": true
                            }
                        ]
                    },
                    {
                        "sbg:altPrefix": "-c",
                        "sbg:category": "Configuration",
                        "id": "columns",
                        "type": "string[]?",
                        "inputBinding": {
                            "prefix": "--columns",
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 8
                        },
                        "label": "Columns list",
                        "doc": "Comma-separated list of columns or tags to carry over from the annotation file. If the annotation file is not a VCF/BCF, list describes the columns of the annotation file and must include CHROM, POS (or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which should be ignored can be indicated by \"-\". If the annotation file is a VCF/BCF, only the edited columns/tags must be present and their order does not matter. The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags can be written both as \"INFO/TAG\" or simply \"TAG\", and FORMAT tags can be written as \"FORMAT/TAG\" or \"FMT/TAG\". To carry over all INFO annotations, use \"INFO\". To add all INFO annotations except \"TAG\", use \"^INFO/TAG\". By default, existing values are replaced. To add annotations without overwriting existing values (that is, to add missing tags or add values to existing tags with missing values), use \"+TAG\" instead of \"TAG\". To append to existing values (rather than replacing or leaving untouched), use \"=TAG\" (instead of \"TAG\" or \"+TAG\"). To replace only existing values without modifying missing annotations, use \"-TAG\". If the annotation file is not a VCF/BCF, all new annotations must be defined via -h, --header-lines."
                    },
                    {
                        "sbg:altPrefix": "-h",
                        "sbg:category": "Configuration",
                        "id": "header_lines",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--header-lines",
                            "shellQuote": false,
                            "position": 14
                        },
                        "label": "Header lines",
                        "doc": "Lines which should be appended to the VCF header."
                    },
                    {
                        "sbg:altPrefix": "-I",
                        "sbg:category": "Configuration",
                        "id": "set_id",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--set-id",
                            "shellQuote": false,
                            "position": 15
                        },
                        "label": "Assign set ID (format)",
                        "doc": "Assign ID on the fly. The format is the same as in the query command (see below). By default all existing IDs are replaced. If the format string is preceded by \"+\", only missing IDs will be set."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "rename_chrs",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--rename-chrs",
                            "shellQuote": false,
                            "position": 22
                        },
                        "label": "Rename chromosomes",
                        "doc": "Rename chromosomes according to the map in file, with OLD_NAME NEW_NAME pairs separated by whitespaces, each on a separate line.",
                        "sbg:fileTypes": "TXT"
                    },
                    {
                        "sbg:altPrefix": "-x",
                        "sbg:category": "Configuration",
                        "id": "remove_annotations",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--remove",
                            "shellQuote": false,
                            "position": 26
                        },
                        "label": "Remove annotations",
                        "doc": "List of annotations to remove. Use \"FILTER\" to remove all filters or \"FILTER/SomeFilter\" to remove a specific filter. Similarly, \"INFO\" can be used to remove all INFO tags and \"FORMAT\" to remove all FORMAT tags except GT. To remove all INFO tags except \"FOO\" and \"BAR\", use \"^INFO/FOO,INFO/BAR\" (and similarly for FORMAT and FILTER). \"INFO\" can be abbreviated to \"INF\" and \"FORMAT\" to \"FMT\"."
                    },
                    {
                        "sbg:toolDefaultValue": "1",
                        "sbg:category": "Execution",
                        "id": "cpu",
                        "type": "int?",
                        "label": "Number of CPUs",
                        "doc": "Number of CPUs. Appropriate instance will be chosen based on this parameter."
                    },
                    {
                        "sbg:toolDefaultValue": "1000",
                        "sbg:category": "Execution",
                        "id": "memory",
                        "type": "int?",
                        "label": "Memory in MB",
                        "doc": "Memory in MB. Appropriate instance will be chosen based on this parameter."
                    },
                    {
                        "sbg:altPrefix": "-k",
                        "sbg:category": "Configuration",
                        "id": "keep_sites",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--keep-sites",
                            "shellQuote": false,
                            "position": 27
                        },
                        "label": "Keep sites",
                        "doc": "Keep sites which do not pass -i and -e expressions instead of discarding them."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "collapse",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "snps",
                                    "indels",
                                    "both",
                                    "all",
                                    "some",
                                    "none",
                                    "id"
                                ],
                                "name": "collapse"
                            }
                        ],
                        "label": "Collapse duplicate positions",
                        "doc": "Controls how to treat records with duplicate positions and defines compatible records across multiple input files. Here by \"compatible\" we mean records which should be considered as identical by the tools. For example, when performing line intersections, the desire may be to consider as identical all sites with matching positions (bcftools isec -c all), or only sites with matching variant type (bcftools isec -c snps  -c indels), or only sites with all alleles identical (bcftools isec -c none)."
                    }
                ],
                "outputs": [
                    {
                        "id": "output_file",
                        "doc": "Annotated file with fields from annotation file.",
                        "label": "Annotated output file",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "${\n    var files_array = [].concat(inputs.input_file)\n    var fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() == 'gz') {\n        fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n    }\n\n    if (inputs.output_name) {\n        var out = inputs.output_name\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".annotated.vcf\"\n        } else if (inputs.output_type == 'CompressedVCF') {\n            out += \".annotated.vcf.gz\"\n        } else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".annotated.bcf\"\n        } else if (inputs.output_type == 'CompressedBCF') {\n            out += \".annotated.bcf.gz\"\n        } else {\n            out += \".annotated.vcf\"\n        }\n    } else if (inputs.output_type == 'UncompressedVCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.vcf'\n    } else if (inputs.output_type == 'CompressedVCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.vcf.gz'\n    } else if (inputs.output_type == 'UncompressedBCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.bcf'\n    } else if (inputs.output_type == 'CompressedBCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.bcf.gz'\n    } else out = fname.split('.vcf')[0] + '.annotated.vcf'\n\n    return out\n}",
                            "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
                        },
                        "secondaryFiles": [
                            {
                                "pattern": ".tbi",
                                "required": false
                            }
                        ],
                        "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ"
                    }
                ],
                "doc": "**BCFtools Annotate**: Add or remove annotations. The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, added or removed.\n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n\n### Common Use Cases\n\n* Remove fields: Provide comma separated strings you want to remove to the **Remove annotations** (`--remove`).\n```\n$bcftools annotate -x ID,INFO/DP,FORMAT/DP file.vcf.gz\n```\n\n* Add ID, QUAL and INFO/TAG, not replacing TAG if already present: Add ID,QUAL,+TAG in **Columns list** (`--columns`) input.\n```\n$bcftools annotate -a src.bcf -c ID,QUAL,+TAG dst.bcf\n```\n\n* Annotate from a BED file: Add a BED file on **Annotation file** (`--annotations`) input, add header file on the **Header lines** (`--header-lines`) input and list of columns to be annotated in **Columns list** (`--columns`) input.\n```\n$bcftools annotate -a annots.bed.gz -h annots.hdr -c CHROM,FROM,TO,TAG input.vcf\n```\n\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands if a VCF file is provided on input. If VCF.GZ is given on input only indexing will be done. Output type can still be chosen with the `output type` command.\n\n* If BED file is given on input for annotation, we added subsequent `bgzip` and `tabix` commands if a BED file is provided on input. If VCF.GZ is given on input only indexing will be done.\n\n### Common Issues and Important Notes\n\n * It is expected for annotation file to be already gzipped and indexed. Otherwise, tool will fail. \n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using an input of 7 MB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
                "label": "Bcftools Annotate",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    var files_array = [].concat(inputs.input_file)\n    var fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() == 'gz') {\n        fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n        return \"bcftools index  -f -t \" + fname + \".gz &&\"\n    } else {\n\n        return \"bgzip -c -f \" + fname + \" > \" + fname + \".gz\" + \" && bcftools index -f -t \" + fname + \".gz &&\"\n\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "bcftools"
                    },
                    {
                        "shellQuote": false,
                        "position": 2,
                        "valueFrom": "annotate"
                    },
                    {
                        "prefix": "--output",
                        "shellQuote": false,
                        "position": 5,
                        "valueFrom": "${\n    var files_array = [].concat(inputs.input_file)\n    var fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() == 'gz') {\n        fname = files_array[0].path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n    }\n\n    if (inputs.output_name) {\n        var out = inputs.output_name\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".annotated.vcf\"\n        } else if (inputs.output_type == 'CompressedVCF') {\n            out += \".annotated.vcf.gz\"\n        } else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".annotated.bcf\"\n        } else if (inputs.output_type == 'CompressedBCF') {\n            out += \".annotated.bcf.gz\"\n        } else {\n            out += \".annotated.vcf\"\n        }\n    } else if (inputs.output_type == 'UncompressedVCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.vcf'\n    } else if (inputs.output_type == 'CompressedVCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.vcf.gz'\n    } else if (inputs.output_type == 'UncompressedBCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.bcf'\n    } else if (inputs.output_type == 'CompressedBCF') {\n        var fname_list = fname.split('.')\n        fname_list.pop() // Remove extension\n        out = fname_list.join('.') + '.annotated' + '.bcf.gz'\n    } else out = fname.split('.vcf')[0] + '.annotated.vcf'\n\n    return out\n}"
                    },
                    {
                        "prefix": "--output-type",
                        "shellQuote": false,
                        "position": 19,
                        "valueFrom": "${\n    if (inputs.output_type === 'CompressedBCF') return 'b'\n    if (inputs.output_type === 'UncompressedBCF') return 'u'\n    if (inputs.output_type === 'CompressedVCF') return 'z'\n    if (inputs.output_type === 'UncompressedVCF') return 'v'\n}"
                    },
                    {
                        "prefix": "--collapse",
                        "shellQuote": false,
                        "position": 6,
                        "valueFrom": "${\n\n    if (inputs.collapse == 'snps') {\n        return 'snps'\n    }\n    if (inputs.collapse == 'indels') {\n        return 'indels'\n    }\n    if (inputs.collapse == 'both') {\n        return 'both'\n    }\n    if (inputs.collapse == 'all') {\n        return 'all'\n    }\n    if (inputs.collapse == 'some') {\n        return 'some'\n    }\n    if (inputs.collapse == 'none') {\n        return 'none'\n    }\n    if (inputs.collapse == 'id') {\n        return 'id'\n    }\n\n\n\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n    if (inputs.memory) {\n        return inputs.memory\n    } else {\n        return 1000\n    }\n}",
                        "coresMin": "${\n    if (inputs.cpu) {\n        return inputs.cpu\n    } else {\n        return 1\n    }\n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools:1.9"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entry": "$(inputs.input_file)",
                                "writable": false
                            },
                            {
                                "entry": "$(inputs.annotations)",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "successCodes": [
                    0
                ],
                "temporaryFailCodes": [
                    1
                ],
                "sbg:toolkitVersion": "1.9",
                "abg:revisionNotes": "Added -f option for indexing",
                "sbg:image_url": null,
                "sbg:license": "MIT License",
                "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
                "sbg:categories": [
                    "VCF-Processing"
                ],
                "sbg:toolkit": "bcftools",
                "sbg:cmdPreview": "bgzip -c -f annotated_input_file.vcf > annotated_input_file.vcf.gz && bcftools index -f -t annotated_input_file.vcf.gz && bcftools annotate --output annotated_input_file.annotated.vcf",
                "sbg:links": [
                    {
                        "id": "http://samtools.github.io/bcftools/",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools",
                        "label": "Source code"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/wiki",
                        "label": "Wiki"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/archive/1.9.zip",
                        "label": "Download"
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105870,
                        "sbg:revisionNotes": null
                    }
                ],
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/bcftools-annotate-cwl1/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": null,
                "sbg:modifiedOn": 1637105870,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105870,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a86e4fb9f4cce7a09cd2063073a8e699f16c1d19d267d2b541e2468d635f14e0e"
            },
            "label": "Bcftools Annotate",
            "sbg:x": -938.9439697265625,
            "sbg:y": -202.12379455566406
        },
        {
            "id": "sbg_merge_and_filter_genesis_results_cwl2",
            "in": [
                {
                    "id": "rdata_files",
                    "source": [
                        "single_variant_association_testing/assoc_combined"
                    ]
                },
                {
                    "id": "outname",
                    "source": "outname"
                },
                {
                    "id": "threshold",
                    "source": "threshold"
                },
                {
                    "id": "lz_threshold",
                    "default": 5e-8,
                    "source": "lz_threshold"
                },
                {
                    "id": "test",
                    "default": "single"
                }
            ],
            "out": [
                {
                    "id": "merged"
                },
                {
                    "id": "out_lz"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/sbg-merge-and-filter-genesis-results-cwl1/2",
                "baseCommand": [
                    "Rscript",
                    "merge.R"
                ],
                "inputs": [
                    {
                        "sbg:category": "Input files",
                        "id": "rdata_files",
                        "type": {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--rdata=",
                                "separate": false,
                                "shellQuote": false,
                                "position": 0
                            }
                        },
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 1,
                            "valueFrom": "${\n    var out = \"\"\n    for (var i = 0; i < [].concat(self).length; i++ ){\n        out += \" --rdata=\" + [].concat(self)[i].path\n    }    \n    return out\n}"
                        },
                        "label": "RDATA files",
                        "doc": "RDATA files to be merged",
                        "sbg:fileTypes": "RDATA"
                    },
                    {
                        "sbg:category": "Options",
                        "id": "outname",
                        "type": "string",
                        "inputBinding": {
                            "prefix": "--outname=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 2
                        },
                        "label": "Output name",
                        "doc": "Name to be given to the merged CSV file."
                    },
                    {
                        "id": "threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--threshold=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 3
                        },
                        "label": "Threshold for results filtering",
                        "doc": "Only the results with p values below this threshold will be kept in the final output file."
                    },
                    {
                        "id": "lz_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--lz_threshold=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 4
                        }
                    },
                    {
                        "id": "test",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "single",
                                    "burden",
                                    "skat"
                                ],
                                "name": "test"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--test=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 5
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "merged",
                        "label": "Merged CSV",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*.csv",
                            "outputEval": "${\n    return inheritMetadata(self, inputs.rdata_files)\n\n}"
                        },
                        "sbg:fileTypes": "CSV"
                    },
                    {
                        "id": "out_lz",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*txt"
                        }
                    }
                ],
                "label": "Merge and filter GENESIS results",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 1000,
                        "coresMin": 1
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "rocker/r-base"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "merge.R",
                                "entry": "## Collect arguments\nargs <- commandArgs(TRUE)\n\n## Parse arguments (we expect the form --arg=value)\nparseArgs <- function(x) strsplit(sub(\"^--\", \"\", x), \"=\")\nargsDF <- as.data.frame(do.call(\"rbind\", parseArgs(args)))\narg <- as.list(as.character(argsDF$V2))\nnames(arg) <- argsDF$V1\n\n\ntest <- arg$test\nif (!(test %in% c(\"single\", \"skat\", \"burden\"))) stop(\"This test is not supported!\")\nind <- which(names(arg) == \"rdata\")\nrdata_files <- unlist(arg[ind])\n\nmultmerge <- function(filenames, test) {\n  datalist <- lapply(filenames, function(x) {\n    load(x)\n    if (test == \"single\") {\n      assoc[!is.na(assoc$SPA.pval) & assoc$SPA.pval<as.numeric(arg$threshold),]\n    } else if (test == \"skat\") {\n      assoc$results[!is.na(assoc$results$Score.pval) & assoc$results$Score.pval<as.numeric(arg$threshold),]\n    } else if (test ==\"burden\") {\n      assoc$results[!is.na(assoc$results$pval) & assoc$results$pval<as.numeric(arg$threshold),]\n    }\n    \n  })\n  Reduce(function(x,y) {rbind(x,y)}, datalist)\n}\n\nres <- multmerge(rdata_files, test)\n\nif (test %in% c(\"skat\", \"burden\")) {\n  res$start <- as.integer(res$start)\n  res <- res[with(res, order(chr, start)), ]\n} else if (test == \"single\") {\n  res$pos <- as.integer(res$pos)\n  res <- res[with(res, order(chr, pos)), ]\n}\n\nif (!is.null(arg$lz_threshold)) {\n    if (test == \"single\") {\n        lz_res <- res[res$SPA.pval < as.numeric(arg$lz_threshold),]\n        lz_res <- lz_res[,c(\"variant.id\", \"chr\")]\n    } else if (test == \"skat\") {\n        lz_res <- res[res$Score.pval < as.numeric(arg$lz_threshold),]\n        lz_res <- lz_res[,c(\"chr\", \"start\", \"end\")]\n    } else if (test ==\"burden\") {\n        lz_res <- res[res$pval < as.numeric(arg$lz_threshold),]\n        lz_res <- lz_res[,c(\"chr\", \"start\", \"end\")]\n    }\n    lz_res$pop <- \"EUR\"\n    write.table(lz_res, file = paste0(test,\"_lz.txt\"), quote = FALSE, sep = \" \", row.names = FALSE)\n}\n\nwrite.csv(res, file = paste0(arg$outname, \".csv\"), quote = FALSE, row.names = FALSE)\n",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "sbg:image_url": null,
                "sbg:cmdPreview": "Rscript csvmerge.R --csv=/path/to/csv_files-1.ext --csv=/path/to/csv_files-2.ext --outname=outname-string-value",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105871,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1639043382,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1646833934,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/sbg-merge-and-filter-genesis-results-cwl1/2",
                "sbg:revision": 2,
                "sbg:revisionNotes": "",
                "sbg:modifiedOn": 1646833934,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105871,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 2,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "ad97ed210c5ca36b40dbe24476515128a59fda77d1e9a66f959e86ee9a8f7e2b1",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "SBG Merge and filter GENESIS results",
            "sbg:x": 1920.607177734375,
            "sbg:y": -143.9972686767578
        },
        {
            "id": "vcf_to_gds_1",
            "in": [
                {
                    "id": "vcf_file",
                    "linkMerge": "merge_flattened",
                    "source": [
                        "snpeff_4_3_cwl1/split_vcfs"
                    ]
                },
                {
                    "id": "memory_gb",
                    "default": 16
                },
                {
                    "id": "cpu",
                    "default": 8
                }
            ],
            "out": [
                {
                    "id": "unique_variant_id_gds_per_chr"
                }
            ],
            "run": {
                "class": "Workflow",
                "cwlVersion": "v1.2",
                "id": "markoz/hgi/vcf-to-gds-1/0",
                "doc": "**VCF to GDS** workflow converts VCF or BCF files into Genomic Data Structure (GDS) format. GDS files are required by all workflows utilizing the GENESIS or SNPRelate R packages.\n\nStep 1 (*VCF to GDS*) converts  VCF or BCF files (one per chromosome) into GDS files, with option to keep a subset of **Format** fields (by default only *GT* field). (BCF files may be used instead of  VCF.) \n\nStep 2 (Unique Variant IDs) ensures that each variant has a unique integer ID across the genome. \n\nStep 3 (Check GDS) ensures that no important information is lost during conversion. If Check GDS fails, it is likely that there was an issue during the conversion.\n**Important note:** This step can be time consuming and therefore expensive. Also, failure of this tool is rare. For that reason we allows this step to be optional and it's turned off by default. To turn it on check yes at *check GDS* port. For information on differences in execution time and cost of the same task with and without check GDS  please refer to Benchmarking section.  \n\n### Common use cases\nThis workflow is used for converting VCF files to GDS files.\n\n### Common issues and important notes\n* This pipeline expects that input **Variant files** are separated to be per chromosome and that files are properly named.  It is expected that chromosome is included in the file name in following format:  chr## , where ## is the name of the chromosome (1-24 or X, Y). Chromosome can be included at any part of the filename. Inputs can be in vcf, vcf.gz and bcf format.  Examples: Data_subset_chr1.vcf, Data_subset_chr1.vcf.gz, Data_chr1_subset.vcf, Data_subset_chr1.bcf. \n\n\n\n* **Number of CPUs** parameter should only be used when working with VCF files. The workflow is unable to utilize more than one thread when working with BCF files, and will fail if number of threads is set for BCF conversion.\n\n* **Note:** Variant IDs in output workflow might be different than expected. Unique variants are assigned for one chromosome at a time, in ascending, natural order (1,2,..,24 or X,Y). Variant IDs are integer IDs unique to your data and do not map to *rsIDs* or any other standard identifier. Be sure to use *variant_id* file for down the line workflows generated based on GDS files created by this workflow.\n\n* **Note** This workflow has not been tested on datasets containing more than 62k samples. Since *check_gds* step is very both ram and storage memory demanding, increasing sample count might lead to task failure. In case of any task failure, feel free to contact our support.\n\n* **Note** **Memory GB** should be set when working with larger number of samples (more than 10k). During benchmarking, 4GB of memory were enough when working with 50k samples. This parameter is related to *VCF to GDS* step, different amount of memory is used in other steps.\n\n\n### Changes introduced by Seven Bridges\nFinal step of the workflow is writing checking status to stderr, and therefore it is stored in *job_err.log*, available in the platform *task stats and logs page*. If workflow completes successfully, it means that validation has passed, if workflow fails on *check_gds* step, it means that validation failed.\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running time and cost. \n      \n\n| Sample Count | Total sample size (GB) | Duration  | Cost - spot ($) |  Instance (AWS)  |\n|-------------------|-------------------|------------------|----------|-------------|------------|------------|\n| 1k samples | 0.2                    | 6m                 | 0.34 |  1x c5.18xlarge |\n| 50k simulated samples (VCF.GZ) | 200        | 1d4h                 | 134     |   4x c5.18xlarge |\n| 62k real samples (BCF) | 400                    | 2d11h                 | 139     | 1x c5.18xlarge |\n\n\n\n*Cost shown here were obtained with **spot instances** enabled. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*      \n\n**Note** Both at 50k samples, and 62k samples, termination of spot instance occurred, leading to higher duration and final cost. These results are not removed from benchmark as this behavior is usual and expected, and should be taken into account when using spot instances.\n\n\n### API Python Implementation\n\nThe app's draft task can also be submitted via the **API**. In order to learn how to get your **Authentication token** and **API endpoint** for the corresponding Platform visit our [documentation](https://github.com/sbg/sevenbridges-python#authentication-and-configuration).\n\n```python\nfrom sevenbridges import Api\n\nauthentication_token, api_endpoint = \"enter_your_token\", \"enter_api_endpoint\"\napi = Api(token=authentication_token, url=api_endpoint)\n# Get project_id/app_id from your address bar. Example: https://f4c.sbgenomics.com/u/your_username/project/app\nproject_id, app_id = \"your_username/project\", \"your_username/project/app\"\n# Get file names from files in your project. Example: Names are taken from Data/Public Reference Files.\ninputs = {\n    \"input_gds_files\": api.files.query(project=project_id, names=[\"basename_chr1.gds\", \"basename_chr2.gds\", ..]),\n    \"phenotype_file\": api.files.query(project=project_id, names=[\"name_of_phenotype_file\"])[0],\n    \"null_model_file\": api.files.query(project=project_id, names=[\"name_of_null_model_file\"])[0]\n}\ntask = api.tasks.create(name='Single Variant Association Testing - API Run', project=project_id, app=app_id, inputs=inputs, run=False)\n```\nInstructions for installing and configuring the API Python client, are provided on [github](https://github.com/sbg/sevenbridges-python#installation). For more information about using the API Python client, consult [sevenbridges-python documentation](http://sevenbridges-python.readthedocs.io/en/latest/). **More examples** are available [here](https://github.com/sbg/okAPI).\n\nAdditionally, [API R](https://github.com/sbg/sevenbridges-r) and [API Java](https://github.com/sbg/sevenbridges-java) clients are available. To learn more about using these API clients please refer to the [API R client documentation](https://sbg.github.io/sevenbridges-r/), and [API Java client documentation](https://docs.sevenbridges.com/docs/java-library-quickstart).",
                "label": "VCF to GDS converter",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "inputs": [
                    {
                        "id": "vcf_file",
                        "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
                        "type": "File[]",
                        "label": "Variants Files",
                        "doc": "Input Variants Files.",
                        "sbg:x": -301.39886474609375,
                        "sbg:y": -11
                    },
                    {
                        "id": "memory_gb",
                        "type": "float?",
                        "label": "memory GB",
                        "doc": "Memory to allocate per job. For low number of samples (up to 10k), 1GB is usually enough. For larger number of samples, value should be set higher (50k samples ~ 4GB). Default: 4.",
                        "sbg:x": -363,
                        "sbg:y": 124
                    },
                    {
                        "id": "format",
                        "type": "string[]?",
                        "label": "Format",
                        "doc": "VCF Format fields to keep in GDS file.",
                        "sbg:toolDefaultValue": "GT",
                        "sbg:x": -351,
                        "sbg:y": 305
                    },
                    {
                        "id": "cpu",
                        "type": "int?",
                        "label": "Number of CPUs",
                        "doc": "Number of CPUs to use per job. Default value: 1.",
                        "sbg:x": -295,
                        "sbg:y": 437
                    },
                    {
                        "id": "check_gds_1",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "Yes",
                                    "No"
                                ],
                                "name": "check_gds"
                            }
                        ],
                        "label": "check GDS",
                        "doc": "Choose “Yes” to check for conversion errors by comparing all values in the output GDS file against the input files. The total cost of the job will be ~10x higher if check GDS is “Yes.”",
                        "sbg:toolDefaultValue": "No",
                        "sbg:x": -71.421875,
                        "sbg:y": 406.5
                    }
                ],
                "outputs": [
                    {
                        "id": "unique_variant_id_gds_per_chr",
                        "outputSource": [
                            "unique_variant_id/unique_variant_id_gds_per_chr"
                        ],
                        "sbg:fileTypes": "GDS",
                        "type": "File[]?",
                        "label": "Unique variant ID corrected GDS files per chromosome",
                        "doc": "GDS files in which each variant has a unique identifier across the entire genome. For example, if chromosome 1 has 100 variants and chromosome 2 has 100 variants, the variant.id field will contain 1:100 in the chromosome 1 file and 101:200 in the chromosome 2 file.",
                        "sbg:x": 385.79351806640625,
                        "sbg:y": 44.093116760253906
                    }
                ],
                "steps": [
                    {
                        "id": "vcf2gds",
                        "in": [
                            {
                                "id": "vcf_file",
                                "source": "vcf_file"
                            },
                            {
                                "id": "memory_gb",
                                "source": "memory_gb"
                            },
                            {
                                "id": "cpu",
                                "source": "cpu"
                            },
                            {
                                "id": "format",
                                "source": [
                                    "format"
                                ]
                            }
                        ],
                        "out": [
                            {
                                "id": "gds_output"
                            },
                            {
                                "id": "config_file"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/vcf2gds/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Input Files",
                                    "id": "vcf_file",
                                    "type": "File",
                                    "label": "Variants Files",
                                    "doc": "Input Variants Files.",
                                    "sbg:fileTypes": "VCF, VCF.GZ, BCF"
                                },
                                {
                                    "id": "gds_file_name",
                                    "type": "string?",
                                    "label": "GDS File",
                                    "doc": "Output GDS file."
                                },
                                {
                                    "sbg:category": "Input options",
                                    "sbg:toolDefaultValue": "4",
                                    "id": "memory_gb",
                                    "type": "float?",
                                    "label": "Memory GB",
                                    "doc": "Memory in GB for each job."
                                },
                                {
                                    "sbg:category": "Input Options",
                                    "sbg:toolDefaultValue": "1",
                                    "id": "cpu",
                                    "type": "int?",
                                    "label": "Cpu",
                                    "doc": "Number of CPUs for each tool job."
                                },
                                {
                                    "sbg:category": "General",
                                    "sbg:toolDefaultValue": "GT",
                                    "id": "format",
                                    "type": "string[]?",
                                    "label": "Format",
                                    "doc": "VCF Format fields to keep in GDS file."
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "gds_output",
                                    "doc": "GDS Output File.",
                                    "label": "GDS Output File",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.gds"
                                    },
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "id": "config_file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.config"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                }
                            ],
                            "label": "vcf2gds",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 5,
                                    "valueFrom": "${\n    return \"Rscript /usr/local/analysis_pipeline/R/vcf2gds.R vcf2gds.config\"\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 1,
                                    "valueFrom": "${\n    if (inputs.cpu)\n        return 'export NSLOTS=' + inputs.cpu + ' &&'\n    else\n        return ''\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return \" >> job.out.log\"\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": "${\n    if(inputs.memory_gb)\n        return parseFloat(inputs.memory_gb * 1024)\n    else\n        return 4*1024.0\n}",
                                    "coresMin": "${ if(inputs.cpu)\n        return inputs.cpu \n    else \n        return 1\n}"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.12.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "vcf2gds.config",
                                            "entry": "${  \n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    var config = \"\";\n    config += \"vcf_file \\\"\" + inputs.vcf_file.path + \"\\\"\\n\"\n    if(inputs.format)\n    {\n        config += \"format \\\"\" + inputs.format.join(' ') + \"\\\"\\n\"\n    }\n    else\n    {\n        config += \"format \\\"GT\\\"\\n\"\n    }\n    if(inputs.gds_file_name)\n        config += \"gds_file \\\"\" + inputs.gds_file_name + \"\\\"\\n\"\n    else\n    {    \n        \n        if (inputs.vcf_file.basename.indexOf('.bcf') == -1){\n        \n         var chromosome = \"chr\" + find_chromosome(inputs.vcf_file.path);\n         config += \"gds_file \\\"\" + inputs.vcf_file.path.split('/').pop().split(chromosome)[0] + chromosome + inputs.vcf_file.path.split('/').pop().split(chromosome)[1].split('.vcf')[0] +\".gds\\\"\";\n        } else{\n            \n             var chromosome = \"chr\" + find_chromosome(inputs.vcf_file.path);\n             config += \"gds_file \\\"\" + inputs.vcf_file.path.split('/').pop().split(chromosome)[0] + chromosome + inputs.vcf_file.path.split('/').pop().split(chromosome)[1].split('.bcf')[0] +\".gds\\\"\";    \n            \n               }\n        \n        \n    }\n    return config\n}\n   ",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105873,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/vcf2gds/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105873,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105873,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "af6578aee66a537003d697a8d54c5ddd84129526a6e999f7a384fa19211471415"
                        },
                        "label": "vcf2gds",
                        "scatter": [
                            "vcf_file"
                        ],
                        "sbg:x": -71,
                        "sbg:y": 184
                    },
                    {
                        "id": "unique_variant_id",
                        "in": [
                            {
                                "id": "gds_file",
                                "source": [
                                    "vcf2gds/gds_output"
                                ]
                            }
                        ],
                        "out": [
                            {
                                "id": "unique_variant_id_gds_per_chr"
                            },
                            {
                                "id": "config"
                            }
                        ],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/unique-variant-id/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "id": "gds_file",
                                    "type": "File[]",
                                    "label": "GDS file",
                                    "doc": "List of GDS files produced by VCF2GDS tool.",
                                    "sbg:fileTypes": "GDS"
                                }
                            ],
                            "outputs": [
                                {
                                    "id": "unique_variant_id_gds_per_chr",
                                    "doc": "Corrected GDS files per chromosome.",
                                    "label": "Unique variant ID corrected GDS files per chromosome",
                                    "type": "File[]?",
                                    "outputBinding": {
                                        "glob": "*chr*",
                                        "outputEval": "$(inheritMetadata(self, inputs.gds_file))"
                                    },
                                    "sbg:fileTypes": "GDS"
                                },
                                {
                                    "id": "config",
                                    "doc": "Config file for running the R script.",
                                    "label": "Config file",
                                    "type": "File?",
                                    "outputBinding": {
                                        "glob": "*.config"
                                    },
                                    "sbg:fileTypes": "CONFIG"
                                }
                            ],
                            "label": "unique_variant_id",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 5,
                                    "valueFrom": "${\n    var cmd_line = \"cp \";\n    for (var i=0; i<inputs.gds_file.length; i++)\n        cmd_line += inputs.gds_file[i].path + \" \"\n    cmd_line += \". && \"\n    if(inputs.merged_gds_file)\n    {\n        cmd_line += \"cp \" + inputs.merged_gds_file.path + \" . && \"\n    }\n    return cmd_line\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 10,
                                    "valueFrom": "${\n    return \" Rscript /usr/local/analysis_pipeline/R/unique_variant_ids.R unique_variant_ids.config\"\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${\n    return ' >> job.out.log'\n}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.12.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "unique_variant_ids.config",
                                            "entry": "${\n    function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n    \n    function compareNatural(a,b){\n        return a.localeCompare(b, 'en', {numeric: true, ignorePunctuation: true})\n    }\n    \n    function find_chromosome(file){\n        var chr_array = [];\n        var chrom_num = file.split(\"/\").pop();\n        chrom_num = chrom_num.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n        return chr_array.toString()\n    }\n    var chr = find_chromosome(inputs.gds_file[0].path);\n    var gds = inputs.gds_file[0].path.split('/').pop();\n    \n    \n    \n    var chr_array = [];\n    for (var i = 0; i < inputs.gds_file.length; i++) \n    {\n        var chrom_num = inputs.gds_file[i].path.split(\"/\").pop();\n        chrom_num = chrom_num.split(\"chr\")[1];\n        \n        if(isNumeric(chrom_num.charAt(1)))\n        {\n            chr_array.push(chrom_num.substr(0,2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(0,1))\n        }\n    }\n    \n    \n    \n    chr_array = chr_array.sort(compareNatural);\n    var chrs = \"\";\n    for (var i = 0; i < chr_array.length; i++) \n    {\n        chrs += chr_array[i] + \" \"\n    }\n    \n    \n    var ind_X = chrs.includes(\"X\")\n    var ind_Y = chrs.includes(\"Y\")\n    var ind_M = chrs.includes(\"M\")\n    \n    var chr_array = chrs.split(\" \")\n    var chr_order = [];\n    var chr_result = \"\"\n    \n    for(i=0; i<chr_array.length; i++){\n        \n    if(!isNaN(chr_array[i]) && chr_array[i]!= \"\") {chr_order.push(parseInt(chr_array[i]))}    \n        \n        \n    }\n    \n    \n    chr_order = chr_order.sort(function(a, b){return a-b})\n    for(i=0; i<chr_order.length; i++){\n        \n        chr_result += chr_order[i].toString() + \" \"\n    }\n    \n    if(ind_X) chr_result += \"X \"\n    if(ind_Y) chr_result += \"Y \"\n    if(ind_M) chr_result += \"M \"\n    \n    chrs = chr_result\n\n    \n    \n    var return_arguments = [];\n    return_arguments.push('gds_file \"' + gds.split(\"chr\")[0] + \"chr \"+gds.split(\"chr\"+chr)[1] +'\"');\n    return_arguments.push('chromosomes \"' + chrs + '\"')\n    \n    return return_arguments.join('\\n') + \"\\n\"\n}\n",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement",
                                    "expressionLib": [
                                        "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                                    ]
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105874,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/unique-variant-id/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105874,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105874,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "aa785d15e5a39e2186acbedd39f42410223352dd0099a8d78410d83193b659276"
                        },
                        "label": "Unique Variant ID",
                        "sbg:x": 138,
                        "sbg:y": 97
                    },
                    {
                        "id": "check_gds",
                        "in": [
                            {
                                "id": "vcf_file",
                                "source": [
                                    "vcf_file"
                                ]
                            },
                            {
                                "id": "gds_file",
                                "source": "unique_variant_id/unique_variant_id_gds_per_chr"
                            },
                            {
                                "id": "check_gds",
                                "source": "check_gds_1"
                            }
                        ],
                        "out": [],
                        "run": {
                            "class": "CommandLineTool",
                            "cwlVersion": "v1.2",
                            "$namespaces": {
                                "sbg": "https://sevenbridges.com"
                            },
                            "id": "markoz/hgi/check-gds/0",
                            "baseCommand": [],
                            "inputs": [
                                {
                                    "sbg:category": "Inputs",
                                    "id": "vcf_file",
                                    "type": "File[]",
                                    "label": "Variants file",
                                    "doc": "VCF or BCF files can have two parts split by chromosome identifier.",
                                    "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ"
                                },
                                {
                                    "sbg:category": "Input",
                                    "id": "gds_file",
                                    "type": "File",
                                    "label": "GDS File",
                                    "doc": "GDS file produced by conversion.",
                                    "sbg:fileTypes": "gds"
                                },
                                {
                                    "sbg:category": "Inputs",
                                    "id": "sample_file",
                                    "type": "File?",
                                    "label": "Sample file",
                                    "doc": "Sample file"
                                },
                                {
                                    "sbg:toolDefaultValue": "No",
                                    "sbg:category": "Inputs",
                                    "id": "check_gds",
                                    "type": [
                                        "null",
                                        {
                                            "type": "enum",
                                            "symbols": [
                                                "Yes",
                                                "No"
                                            ],
                                            "name": "check_gds"
                                        }
                                    ],
                                    "label": "check GDS",
                                    "doc": "Choose “Yes” to check for conversion errors by comparing all values in the output GDS file against the input files. The total cost of the job will be ~10x higher if check GDS is “Yes”."
                                }
                            ],
                            "outputs": [],
                            "label": "check_gds",
                            "arguments": [
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 0,
                                    "valueFrom": "${\n   \n   function find_chromosome(file){\n        var chr_array = []\n        var chrom_num = file.split(\"/\").pop()\n        chrom_num = chrom_num.substr(0,chrom_num.lastIndexOf(\".\")).split('_').slice(0,-1).join('_')\n        if(isNumeric(chrom_num.charAt(chrom_num.length-2)))\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 2))\n        }\n        else\n        {\n            chr_array.push(chrom_num.substr(chrom_num.length - 1))\n        }\n        return chr_array.toString()\n    }\n   \n   function isNumeric(s) {\n        return !isNaN(s - parseFloat(s));\n    }\n   \n   var chr = inputs.gds_file.path.split('chr')[1].split('.')[0];\n\n   \n   if(inputs.check_gds == 'Yes'){\n    return \"Rscript /usr/local/analysis_pipeline/R/check_gds.R check_gds.config \" +\"--chromosome \" + chr}\n    else{\n         return 'echo Check GDS step skipped.'\n    }\n}"
                                },
                                {
                                    "prefix": "",
                                    "shellQuote": false,
                                    "position": 100,
                                    "valueFrom": "${return \" >> job.out.log\"}"
                                }
                            ],
                            "requirements": [
                                {
                                    "class": "ShellCommandRequirement"
                                },
                                {
                                    "class": "ResourceRequirement",
                                    "ramMin": 10000
                                },
                                {
                                    "class": "DockerRequirement",
                                    "dockerPull": "uwgac/topmed-master:2.12.0"
                                },
                                {
                                    "class": "InitialWorkDirRequirement",
                                    "listing": [
                                        {
                                            "entryname": "check_gds.config",
                                            "entry": "${  \n    var config = \"\";\n    var vcf_array = [].concat(inputs.vcf_file);\n    var vcf_first_part = vcf_array[0].path.split('chr')[0];\n    var vcf_second_part = vcf_array[0].path.split('chr')[1].split('.');\n    vcf_second_part.shift();\n    config += \"vcf_file \\\"\" + vcf_first_part + 'chr .' + vcf_second_part.join('.') + \"\\\"\\n\";\n    var gds_first_part = inputs.gds_file.path.split('chr')[0];\n    var gds_second_part = inputs.gds_file.path.split('chr')[1].split('.');\n    gds_second_part.shift();\n    config += \"gds_file \\\"\" + gds_first_part + 'chr .' + gds_second_part.join('.') + \"\\\"\\n\";\n    if(inputs.sample_file)\n        config += \"sample_file \\\"\" + inputs.sample_file.path + \"\\\"\\n\"\n    return config\n}\n   ",
                                            "writable": false
                                        }
                                    ]
                                },
                                {
                                    "class": "InlineJavascriptRequirement"
                                }
                            ],
                            "hints": [
                                {
                                    "class": "sbg:SaveLogs",
                                    "value": "job.out.log"
                                }
                            ],
                            "sbg:projectName": "HGI",
                            "sbg:revisionsInfo": [
                                {
                                    "sbg:revision": 0,
                                    "sbg:modifiedBy": "marko_zecevic",
                                    "sbg:modifiedOn": 1637105875,
                                    "sbg:revisionNotes": null
                                }
                            ],
                            "sbg:image_url": null,
                            "sbg:appVersion": [
                                "v1.2"
                            ],
                            "sbg:id": "markoz/hgi/check-gds/0",
                            "sbg:revision": 0,
                            "sbg:revisionNotes": null,
                            "sbg:modifiedOn": 1637105875,
                            "sbg:modifiedBy": "marko_zecevic",
                            "sbg:createdOn": 1637105875,
                            "sbg:createdBy": "marko_zecevic",
                            "sbg:project": "markoz/hgi",
                            "sbg:sbgMaintained": false,
                            "sbg:validationErrors": [],
                            "sbg:contributors": [
                                "marko_zecevic"
                            ],
                            "sbg:latestRevision": 0,
                            "sbg:publisher": "sbg",
                            "sbg:content_hash": "a09ff33348281e661a17f1324bc60c66f3f62187f71732d1a08bce7749a623281"
                        },
                        "label": "Check GDS",
                        "scatter": [
                            "gds_file"
                        ],
                        "sbg:x": 374.6356201171875,
                        "sbg:y": 303.9109191894531
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:AWSInstanceType",
                        "value": "c5.18xlarge;ebs-gp2;700"
                    },
                    {
                        "class": "sbg:maxNumberOfParallelInstances",
                        "value": "5"
                    },
                    {
                        "class": "sbg:AzureInstanceType",
                        "value": "Standard_F64s_v2;PremiumSSD;1024"
                    }
                ],
                "requirements": [
                    {
                        "class": "ScatterFeatureRequirement"
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    },
                    {
                        "class": "StepInputExpressionRequirement"
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105875,
                        "sbg:revisionNotes": "Workflow decomposed"
                    }
                ],
                "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/vcf-to-gds-1/0.png",
                "sbg:license": "MIT",
                "sbg:toolAuthor": "TOPMed DCC",
                "sbg:links": [
                    {
                        "id": "https://github.com/UW-GAC/analysis_pipeline",
                        "label": "Source Code, Download"
                    },
                    {
                        "id": "https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz567/5536872?redirectedFrom=fulltext",
                        "label": "Publication"
                    },
                    {
                        "id": "https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test.html",
                        "label": "Home Page"
                    },
                    {
                        "id": "https://bioconductor.org/packages/devel/bioc/manuals/GENESIS/man/GENESIS.pdf",
                        "label": "Documentation"
                    }
                ],
                "sbg:categories": [
                    "GWAS",
                    "VCF Processing",
                    "CWL1.0"
                ],
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/vcf-to-gds-1/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": "Workflow decomposed",
                "sbg:modifiedOn": 1637105875,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105875,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a30aeaa0b8f6e18078799490e085a41d511b3cd8540292574f37eaf670e6bcb14"
            },
            "label": "VCF to GDS converter",
            "sbg:x": -538,
            "sbg:y": -220.33895874023438
        },
        {
            "id": "hwe_filter_cwl1",
            "in": [
                {
                    "id": "gdsfile",
                    "source": "vcf_to_gds_1/unique_variant_id_gds_per_chr"
                },
                {
                    "id": "phenofile",
                    "source": "phenotype_file"
                },
                {
                    "id": "field",
                    "source": "field"
                },
                {
                    "id": "cases_threshold",
                    "source": "cases_threshold"
                },
                {
                    "id": "controls_threshold",
                    "source": "controls_threshold"
                }
            ],
            "out": [
                {
                    "id": "variants_kept"
                },
                {
                    "id": "filtered_out_variants"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/hwe-filter-cwl1/2",
                "baseCommand": [
                    "Rscript",
                    "hwe.R"
                ],
                "inputs": [
                    {
                        "sbg:category": "Input files",
                        "id": "gdsfile",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "--gdsfile=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 1
                        },
                        "label": "GDS file",
                        "doc": "GDS file with variants",
                        "sbg:fileTypes": "GDS"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "phenofile",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "--phenofile=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 2
                        },
                        "label": "Phenotype file",
                        "doc": "An RData file containing phenotype information.",
                        "sbg:fileTypes": "RDATA"
                    },
                    {
                        "id": "field",
                        "type": "string",
                        "inputBinding": {
                            "prefix": "--field=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 3
                        },
                        "label": "Case/control indicator variable",
                        "doc": "A variable in the phenotype which indicates samples belonging to either cases or control group."
                    },
                    {
                        "sbg:toolDefaultValue": "1e-10",
                        "id": "cases_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--cases_threshold=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "HWE threshold for cases",
                        "doc": "HWE test p.value threshold for cases"
                    },
                    {
                        "sbg:toolDefaultValue": "1e-6",
                        "id": "controls_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--controls_threshold=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 5
                        },
                        "label": "HWE test p.value threshold for controls",
                        "doc": "HWE test p.value threshold for controls"
                    }
                ],
                "outputs": [
                    {
                        "id": "variants_kept",
                        "label": "Variants which passed filtering",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*rds",
                            "outputEval": "${\n    return inheritMetadata(self, inputs.gdsfile)\n\n}"
                        },
                        "sbg:fileTypes": "RDS"
                    },
                    {
                        "id": "filtered_out_variants",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "*csv"
                        },
                        "sbg:fileTypes": "CSV"
                    }
                ],
                "label": "HWE filter",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 4000,
                        "coresMin": 2
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/marko_zecevic/seqvartools:1.0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "hwe.R",
                                "entry": "## Collect arguments\nargs <- commandArgs(TRUE)\n\n## Parse arguments (we expect the form --arg=value)\nparseArgs <- function(x) strsplit(sub(\"^--\", \"\", x), \"=\")\nargsDF <- as.data.frame(do.call(\"rbind\", parseArgs(args)))\narg <- as.list(as.character(argsDF$V2))\nnames(arg) <- argsDF$V1\n\n\nlibrary(SeqVarTools)\n\n\n# gdsfile <- \"/sbgenomics/project-files/_2_imgge.covid.final.concatenated.analyzed.recode.sorted.annotated.chr9.gds\"\ngdsfile <- arg$gdsfile\n# phenofile <- \"/sbgenomics/project-files/phenotype_pca.RData\"\nphenofile <- arg$phenofile\n# field <- \"hgi\"\nfield <- arg$field\n# cases_threshold <- 1e-10\nif (is.null(arg$cases_threshold)) {\n  cases_threshold <- 1e-10\n} else cases_threshold <- arg$cases_threshold\n# controls_threshold <- 1e-6\nif (is.null(arg$controls_threshold)) {\n  controls_threshold <- 1e-6\n} else controls_threshold <- arg$controls_threshold\n\n\ngds <- seqOpen(gdsfile)\nload(phenofile)\ncases <- annot$sample.id[annot[[field]] == 1]\ncontrols <- annot$sample.id[annot[[field]] == 0]\nlowp_ids <- integer(0)\n\nif (seqGetData(gds, \"chromosome\")[1] == \"X\") {\n  female <- annot$sample.id[annot$sex == \"F\"]\n  cases <- cases[cases %in% female]\n  controls <- controls[controls %in% female]\n}\n\n## CASES\nseqSetFilter(gds, sample.id=cases)\n# run HWE test\nhwe.res <- hwe(gds)\n# identify variants with small p-values\nlowp_cases <- !is.na(hwe.res$p) & hwe.res$p < cases_threshold\nif (sum(lowp_cases) > 0) {\n  df_cases <- hwe.res[lowp_cases,]\n  lowp_ids <- c(lowp_ids, df_cases$variant.id)\n  write.csv(df_cases, file = paste0(tools::file_path_sans_ext(basename(gdsfile)),\"lowp_cases.csv\"), quote = FALSE, row.names = FALSE)\n}\nseqResetFilter(gds)\n\n## CONTROLS\nseqSetFilter(gds, sample.id=controls)\n# run HWE test\nhwe.res <- hwe(gds)\n# identify variants with small p-values\nlowp_controls <- !is.na(hwe.res$p) & hwe.res$p < controls_threshold\nif (sum(lowp_controls) > 0) {\n  df_controls <- hwe.res[lowp_controls,]\n  lowp_ids <- c(lowp_ids, df_controls$variant.id)\n  write.csv(df_controls, file = paste0(tools::file_path_sans_ext(basename(gdsfile)),\".lowp_controls.csv\"), quote = FALSE, row.names = FALSE)\n}\nseqResetFilter(gds)\n\n\nvariant.id <- seqGetData(gds, \"variant.id\")\nkeep <- variant.id[!variant.id %in% lowp_ids]\n\nsaveRDS(keep, file = paste0(tools::file_path_sans_ext(basename(gdsfile)), \".variants.hwe_filtered.rds\"))\n\nseqClose(gds)\n\n\n",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "sbg:image_url": null,
                "sbg:cmdPreview": "Rscript csvmerge.R --csv=/path/to/csv_files-1.ext --csv=/path/to/csv_files-2.ext --outname=outname-string-value",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105877,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1638920179,
                        "sbg:revisionNotes": "test X on females only"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1638974556,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/hwe-filter-cwl1/2",
                "sbg:revision": 2,
                "sbg:revisionNotes": "",
                "sbg:modifiedOn": 1638974556,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105877,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 2,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a1cc9b239673a69c1486cb17a043dfbcc79333dd0ecfd32ad48207fcf83ff86b5"
            },
            "label": "HWE filter",
            "scatter": [
                "gdsfile"
            ],
            "sbg:x": -307.62152099609375,
            "sbg:y": -119.7964096069336
        },
        {
            "id": "merge_rds_arrays",
            "in": [
                {
                    "id": "rds_files",
                    "source": [
                        "hwe_filter_cwl1/variants_kept"
                    ]
                }
            ],
            "out": [
                {
                    "id": "output"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/merge-rds-arrays/0",
                "baseCommand": [
                    "Rscript",
                    "merge.R"
                ],
                "inputs": [
                    {
                        "sbg:category": "Input files",
                        "id": "rds_files",
                        "type": {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--rds=",
                                "separate": false,
                                "shellQuote": false,
                                "position": 0
                            }
                        },
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 1
                        },
                        "label": "RDS files",
                        "doc": "RDS files with ids of variants that passed HWE test.",
                        "sbg:fileTypes": "RDS"
                    }
                ],
                "outputs": [
                    {
                        "id": "output",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*.rds"
                        }
                    }
                ],
                "label": "Merge RDS arrays",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 8000,
                        "coresMin": 4
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "rocker/r-base"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "merge.R",
                                "entry": "## Collect arguments\nargs <- commandArgs(TRUE)\n\n## Parse arguments (we expect the form --arg=value)\nparseArgs <- function(x) strsplit(sub(\"^--\", \"\", x), \"=\")\nargsDF <- as.data.frame(do.call(\"rbind\", parseArgs(args)))\narg <- as.list(as.character(argsDF$V2))\nnames(arg) <- argsDF$V1\n\n\nrds_files <- unlist(arg, use.names = FALSE)\n\nvars <- sort(unlist(sapply(rds_files, readRDS), use.names = FALSE))\n\nsaveRDS(vars, \"variant_ids.rds\")",
                                "writable": false
                            }
                        ]
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105877,
                        "sbg:revisionNotes": null
                    }
                ],
                "sbg:image_url": null,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/merge-rds-arrays/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": null,
                "sbg:modifiedOn": 1637105877,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1637105877,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "aed65db21b3156bd3ef319bf77ecfc5a0a61eb1f220dc95f39ee6b18f35caa17f"
            },
            "label": "Merge RDS arrays",
            "sbg:x": -112.14628601074219,
            "sbg:y": -54.788116455078125
        },
        {
            "id": "bcftools_view_1_10_1",
            "in": [
                {
                    "id": "in_variants",
                    "source": "in_variants"
                },
                {
                    "id": "samples_file",
                    "source": "samples_file"
                }
            ],
            "out": [
                {
                    "id": "out_variants"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.0",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/bcftools-view-1-10-1/0",
                "baseCommand": [],
                "inputs": [
                    {
                        "sbg:category": "File inputs",
                        "id": "in_variants",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 48,
                            "valueFrom": "${\n    var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        return froot + \".gz\"} \n    else {\n        if(fname.split('.').pop() == 'bcf'){\n            var index_csi_file = fname + '.csi'\n            if (in_file.secondaryFiles[0]){\n                var secondary_given = in_file.secondaryFiles[0].path.replace(/^.*[\\\\\\/]/, '');}\n                if(secondary_given == index_csi_file){\n                    return fname;}\n            \n        }\n        return fname + \".gz\"\n    }\n}"
                        },
                        "label": "Input variants file",
                        "doc": "Input variants file.",
                        "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
                        "secondaryFiles": [
                            "${ \n    if(self.basename.split('.').pop() == 'gz'){\n    if(self.nameroot.split('.').pop() == 'bcf'){\n        return self.nameroot + \".gz.csi\"}\n    else{\n        return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n    if(self.basename.split('.').pop() == 'bcf'){\n        return self.basename + \".csi\"\n    }\n    else{\n    return self.basename + \".tbi\"}\n}\n\n}"
                        ]
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-i",
                        "id": "include_expression",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--include",
                            "shellQuote": false,
                            "position": 24
                        },
                        "label": "Include expression",
                        "doc": "Include only sites for which the expression is true."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-e",
                        "id": "exclude_expression",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--exclude",
                            "shellQuote": false,
                            "position": 23
                        },
                        "label": "Exclude expression",
                        "doc": "Exclude sites for which the expression is true."
                    },
                    {
                        "sbg:category": "Output options",
                        "id": "output_name",
                        "type": "string?",
                        "label": "Output file name",
                        "doc": "Name of the output file."
                    },
                    {
                        "sbg:category": "Config inputs",
                        "sbg:altPrefix": "-O",
                        "id": "output_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "CompressedBCF",
                                    "UncompressedBCF",
                                    "CompressedVCF",
                                    "UncompressedVCF"
                                ],
                                "name": "output_type"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--output-type",
                            "shellQuote": false,
                            "position": 9,
                            "valueFrom": "${\n    if (self == 0) {\n        self = null;\n        inputs.output_type = null\n    };\n\n\n    if (inputs.output_type === 'CompressedBCF') return 'b'\n    if (inputs.output_type === 'UncompressedBCF') return 'u'\n    if (inputs.output_type === 'CompressedVCF') return 'z'\n    if (inputs.output_type === 'UncompressedVCF') return 'v'\n}"
                        },
                        "label": "Output type",
                        "doc": "b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v].",
                        "default": 0
                    },
                    {
                        "sbg:category": "Config inputs",
                        "sbg:altPrefix": "-r",
                        "id": "regions",
                        "type": "string[]?",
                        "inputBinding": {
                            "prefix": "--regions",
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 10
                        },
                        "label": "Regions for processing",
                        "doc": "Restrict to comma-separated list of regions (e.g. chr|chr:pos|chr:from-to|chr:from-[,…])."
                    },
                    {
                        "sbg:altPrefix": "-R",
                        "sbg:category": "File inputs",
                        "id": "regions_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--regions-file",
                            "shellQuote": false,
                            "position": 11
                        },
                        "label": "Regions from file",
                        "doc": "Regions listed in a file.",
                        "sbg:fileTypes": "BED, TXT"
                    },
                    {
                        "sbg:category": "Config inputs",
                        "sbg:altPrefix": "-t",
                        "id": "targets",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--targets",
                            "shellQuote": false,
                            "position": 12
                        },
                        "label": "Targets",
                        "doc": "Similar to regions option but streams rather than index-jumps."
                    },
                    {
                        "sbg:category": "File inputs",
                        "sbg:altPrefix": "-T",
                        "id": "targets_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--targets-file",
                            "shellQuote": false,
                            "position": 13
                        },
                        "label": "Targets file",
                        "doc": "Similar to regions file option but streams rather than index-jumps.",
                        "sbg:fileTypes": "BED, TXT"
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:toolDefaultValue": "0",
                        "id": "threads",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--threads",
                            "shellQuote": false,
                            "position": 38
                        },
                        "label": "Threads",
                        "doc": "Number of output compression threads to use in addition to main thread. Only used when output type is compressed BCF compressed VCF."
                    },
                    {
                        "sbg:category": "Config inputs",
                        "id": "no_version",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--no-version",
                            "shellQuote": false,
                            "position": 8
                        },
                        "label": "Don't append version to header",
                        "doc": "Do not append version and command line to the header."
                    },
                    {
                        "sbg:toolDefaultValue": "1",
                        "sbg:category": "Execution",
                        "id": "cpu_per_job",
                        "type": "int?",
                        "label": "CPU per job",
                        "doc": "Number of CPUs per job. Appropriate instance will be chosen based on this parameter."
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:toolDefaultValue": "1000",
                        "id": "mem_per_job",
                        "type": "int?",
                        "label": "Memory per job",
                        "doc": "Memory per job in MB. Appropriate instance will be chosen based on this parameter."
                    },
                    {
                        "sbg:category": "Output options",
                        "sbg:altPrefix": "-G",
                        "id": "drop_genotypes",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--drop-genotypes",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Drop genotypes",
                        "doc": "Drop individual genotype information (after subsetting if -s option set)."
                    },
                    {
                        "sbg:category": "Output options",
                        "sbg:altPrefix": "-h",
                        "id": "header_only",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--header-only",
                            "shellQuote": false,
                            "position": 5
                        },
                        "label": "Header only",
                        "doc": "Output the VCF header only."
                    },
                    {
                        "sbg:category": "Output options",
                        "sbg:altPrefix": "-H",
                        "id": "suppress_header",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--no-header",
                            "shellQuote": false,
                            "position": 6
                        },
                        "label": "Suppress header",
                        "doc": "Suppress the header in VCF output."
                    },
                    {
                        "sbg:category": "Output options",
                        "sbg:altPrefix": "-l",
                        "id": "compression_level",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--compression-level",
                            "shellQuote": false,
                            "position": 7
                        },
                        "label": "Compression level",
                        "doc": "Compression level: 0 uncompressed, 1 best speed, 9 best compression."
                    },
                    {
                        "sbg:category": "Subset options",
                        "sbg:altPrefix": "-a",
                        "id": "trim_alt_alleles",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--trim-alt-alleles",
                            "shellQuote": false,
                            "position": 14
                        },
                        "label": "Trim alt alleles",
                        "doc": "Trim alternate alleles not seen in the subset."
                    },
                    {
                        "sbg:category": "Subset options",
                        "sbg:altPrefix": "-I",
                        "id": "no_update",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--no-update",
                            "shellQuote": false,
                            "position": 15
                        },
                        "label": "Do not update INFO fields",
                        "doc": "Do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)."
                    },
                    {
                        "sbg:category": "Subset options",
                        "id": "force_sample",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--force-samples",
                            "shellQuote": false,
                            "position": 16
                        },
                        "label": "Force samples",
                        "doc": "Only warn about unknown subset samples."
                    },
                    {
                        "sbg:category": "Subset options",
                        "sbg:altPrefix": "-s",
                        "id": "samples",
                        "type": "string[]?",
                        "inputBinding": {
                            "prefix": "--samples",
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 17
                        },
                        "label": "Samples list",
                        "doc": "Comma separated list of samples to include (or exclude with \"^\" prefix)."
                    },
                    {
                        "sbg:category": "Subset options",
                        "sbg:altPrefix": "-S",
                        "id": "samples_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--samples-file",
                            "shellQuote": false,
                            "position": 18
                        },
                        "label": "Samples file",
                        "doc": "File of samples to include (or exclude with \"^\" prefix).",
                        "sbg:fileTypes": "TXT"
                    },
                    {
                        "sbg:category": "Filter options",
                        "id": "min_count",
                        "type": [
                            "null",
                            {
                                "type": "record",
                                "fields": [
                                    {
                                        "sbg:stageInput": null,
                                        "sbg:category": "Filter options",
                                        "sbg:altPrefix": "-m",
                                        "name": "minimum_count",
                                        "type": "int?",
                                        "label": "Minimum allele count",
                                        "doc": "Minimum allele count."
                                    },
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:altPrefix": "-m",
                                        "name": "min_count_type",
                                        "type": [
                                            "null",
                                            {
                                                "type": "enum",
                                                "symbols": [
                                                    ":nref",
                                                    ":alt1",
                                                    ":minor",
                                                    ":major",
                                                    ":nonmajor"
                                                ],
                                                "name": "min_count_type"
                                            }
                                        ],
                                        "label": "Minimum allele count type",
                                        "doc": "Type of alleles."
                                    }
                                ],
                                "name": "min_count"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 20
                        },
                        "label": "Minimum allele count",
                        "doc": "Minimum allele count (INFO/AC) of sites to be printed. Specifying the type of allele is optional."
                    },
                    {
                        "sbg:category": "Filter options",
                        "id": "max_count",
                        "type": [
                            "null",
                            {
                                "type": "record",
                                "fields": [
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "maximum_count",
                                        "type": "int?",
                                        "label": "Maximum allele count",
                                        "doc": "Maximum allele count."
                                    },
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "max_count_type",
                                        "type": [
                                            "null",
                                            {
                                                "type": "enum",
                                                "symbols": [
                                                    ":nref",
                                                    ":alt1",
                                                    ":minor",
                                                    "major",
                                                    "nonmajor"
                                                ],
                                                "name": "max_count_type"
                                            }
                                        ],
                                        "label": "Maximum allele count type",
                                        "doc": "Type of alleles."
                                    }
                                ],
                                "name": "max_count"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 20
                        },
                        "label": "Maximum allele count",
                        "doc": "Maximum allele count (INFO/AC) of sites to be printed. Specifying the type of allele is optional and can be set to non-reference (nref, the default), 1st alternate (alt1), the least frequent (minor), the most frequent (major) or sum of all but the most frequent (nonmajor) alleles."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-f",
                        "id": "apply_filters",
                        "type": "string[]?",
                        "inputBinding": {
                            "prefix": "--apply-filters",
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 21
                        },
                        "label": "Apply filters",
                        "doc": "Require at least one of the listed FILTER strings (e.g. \"PASS,.\")."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-g",
                        "id": "genotype",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "hom",
                                    "^hom",
                                    "het",
                                    "^het",
                                    "miss",
                                    "^miss"
                                ],
                                "name": "genotype"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--genotype",
                            "shellQuote": false,
                            "position": 22
                        },
                        "label": "Genotype",
                        "doc": "Require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-k",
                        "id": "known",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--known",
                            "shellQuote": false,
                            "position": 25
                        },
                        "label": "Known sites",
                        "doc": "Print known sites only (ID column is not \".\")."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-n",
                        "id": "novel",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--novel",
                            "shellQuote": false,
                            "position": 26
                        },
                        "label": "Novel sites",
                        "doc": "Print novel sites only (ID column is \".\")."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-m",
                        "id": "min_alleles",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-alleles",
                            "shellQuote": false,
                            "position": 27
                        },
                        "label": "Minimum alleles",
                        "doc": "Print sites with at least INT alleles listed in REF and ALT columns."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-M",
                        "id": "max_alleles",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-alleles",
                            "shellQuote": false,
                            "position": 28
                        },
                        "label": "Maximum alleles",
                        "doc": "Print sites with at most INT alleles listed in REF and ALT columns."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-p",
                        "id": "phased",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--phased",
                            "shellQuote": false,
                            "position": 29
                        },
                        "label": "Phased",
                        "doc": "Print sites where all samples are phased. Haploid genotypes are considered phased. Missing genotypes considered unphased unless the phased bit is set."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-P",
                        "id": "exclude_phased",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--exclude-phased",
                            "shellQuote": false,
                            "position": 30
                        },
                        "label": "Exclude phased",
                        "doc": "Exclude sites where all samples are phased."
                    },
                    {
                        "sbg:category": "Filter options",
                        "id": "min_freq",
                        "type": [
                            "null",
                            {
                                "type": "record",
                                "fields": [
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "minimum_freq",
                                        "type": "float?",
                                        "label": "Minimum frequency",
                                        "doc": "Minimum allele frequency."
                                    },
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "min_freq_type",
                                        "type": [
                                            "null",
                                            {
                                                "type": "enum",
                                                "symbols": [
                                                    ":nref",
                                                    ":alt1",
                                                    ":minor",
                                                    ":major",
                                                    ":nonmajor"
                                                ],
                                                "name": "min_freq_type"
                                            }
                                        ],
                                        "label": "Minimum frequency type",
                                        "doc": "Frequency type."
                                    }
                                ],
                                "name": "min_freq"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 31
                        },
                        "label": "Minimum frequency",
                        "doc": "Minimum allele frequency (INFO/AC / INFO/AN) of sites to be printed. Specifying the type of allele is optional and can be set to non-reference (nref, the default), 1st alternate (alt1), the least frequent (minor), the most frequent (major) or sum of all but the most frequent (nonmajor) alleles."
                    },
                    {
                        "sbg:category": "Filter options",
                        "id": "max_freq",
                        "type": [
                            "null",
                            {
                                "type": "record",
                                "fields": [
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "maximum_freq",
                                        "type": "float?",
                                        "label": "Maximum frequency",
                                        "doc": "Maximum allele frequency."
                                    },
                                    {
                                        "sbg:category": "Filter options",
                                        "sbg:stageInput": null,
                                        "name": "max_freq_type",
                                        "type": [
                                            "null",
                                            {
                                                "type": "enum",
                                                "symbols": [
                                                    ":nref",
                                                    ":alt1",
                                                    ":minor",
                                                    ":major",
                                                    ":nonmajor"
                                                ],
                                                "name": "max_freq_type"
                                            }
                                        ],
                                        "label": "Maximum frequency type",
                                        "doc": "Maximum frequency type."
                                    }
                                ],
                                "name": "max_freq"
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 32
                        },
                        "label": "Maximum frequency",
                        "doc": "Maximum allele frequency (INFO/AC / INFO/AN) of sites to be printed. Specifying the type of allele is optional and can be set to non-reference (nref, the default), 1st alternate (alt1), the least frequent (minor), the most frequent (major) or sum of all but the most frequent (nonmajor) alleles."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-u",
                        "id": "uncalled",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--uncalled",
                            "shellQuote": false,
                            "position": 33
                        },
                        "label": "Uncalled genotype",
                        "doc": "Print sites without a called genotype."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-U",
                        "id": "exclude_uncalled",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--exclude-uncalled",
                            "shellQuote": false,
                            "position": 34
                        },
                        "label": "Exclude uncalled genotype",
                        "doc": "Exclude sites without a called genotype."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-v",
                        "id": "types",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "snps",
                                    "indels",
                                    "mnps",
                                    "ref",
                                    "bnd",
                                    "other"
                                ],
                                "name": "types"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--types",
                            "shellQuote": false,
                            "position": 35
                        },
                        "label": "Variant types",
                        "doc": "Comma-separated list of variant types to select. Site is selected if any of the ALT alleles is of the type requested. Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags like INFO/INDEL or INFO/VT."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-V",
                        "id": "exclude_types",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "snps",
                                    "indels",
                                    "mnps",
                                    "ref",
                                    "bnd",
                                    "other"
                                ],
                                "name": "exclude_types"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--exclude-types",
                            "shellQuote": false,
                            "position": 37
                        },
                        "label": "Exclude variant types",
                        "doc": "Comma-separated list of variant types to exclude. Site is excluded if any of the ALT alleles is of the type requested. Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags like INFO/INDEL or INFO/VT."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-x",
                        "id": "private",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--private",
                            "shellQuote": false,
                            "position": 39
                        },
                        "label": "Private to subset samples",
                        "doc": "Print sites where only the subset samples carry an non-reference allele. Requires --samples or --samples-file."
                    },
                    {
                        "sbg:category": "Filter options",
                        "sbg:altPrefix": "-X",
                        "id": "exclude_private",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--exclude-private",
                            "shellQuote": false,
                            "position": 40
                        },
                        "label": "Exclude private sites to subset samples",
                        "doc": "Exclude sites where only the subset samples carry an non-reference allele."
                    }
                ],
                "outputs": [
                    {
                        "id": "out_variants",
                        "doc": "Output subset VCF file.",
                        "label": "Output VCF file",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "${  var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.vcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id}\n\n    if (inputs.output_name) {\n        var out = inputs.output_name\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".subset.vcf\"} \n        else if (inputs.output_type == 'CompressedVCF') {\n            out += \".subset.vcf.gz\"} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".subset.bcf\"} \n        else if (inputs.output_type == 'CompressedBCF') {\n            out += \".subset.bcf.gz\"} \n        else {\n            out += \".subset.vcf\"}\n    }\n    \n    else if (inputs.output_type == 'UncompressedVCF') {\n        var out = froot + '.subset' + '.vcf'} \n    else if (inputs.output_type == 'CompressedVCF') {\n        var out = froot + '.subset' + '.vcf.gz'} \n    else if (inputs.output_type == 'UncompressedBCF') {\n        var out = froot + '.subset' + '.bcf'} \n    else if (inputs.output_type == 'CompressedBCF') {\n        var out = froot + '.subset' + '.bcf.gz'} \n    else var out = froot + '.subset.vcf'\n\n    return out\n}",
                            "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
                        },
                        "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ"
                    }
                ],
                "doc": "**BCFtools View**: View, subset and filter VCF or BCF files by position and filtering expression (Former BCFtools Subset).\n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding docs can be found at the bottom of the page.\n\n### Common Use Cases\n\n* Convert between VCF and BCF using the **Output type** (`--output-type`) option\n```\n$bcftools view --output-type v input.bcf.gz\n```\n\n* Output only header of VCF using the **Header only** (`--header-only`) option\n```\n$bcftools view --header-only input.vcf.gz\n```\n\n* Subset VCF to output SNPs only using the **Variant types** (`--types`) option\n```\n$bcftools view --types snps input.vcf.gz\n```\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands if a VCF file is provided on input. If VCF.GZ is given on input, only indexing will be done. Index file `.tbi` is added as secondary file of `in_variants` input which means if VCF.GZ is provided on input, tool will look for index file in project or previous tool (in case of usage in workflow) and if present, it will not perform nor compressing nor indexing. Output type can still be chosen with the `output type` command.\n\n### Common Issues and Important Notes\n\n * No common issues specific to the tool's execution on the Seven Bridges platform have been detected.\n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using an input of 7 MB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
                "label": "Bcftools View",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    var in_files_array = [].concat(inputs.in_variants);\n    var in_file = in_files_array[0];\n    var fname = in_file.basename;\n    var fname_ext = in_file.nameext;\n    var froot = in_file.nameroot;\n    if (fname_ext == '.gz') {\n        var index_tbi_file = fname + '.tbi';\n        var index_csi_file = fname + '.csi'\n        if (in_file.secondaryFiles[0]){\n        var secondary_given = in_file.secondaryFiles[0].path.replace(/^.*[\\\\\\/]/, '');}\n            if(secondary_given == index_tbi_file || secondary_given == index_csi_file){\n                return \"\";\n        }\n        else {\n            if(froot.split('.').pop() == 'bcf'){\n                return \"bcftools index  -f -c \" + froot + \".gz &&\";\n            }\n            else{\n            return \"bcftools index  -f -t \" + froot + \".gz &&\";}\n        }\n    } else {\n        if(fname.split('.').pop() == 'bcf'){\n            var index_csi_file = fname + '.csi'\n            if (in_file.secondaryFiles[0]){\n                var secondary_given = in_file.secondaryFiles[0].path.replace(/^.*[\\\\\\/]/, '');}\n                if(secondary_given == index_csi_file){\n                return \"\";}\n            else{\n                return \"bgzip -c -f \" + fname + \" > \" + fname + \".gz\" + \" && bcftools index -f -c \" + fname + \".gz &&\";}\n            }    \n        return \"bgzip -c -f \" + fname + \" > \" + fname + \".gz\" + \" && bcftools index -f -t \" + fname + \".gz &&\";\n\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "bcftools"
                    },
                    {
                        "shellQuote": false,
                        "position": 2,
                        "valueFrom": "view"
                    },
                    {
                        "prefix": "--output-file",
                        "shellQuote": false,
                        "position": 6,
                        "valueFrom": "${  var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.vcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id}\n\n    if (inputs.output_name) {\n        var out = inputs.output_name\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".subset.vcf\"} \n        else if (inputs.output_type == 'CompressedVCF') {\n            out += \".subset.vcf.gz\"} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".subset.bcf\"} \n        else if (inputs.output_type == 'CompressedBCF') {\n            out += \".subset.bcf.gz\"} \n        else {\n            out += \".subset.vcf\"}\n    }\n    \n    else if (inputs.output_type == 'UncompressedVCF') {\n        var out = froot + '.subset' + '.vcf'} \n    else if (inputs.output_type == 'CompressedVCF') {\n        var out = froot + '.subset' + '.vcf.gz'} \n    else if (inputs.output_type == 'UncompressedBCF') {\n        var out = froot + '.subset' + '.bcf'} \n    else if (inputs.output_type == 'CompressedBCF') {\n        var out = froot + '.subset' + '.bcf.gz'} \n    else var out = froot + '.subset.vcf'\n\n    return out\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 19,
                        "valueFrom": "${\n    if (inputs.min_count && inputs.min_count.minimum_count) {\n        if (inputs.min_count.min_count_type) {\n            return \"--min-ac \" + inputs.min_count.minimum_count + inputs.min_count.min_count_type\n        } else {\n            return \"--min-ac \" + inputs.min_count.minimum_count\n        }\n    } else {\n        return \"\"\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 20,
                        "valueFrom": "${\n    if (inputs.max_count && inputs.max_count.maximum_count) {\n        if (inputs.max_count.max_count_type) {\n            return \"--max-ac \" + inputs.max_count.maximum_count + inputs.max_count.max_count_type\n        } else {\n            return \"--max-ac \" + inputs.max_count.maximum_count\n        }\n    } else {\n        return \"\"\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 31,
                        "valueFrom": "${\n    if (inputs.min_freq && inputs.min_freq.minimum_freq) {\n        if (inputs.min_freq.min_freq_type) {\n            return \"--min-af \" + inputs.min_freq.minimum_freq + inputs.min_freq.min_freq_type\n        } else {\n            return \"--min-af \" + inputs.min_freq.minimum_freq\n        }\n    } else {\n        return \"\"\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 32,
                        "valueFrom": "${\n    if (inputs.max_freq && inputs.max_freq.maximum_freq) {\n        if (inputs.max_freq.max_freq_type) {\n            return \"--max-af \" + inputs.max_freq.maximum_freq + inputs.max_freq.max_freq_type\n        } else {\n            return \"--max-af \" + inputs.max_freq.maximum_freq\n        }\n    } else {\n        return \"\"\n    }\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n    if (inputs.mem_per_job) {\n        return inputs.mem_per_job;\n    } else {\n        return 1000;\n    }\n\n}",
                        "coresMin": "${\n    if (inputs.cpu_per_job) {\n        return inputs.cpu_per_job;\n    } else {\n        return 1;\n    }\n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerImageId": "21caaa02f72e",
                        "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools-1.10.1:0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            "$(inputs.in_variants)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};",
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "successCodes": [
                    0
                ],
                "temporaryFailCodes": [
                    1
                ],
                "sbg:toolkitVersion": "1.10.1",
                "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
                "sbg:categories": [
                    "VCF-Processing"
                ],
                "sbg:links": [
                    {
                        "id": "http://samtools.github.io/bcftools/",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools",
                        "label": "Source code"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/wiki",
                        "label": "Wiki"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/archive/1.9.zip",
                        "label": "Download"
                    }
                ],
                "sbg:cmdPreview": "bcftools index -t -f input_file.some.vcf.gz && bcftools view --output-file test.subset.vcf          input_file.some.vcf.gz",
                "sbg:toolkit": "bcftools",
                "sbg:image_url": null,
                "sbg:license": "MIT License",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1636982060,
                        "sbg:revisionNotes": null
                    }
                ],
                "sbg:appVersion": [
                    "v1.0"
                ],
                "sbg:id": "markoz/hgi/bcftools-view-1-10-1/0",
                "sbg:revision": 0,
                "sbg:revisionNotes": null,
                "sbg:modifiedOn": 1636982060,
                "sbg:modifiedBy": "marko_zecevic",
                "sbg:createdOn": 1636982060,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 0,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a716be8a87e57bbbaab7d6f94b64a574343c39c6937313e6bc9b9ee390752a45d"
            },
            "label": "Bcftools View",
            "sbg:x": -1148.3094482421875,
            "sbg:y": -48.0619010925293
        },
        {
            "id": "snpeff_4_3_cwl1",
            "in": [
                {
                    "id": "input_vcfs",
                    "source": [
                        "bcftools_annotate_cwl1/output_file"
                    ]
                }
            ],
            "out": [
                {
                    "id": "split_vcfs"
                },
                {
                    "id": "joined_vcf"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/snpeff-4-3-cwl1/1",
                "baseCommand": [
                    "java"
                ],
                "inputs": [
                    {
                        "sbg:category": "Input",
                        "id": "input_vcfs",
                        "type": "File[]",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 9
                        },
                        "label": "VCF files to split or join",
                        "doc": "VCF files to split or join together",
                        "sbg:fileTypes": "VCF, VCF.GZ"
                    },
                    {
                        "sbg:toolDefaultValue": "8192",
                        "id": "mem_mb",
                        "type": "int?",
                        "label": "Memory to use for the task [MB]",
                        "doc": "Memory to use for the task [MB]."
                    },
                    {
                        "id": "join_flag",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-j",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Join input VCF files",
                        "doc": "Join input VCF files."
                    },
                    {
                        "id": "num_lines_split",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-l",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Split by <num> lines",
                        "doc": "Split by 'num' lines."
                    }
                ],
                "outputs": [
                    {
                        "id": "split_vcfs",
                        "doc": "Split VCFs.",
                        "label": "Split VCFs",
                        "type": "File[]?",
                        "outputBinding": {
                            "glob": "${\n    return [\"*chrX*\", \"*chrY*\", \"*chr[1-9]*\"]\n}",
                            "outputEval": "${\n    return inheritMetadata(self, inputs.input_vcfs)\n\n}"
                        },
                        "sbg:fileTypes": "VCF"
                    },
                    {
                        "id": "joined_vcf",
                        "doc": "Joined VCF file.",
                        "label": "Joined VCF file",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*joined.vcf"
                        },
                        "sbg:fileTypes": "VCF"
                    }
                ],
                "doc": "**SnpSift split** tool splits VCF files by chromosome (default) or specified number of lines [1].\n\n### Common Use Cases\nSplitting a VCF file for post-processing in parallel.\n\n### Changes Introduced by Seven Bridges\n\nNo significant changes were introduced.\n\n### Common Issues and Important Notes\n\n* Input **Input VCF files to split or join** is required.\n\n### Performance Benchmarking\n\nTypical runs take <5 minutes and cost <$0.05.\n\n### References\n\n[1] [SnpSift split documentation](http://snpeff.sourceforge.net/SnpSift.html#Split)",
                "label": "SnpSift Split 4.3k",
                "arguments": [
                    {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    if (inputs.mem_mb) {\n        return \"-Xmx\" + (inputs.mem_mb).toString() + \"m\"\n    } else {\n        return \"-Xmx8192m\"\n    }\n}"
                    },
                    {
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "-jar"
                    },
                    {
                        "shellQuote": false,
                        "position": 2,
                        "valueFrom": "/opt/snpEff/SnpSift.jar"
                    },
                    {
                        "shellQuote": false,
                        "position": 3,
                        "valueFrom": "split"
                    },
                    {
                        "shellQuote": false,
                        "position": 104,
                        "valueFrom": "> debug.log 2>&1"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 1000,
                        "coresMin": 1
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/jrandjelovic/snpeff:v4.3k"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entry": "$(inputs.input_vcfs)",
                                "writable": false
                            }
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "stdout": "${\n    if (inputs.join_flag) {\n        tempout = ''\n        filename = inputs.input_vcfs[0].path.split('/').pop()\n        tempout = tempout.concat(filename.split('.')[0], '_')\n        ///for (i=0; i < inputs.input_vcfs.length; i++) {\n        ///  filename = inputs.input_vcfs[i].path.split('/').pop()\n        //  tempout = tempout.concat(filename.split('.')[0], '_')\n        ///}\n        tempout = tempout.concat(\"joined.vcf\").replace(/^.*[\\\\\\/]/, '')\n        return tempout\n    }\n}",
                "sbg:image_url": null,
                "sbg:cmdPreview": "java -Xmx1m -jar /opt/snpEff/SnpSift.jar split  /path/to/input_vcf-1.ext  > debug.log 2>&1",
                "sbg:toolkitVersion": "4.3k",
                "sbg:categories": [
                    "VCF-Processing"
                ],
                "sbg:license": "GNU Lesser General Public License v3.0 only",
                "sbg:links": [
                    {
                        "label": "Homepage",
                        "id": "http://snpeff.sourceforge.net/SnpSift.html"
                    },
                    {
                        "label": "Source Code",
                        "id": "https://github.com/pcingola/SnpEff"
                    },
                    {
                        "label": "Documentation",
                        "id": "http://snpeff.sourceforge.net/SnpSift.html"
                    },
                    {
                        "label": "Download",
                        "id": "https://sourceforge.net/projects/snpeff/files/"
                    }
                ],
                "sbg:toolkit": "SnpEff",
                "sbg:toolAuthor": "Pablo Cingolani/Broad Institute",
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "marko_zecevic",
                        "sbg:modifiedOn": 1637105871,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1638920038,
                        "sbg:revisionNotes": ""
                    }
                ],
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/snpeff-4-3-cwl1/1",
                "sbg:revision": 1,
                "sbg:revisionNotes": "",
                "sbg:modifiedOn": 1638920038,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1637105871,
                "sbg:createdBy": "marko_zecevic",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz",
                    "marko_zecevic"
                ],
                "sbg:latestRevision": 1,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a8186b378f218e5990623939255cd95aee9898d57c7e399649687266642584a10"
            },
            "label": "SnpSift Split 4.3k",
            "sbg:x": -752.081787109375,
            "sbg:y": -67.63298034667969
        },
        {
            "id": "remap",
            "in": [
                {
                    "id": "results",
                    "source": "sbg_merge_and_filter_genesis_results_cwl2/merged"
                },
                {
                    "id": "chain",
                    "source": "chain"
                },
                {
                    "id": "custom_input",
                    "source": "custom_input_1"
                }
            ],
            "out": [
                {
                    "id": "hg19_out"
                },
                {
                    "id": "hg38_out"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "markoz/hgi/remap/8",
                "baseCommand": [
                    "Rscript",
                    "remap.R"
                ],
                "inputs": [
                    {
                        "id": "results",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--results=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 1
                        }
                    },
                    {
                        "id": "chain",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--chain=",
                            "separate": false,
                            "shellQuote": false,
                            "position": 2
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "hg19_out",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*hg19.csv"
                        },
                        "sbg:fileTypes": "CSV"
                    },
                    {
                        "id": "hg38_out",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*hg38.csv"
                        },
                        "sbg:fileTypes": "CSV"
                    }
                ],
                "label": "Remap",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/marko_zecevic/rtracklayer:1.54.0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            {
                                "entryname": "remap.R",
                                "entry": "## Collect arguments\nargs <- commandArgs(TRUE)\n\n## Parse arguments (we expect the form --arg=value)\nparseArgs <- function(x) strsplit(sub(\"^--\", \"\", x), \"=\")\nargsDF <- as.data.frame(do.call(\"rbind\", parseArgs(args)))\narg <- as.list(as.character(argsDF$V2))\nnames(arg) <- argsDF$V1\n\nres <- read.csv(arg$results)\nres <- res[order(res$chr, res$SPA.pval),]\n\nfile_prefix <- tools::file_path_sans_ext(basename(arg$results))\n\n\nlibrary(rtracklayer)\n# download chain file here: http://hgdownload.soe.ucsc.edu/downloads.html >> Feb. 2009 (GRCh37/hg19) >> LiftOver files >> hg19ToHg38.over.chain.gz\n\n# specify coordinates to liftover\ngrObject <- GRanges(seqnames=paste0(\"chr\", res$chr), ranges=IRanges(start=res$pos, end=res$pos))\n\n# import the chain file\nchainObject <- import.chain(arg$chain)\n\n# run liftOver\nresults <- as.data.frame(liftOver(grObject, chainObject))\n# res <- res[,c(\"chr\", \"pos\", \"n.obs\", \"SPA.pval\", \"Est\", \"Est.SE\")] # these columns are required for running FUMA\nwrite.csv(res, paste0(file_prefix, \"_hg19.csv\"), quote = FALSE)\n\nres$pos <- NA\nres$pos[results$group] <- results$start\n\nwrite.csv(res, paste0(file_prefix, \"_hg38.csv\"), quote = FALSE)",
                                "writable": false
                            }
                        ]
                    }
                ],
                "sbg:projectName": "HGI",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1644679024,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1644680496,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1644702579,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1644702910,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1644741167,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1645379758,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1647787819,
                        "sbg:revisionNotes": "do not output all columns"
                    },
                    {
                        "sbg:revision": 7,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1648132941,
                        "sbg:revisionNotes": "added Est and Est.SE"
                    },
                    {
                        "sbg:revision": 8,
                        "sbg:modifiedBy": "markoz",
                        "sbg:modifiedOn": 1653663478,
                        "sbg:revisionNotes": "output all columns"
                    }
                ],
                "sbg:image_url": null,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "markoz/hgi/remap/8",
                "sbg:revision": 8,
                "sbg:revisionNotes": "output all columns",
                "sbg:modifiedOn": 1653663478,
                "sbg:modifiedBy": "markoz",
                "sbg:createdOn": 1644679024,
                "sbg:createdBy": "markoz",
                "sbg:project": "markoz/hgi",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "markoz"
                ],
                "sbg:latestRevision": 8,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a57b77bfd153a8c0490749909959c8fe393a8f2291f07ecfcc1a8642f9dd357a8",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "Remap",
            "sbg:x": 2164,
            "sbg:y": -19.246591567993164,
            "when": "$(inputs.custom_input)"
        }
    ],
    "hints": [
        {
            "class": "sbg:maxNumberOfParallelInstances",
            "value": "5"
        }
    ],
    "requirements": [
        {
            "class": "SubworkflowFeatureRequirement"
        },
        {
            "class": "ScatterFeatureRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "sbg:projectName": "HGI",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637105387,
            "sbg:revisionNotes": null
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637105400,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 2,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1637105879,
            "sbg:revisionNotes": "Workflow decomposed"
        },
        {
            "sbg:revision": 3,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637106191,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 4,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637180569,
            "sbg:revisionNotes": "added BCFtools view"
        },
        {
            "sbg:revision": 5,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637330628,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 6,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637528748,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 7,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637528788,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 8,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637533934,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 9,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637537447,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 10,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637538327,
            "sbg:revisionNotes": "test"
        },
        {
            "sbg:revision": 11,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637538379,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 12,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637541089,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 13,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637658586,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 14,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637658656,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 15,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637660472,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 16,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637683342,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 17,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637690714,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 18,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637690730,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 19,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637874822,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 20,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637880749,
            "sbg:revisionNotes": "lowered the RAM for sliding window testing"
        },
        {
            "sbg:revision": 21,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637881427,
            "sbg:revisionNotes": "Sliding window teting made a conditional step"
        },
        {
            "sbg:revision": 22,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637914065,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 23,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638909303,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 24,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638920099,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 25,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638920205,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 26,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638920306,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 27,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638974597,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 28,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638976411,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 29,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638993995,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 30,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1638999485,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 31,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639043418,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 32,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639258968,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 33,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639912769,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 34,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639934716,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 35,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1640203255,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 36,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1641862228,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 37,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1643552684,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 38,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1643568144,
            "sbg:revisionNotes": "hwe thresholds exposed"
        },
        {
            "sbg:revision": 39,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1644760951,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 40,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1644795034,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 41,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1644870813,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 42,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645099615,
            "sbg:revisionNotes": "conditional variant file exposed"
        },
        {
            "sbg:revision": 43,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645100677,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 44,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645167494,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 45,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645168276,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 46,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645179610,
            "sbg:revisionNotes": "exposed single snp test inputs"
        },
        {
            "sbg:revision": 47,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645195136,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 48,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645311371,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 49,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645311393,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 50,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1645379801,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 51,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1646834037,
            "sbg:revisionNotes": "Removed sliding window testing"
        },
        {
            "sbg:revision": 52,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647079104,
            "sbg:revisionNotes": "lz update1"
        },
        {
            "sbg:revision": 53,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647085875,
            "sbg:revisionNotes": "locuszoom edit2"
        },
        {
            "sbg:revision": 54,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647087399,
            "sbg:revisionNotes": "locuszoom update3"
        },
        {
            "sbg:revision": 55,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647089016,
            "sbg:revisionNotes": "locuszoom update 4"
        },
        {
            "sbg:revision": 56,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647090777,
            "sbg:revisionNotes": "locuszoom update5"
        },
        {
            "sbg:revision": 57,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647092472,
            "sbg:revisionNotes": "reverted to original locuszoom"
        },
        {
            "sbg:revision": 58,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647108098,
            "sbg:revisionNotes": "locuszoom update 11"
        },
        {
            "sbg:revision": 59,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647162219,
            "sbg:revisionNotes": "locuszoom update"
        },
        {
            "sbg:revision": 60,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647167756,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 61,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647174693,
            "sbg:revisionNotes": "locuszoom edit"
        },
        {
            "sbg:revision": 62,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647178967,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 63,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647284579,
            "sbg:revisionNotes": "locuszoom update"
        },
        {
            "sbg:revision": 64,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647287209,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 65,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1647787853,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 66,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648133108,
            "sbg:revisionNotes": "output est and est.se for FUMA"
        },
        {
            "sbg:revision": 67,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648482232,
            "sbg:revisionNotes": "rfrows=6"
        },
        {
            "sbg:revision": 68,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648492499,
            "sbg:revisionNotes": "locuszoom updated"
        },
        {
            "sbg:revision": 69,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648706093,
            "sbg:revisionNotes": "ld pruning outputs exposed"
        },
        {
            "sbg:revision": 70,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648761069,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 71,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648763969,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 72,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648766028,
            "sbg:revisionNotes": "10 11"
        },
        {
            "sbg:revision": 73,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648835506,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 74,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1648835649,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 75,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1653663521,
            "sbg:revisionNotes": "output all columns of the summary file"
        },
        {
            "sbg:revision": 76,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1653841629,
            "sbg:revisionNotes": "remove kinship matrix connection between PC-Relate and GENESIS Null Model"
        },
        {
            "sbg:revision": 77,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1653845763,
            "sbg:revisionNotes": "score.spa -> score"
        },
        {
            "sbg:revision": 78,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1653846846,
            "sbg:revisionNotes": "nothing after single variant association testing"
        },
        {
            "sbg:revision": 79,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1653923099,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 80,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1654734016,
            "sbg:revisionNotes": "manhattan back to original state"
        },
        {
            "sbg:revision": 81,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1654734053,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 82,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1654739983,
            "sbg:revisionNotes": "no GRM in null model fitting"
        },
        {
            "sbg:revision": 83,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1654767167,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 84,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1657779021,
            "sbg:revisionNotes": ""
        }
    ],
    "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/genesis-pipeline/84.png",
    "sbg:appVersion": [
        "v1.2",
        "v1.0",
        "v1.1"
    ],
    "id": "https://cgc-api.sbgenomics.com/v2/apps/markoz/hgi/genesis-pipeline/84/raw/",
    "sbg:id": "markoz/hgi/genesis-pipeline/84",
    "sbg:revision": 84,
    "sbg:revisionNotes": "",
    "sbg:modifiedOn": 1657779021,
    "sbg:modifiedBy": "markoz",
    "sbg:createdOn": 1637105387,
    "sbg:createdBy": "markoz",
    "sbg:project": "markoz/hgi",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:contributors": [
        "marko_zecevic",
        "markoz"
    ],
    "sbg:latestRevision": 84,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "aa1e78153e28a0c14021c45807a9f955ea0db989604795f870b04d6d0542191d0",
    "sbg:workflowLanguage": "CWL"
}