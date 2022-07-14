{
  "class": "Workflow",
  "cwlVersion": "v1.2",
  "label": "Infinium array genotyping",
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "inputs": [
    {
      "id": "idat_archive",
      "sbg:fileTypes": "ZIP",
      "type": "File?",
      "sbg:x": -1137.9208984375,
      "sbg:y": -391.8548889160156
    },
    {
      "id": "egt_cluster_file",
      "sbg:fileTypes": "EGT",
      "type": "File?",
      "sbg:x": -1193.1630859375,
      "sbg:y": -543.225830078125
    },
    {
      "id": "bpm_manifest_file",
      "sbg:fileTypes": "BPM",
      "type": "File?",
      "sbg:x": -1045.5985107421875,
      "sbg:y": -290
    },
    {
      "id": "reference_fasta",
      "sbg:fileTypes": "FASTA, FA",
      "type": "File?",
      "sbg:x": -1023.0502319335938,
      "sbg:y": -643.8870849609375
    },
    {
      "id": "in_map",
      "sbg:fileTypes": "GMAP.GZ",
      "type": "File[]",
      "label": "Genetic map",
      "doc": "Genetic map.",
      "sbg:x": 1159.4847412109375,
      "sbg:y": -502.92156982421875
    },
    {
      "id": "output_name",
      "type": "string?",
      "label": "Output file name",
      "doc": "Name of the output file.",
      "sbg:exposed": true
    },
    {
      "id": "out_name",
      "type": "string?",
      "sbg:exposed": true
    },
    {
      "id": "in_refhaps",
      "sbg:fileTypes": "M3VCF, M3VCF.GZ",
      "type": "File[]",
      "label": "Reference haplotypes",
      "doc": "M3VCF file containing haplotype data for reference panel.",
      "sbg:x": 1704.4554443359375,
      "sbg:y": -660.6348266601562
    },
    {
      "id": "samples_file_1",
      "sbg:fileTypes": "TXT",
      "type": "File?",
      "label": "New sample names",
      "doc": "New sample names, one name per line, in the same order as they appear in the VCF file. Alternatively, only samples which need to be renamed can be listed as \"old_name new_name\\n\" pairs separated by whitespaces, each on a separate line. If a sample name contains spaces, the spaces can be escaped using the backslash character, for example \"Not\\ a\\ good\\ sample\\ name\".",
      "sbg:x": -268.7236633300781,
      "sbg:y": -275.28900146484375
    },
    {
      "id": "output_name_1",
      "type": "string?",
      "label": "Output file name",
      "doc": "Name of the output file.",
      "sbg:exposed": true
    },
    {
      "id": "max-missing",
      "type": "float",
      "label": "Threshold variant call rate",
      "doc": "Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).",
      "sbg:exposed": true
    },
    {
      "id": "keep_ids_file",
      "type": "File?",
      "label": "Keep variants",
      "doc": "File containing IDs of variants to keep",
      "sbg:x": -92.08368682861328,
      "sbg:y": -523.9205322265625
    },
    {
      "id": "include_expression",
      "type": "string?",
      "label": "Include expression",
      "doc": "Include only sites for which the expression is true.",
      "sbg:exposed": true
    },
    {
      "id": "chr",
      "type": "string[]?",
      "label": "Chromosomes to analyze",
      "doc": "Chromosomes to analyze.",
      "sbg:exposed": true
    },
    {
      "id": "hwe",
      "type": "float",
      "label": "Threshold p-value for HW equilibrium",
      "doc": "p-value for elimination via the Hardy-Weinberg Equilibrium filtering.",
      "sbg:exposed": true
    },
    {
      "id": "not_chr",
      "type": "string[]?",
      "label": "Chromosomes to ommit",
      "doc": "Chromosomes to omit.",
      "sbg:exposed": true
    },
    {
      "id": "samples_file",
      "sbg:fileTypes": "TXT",
      "type": "File?",
      "label": "Keep samples",
      "doc": "File of samples to include (or exclude with \"^\" prefix).",
      "sbg:x": -26.66947364807129,
      "sbg:y": -654.259521484375
    },
    {
      "id": "include_expression_1",
      "type": "string?",
      "label": "Include expression",
      "doc": "Include only sites for which the expression is true.",
      "sbg:exposed": true
    }
  ],
  "outputs": [
    {
      "id": "output_file",
      "outputSource": [
        "bcftools_view_1/out_variants"
      ],
      "sbg:fileTypes": "VCF, VCF.GZ",
      "type": "File?",
      "label": "Final VCF",
      "sbg:x": 2519.125732421875,
      "sbg:y": -575.370849609375
    },
    {
      "id": "het_file",
      "outputSource": [
        "vcftools_het/het_file"
      ],
      "sbg:fileTypes": "HET",
      "type": "File?",
      "label": "HET file pre",
      "doc": "HET file",
      "sbg:x": 1625.43701171875,
      "sbg:y": -289.9977111816406
    }
  ],
  "steps": [
    {
      "id": "iaap_cli_gencall_1_1",
      "in": [
        {
          "id": "bpm_manifest_file",
          "source": "bpm_manifest_file"
        },
        {
          "id": "egt_cluster_file",
          "source": "egt_cluster_file"
        },
        {
          "id": "idat_archive",
          "source": "idat_archive"
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
        "id": "markoz/hgi/iaap-cli-gencall-1-1/0",
        "baseCommand": [
          "unzip"
        ],
        "inputs": [
          {
            "id": "bpm_manifest_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 4
            },
            "sbg:fileTypes": "BPM"
          },
          {
            "id": "egt_cluster_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 5
            },
            "sbg:fileTypes": "EGT"
          },
          {
            "id": "idat_archive",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            },
            "sbg:fileTypes": "ZIP"
          }
        ],
        "outputs": [
          {
            "id": "output",
            "type": "File?",
            "outputBinding": {
              "glob": "*.zip",
              "outputEval": "$(inheritMetadata(self, inputs.idat_archive))"
            }
          }
        ],
        "label": "iaap-cli gencall 1.1",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 6,
            "valueFrom": "output_folder"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 8,
            "valueFrom": "--output-gtc --gender-estimate-call-rate-threshold -0.1"
          },
          {
            "prefix": "--idat-folder",
            "shellQuote": false,
            "position": 7,
            "valueFrom": "idat_files"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 3,
            "valueFrom": "-d idat_files && LANG=\"en_US.UTF-8\" && /bin/iaap-cli/iaap-cli gencall"
          },
          {
            "prefix": "> out.log && cd output_folder; zip -r",
            "shellQuote": false,
            "position": 10,
            "valueFrom": "${\n    return \"../\" + inputs.idat_archive.metadata.sample_id + \".zip *\"\n\n}"
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
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/marko_zecevic/gencall:1.1"
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982049,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/iaap-cli-gencall-1-1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982049,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982049,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a60586dd7b1ac27186e410ac4dea413098c015f92ff87f500db4a9171246275f8"
      },
      "label": "iaap-cli gencall 1.1",
      "sbg:x": -825.8085327148438,
      "sbg:y": -399.1614990234375
    },
    {
      "id": "bcftools_gtc2vcf_1_10",
      "in": [
        {
          "id": "bpm_manifest_file",
          "source": "bpm_manifest_file"
        },
        {
          "id": "egt_cluster_file",
          "source": "egt_cluster_file"
        },
        {
          "id": "reference_fasta",
          "source": "reference_fasta"
        },
        {
          "id": "out_name",
          "source": "out_name"
        },
        {
          "id": "gtcs",
          "source": "iaap_cli_gencall_1_1/output"
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
        "id": "markoz/hgi/bcftools-gtc2vcf-1-10/0",
        "baseCommand": [
          "unzip"
        ],
        "inputs": [
          {
            "id": "bpm_manifest_file",
            "type": "File?",
            "inputBinding": {
              "prefix": "--bpm",
              "shellQuote": false,
              "position": 13
            },
            "sbg:fileTypes": "BPM"
          },
          {
            "id": "egt_cluster_file",
            "type": "File?",
            "inputBinding": {
              "prefix": "--egt",
              "shellQuote": false,
              "position": 14
            },
            "sbg:fileTypes": "EGT"
          },
          {
            "id": "reference_fasta",
            "type": "File?",
            "inputBinding": {
              "prefix": "--fasta-ref",
              "shellQuote": false,
              "position": 16
            },
            "sbg:fileTypes": "FASTA, FA"
          },
          {
            "id": "out_name",
            "type": "string?",
            "inputBinding": {
              "prefix": "-o",
              "shellQuote": false,
              "position": 20,
              "valueFrom": "${\n    return inputs.out_name + \".vcf\"\n}"
            }
          },
          {
            "id": "gtcs",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            },
            "sbg:fileTypes": "ZIP"
          }
        ],
        "outputs": [
          {
            "id": "output",
            "type": "File?",
            "outputBinding": {
              "glob": "*.vcf",
              "outputEval": "$(inheritMetadata(self, inputs.gtcs))"
            }
          }
        ],
        "label": "bcftools gtc2vcf 1.10",
        "arguments": [
          {
            "prefix": "--gtcs",
            "shellQuote": false,
            "position": 15,
            "valueFrom": "gtcs/"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 4,
            "valueFrom": "-d . && mkdir gtcs && cp -r */*.gtc gtcs/ && bcftools +gtc2vcf --no-version -Ov"
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
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/marko_zecevic/gtc2vcf:1.1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entry": "$(inputs.gtcs)",
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
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982050,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/bcftools-gtc2vcf-1-10/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982050,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982050,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a24c2a83021840bbbfeaab41549d96d116c9cbbb0520b7162778e00aab78ab7e6"
      },
      "label": "bcftools gtc2vcf 1.10",
      "sbg:x": -660.80908203125,
      "sbg:y": -478.5002746582031
    },
    {
      "id": "tabix_bgzip_1_9_cwl1_0",
      "in": [
        {
          "id": "input_file",
          "source": "ebi_vcf_debugulator_cwl1/corrected_vcf_file"
        }
      ],
      "out": [
        {
          "id": "output_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/tabix-bgzip-1-9-cwl1-0/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:category": "File inputs",
            "id": "input_file",
            "type": "File",
            "label": "Input file",
            "doc": "Input file to be compressed/decompressed.",
            "sbg:fileTypes": "VCF, VCF.GZ, GFF, GFF.GZ, BED, BED.GZ, SAM, SAM.GZ, PSLTAB, PSLTAB.GZ"
          },
          {
            "sbg:category": "Config Inputs",
            "sbg:toolDefaultValue": "False",
            "id": "decompress",
            "type": "boolean?",
            "label": "Decompress input file",
            "doc": "Decompress input file."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "2048",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job",
            "doc": "Memory per job."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "1",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "CPU per job",
            "doc": "CPU per job."
          },
          {
            "sbg:category": "Config Inputs",
            "sbg:toolDefaultValue": "-1",
            "id": "compress_level",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "-1",
                  "0",
                  "1",
                  "2",
                  "3",
                  "4",
                  "5",
                  "6",
                  "7",
                  "8",
                  "9"
                ],
                "name": "compress_level"
              }
            ],
            "inputBinding": {
              "prefix": "",
              "shellQuote": false,
              "position": 0,
              "valueFrom": "${\n    if(inputs.compress_level != undefined && (inputs.decompress==undefined || inputs.decompress == false))\n        return ' -I ' + inputs.compress_level\n    \n    return ''\n}"
            },
            "label": "Compress level",
            "doc": "Compression level to use when compressing; 0 to 9, or -1 for default."
          },
          {
            "sbg:category": "Config Inputs",
            "sbg:toolDefaultValue": "1",
            "id": "threads",
            "type": "int?",
            "inputBinding": {
              "prefix": "-@",
              "shellQuote": false,
              "position": 0
            },
            "label": "Number of threads to be used",
            "doc": "Number of threads to be used."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "doc": "Compressed or decompressed file.",
            "label": "Compressed or decompressed file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    filename = ''\n    if (inputs.input_file instanceof Array)\n        filename = inputs.input_file[0].path\n    else\n        filename = inputs.input_file.path\n\n    paths = filename.split('/')\n    names = filename.split('/')[paths.length - 1].split('.')\n\n    nn = ''\n    if (inputs.decompress) {\n        lind = names.length - 1\n    } else if ((filename.charAt(filename.length - 3) + filename.charAt(filename.length - 2) + filename.charAt(filename.length - 1)) ==\n        '.gz' && (inputs.decompress == false || inputs.decompress == undefined))\n\n    {\n        lind = names.length - 1\n    } else {\n        lind = names.length\n    }\n    for (i = 0; i < lind; i++) {\n        if (i != 0)\n            nn += '.'\n        nn += names[i]\n        if ((filename.charAt(filename.length - 3) + filename.charAt(filename.length - 2) + filename.charAt(filename.length - 1)) ==\n            '.gz' && (inputs.decompress == false || inputs.decompress == undefined) && i == 0)\n            nn += '.tab'\n        else if ((inputs.decompress == false || inputs.decompress == undefined) && i == 0 && inputs.suffix_append == true)\n            nn += '.tab'\n    }\n\n    if (inputs.decompress) {\n        return nn\n    } else {\n        return nn + '.gz'\n    }\n}",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
            },
            "sbg:fileTypes": "VCF.GZ, VCF, BED.GZ, BED, GFF.GZ, GFF, SAM.GZ, SAM, PSLTAB.GZ, PSLTAB"
          }
        ],
        "doc": "**Tabix BGZIP** is used for compressing/decompressing (BAM, VCF, BED, ...) any file in BGZF and from BGZF format.\n\nA list of all inputs and parameters with corresponding descriptions can be found at the bottom of this page.\n\n###Common Use Cases\n\n**Tabix BGZIP** is used in cases where tools expect inputs which are \u2018BGZF\u2019 formatted. In most of these tools, it is expected that \u2018BGZF\u2019 file is indexed, so **Tabix Index** is called after **Tabix BGZIP**.\n\nThere are three modes for **Tabix BGZIP**:\n 1. it will compress the input if the **decompress** parameter is not set at all or set to *False*,\n 2. it will decompress the input if the **decompress** parameter is set to *True*,\n 3. if the input file has suffix \u2018.gz\u2019 and **decompress** parameter is not set at all or set to *False*, it will decompress the input using the gzip command and then compress it using **Tabix BGZIP**.\n\n###Changes Introduced by Seven Bridges\n**Tabix BGZIP** is extended so that it will work with the given compressed input and the **decompress** parameter not set or set to *False* in the way described in the **Common Use Cases** section, case 3.\nAlso, the wrapper around the tool doesn\u2019t support following parameters:  '-- offset', '--help', '--stdout', '--size' as the tool used with these options outputs the result to stdout. Options '--index', '--index-name', '--reindex', '--rebgzip' are not part of the wrapper as these functions are redundant either with wrapped options of **Tabix BGZIP** or they can be achieved using **Tabix Index**. Option '--force' isn\u2019t a part of the wrapper as it is used by default on the command line.\n\n###Common Issues and Important Notes\nThere aren't any common issues.\n\n###Performance Benchmarking\n**Tabix BGZIP** is not CPU/Memory intensive. The default c4.2 AWS instance can be used.\nCost can be significantly reduced by using spot instances. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.",
        "label": "Tabix BGZIP CWL1.0",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    com = ''\n    filename = ''\n    if (inputs.input_file instanceof Array)\n        filename = inputs.input_file[0].path\n    else\n        filename = inputs.input_file.path\n\n    paths = filename.split('/')\n    names = filename.split('/')[paths.length - 1].split('.')\n\n    name = filename\n    if ((filename.charAt(filename.length - 3) + filename.charAt(filename.length - 2) + filename.charAt(filename.length - 1)) ==\n        '.gz' && (inputs.decompress == false || inputs.decompress == undefined)) {\n        com += 'gzip -d -c '\n        com += filename\n        com += '>'\n        com += filename.split('/')[paths.length - 1].substring(0, (filename.split('/')[paths.length - 1].length) - 3)\n        com += ' ; '\n        name = filename.split('/')[paths.length - 1].substring(0, (filename.split('/')[paths.length - 1].length) - 3)\n    }\n    com += '/opt/htslib-1.9/bgzip '\n    com += ' -c -f '\n    if (inputs.decompress == true) {\n        com += '-d '\n    }\n    com += name\n    return com\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n  if (inputs.mem_per_job)\n  {\n    return inputs.mem_per_job\n  }\n  else\n  {\n    return 2048\n  }\n}",
            "coresMin": "${\n  if (inputs.cpu_per_job)\n  {\n    return inputs.cpu_per_job\n  }\n  else\n  {\n    return 1\n  }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/nevenam/htslib-1-9:0"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};"
            ]
          }
        ],
        "stdout": "${\n    filename = ''\n    if (inputs.input_file instanceof Array)\n        filename = inputs.input_file[0].path\n    else\n        filename = inputs.input_file.path\n\n    paths = filename.split('/')\n    names = filename.split('/')[paths.length - 1].split('.')\n\n    nn = ''\n    if (inputs.decompress) {\n        lind = names.length - 1\n    } else if ((filename.charAt(filename.length - 3) + filename.charAt(filename.length - 2) + filename.charAt(filename.length - 1)) ==\n        '.gz' && (inputs.decompress == false || inputs.decompress == undefined))\n\n    {\n        lind = names.length - 1\n    } else {\n        lind = names.length\n    }\n    for (i = 0; i < lind; i++) {\n        if (i != 0)\n            nn += '.'\n        nn += names[i]\n        if ((filename.charAt(filename.length - 3) + filename.charAt(filename.length - 2) + filename.charAt(filename.length - 1)) ==\n            '.gz' && (inputs.decompress == false || inputs.decompress == undefined) && i == 0)\n            nn += '.tab'\n        else if ((inputs.decompress == false || inputs.decompress == undefined) && i == 0 && inputs.suffix_append == true)\n            nn += '.tab'\n    }\n\n    if (inputs.decompress) {\n        return nn\n    } else {\n        return nn + '.gz'\n    }\n}",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982054,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "/opt/samtools-1.3/tabix-0.2.6/bgzip  -c -f /path/to/input_vcf_file.vcf  > input_vcf_file.vcf.gz",
        "sbg:license": "The MIT/Expat License",
        "sbg:links": [
          {
            "id": "http://www.htslib.org/",
            "label": "Homepage"
          },
          {
            "id": "https://github.com/samtools/htslib/tree/master",
            "label": "Source Code"
          },
          {
            "id": "http://www.htslib.org/doc/#manual-pages",
            "label": "Wiki"
          },
          {
            "id": "http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.5.tar.bz2/download",
            "label": "Download"
          },
          {
            "id": "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/",
            "label": "Publication"
          },
          {
            "id": "http://www.htslib.org/doc/#manual-pages",
            "label": "Documentation"
          }
        ],
        "sbg:toolAuthor": "Heng Li -  Broad Institue",
        "sbg:toolkit": "Tabix",
        "sbg:toolkitVersion": "1.9.0",
        "sbg:projectName": "HGI",
        "sbg:categories": [
          "CWL1.0",
          "Utilities",
          "VCF Processing"
        ],
        "sbg:appVersion": [
          "v1.0"
        ],
        "sbg:id": "markoz/hgi/tabix-bgzip-1-9-cwl1-0/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982054,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982054,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "ae81d79dd2b45e4301dad6484a4f7160faf83f9ca36b5cb4f756a8867cbe6969b"
      },
      "label": "Tabix BGZIP CWL1.0",
      "sbg:x": 967.47412109375,
      "sbg:y": -392.68658447265625
    },
    {
      "id": "shapeit_4_4_2_1_cwl1_1",
      "in": [
        {
          "id": "in_genotypes",
          "source": "bcftools_index_cwl1/output_file"
        },
        {
          "id": "in_map",
          "source": "in_map"
        },
        {
          "id": "suffix",
          "default": "vcf.gz"
        }
      ],
      "out": [
        {
          "id": "out_phased_haplotypes"
        },
        {
          "id": "out_bin_phased_haplotypes"
        },
        {
          "id": "out_log"
        },
        {
          "id": "out_phased_haplotypes_par"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/shapeit-4-4-2-1-cwl1-1/0",
        "baseCommand": [
          "/opt/shapeit4-4.2.1/bin/shapeit4.2"
        ],
        "inputs": [
          {
            "sbg:category": "Basic Options",
            "sbg:toolDefaultValue": "15052011",
            "id": "seed",
            "type": "int?",
            "inputBinding": {
              "prefix": "--seed",
              "shellQuote": false,
              "position": 1
            },
            "label": "Random number generator seed",
            "doc": "Seed of the random number generator. Default: 15052011."
          },
          {
            "sbg:altPrefix": "-T",
            "sbg:category": "Basic Options",
            "sbg:toolDefaultValue": "1",
            "id": "thread",
            "type": "int?",
            "inputBinding": {
              "prefix": "--thread",
              "shellQuote": false,
              "position": 2
            },
            "label": "Number of threads",
            "doc": "Number of threads used. Default: 1."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "1",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "CPUs per job",
            "doc": "CPUs per job. Default: 1."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "1000",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job [MB]",
            "doc": "Memory per job [MB]."
          },
          {
            "sbg:altPrefix": "-I",
            "sbg:category": "Input Files",
            "id": "in_genotypes",
            "type": "File",
            "inputBinding": {
              "prefix": "--input",
              "shellQuote": false,
              "position": 3
            },
            "label": "Genotypes to phase",
            "doc": "Genotypes to be phased in VCF/BCF format. The file should be indexed (bcftools index unphased.vcf.gz).",
            "sbg:fileTypes": "VCF.GZ, BCF.GZ, BCF, VCF",
            "secondaryFiles": [
              {
                "pattern": ".csi",
                "required": false
              },
              {
                "pattern": ".tbi",
                "required": false
              }
            ]
          },
          {
            "sbg:altPrefix": "-H",
            "sbg:category": "Input Files",
            "id": "in_reference_panel",
            "type": "File?",
            "inputBinding": {
              "prefix": "--reference",
              "shellQuote": false,
              "position": 4
            },
            "label": "Reference panel",
            "doc": "Reference panel of haplotypes in  VCF/BCF format. The file should be indexed (bcftools index).",
            "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ",
            "secondaryFiles": [
              {
                "pattern": ".csi",
                "required": false
              },
              {
                "pattern": ".tbi",
                "required": false
              }
            ]
          },
          {
            "sbg:altPrefix": "-S",
            "sbg:category": "Input Files",
            "id": "in_scaffold",
            "type": "File?",
            "inputBinding": {
              "prefix": "--scaffold",
              "shellQuote": false,
              "position": 5
            },
            "label": "Scaffold of haplotypes",
            "doc": "Scaffold of haplotypes in VCF/BCF format.",
            "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ",
            "secondaryFiles": [
              {
                "pattern": ".csi",
                "required": false
              },
              {
                "pattern": ".tbi",
                "required": false
              }
            ]
          },
          {
            "sbg:altPrefix": "-M",
            "sbg:category": "Input Files",
            "id": "in_map",
            "type": "File",
            "inputBinding": {
              "prefix": "--map",
              "shellQuote": false,
              "position": 6
            },
            "label": "Genetic map",
            "doc": "Genetic map.",
            "sbg:fileTypes": "GMAP.GZ"
          },
          {
            "sbg:category": "Config Inputs",
            "id": "use_ps",
            "type": "string?",
            "inputBinding": {
              "prefix": "--use-PS",
              "shellQuote": false,
              "position": 8
            },
            "label": "PS field to use",
            "doc": "Informs phasing using PS field from read based phasing."
          },
          {
            "sbg:category": "Config Inputs",
            "sbg:toolDefaultValue": "False",
            "id": "sequencing",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--sequencing",
              "shellQuote": false,
              "position": 9
            },
            "label": "Parameter settings for sequencing data",
            "doc": "Default parameter setting for sequencing data (this divides by 50 the default value of --pbwt-modulo). Default: False."
          },
          {
            "sbg:toolDefaultValue": "5b,1p,1b,1p,1b,1p,5m",
            "sbg:category": "MCMC Parameters",
            "id": "mcmc_iterations",
            "type": "string?",
            "inputBinding": {
              "prefix": "--mcmc-iterations",
              "shellQuote": false,
              "position": 10
            },
            "label": "Iteration scheme of the MCMC",
            "doc": "Iteration scheme of the MCMC. Default: 5b,1p,1b,1p,1b,1p,5m."
          },
          {
            "sbg:category": "MCMC Parameters",
            "sbg:toolDefaultValue": "0.999",
            "id": "mcmc_prune",
            "type": "float?",
            "inputBinding": {
              "prefix": "--mcmc-prune",
              "shellQuote": false,
              "position": 11
            },
            "label": "Pruning threshold for genotype graphs",
            "doc": "Pruning threshold for genotype graphs. Default: 0.999."
          },
          {
            "sbg:category": "PBWT Parameters",
            "sbg:toolDefaultValue": "0.02",
            "id": "pbwt_modulo",
            "type": "float?",
            "inputBinding": {
              "prefix": "--pbwt-modulo",
              "shellQuote": false,
              "position": 12
            },
            "label": "PBWT modulo",
            "doc": "Storage frequency of PBWT indexes in cM (i.e. storage every 0.02 cM by default)."
          },
          {
            "sbg:category": "PBWT Parameters",
            "sbg:toolDefaultValue": "4",
            "id": "pbwt_depth",
            "type": "int?",
            "inputBinding": {
              "prefix": "--pbwt-depth",
              "shellQuote": false,
              "position": 13
            },
            "label": "PBWT depth",
            "doc": "Depth of PBWT indexes to condition on. Default: 4."
          },
          {
            "sbg:category": "PBWT Parameters",
            "sbg:toolDefaultValue": "2",
            "id": "pbwt_mac",
            "type": "int?",
            "inputBinding": {
              "prefix": "--pbwt-mac",
              "shellQuote": false,
              "position": 14
            },
            "label": "PBWT MAC",
            "doc": "Minimal Minor Allele Count at which PBWT is evaluated. Default: 2."
          },
          {
            "sbg:category": "PBWT Parameters",
            "sbg:toolDefaultValue": "0.5",
            "id": "pbwt_mdr",
            "type": "float?",
            "inputBinding": {
              "prefix": "--pbwt-mdr",
              "shellQuote": false,
              "position": 15
            },
            "label": "PBWT MDR",
            "doc": "Maximal Missing Data Rate at which PBWT is evaluated. Default: 0.5."
          },
          {
            "sbg:category": "PBWT Parameters",
            "id": "pbwt_disable_init",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--pbwt-disable-init",
              "shellQuote": false,
              "position": 16
            },
            "label": "Disable initialization by PBWT sweep",
            "doc": "Disable initialization by PBWT sweep."
          },
          {
            "sbg:altPrefix": "-W",
            "sbg:category": "HMM Parameters",
            "sbg:toolDefaultValue": "2.5",
            "id": "window",
            "type": "float?",
            "inputBinding": {
              "prefix": "--window",
              "shellQuote": false,
              "position": 17
            },
            "label": "Minimal size of the phasing window in cM",
            "doc": "Minimal size of the phasing window in cM. Default: 2.5."
          },
          {
            "sbg:category": "HMM Parameters",
            "sbg:toolDefaultValue": "15000",
            "id": "effective_size",
            "type": "int?",
            "inputBinding": {
              "prefix": "--effective-size",
              "shellQuote": false,
              "position": 18
            },
            "label": "Effective size of the population",
            "doc": "Effective size of the population. Default: 15000."
          },
          {
            "sbg:category": "Output Files",
            "sbg:altPrefix": "-O",
            "id": "prefix",
            "type": "string?",
            "label": "Output file name prefix",
            "doc": "Output file name for the phased haplotypes."
          },
          {
            "sbg:category": "Config Inputs",
            "id": "suffix",
            "type": {
              "type": "enum",
              "symbols": [
                "vcf",
                "vcf.gz",
                "bcf",
                "bcf.gz"
              ],
              "name": "suffix"
            },
            "label": "Phased haplotypes output file type",
            "doc": "Phased haplotypes output file type."
          },
          {
            "sbg:category": "Output Files",
            "id": "bingraph",
            "type": "string?",
            "inputBinding": {
              "prefix": "--bingraph",
              "shellQuote": false,
              "position": 20
            },
            "label": "File name for phased haplotypes in BIN format",
            "doc": "Phased haplotypes in BIN format. Useful to sample multiple likely haplotype configurations per sample."
          },
          {
            "sbg:category": "Output Files",
            "id": "log",
            "type": "string?",
            "inputBinding": {
              "prefix": "--log",
              "shellQuote": false,
              "position": 21
            },
            "label": "Log file name",
            "doc": "Log file name."
          }
        ],
        "outputs": [
          {
            "id": "out_phased_haplotypes",
            "doc": "Phased haplotypes.",
            "label": "Phased haplotypes",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    return [\"*_X.phased*\", \"*_[1-9].phased*\", \"*_1[0-9].phased*\", \"*_2[0-2].phased*\"]\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_genotypes))"
            },
            "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ"
          },
          {
            "id": "out_bin_phased_haplotypes",
            "doc": "Optional BIN phased haplotypes file.",
            "label": "Optional BIN phased haplotypes file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    if (inputs.bingraph)\n    {\n        return inputs.bingraph\n    }\n    else\n    {\n        return ''\n    }\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_genotypes))"
            },
            "sbg:fileTypes": "BIN"
          },
          {
            "id": "out_log",
            "doc": "Optional output log file.",
            "label": "Optional output log file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    if (inputs.log)\n    {\n        return inputs.log\n    }\n    else\n    {\n        return ''\n    }\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_genotypes))"
            },
            "sbg:fileTypes": "LOG, TXT"
          },
          {
            "id": "out_phased_haplotypes_par",
            "type": "File[]?",
            "outputBinding": {
              "glob": "*_par*",
              "outputEval": "$(inheritMetadata(self, inputs.in_genotypes))"
            }
          }
        ],
        "doc": "**SHAPEIT 4** is a phasing tool for sequencing and SNP array data  [1,2].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n_**Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.**_\n\n### Common Use Cases\n\n**SHAPEIT 4** can be used to estimate haplotypes starting from SNP array or sequencing data (**Genotypes to phase**). The tool phases a specified genomic region (**Target region**) and requires a genomic map (**Genomic map**). Genomic maps for humans are [available at the tool GitHub repository](https://github.com/odelaneau/shapeit4/tree/master/maps). \n\n### Changes Introduced by Seven Bridges\n\n* Parameter `--help` was omitted from the wrapper.\n* Parameter `--output` is no longer required, but the output file type for the phased haplotypes must be provided instead (**Phased haplotypes output file type**).\n\n### Common Issues and Important Notes\n\n* Inputs **Genotypes to phase**, **Genetic map**, **Target region** and **Phased haplotypes output file type**  are required.\n* **Genotypes to phase** input should be in BCF or VCF format, with the accompanying index (CSI or TBI). The authors recommend indexing this file with **BCFtools Index**.\n* If used, the **Reference panel** or **Scaffold of haplotypes** inputs should be in BCF/VCF format with the accompanying index (CSI or TBI).\n* Genetic maps for humans are available [at the tool GitHub repository](https://github.com/odelaneau/shapeit4/tree/master/maps). Please untar the archives before use and provide the chromosome file matching the other inputs.\n\n### Performance Benchmarking\n\nChromosome 20 of the 1000 genomes dataset (GRCh38 coordinates, 1647102 variants, 3202 individuals) was used for testing (chr20:10000000-30000000 and the entire chromosome) with the default tool parameters, except for `--sequencing` and `--thread`.\n\n| Experiment type  | Duration | Cost | Instance (AWS on-demand)|\n|-----------------------|-----------------|------------|-----------------|-------------|--------------|------------------|-------------|---------------|\n| 20 Mb - 15 cores  | 272 min | $3.08 + $0.06 |  c5.4xlarge 100 GB EBS |\n| 20 Mb - 30 cores  | 181 min | $4.61 + $0.08 |  c5.9xlarge 200 GB EBS |\n| chr20 - 15 cores  | 838 min | $10.73 + $0.39 |  m5.4xlarge 200 GB EBS |\n| chr20 - 30 cores  | 485 min | $12.37 + $0.22 |  c5.9xlarge 200 GB EBS |\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.* \n\n### References\n\n[1] [SHAPEIT 4 publication](https://www.nature.com/articles/s41467-019-13225-y)\n\n[2] [SHAPEIT 4 documentation](https://odelaneau.github.io/shapeit4/)",
        "label": "SHAPEIT 4",
        "arguments": [
          {
            "prefix": "--output",
            "shellQuote": false,
            "position": 19,
            "valueFrom": "${\n    if (inputs.prefix)\n    {\n        return inputs.prefix.concat('.phased.', inputs.suffix)\n    }\n    else\n    {\n        var pref = [].concat(inputs.in_genotypes)[0].nameroot.split('.vcf')[0].split('.bcf')[0] + \"_\" + inputs.in_map.metadata.region\n        if (inputs.in_map.path.includes('par1')) {\n                pref = pref + \"_par1\"\n            } else if  (inputs.in_map.path.includes('par2')) {\n                pref = pref + \"_par2\"\n            }\n        return pref.concat('.phased.', inputs.suffix)\n    }\n}"
          },
          {
            "prefix": "--region",
            "shellQuote": false,
            "position": 7,
            "valueFrom": "${\n    return inputs.in_map.metadata.region\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_per_job)\n    {\n        return inputs.mem_per_job\n    }\n    else\n    {\n        return 1000\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu_per_job)\n    {\n        return inputs.cpu_per_job\n    }\n    else if (inputs.thread)\n    {\n        return inputs.thread\n    }\n    else\n    {\n        return 1\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/jrandjelovic/shapeit4-4-2-1:0"
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "stdout": "stdout.log",
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982056,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:toolkit": "SHAPEIT",
        "sbg:toolkitVersion": "4.2.1 (5d26190)",
        "sbg:license": "MIT",
        "sbg:toolAuthor": "Olivier Delaneau",
        "sbg:links": [
          {
            "id": "https://odelaneau.github.io/shapeit4/",
            "label": "Homepage"
          },
          {
            "id": "https://github.com/odelaneau/shapeit4",
            "label": "Source Code"
          },
          {
            "id": "https://github.com/odelaneau/shapeit4/releases/tag/v4.2.1",
            "label": "Downloads"
          },
          {
            "id": "https://www.nature.com/articles/s41467-019-13225-y",
            "label": "Publication"
          },
          {
            "id": "https://odelaneau.github.io/shapeit4/#documentation",
            "label": "Documentation"
          }
        ],
        "sbg:categories": [
          "CWL1.1",
          "Phasing",
          "VCF Processing",
          "Utilities"
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/shapeit-4-4-2-1-cwl1-1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982056,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982056,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a7cea539f8f317662cd2be50c93679fb76f4b8a9fd44e9bf6c671fd52f5d6f35a"
      },
      "label": "SHAPEIT 4",
      "scatter": [
        "in_map"
      ],
      "scatterMethod": "dotproduct",
      "sbg:x": 1407.3214111328125,
      "sbg:y": -382.7875061035156
    },
    {
      "id": "bcftools_concat_1_10_1",
      "in": [
        {
          "id": "in_variants",
          "linkMerge": "merge_flattened",
          "source": [
            "shapeit_4_4_2_1_cwl1_1/out_phased_haplotypes_par"
          ],
          "pickValue": "all_non_null"
        },
        {
          "id": "output_name",
          "source": "output_name"
        },
        {
          "id": "output_type",
          "default": "CompressedVCF"
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
        "id": "markoz/hgi/bcftools-concat-1-10-1/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:category": "Input",
            "id": "in_variants",
            "type": "File[]",
            "label": "Input variants file",
            "doc": "Input files which will be concatenated.",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
            "secondaryFiles": [
              "${ \n    if(self.basename.split('.').pop() == 'gz'){\n    if(self.nameroot.split('.').pop() == 'bcf'){\n        return self.nameroot + \".gz.csi\"}\n    else{\n        return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n    if(self.basename.split('.').pop() == 'bcf'){\n        return self.basename + \".csi\"\n    }\n    else{\n    return self.basename + \".tbi\"}\n}\n\n}"
            ]
          },
          {
            "sbg:category": "File format options",
            "id": "output_name",
            "type": "string?",
            "label": "Output file name",
            "doc": "Name of the output file."
          },
          {
            "sbg:category": "File format options",
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
              "position": 5,
              "valueFrom": "${\n    if (self == 0) {\n        self = null;\n        inputs.output_type = null\n    };\n\n\n    if (inputs.output_type === 'CompressedBCF') return 'b'\n    if (inputs.output_type === 'UncompressedBCF') return 'u'\n    if (inputs.output_type === 'CompressedVCF') return 'z'\n    if (inputs.output_type === 'UncompressedVCF') return 'v'\n}"
            },
            "label": "Output type",
            "doc": "Output types: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v].",
            "default": 0
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "0",
            "id": "threads",
            "type": "int?",
            "inputBinding": {
              "prefix": "--threads",
              "shellQuote": false,
              "position": 41
            },
            "label": "Threads",
            "doc": "Number of output compression threads to use in addition to main thread. Only used when output type is CompressedBCF or CompressedVCF."
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
            "sbg:toolDefaultValue": "1000",
            "sbg:category": "Execution",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job",
            "doc": "Memory per job in MB. Appropriate instance will be chosen based on this parameter."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-a",
            "id": "allow_overlaps",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--allow-overlaps",
              "shellQuote": false,
              "position": 6
            },
            "label": "Allow overlaps",
            "doc": "First coordinate of the next file can precede last record of the current file."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-c",
            "id": "compact_ps",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--compact-PS",
              "shellQuote": false,
              "position": 7
            },
            "label": "Compact phase set",
            "doc": "Do not output PS tag at each site, only at the start of a new phase set block."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-d",
            "id": "remove_duplicates",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "snps",
                  "indels",
                  "both",
                  "all",
                  "none"
                ],
                "name": "remove_duplicates"
              }
            ],
            "inputBinding": {
              "prefix": "--rm-dups",
              "shellQuote": false,
              "position": 8,
              "valueFrom": "${\n    if (self == 0) {\n        self = null;\n        inputs.remove_duplicates = null\n    };\n\n\n    if (inputs.remove_duplicates === 'snps') return 'snps'\n    if (inputs.remove_duplicates === 'indels') return 'indels'\n    if (inputs.remove_duplicates === 'both') return 'both'\n    if (inputs.remove_duplicates === 'all') return 'all'\n    if (inputs.remove_duplicates === 'none') return 'none'\n}"
            },
            "label": "Remove duplicates",
            "doc": "Output duplicate records present in multiple files only once: <snps|indels|both|all|none>. Requires -a, --allow-overlaps.",
            "default": 0
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-D",
            "id": "remove_all_duplicates",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--remove-duplicates",
              "shellQuote": false,
              "position": 9
            },
            "label": "Remove all duplicates",
            "doc": "Alias for -d none."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-l",
            "id": "ligate",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--ligate",
              "shellQuote": false,
              "position": 10
            },
            "label": "Ligate phased VCF",
            "doc": "Ligate phased VCFs by matching phase at overlapping haplotypes."
          },
          {
            "sbg:category": "General Options",
            "id": "no_version",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--no-version",
              "shellQuote": false,
              "position": 13
            },
            "label": "No version",
            "doc": "Do not append version and command line to the header."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-n",
            "id": "naive_concat",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--naive",
              "shellQuote": false,
              "position": 12
            },
            "label": "Naive concat",
            "doc": "Concatenate files without recompression (dangerous, use with caution)."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-q",
            "id": "min_pq",
            "type": "float?",
            "inputBinding": {
              "prefix": "--min-PQ",
              "shellQuote": false,
              "position": 14
            },
            "label": "Min phase quality",
            "doc": "Break phase set if phasing quality is lower than."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-r",
            "id": "regions",
            "type": "string[]?",
            "inputBinding": {
              "prefix": "--regions",
              "shellQuote": false,
              "position": 15
            },
            "label": "Regions",
            "doc": "Restrict to comma-separated list of regions (e.g. chr|chr:pos|chr:from-to|chr:from-[,\u2026])."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-R",
            "id": "regions_file",
            "type": "File?",
            "inputBinding": {
              "prefix": "--regions-file",
              "shellQuote": false,
              "position": 16
            },
            "label": "Regions file",
            "doc": "Restrict to regions listed in a file.",
            "sbg:fileTypes": "BED, TXT"
          }
        ],
        "outputs": [
          {
            "id": "out_variants",
            "doc": "Concatenated VCF from multiple inputs.",
            "label": "Concatenated VCF file",
            "type": "File?",
            "outputBinding": {
              "glob": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    var array_length = files_array.length;\n    if (array_length != 1){\n        if (fext == '.gz') {\n            if (froot.split('.').pop() == 'vcf'){\n                froot = froot.split('.vcf')[0]}\n            else if (froot.split('.').pop() == 'bcf'){\n                froot = froot.split('.bcf')[0]}\n            \n        }\n    \n        if(in_file.metadata.sample_id){\n            var froot = in_file.metadata.sample_id;}\n    \n        if (inputs.output_name) {\n            var out = inputs.output_name\n            if (inputs.output_type == 'UncompressedVCF') {\n                out += \".concatenated.vcf\";} \n            else if (inputs.output_type == 'CompressedVCF') {\n                out += \".concatenated.vcf.gz\";} \n            else if (inputs.output_type == 'UncompressedBCF') {\n                out += \".concatenated.bcf\";} \n            else if (inputs.output_type == 'CompressedBCF') {\n                out += \".concatenated.bcf.gz\";} \n            else {\n                out += \".concatenated.vcf\";}\n        }\n        \n        else if (inputs.output_type == 'UncompressedVCF') {\n            var out = froot + '.concatenated' + '.vcf';} \n        else if (inputs.output_type == 'CompressedVCF') {\n            var out = froot + '.concatenated' + '.vcf.gz';} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            var out = froot + '.concatenated' + '.bcf';} \n        else if (inputs.output_type == 'CompressedBCF') {\n            var out = froot + '.concatenated' + '.bcf.gz';} \n        else var out = froot + '.concatenated.vcf';\n    \n        return out;}\n    \n    else{return fname;}\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
            },
            "secondaryFiles": [
              ".tbi"
            ],
            "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ"
          }
        ],
        "doc": "**BCFtools Concat**: Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. \n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n\n### Common Use Cases\n\nCan be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an INDEL VCF into one. The input files must be sorted by chr and position. The files must be given in the correct order to produce sorted VCF on the output unless the **Allow overlaps** (`--allow-overlaps`) option is specified. \n\nWith the **Naive option** (`--naive`) set to True, the files are concatenated without being recompressed, which is very fast but dangerous if the BCF headers differ. \n```\n$bcftools concat --allow-overlaps --naive vcf_file1.vcf.gz vcf_file2.vcf.gz\n```\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands afterwards if a VCF file is provided on input. If VCF.GZ is given on input only indexing will be done. Output type can still be chosen with the `output type` command.\n \n\n### Common Issues and Important Notes\n\n * All VCF files must have same sample names, otherwise tool will fail.\n * For option `--allow-overlaps` files must be compressed. \n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using two inputs sized 12.4 MB and 56 KB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
        "label": "Bcftools Concat",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    if (inputs.in_variants.length == 1) {\n        return \"python copy_files.py && echo pass\";\n    } else {\n        return \"python copy_files.py && bcftools concat \";\n    }\n\n}"
          },
          {
            "prefix": "--output",
            "shellQuote": false,
            "position": 6,
            "valueFrom": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.bcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id;}\n\n    if (inputs.output_name) {\n        var out = inputs.output_name;\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".concatenated.vcf\";} \n        else if (inputs.output_type == 'CompressedVCF') {\n            out += \".concatenated.vcf.gz\";} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".concatenated.bcf\";} \n        else if (inputs.output_type == 'CompressedBCF') {\n            out += \".concatenated.bcf.gz\";} \n        else {\n            out += \".concatenated.vcf\";}\n    }\n    \n    else if (inputs.output_type == 'UncompressedVCF') {\n        var out = froot + '.concatenated' + '.vcf';} \n    else if (inputs.output_type == 'CompressedVCF') {\n        var out = froot + '.concatenated' + '.vcf.gz';} \n    else if (inputs.output_type == 'UncompressedBCF') {\n        var out = froot + '.concatenated' + '.bcf';} \n    else if (inputs.output_type == 'CompressedBCF') {\n        var out = froot + '.concatenated' + '.bcf.gz';} \n    else var out = froot + '.concatenated.vcf';\n\n    return out;\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 100,
            "valueFrom": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    if(fext == '.vcf'){\n        return \"./input_files/*.vcf.gz\"\n    }\n    if(fext == '.bcf'){\n        return \"./input_files/*.bcf\"\n    }\n    else if(fext == '.gz'){\n        if (froot.split('.').pop() == 'bcf'){\n        return \"./input_files/*.bcf.gz\"}\n        else if (froot.split('.').pop() == 'vcf'){\n        return \"./input_files/*.vcf.gz\"}\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.memory) {\n\n        return inputs.memory\n\n    } else {\n        return 1000\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu) {\n        return inputs.cpu\n    } else {\n        return 1\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "21caaa02f72e",
            "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools-1.10.1:1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "file_paths.txt",
                "entry": "${  var i;\n    var text = ''\n    var file_array = [].concat(inputs.in_variants);\n    for (i = 0; i < file_array.length; i++) {\n       text += file_array[i].path + '\\n';}\n    \n    return text;\n}",
                "writable": false
              },
              {
                "entryname": "copy_files.py",
                "entry": "import os\n\nos.system('mkdir input_files')\n\nwith open('./file_paths.txt') as file:\n    for counter, line in enumerate(file):\n        mv_line = line.strip('\\n')\n        basename = line.split('/')[-1].strip('\\n')\nfile.close()\n\n\nif (counter > 0):\n    with open('./file_paths.txt') as txt:\n        for counter, line in enumerate(txt):\n            mv_line = line.strip('\\n')\n            basename = line.split('/')[-1].strip('\\n')\n            if basename.endswith('.vcf.gz'):\n                nameroot = basename.split('.vcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf.gz'\n            elif basename.endswith('.vcf'):\n                nameroot = basename.split('.vcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf'\n            elif basename.endswith('.bcf.gz'):\n                nameroot = basename.split('.bcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf.gz'\n            elif basename.endswith('.bcf'):\n                nameroot = basename.split('.bcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf'\n\n            cmd = \"cp \" + mv_line + \" ./input_files/\" + new_name\n            status = os.system(cmd)\n\n    bash_cmd = \"bash /opt/gzip_vcf.sh\"\n    os.system(bash_cmd)\n\nelse:\n    cmd = \"cp \" + mv_line + \" ./\" + basename\n    print(cmd)\n    os.system(cmd)\n            ",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};",
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:toolkitVersion": "1.10.1",
        "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
        "sbg:license": "MIT License",
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
        "sbg:cmdPreview": "bash gzip_vcf.sh && bcftools concat --output input_file-1.concated.vcf  input_file-1.vcf.gz input_file-2.vcf.gz",
        "sbg:toolkit": "bcftools",
        "sbg:image_url": null,
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982057,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:appVersion": [
          "v1.0"
        ],
        "sbg:id": "markoz/hgi/bcftools-concat-1-10-1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982057,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982057,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "af5855b119c26e49f6aeca0a140859deda105ad0adc92d1ee3fd0eb2d0b2a0731"
      },
      "label": "Bcftools Concat",
      "sbg:x": 1619.142578125,
      "sbg:y": -514.5716552734375
    },
    {
      "id": "sbg_reorder",
      "in": [
        {
          "id": "input",
          "linkMerge": "merge_flattened",
          "source": [
            "bcftools_concat_1_10_1/out_variants",
            "shapeit_4_4_2_1_cwl1_1/out_phased_haplotypes"
          ],
          "pickValue": "all_non_null"
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
        "id": "markoz/hgi/sbg-reorder/0",
        "baseCommand": [
          "echo"
        ],
        "inputs": [
          {
            "id": "input",
            "type": "File[]?",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "output",
            "type": "File[]?",
            "outputBinding": {
              "outputEval": "${\nfunction compare( a, b ) {\n  if ( a.nameroot < b.metadata.nameroot ){\n    return -1;\n  }\n  if ( a.metadata.nameroot > b.metadata.nameroot ){\n    return 1;\n  }\n  return 0;\n}\n\n    return inputs.input.sort(compare)\n}"
            }
          }
        ],
        "label": "reorder",
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entry": "$(inputs.input)",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982058,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/sbg-reorder/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982058,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982058,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a7fc0d07691e31df34cd2aaf1eb1df30ab80aa61211e83af5e4199098fbe01c49"
      },
      "label": "reorder",
      "sbg:x": 1761.57763671875,
      "sbg:y": -393.945556640625
    },
    {
      "id": "minimac4_1_0_2_cwl1_2",
      "in": [
        {
          "id": "in_refhaps",
          "source": "in_refhaps"
        },
        {
          "id": "in_haps",
          "source": "sbg_reorder/output"
        },
        {
          "id": "cpu_per_job",
          "default": 2
        },
        {
          "id": "mem_per_job",
          "default": 4096
        }
      ],
      "out": [
        {
          "id": "out_log"
        },
        {
          "id": "out_dose"
        },
        {
          "id": "out_info"
        },
        {
          "id": "out_empirical_dose"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/minimac4-1-0-2-cwl1-2/0",
        "baseCommand": [
          "/opt/minimac4-1.0.2/bin/minimac4",
          "--noPhoneHome",
          "--ignoreDuplicates"
        ],
        "inputs": [
          {
            "sbg:category": "Reference Haplotypes",
            "id": "in_refhaps",
            "type": "File",
            "inputBinding": {
              "prefix": "--refHaps",
              "shellQuote": false,
              "position": 1
            },
            "label": "Reference haplotypes",
            "doc": "M3VCF file containing haplotype data for reference panel.",
            "sbg:fileTypes": "M3VCF, M3VCF.GZ"
          },
          {
            "sbg:category": "Reference Haplotypes",
            "sbg:toolDefaultValue": "False",
            "id": "pass_only",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--passOnly",
              "shellQuote": false,
              "position": 2
            },
            "label": "Import only PASS variants",
            "doc": "This option only imports variants with FILTER = PASS."
          },
          {
            "sbg:category": "Reference Haplotypes",
            "sbg:toolDefaultValue": "False",
            "id": "rsid",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--rsid",
              "shellQuote": false,
              "position": 3
            },
            "label": "Import only RS IDs",
            "doc": "This option only imports RS ID of variants from ID column (if available)."
          },
          {
            "sbg:category": "GWAS Haplotypes",
            "id": "in_haps",
            "type": "File",
            "inputBinding": {
              "prefix": "--haps",
              "shellQuote": false,
              "position": 4
            },
            "label": "Pre-phased target genotype data",
            "doc": "File containing haplotype data for target (GWAS) samples in VCF format.",
            "sbg:fileTypes": "VCF, VCF.GZ",
            "secondaryFiles": [
              {
                "pattern": ".tbi",
                "required": false
              }
            ]
          },
          {
            "sbg:category": "Output Parameters",
            "id": "prefix",
            "type": "string?",
            "label": "Output file names prefix",
            "doc": "Prefix for all output files generated."
          },
          {
            "sbg:category": "Output Parameters",
            "sbg:toolDefaultValue": "False",
            "id": "nobgzip",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--nobgzip",
              "shellQuote": false,
              "position": 6
            },
            "label": "Do not GZIP outputs",
            "doc": "If ON, output files will NOT be gzipped."
          },
          {
            "sbg:category": "Output Parameters",
            "sbg:toolDefaultValue": "GT,DS",
            "id": "format",
            "type": [
              "null",
              {
                "type": "array",
                "items": {
                  "type": "enum",
                  "name": "format",
                  "symbols": [
                    "GT",
                    "DS",
                    "HDS",
                    "GP",
                    "SD"
                  ]
                }
              }
            ],
            "inputBinding": {
              "prefix": "--format",
              "shellQuote": false,
              "position": 7
            },
            "label": "FORMAT fields",
            "doc": "Specifies which fields to output for the FORMAT field in output VCF file. Default: GT,DS.\nGT - Estimated most likely genotype.\nDS - Estimated alternate allele dosage [P(0/1)+2*P(1/1)].\nHDS - Estimated phased haploid alternate allele dosage.\nGP - Estimated Posterior Genotype Probabilities P(0/0), P(0/1) and P(1/1).\nSD - Estimated Variance of Posterior Genotype Probabilities."
          },
          {
            "sbg:category": "Output Parameters",
            "sbg:toolDefaultValue": "False",
            "id": "all_typed_sites",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--allTypedSites",
              "shellQuote": false,
              "position": 8
            },
            "label": "Include all genotyped sites",
            "doc": "If ON, Minimac4 will also include variants that were genotyped but NOT in the reference panel in the output files (and imputes any missing data in such variants to the major allele frequency)."
          },
          {
            "sbg:category": "Subset Parameters",
            "id": "chr",
            "type": "string?",
            "inputBinding": {
              "prefix": "--chr",
              "shellQuote": false,
              "position": 9
            },
            "label": "Chromosome for imputation",
            "doc": "Chromosome number for which to carry out imputation."
          },
          {
            "sbg:category": "Subset Parameters",
            "id": "start",
            "type": "int?",
            "inputBinding": {
              "prefix": "--start",
              "shellQuote": false,
              "position": 11
            },
            "label": "Start position for imputation",
            "doc": "Start position for imputation by chunking."
          },
          {
            "sbg:category": "Subset Parameters",
            "id": "end",
            "type": "int?",
            "inputBinding": {
              "prefix": "--end",
              "shellQuote": false,
              "position": 12
            },
            "label": "End position for imputation",
            "doc": "End position for imputation by chunking."
          },
          {
            "sbg:category": "Subset Parameters",
            "sbg:toolDefaultValue": "500000",
            "id": "window",
            "type": "int?",
            "inputBinding": {
              "prefix": "--window",
              "shellQuote": false,
              "position": 13
            },
            "label": "Buffer region size",
            "doc": "Length of buffer region on either side of --start and --end."
          },
          {
            "sbg:category": "Other Parameters",
            "sbg:toolDefaultValue": "False",
            "id": "log",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--log",
              "shellQuote": false,
              "position": 14
            },
            "label": "Write log file",
            "doc": "If ON, log will be written to $prefix.logfile."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "5",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "CPUs per job",
            "doc": "CPUs per job."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "1000",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job [MB]",
            "doc": "Memory per job [MB]."
          },
          {
            "sbg:category": "Other Parameters",
            "sbg:toolDefaultValue": "5",
            "id": "cpus",
            "type": "int?",
            "inputBinding": {
              "prefix": "--cpus",
              "shellQuote": false,
              "position": 15
            },
            "label": "Number of CPUs for parallel computing",
            "doc": "Number of CPUs for parallel computing."
          },
          {
            "sbg:toolDefaultValue": "False",
            "sbg:category": "Output Parameters",
            "id": "meta",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--meta",
              "shellQuote": false,
              "position": 16
            },
            "label": "Create a file for meta-imputation",
            "doc": "OFF by default. If ON, Minimac4 will generate a separate file that can be used by MetaMinimac2 for meta-imputation."
          },
          {
            "sbg:toolDefaultValue": "200",
            "sbg:category": "Output Parameters",
            "id": "vcfbuffer",
            "type": "int?",
            "inputBinding": {
              "prefix": "--vcfBuffer",
              "shellQuote": false,
              "position": 17
            },
            "label": "VCF buffer",
            "doc": "This option defines the maximum number of samples in the target genotype data to be imputed at a time. By default, it is set as 200, or the total number of samples, whichever is smaller. Note that the larger the value is, the more memory Minimac4 will consume."
          },
          {
            "sbg:category": "Output Parameters",
            "sbg:toolDefaultValue": "False",
            "id": "mem_usage",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--memUsage",
              "shellQuote": false,
              "position": 18
            },
            "label": "Estimate memory usage only",
            "doc": "OFF by default. If ON, Minimac4 will not perform imputation. Instead, it will estimate memory that imputation would consume based on a single chunk."
          },
          {
            "sbg:category": "Chunking Parameters",
            "sbg:toolDefaultValue": "20",
            "id": "chunk_length_mb",
            "type": "float?",
            "inputBinding": {
              "prefix": "--ChunkLengthMb",
              "shellQuote": false,
              "position": 19
            },
            "label": "Chunk length [Mb]",
            "doc": "This option defines the average length of chunks in units of million base pairs (Mb). The input value should be within (0.001, 300]. The default setting is 20."
          },
          {
            "sbg:category": "Chunking Parameters",
            "sbg:toolDefaultValue": "3",
            "id": "chunk_overlap_mb",
            "type": "float?",
            "inputBinding": {
              "prefix": "--ChunkOverlapMb",
              "shellQuote": false,
              "position": 20
            },
            "label": "Chunk overlap [Mb]",
            "doc": "This option defines the length of overlap between chunks in units of Mb, 3Mb by default. The valid input value should be within (0.001, 300].  The overlap length should be at most 1/3 of the chunk length, if larger, Minimac4 will automatically reduce it to 1/3 of the chunk length."
          },
          {
            "sbg:category": "Approximation Parameters",
            "sbg:toolDefaultValue": "False",
            "id": "minimac3",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--minimac3",
              "shellQuote": false,
              "position": 21
            },
            "label": "Use Minimac3 algorithm for imputation",
            "doc": "OFF by default. If ON, Minimac3 algorithm will be used for imputation."
          },
          {
            "sbg:category": "Approximation Parameters",
            "sbg:toolDefaultValue": "0.01",
            "id": "prob_threshold",
            "type": "float?",
            "inputBinding": {
              "prefix": "--probThreshold",
              "shellQuote": false,
              "position": 22
            },
            "label": "Approximation levels - probThreshold",
            "doc": "Approximation levels - probThreshold."
          },
          {
            "sbg:category": "Approximation Parameters",
            "sbg:toolDefaultValue": "0.01",
            "id": "diff_threshold",
            "type": "float?",
            "inputBinding": {
              "prefix": "--diffThreshold",
              "shellQuote": false,
              "position": 23
            },
            "label": "Approximation levels - diffThreshold",
            "doc": "Approximation levels - diffThreshold."
          },
          {
            "sbg:category": "Approximation Parameters",
            "sbg:toolDefaultValue": "0.01",
            "id": "top_threshold",
            "type": "float?",
            "inputBinding": {
              "prefix": "--topThreshold",
              "shellQuote": false,
              "position": 24
            },
            "label": "Approximation levels - topThreshold",
            "doc": "Approximation levels - topThreshold."
          },
          {
            "sbg:category": "Platform Options",
            "sbg:toolDefaultValue": "False",
            "id": "debug_log",
            "type": "boolean?",
            "label": "Create a log file for debugging",
            "doc": "Create a log file for debugging."
          }
        ],
        "outputs": [
          {
            "id": "out_log",
            "doc": "Optional Minimac4 log file.",
            "label": "Optional Minimac4 log file",
            "type": "File?",
            "outputBinding": {
              "glob": "*.logfile",
              "outputEval": "$(inheritMetadata(self, inputs.in_haps))"
            },
            "sbg:fileTypes": "LOGFILE"
          },
          {
            "id": "out_dose",
            "doc": "Minimap4 imputation results.",
            "label": "Minimap4 imputation results",
            "type": "File?",
            "outputBinding": {
              "glob": "*.dose.vc*",
              "outputEval": "$(inheritMetadata(self, inputs.in_haps))"
            },
            "sbg:fileTypes": "VCF.GZ, VCF"
          },
          {
            "id": "out_info",
            "doc": "Minimac4 info file.",
            "label": "Minimac4 info file",
            "type": "File?",
            "outputBinding": {
              "glob": "*.info",
              "outputEval": "$(inheritMetadata(self, inputs.in_haps))"
            },
            "sbg:fileTypes": "INFO"
          },
          {
            "id": "out_empirical_dose",
            "doc": "Optional empirical dosage file.",
            "label": "Optional empirical dosage file",
            "type": "File?",
            "outputBinding": {
              "glob": "*.empiricalDose.vc*",
              "outputEval": "$(inheritMetadata(self, inputs.in_haps))"
            },
            "sbg:fileTypes": "VCF, VCF.GZ"
          }
        ],
        "doc": "**Minimac4** is a genetic imputation algorithm  [1,2].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n_**Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.**_\n\n### Common Use Cases\n\n**Minimac4** can be used to impute genotypes in a genomic region starting from a reference panel in M3VCF format (**Reference haplotypes**) and pre-phased target GWAS haplotypes (**Pre-phased target genotype data**). By default, **Minimac4** automates the chunking process.\n\n### Changes Introduced by Seven Bridges\n\n* Parameter `--processReference` was omitted from the wrapper as the option is currently deactivated in the tool.\n* Parameter `--help` was omitted from the wrapper.\n* Parameter `--noPhoneHome` is hardcoded in the wrapper, consequently `--phoneHomeThinning` parameter was omitted from the wrapper.\n* Parameters `--meta`, `--vcfBuffer`, `--memUsage`, `--ChunkLengthMb`, `--ChunkOverlapMb`, `--minimac3`, `--probThreshold`, `--diffThreshold` and `--topThreshold` were absent from the tool usage print out, but were added to the wrapper based on the tool documentation page.\n* If not provided by the user, the `--prefix` parameter was set to default to the root of the **Pre-phased target genotype data** input file, in order to facilitate batch executions.\n* If the tool is run in **Estimate memory usage** mode, the log file will be automatically created (`--log`) to store the results of the task.\n* **Create a log file for debugging** input was added to help with debugging failed tasks. If set, a debug.log file will be created and kept among task execution logs.\n\n### Common Issues and Important Notes\n\n* Inputs **Reference haplotypes** and **Pre-phased target genotype data** are required.\n* **Reference haplotypes** file must be in M3VCF format. Reference panels in VCF format should be converted to M3VCF (along with parameter estimation) with Minimac3 before use.\n* If **Chromosome for imputation** input is used **Start position for imputation** and **End position for imputation** must be specified as well.\n* When run in **Estimate memory usage** mode, the tool will create dummy outputs containing only headers. These files can be ignored - the memory estimation results are written in the log file.\n\n\n### Performance Benchmarking\n\n| Experiment type  | Duration | Cost | Instance (AWS on-demand)|\n|---------------------------|------------------------|-----------------------|--------------------------------|\n| 1000g chr20 | 102 mins | $0.67 + $0.23 | c4.2xlarge - 1024 GB EBS | \n| 1000g chr20 - 15 cpus | 48 mins | $0.54 + $0.11 | c5.4xlarge - 1024 GB EBS | \n| 1000g chr1 - 30 cpus | 136 mins | $3.47 + $0.31 | c5.9xlarge - 1024 GB EBS | \n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.* \n\n### References\n\n[1] [Minimac2 publication](https://pubmed.ncbi.nlm.nih.gov/25338720/)\n\n[2] [Minimac4 documentation](https://genome.sph.umich.edu/wiki/Minimac4)",
        "label": "Minimac4",
        "arguments": [
          {
            "prefix": "--prefix",
            "shellQuote": false,
            "position": 5,
            "valueFrom": "${\n    if (inputs.prefix)\n    {\n        return inputs.prefix\n    }\n    else\n    {\n        var test = [].concat(inputs.in_haps)[0].nameroot\n        return test.concat('.minimac4')\n    }\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 25,
            "valueFrom": "${\n    if ((inputs.mem_usage) && !(inputs.log))\n    {\n        return ' --log '\n    }\n    else\n    {\n        return ''\n    }\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 100,
            "valueFrom": "${\n    if (inputs.debug_log)\n    {\n        return ' > debug.log 2>&1'\n    }\n    else\n    {\n        return ''\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_per_job)\n    {\n        return inputs.mem_per_job\n    }\n    else\n    {\n        return 1000\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu_per_job)\n    {\n        return inputs.cpu_per_job\n    }\n    else if (inputs.cpus)\n    {\n        return inputs.cpus\n    }\n    else\n    {\n        return 5\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/jrandjelovic/minimac4-1-0-2:1"
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982059,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:image_url": null,
        "sbg:toolkit": "Minimac4",
        "sbg:toolkitVersion": "1.0.2",
        "sbg:toolAuthor": "The Center for Statistical Genetics at the University of Michigan School of Public Health",
        "sbg:license": "GPL-3.0",
        "sbg:categories": [
          "CWL1.2",
          "Imputation"
        ],
        "sbg:links": [
          {
            "id": "https://genome.sph.umich.edu/wiki/Minimac4",
            "label": "Homepage"
          },
          {
            "id": "https://github.com/statgen/Minimac4",
            "label": "Source Code"
          },
          {
            "id": "https://github.com/statgen/Minimac4/releases/tag/v1.0.2",
            "label": "Download"
          },
          {
            "id": "https://pubmed.ncbi.nlm.nih.gov/25338720/",
            "label": "Publication"
          },
          {
            "id": "https://genome.sph.umich.edu/wiki/Minimac4_Documentation",
            "label": "Documentation"
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/minimac4-1-0-2-cwl1-2/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982059,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982059,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a5f26b5302ca017f6e0d864d41c703638743eb7521c8244f76913dd391506b5b6"
      },
      "label": "Minimac4",
      "scatter": [
        "in_refhaps",
        "in_haps"
      ],
      "scatterMethod": "dotproduct",
      "sbg:x": 1928.357666015625,
      "sbg:y": -535.582763671875
    },
    {
      "id": "bcftools_concat_1_10_2",
      "in": [
        {
          "id": "in_variants",
          "source": [
            "minimac4_1_0_2_cwl1_2/out_dose"
          ]
        },
        {
          "id": "output_name",
          "default": "imgge.covid.final",
          "source": "output_name_1"
        },
        {
          "id": "output_type",
          "default": "CompressedVCF"
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
        "id": "markoz/hgi/bcftools-concat-1-10-1/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:category": "Input",
            "id": "in_variants",
            "type": "File[]",
            "label": "Input variants file",
            "doc": "Input files which will be concatenated.",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
            "secondaryFiles": [
              "${ \n    if(self.basename.split('.').pop() == 'gz'){\n    if(self.nameroot.split('.').pop() == 'bcf'){\n        return self.nameroot + \".gz.csi\"}\n    else{\n        return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n    if(self.basename.split('.').pop() == 'bcf'){\n        return self.basename + \".csi\"\n    }\n    else{\n    return self.basename + \".tbi\"}\n}\n\n}"
            ]
          },
          {
            "sbg:category": "File format options",
            "id": "output_name",
            "type": "string?",
            "label": "Output file name",
            "doc": "Name of the output file."
          },
          {
            "sbg:category": "File format options",
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
              "position": 5,
              "valueFrom": "${\n    if (self == 0) {\n        self = null;\n        inputs.output_type = null\n    };\n\n\n    if (inputs.output_type === 'CompressedBCF') return 'b'\n    if (inputs.output_type === 'UncompressedBCF') return 'u'\n    if (inputs.output_type === 'CompressedVCF') return 'z'\n    if (inputs.output_type === 'UncompressedVCF') return 'v'\n}"
            },
            "label": "Output type",
            "doc": "Output types: b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v].",
            "default": 0
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "0",
            "id": "threads",
            "type": "int?",
            "inputBinding": {
              "prefix": "--threads",
              "shellQuote": false,
              "position": 41
            },
            "label": "Threads",
            "doc": "Number of output compression threads to use in addition to main thread. Only used when output type is CompressedBCF or CompressedVCF."
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
            "sbg:toolDefaultValue": "1000",
            "sbg:category": "Execution",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job",
            "doc": "Memory per job in MB. Appropriate instance will be chosen based on this parameter."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-a",
            "id": "allow_overlaps",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--allow-overlaps",
              "shellQuote": false,
              "position": 6
            },
            "label": "Allow overlaps",
            "doc": "First coordinate of the next file can precede last record of the current file."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-c",
            "id": "compact_ps",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--compact-PS",
              "shellQuote": false,
              "position": 7
            },
            "label": "Compact phase set",
            "doc": "Do not output PS tag at each site, only at the start of a new phase set block."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-d",
            "id": "remove_duplicates",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "snps",
                  "indels",
                  "both",
                  "all",
                  "none"
                ],
                "name": "remove_duplicates"
              }
            ],
            "inputBinding": {
              "prefix": "--rm-dups",
              "shellQuote": false,
              "position": 8,
              "valueFrom": "${\n    if (self == 0) {\n        self = null;\n        inputs.remove_duplicates = null\n    };\n\n\n    if (inputs.remove_duplicates === 'snps') return 'snps'\n    if (inputs.remove_duplicates === 'indels') return 'indels'\n    if (inputs.remove_duplicates === 'both') return 'both'\n    if (inputs.remove_duplicates === 'all') return 'all'\n    if (inputs.remove_duplicates === 'none') return 'none'\n}"
            },
            "label": "Remove duplicates",
            "doc": "Output duplicate records present in multiple files only once: <snps|indels|both|all|none>. Requires -a, --allow-overlaps.",
            "default": 0
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-D",
            "id": "remove_all_duplicates",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--remove-duplicates",
              "shellQuote": false,
              "position": 9
            },
            "label": "Remove all duplicates",
            "doc": "Alias for -d none."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-l",
            "id": "ligate",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--ligate",
              "shellQuote": false,
              "position": 10
            },
            "label": "Ligate phased VCF",
            "doc": "Ligate phased VCFs by matching phase at overlapping haplotypes."
          },
          {
            "sbg:category": "General Options",
            "id": "no_version",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--no-version",
              "shellQuote": false,
              "position": 13
            },
            "label": "No version",
            "doc": "Do not append version and command line to the header."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-n",
            "id": "naive_concat",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--naive",
              "shellQuote": false,
              "position": 12
            },
            "label": "Naive concat",
            "doc": "Concatenate files without recompression (dangerous, use with caution)."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-q",
            "id": "min_pq",
            "type": "float?",
            "inputBinding": {
              "prefix": "--min-PQ",
              "shellQuote": false,
              "position": 14
            },
            "label": "Min phase quality",
            "doc": "Break phase set if phasing quality is lower than."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-r",
            "id": "regions",
            "type": "string[]?",
            "inputBinding": {
              "prefix": "--regions",
              "shellQuote": false,
              "position": 15
            },
            "label": "Regions",
            "doc": "Restrict to comma-separated list of regions (e.g. chr|chr:pos|chr:from-to|chr:from-[,\u2026])."
          },
          {
            "sbg:category": "General Options",
            "sbg:altPrefix": "-R",
            "id": "regions_file",
            "type": "File?",
            "inputBinding": {
              "prefix": "--regions-file",
              "shellQuote": false,
              "position": 16
            },
            "label": "Regions file",
            "doc": "Restrict to regions listed in a file.",
            "sbg:fileTypes": "BED, TXT"
          }
        ],
        "outputs": [
          {
            "id": "out_variants",
            "doc": "Concatenated VCF from multiple inputs.",
            "label": "Concatenated VCF file",
            "type": "File?",
            "outputBinding": {
              "glob": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    var array_length = files_array.length;\n    if (array_length != 1){\n        if (fext == '.gz') {\n            if (froot.split('.').pop() == 'vcf'){\n                froot = froot.split('.vcf')[0]}\n            else if (froot.split('.').pop() == 'bcf'){\n                froot = froot.split('.bcf')[0]}\n            \n        }\n    \n        if(in_file.metadata.sample_id){\n            var froot = in_file.metadata.sample_id;}\n    \n        if (inputs.output_name) {\n            var out = inputs.output_name\n            if (inputs.output_type == 'UncompressedVCF') {\n                out += \".concatenated.vcf\";} \n            else if (inputs.output_type == 'CompressedVCF') {\n                out += \".concatenated.vcf.gz\";} \n            else if (inputs.output_type == 'UncompressedBCF') {\n                out += \".concatenated.bcf\";} \n            else if (inputs.output_type == 'CompressedBCF') {\n                out += \".concatenated.bcf.gz\";} \n            else {\n                out += \".concatenated.vcf\";}\n        }\n        \n        else if (inputs.output_type == 'UncompressedVCF') {\n            var out = froot + '.concatenated' + '.vcf';} \n        else if (inputs.output_type == 'CompressedVCF') {\n            var out = froot + '.concatenated' + '.vcf.gz';} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            var out = froot + '.concatenated' + '.bcf';} \n        else if (inputs.output_type == 'CompressedBCF') {\n            var out = froot + '.concatenated' + '.bcf.gz';} \n        else var out = froot + '.concatenated.vcf';\n    \n        return out;}\n    \n    else{return fname;}\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
            },
            "secondaryFiles": [
              ".tbi"
            ],
            "sbg:fileTypes": "VCF, BCF, VCF.GZ, BCF.GZ"
          }
        ],
        "doc": "**BCFtools Concat**: Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. \n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n\n### Common Use Cases\n\nCan be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an INDEL VCF into one. The input files must be sorted by chr and position. The files must be given in the correct order to produce sorted VCF on the output unless the **Allow overlaps** (`--allow-overlaps`) option is specified. \n\nWith the **Naive option** (`--naive`) set to True, the files are concatenated without being recompressed, which is very fast but dangerous if the BCF headers differ. \n```\n$bcftools concat --allow-overlaps --naive vcf_file1.vcf.gz vcf_file2.vcf.gz\n```\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands afterwards if a VCF file is provided on input. If VCF.GZ is given on input only indexing will be done. Output type can still be chosen with the `output type` command.\n \n\n### Common Issues and Important Notes\n\n * All VCF files must have same sample names, otherwise tool will fail.\n * For option `--allow-overlaps` files must be compressed. \n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using two inputs sized 12.4 MB and 56 KB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
        "label": "Bcftools Concat",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    if (inputs.in_variants.length == 1) {\n        return \"python copy_files.py && echo pass\";\n    } else {\n        return \"python copy_files.py && bcftools concat \";\n    }\n\n}"
          },
          {
            "prefix": "--output",
            "shellQuote": false,
            "position": 6,
            "valueFrom": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.bcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id;}\n\n    if (inputs.output_name) {\n        var out = inputs.output_name;\n        if (inputs.output_type == 'UncompressedVCF') {\n            out += \".concatenated.vcf\";} \n        else if (inputs.output_type == 'CompressedVCF') {\n            out += \".concatenated.vcf.gz\";} \n        else if (inputs.output_type == 'UncompressedBCF') {\n            out += \".concatenated.bcf\";} \n        else if (inputs.output_type == 'CompressedBCF') {\n            out += \".concatenated.bcf.gz\";} \n        else {\n            out += \".concatenated.vcf\";}\n    }\n    \n    else if (inputs.output_type == 'UncompressedVCF') {\n        var out = froot + '.concatenated' + '.vcf';} \n    else if (inputs.output_type == 'CompressedVCF') {\n        var out = froot + '.concatenated' + '.vcf.gz';} \n    else if (inputs.output_type == 'UncompressedBCF') {\n        var out = froot + '.concatenated' + '.bcf';} \n    else if (inputs.output_type == 'CompressedBCF') {\n        var out = froot + '.concatenated' + '.bcf.gz';} \n    else var out = froot + '.concatenated.vcf';\n\n    return out;\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 100,
            "valueFrom": "${  var files_array = [].concat(inputs.in_variants);\n    var in_file = files_array[0];\n    var fname = in_file.basename;\n    var fext = in_file.nameext;\n    var froot = in_file.nameroot;\n    if(fext == '.vcf'){\n        return \"./input_files/*.vcf.gz\"\n    }\n    if(fext == '.bcf'){\n        return \"./input_files/*.bcf\"\n    }\n    else if(fext == '.gz'){\n        if (froot.split('.').pop() == 'bcf'){\n        return \"./input_files/*.bcf.gz\"}\n        else if (froot.split('.').pop() == 'vcf'){\n        return \"./input_files/*.vcf.gz\"}\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.memory) {\n\n        return inputs.memory\n\n    } else {\n        return 1000\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu) {\n        return inputs.cpu\n    } else {\n        return 1\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "21caaa02f72e",
            "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools-1.10.1:1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "file_paths.txt",
                "entry": "${  var i;\n    var text = ''\n    var file_array = [].concat(inputs.in_variants);\n    for (i = 0; i < file_array.length; i++) {\n       text += file_array[i].path + '\\n';}\n    \n    return text;\n}",
                "writable": false
              },
              {
                "entryname": "copy_files.py",
                "entry": "import os\n\nos.system('mkdir input_files')\n\nwith open('./file_paths.txt') as file:\n    for counter, line in enumerate(file):\n        mv_line = line.strip('\\n')\n        basename = line.split('/')[-1].strip('\\n')\nfile.close()\n\n\nif (counter > 0):\n    with open('./file_paths.txt') as txt:\n        for counter, line in enumerate(txt):\n            mv_line = line.strip('\\n')\n            basename = line.split('/')[-1].strip('\\n')\n            if basename.endswith('.vcf.gz'):\n                nameroot = basename.split('.vcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf.gz'\n            elif basename.endswith('.vcf'):\n                nameroot = basename.split('.vcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.vcf'\n            elif basename.endswith('.bcf.gz'):\n                nameroot = basename.split('.bcf.gz')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf.gz'\n            elif basename.endswith('.bcf'):\n                nameroot = basename.split('.bcf')\n                new_name = nameroot[0] + \".\" + str(counter) + '.bcf'\n\n            cmd = \"cp \" + mv_line + \" ./input_files/\" + new_name\n            status = os.system(cmd)\n\n    bash_cmd = \"bash /opt/gzip_vcf.sh\"\n    os.system(bash_cmd)\n\nelse:\n    cmd = \"cp \" + mv_line + \" ./\" + basename\n    print(cmd)\n    os.system(cmd)\n            ",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "var updateMetadata = function(file, key, value) {\n    file['metadata'][key] = value;\n    return file;\n};\n\n\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};\n\nvar toArray = function(file) {\n    return [].concat(file);\n};\n\nvar groupBy = function(files, key) {\n    var groupedFiles = [];\n    var tempDict = {};\n    for (var i = 0; i < files.length; i++) {\n        var value = files[i]['metadata'][key];\n        if (value in tempDict)\n            tempDict[value].push(files[i]);\n        else tempDict[value] = [files[i]];\n    }\n    for (var key in tempDict) {\n        groupedFiles.push(tempDict[key]);\n    }\n    return groupedFiles;\n};\n\nvar orderBy = function(files, key, order) {\n    var compareFunction = function(a, b) {\n        if (a['metadata'][key].constructor === Number) {\n            return a['metadata'][key] - b['metadata'][key];\n        } else {\n            var nameA = a['metadata'][key].toUpperCase();\n            var nameB = b['metadata'][key].toUpperCase();\n            if (nameA < nameB) {\n                return -1;\n            }\n            if (nameA > nameB) {\n                return 1;\n            }\n            return 0;\n        }\n    };\n\n    files = files.sort(compareFunction);\n    if (order == undefined || order == \"asc\")\n        return files;\n    else\n        return files.reverse();\n};",
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:toolkitVersion": "1.10.1",
        "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
        "sbg:license": "MIT License",
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
        "sbg:cmdPreview": "bash gzip_vcf.sh && bcftools concat --output input_file-1.concated.vcf  input_file-1.vcf.gz input_file-2.vcf.gz",
        "sbg:toolkit": "bcftools",
        "sbg:image_url": null,
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982057,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:appVersion": [
          "v1.0"
        ],
        "sbg:id": "markoz/hgi/bcftools-concat-1-10-1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982057,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982057,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "af5855b119c26e49f6aeca0a140859deda105ad0adc92d1ee3fd0eb2d0b2a0731"
      },
      "label": "Bcftools Concat",
      "sbg:x": 2072.980224609375,
      "sbg:y": -347
    },
    {
      "id": "bcftools_reheader_1_10_1",
      "in": [
        {
          "id": "samples_file",
          "source": "samples_file_1"
        },
        {
          "id": "in_variants",
          "source": [
            "tabix_index_1_9_cwl1_0/out_variants"
          ]
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
        "id": "markoz/hgi/bcftools-reheader-1-10-1/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "0",
            "id": "threads",
            "type": "int?",
            "inputBinding": {
              "prefix": "--threads",
              "shellQuote": false,
              "position": 3
            },
            "label": "Threads",
            "doc": "Number of extra compression threads (BCF only)."
          },
          {
            "sbg:category": "Configuration",
            "sbg:includeInPorts": true,
            "id": "new_sample_names",
            "type": "string[]?",
            "label": "New Sample Names",
            "doc": "Strings describing changes: NEW_NAME should be given in New Sample Names in the same order as the OLD_NAME pairs."
          },
          {
            "sbg:category": "Configuration",
            "id": "samples_file",
            "type": "File?",
            "label": "Samples file",
            "doc": "New sample names, one name per line, in the same order as they appear in the VCF file. Alternatively, only samples which need to be renamed can be listed as \"old_name new_name\\n\" pairs separated by whitespaces, each on a separate line. If a sample name contains spaces, the spaces can be escaped using the backslash character, for example \"Not\\ a\\ good\\ sample\\ name\".",
            "sbg:fileTypes": "TXT"
          },
          {
            "sbg:toolDefaultValue": "Derived from input",
            "sbg:altPrefix": "-o",
            "sbg:category": "Configuration",
            "id": "output_name",
            "type": "string?",
            "inputBinding": {
              "prefix": "--output",
              "shellQuote": false,
              "position": 8,
              "valueFrom": "${  var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        if (froot.split('.').pop() == 'vcf'){\n        froot = froot.split('.vcf')[0]}\n        else if (froot.split('.').pop() == 'bcf'){\n        froot = froot.split('.bcf')[0]}\n    }\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id}\n\n    if (inputs.output_name) {\n       var out = inputs.output_name}\n    \n    else var out = froot + '.reheadered.vcf.gz'\n\n    return out\n}"
            },
            "label": "Output file name",
            "doc": "Name of the output file.",
            "default": 0
          },
          {
            "sbg:category": "Configuration",
            "id": "old_sample_names",
            "type": "string[]?",
            "label": "Old Sample Names",
            "doc": "Strings describing changes: OLD_NAME should be given in Old Sample Names in the same order as the NEW_NAME pairs."
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
            "sbg:category": "File Input",
            "id": "in_variants",
            "type": {
              "type": "array",
              "items": "File",
              "inputBinding": {
                "separate": true
              }
            },
            "inputBinding": {
              "shellQuote": false,
              "position": 40,
              "valueFrom": "${\n    var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        return froot + \".gz\"} \n    else {\n        if(fname.split('.').pop() == 'bcf'){\n            var index_csi_file = fname + '.csi'\n            if (in_file.secondaryFiles[0]){\n                var secondary_given = in_file.secondaryFiles[0].path.replace(/^.*[\\\\\\/]/, '');}\n                if(secondary_given == index_csi_file){\n                    return fname;}\n            \n        }\n        return fname + \".gz\"\n    }\n}"
            },
            "label": "Input variants file",
            "doc": "Input file which will be reheadered.",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF, BCF.GZ",
            "secondaryFiles": [
              "${return self.basename + \".tbi\"}"
            ]
          },
          {
            "sbg:category": "File Input",
            "sbg:altPrefix": "-h",
            "id": "header_file",
            "type": "File?",
            "inputBinding": {
              "prefix": "--header",
              "shellQuote": false,
              "position": 13
            },
            "label": "Header file",
            "doc": "File with the new header.",
            "sbg:fileTypes": "VCF, TXT"
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
            "sbg:toolDefaultValue": "False",
            "sbg:category": "Configuration",
            "id": "output_index",
            "type": "boolean?",
            "inputBinding": {
              "shellQuote": false,
              "position": 100,
              "valueFrom": "${ var files_array = [].concat(inputs.in_variants)\n  var fname = files_array[0].basename\n  if(inputs.output_index == true){ \n    return \" && bcftools index  -f -t \" + \"./*.reheadered.vcf.gz\"\n  }\n  else{return \"\"}\n}"
            },
            "label": "Output index file",
            "doc": "If set to True, output file will be indexed."
          },
          {
            "sbg:category": "Configuration",
            "id": "fai",
            "type": "File?",
            "inputBinding": {
              "prefix": "--fai",
              "shellQuote": false,
              "position": 9
            },
            "label": "Fai file",
            "doc": "Update sequences and their lengths from the FAI file.",
            "sbg:fileTypes": "FAI"
          }
        ],
        "outputs": [
          {
            "id": "out_variants",
            "doc": "Reheadered output file.",
            "label": "Output VCF file",
            "type": "File?",
            "outputBinding": {
              "glob": "${  var files_array = [].concat(inputs.in_variants)\n    var in_file = files_array[0]\n    var fname = in_file.basename\n    var fext = in_file.nameext\n    var froot = in_file.nameroot\n    if (fext == '.gz') {\n        var froot = froot.split('.vcf')[0]}\n\n    if(in_file.metadata.sample_id){\n        var froot = in_file.metadata.sample_id}\n\n    if (inputs.output_name) {\n       var out = inputs.output_name}\n    \n    else var out = froot + '.reheadered.vcf.gz'\n\n    return out\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
            },
            "secondaryFiles": [
              ".tbi"
            ],
            "sbg:fileTypes": "VCF.GZ"
          }
        ],
        "doc": "**BCFtools Reheader**: Modify header of VCF/BCF files, and change sample names.\n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n\n### Common Use Cases\n\nChange the header of a VCF file using a header file on the **Header file** (`--header`) input\n```\n$bcftools reheader --header header_file.txt input.vcf.gz\n```\nChange sample names in a VCF file, one name per line, in the same order as they appear in the VCF file provided on the **Samples file** input or a list of strings on the **Sample strings** input.\n\n```\n$bcftools reheader --samples samples_file.txt input.vcf.gz\n```\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands if a VCF file is provided on input. Index file `.tbi` is added as secondary file of `in_variants` input which means if VCF.GZ is provided on input, tool will look for index file in project or previous tool (in case of usage in workflow) and if present, it will not perform nor compressing nor indexing.  If VCF.GZ is given on input only indexing will be done.\n\n### Common Issues and Important Notes\n\n* This tool can only output compressed VCF file.\n\n* Tool doesn't work with BCF files. \n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using a VCF input of 7 MB and header input of 8 KB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
        "label": " Bcftools Reheader",
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
            "valueFrom": "reheader"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 23,
            "valueFrom": "${ //Samples file selection\n    if (inputs.new_sample_names) {\n        return '--samples samples.txt';} \n    else {\n        if (inputs.samples_file) {\n            return \"--samples \" + inputs.samples_file.path;} \n        else {return \"\";}\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_per_job) {\n        return inputs.mem_per_job;\n    } else {\n        return 1000;\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu_per_job) {\n        return inputs.cpu_per_job;\n    } else {\n        return 1;\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools-1.10.1:0"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "samples.txt",
                "entry": "${\n    if (inputs.new_sample_names) {\n\n        var content = ''\n        var samples_strings = [].concat(inputs.new_sample_names)\n        var i;\n        \n        for (i = 0; i < samples_strings.length; i++) {\n            if (inputs.old_sample_names) {\n                content += inputs.old_sample_names[i] + ' '\n            }\n            content += samples_strings[i] + '\\n'\n        }\n\n        return content} \n        \n    else {return ''}\n\n    \n}",
                "writable": false
              },
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
        "sbg:image_url": null,
        "sbg:toolkitVersion": "1.10.1",
        "sbg:toolkit": "bcftools",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "http://samtools.github.io/bcftools/"
          },
          {
            "label": "Source code",
            "id": "https://github.com/samtools/bcftools"
          },
          {
            "label": "Wiki",
            "id": "https://github.com/samtools/bcftools/wiki"
          },
          {
            "label": "Download",
            "id": "https://github.com/samtools/bcftools/archive/develop.zip"
          }
        ],
        "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
        "sbg:license": "MIT licence",
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:cmdPreview": "bcftools index  -f -t input_file.vcf.gz && bcftools reheader  --samples samples.txt  input_file.vcf.gz",
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982061,
            "sbg:revisionNotes": null
          }
        ],
        "sbg:appVersion": [
          "v1.0"
        ],
        "sbg:id": "markoz/hgi/bcftools-reheader-1-10-1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": null,
        "sbg:modifiedOn": 1636982061,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982061,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a38b28a5e8caabc15a90924b36d897eaf36b3ca627dbc155c61198a3e18641af4"
      },
      "label": " Bcftools Reheader",
      "sbg:x": -82.13785552978516,
      "sbg:y": -381.5942687988281
    },
    {
      "id": "vcftools_sort_0_1_14_cwl1",
      "in": [
        {
          "id": "input_file",
          "source": "bcftools_gtc2vcf_1_10/output"
        },
        {
          "id": "compressed",
          "default": true
        }
      ],
      "out": [
        {
          "id": "output_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.1",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/vcftools-sort-0-1-14-cwl1/0",
        "baseCommand": [
          "vcf-sort"
        ],
        "inputs": [
          {
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "separate": false,
              "shellQuote": false,
              "position": 1
            },
            "label": "Input file",
            "doc": "Input file (vcf or vcf.gz)",
            "sbg:fileTypes": "VCF, VCF.GZ"
          },
          {
            "id": "chromosomal_order",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "-c",
              "shellQuote": false,
              "position": 0
            },
            "label": "Chromosomal order",
            "doc": "Use natural ordering (1,2,10,MT,X) rather then the default (1,10,2,MT,X). This requires                                      new version of the unix \"sort\" command which supports the --version-sort option."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "1",
            "id": "parallel",
            "type": "int?",
            "inputBinding": {
              "prefix": "-p",
              "shellQuote": false,
              "position": 2
            },
            "label": "Parallel threads",
            "doc": "Change the number of sorts run concurrently to <int>."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "N/A",
            "id": "mem_mb",
            "type": "int?",
            "label": "Memory in MB",
            "doc": "Memory in MB for execution."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "FALSE",
            "id": "compressed",
            "type": "boolean?",
            "label": "Compressed output",
            "doc": "Check to make the output compressed (usually for further processing)."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "label": "Output file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    filepath = inputs.input_file.path\n\n    filename = filepath.split(\"/\").pop();\n\n    file_dot_sep = filename.split(\".\");\n    file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n\n    if (inputs.compressed) {\n        new_filename += \".gz\";\n    }\n\n    return new_filename;\n\n}",
              "outputEval": "${\n\n    for (var i = 0; i < self.length; i++) {\n        var out_metadata = {\n            '__inherit__': 'input_file'\n        };\n        self[i] = setMetadata(self[i], out_metadata)\n    };\n\n    return self\n\n}"
            },
            "sbg:fileTypes": "VCF, VCF.GZ"
          }
        ],
        "doc": "VCFtools sort sorts a VCF file.",
        "label": "VCFtools Sort",
        "arguments": [
          {
            "shellQuote": false,
            "position": 100,
            "valueFrom": "${\n\n\n    if (inputs.compressed) {\n\n        filepath = inputs.input_file.path\n\n        filename = filepath.split(\"/\").pop();\n\n        file_dot_sep = filename.split(\".\");\n        file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n        new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n        return '&& bgzip -c -f ' + new_filename + ' > ' + new_filename + '.gz'\n\n    }\n\n\n\n}"
          },
          {
            "shellQuote": false,
            "position": 50,
            "valueFrom": "${\n    filepath = inputs.input_file.path\n\n    filename = filepath.split(\"/\").pop();\n\n    file_dot_sep = filename.split(\".\");\n    file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n\n    return '> ' + new_filename;\n\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_mb) {\n\n        return inputs.mem_mb\n\n    } else {\n\n        return 1000\n\n    }\n}",
            "coresMin": "${\n    if (inputs.parallel) {\n        return inputs.parallel\n    } else {\n        return 1\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/ognjenm/vcftools:0.1.14"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "vcf-sort sample1.vcf  > sample1.sorted.vcf  && bgzip -c -f sample1.sorted.vcf > sample1.sorted.vcf.gz",
        "sbg:toolkitVersion": "0.1.14",
        "sbg:license": "GNU General Public License version 3.0 (GPLv3)",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "https://vcftools.github.io"
          },
          {
            "label": "Source code",
            "id": "https://github.com/vcftools/vcftools"
          },
          {
            "label": "Publications",
            "id": "http://bioinformatics.oxfordjournals.org/content/27/15/2156"
          }
        ],
        "sbg:toolkit": "VCFtools",
        "sbg:toolAuthor": "Adam Auton, Petr Danecek, Anthony Marcketta",
        "sbg:appVersion": [
          "v1.1"
        ],
        "sbg:id": "markoz/hgi/vcftools-sort-0-1-14-cwl1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-sort-0-1-14",
        "sbg:modifiedOn": 1636982240,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982240,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:projectName": "HGI",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982240,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-sort-0-1-14"
          }
        ],
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a8be3a0977da128b3cf37fc288b3ed90823dc54be3bcc4992f7e5607de593c9a2"
      },
      "label": "VCFtools Sort",
      "sbg:x": -472.26080322265625,
      "sbg:y": -392.45062255859375
    },
    {
      "id": "ebi_vcf_validator_cwl1",
      "in": [
        {
          "id": "input_vcf",
          "source": "vcftools_hwe_0_1_14_cwl1/output_file"
        },
        {
          "id": "report_types",
          "default": [
            "database",
            "summary",
            "text"
          ]
        }
      ],
      "out": [
        {
          "id": "text_out"
        },
        {
          "id": "db_out"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.1",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/ebi-vcf-validator-cwl1/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:altPrefix": "-i",
            "id": "input_vcf",
            "type": "File",
            "label": "Input VCF file",
            "doc": "VCF file to validate.",
            "sbg:fileTypes": "VCF, VCF.GZ"
          },
          {
            "sbg:altPrefix": "-l",
            "sbg:toolDefaultValue": "warning",
            "id": "validation_level",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "error",
                  "warning",
                  "stop"
                ],
                "name": "validation_level"
              }
            ],
            "inputBinding": {
              "prefix": "--level",
              "shellQuote": false,
              "position": 2
            },
            "label": "Validation level",
            "doc": "Validation level can be \"error\" - in which case Validator will display only syntax errors; \"warning\" - when both errors and warnings are shown; or \"stop\" - in which case the app will stop after the first syntax error is found."
          },
          {
            "sbg:altPrefix": "-r",
            "sbg:toolDefaultValue": "summary",
            "id": "report_types",
            "type": "string[]?",
            "inputBinding": {
              "prefix": "--report",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 3
            },
            "label": "Types of reports",
            "doc": "Comma separated values for types of reports (allowed values: \"summary\" for a human-readable summary, \"text\" for a full human-readable report or \"database\" for saving a structured report into a database file)."
          }
        ],
        "outputs": [
          {
            "id": "text_out",
            "doc": "VCF validator reports.",
            "label": "VCF validator reports",
            "type": "File[]?",
            "outputBinding": {
              "glob": "*.txt"
            },
            "sbg:fileTypes": "TXT, DB"
          },
          {
            "id": "db_out",
            "type": "File?",
            "outputBinding": {
              "glob": "*.db"
            }
          }
        ],
        "doc": "**EBI vcf-validator** is a validation suite for VCF files, including checks from vcftools and additional lexical, syntactic, and semantic analysis of VCF files [1].\n\n### Common Use Cases\n\n**EBI vcf-validator** can be used to ensure that a VCF file adheres to VCF specification. The tool outputs reports which can be used to manually address any potential inconsistencies in the input VCF, or as input (if produced as a DB file) to the **EBI vcf-debugulator**  tool.\n\n### Changes Introduced by Seven Bridges\n\n* `--outdir` flag was hardcoded in the tool wrapper.\n\n### Common Issues and Important Notes\n\n* Input **Input VCF file** is required and should be a VCF or VCF.GZ file.\n\n### Performance Benchmarking\n\nTypical runs take <5 minutes and cost <$0.05.\n\n### References\n\n[1] [EBI vcf-validator GitHub page](https://github.com/EBIvariation/vcf-validator)",
        "label": "EBI vcf-validator",
        "arguments": [
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    filename = inputs.input_vcf.path.split('/').pop()\n    ext = filename.split('.').pop()\n    if (ext == 'gz') {\n        return \"zcat \".concat(inputs.input_vcf.path, \" | /opt/vcf-validator-0.7/build/bin/vcf_validator\")\n    } else {\n        return \"/opt/vcf-validator-0.7/build/bin/vcf_validator -i \".concat(inputs.input_vcf.path)\n    }\n}"
          },
          {
            "prefix": "-o",
            "shellQuote": false,
            "position": 4,
            "valueFrom": "./"
          },
          {
            "separate": false,
            "shellQuote": false,
            "position": 102,
            "valueFrom": "${\n    return \" ; echo Validation Done\"\n}"
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
            "dockerPull": "images.sbgenomics.com/jrandjelovic/ebi-vcf-validator:v1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "/opt/vcf-validator-0.7/build/bin/vcf_validator -i /path/to/input_vcf.vcf -o ./  ; echo Validation Done",
        "sbg:toolkitVersion": "0.7",
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:license": "Apache License 2.0",
        "sbg:links": [
          {
            "label": "Source Code",
            "id": "https://github.com/EBIvariation/vcf-validator"
          }
        ],
        "sbg:toolkit": "EBI vcf-validator",
        "sbg:toolAuthor": "The European Bioinformatics Institute",
        "sbg:appVersion": [
          "v1.1"
        ],
        "sbg:id": "markoz/hgi/ebi-vcf-validator-cwl1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/ebi-vcf-validator",
        "sbg:modifiedOn": 1636982206,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982206,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:projectName": "HGI",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982206,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/ebi-vcf-validator"
          }
        ],
        "sbg:publisher": "sbg",
        "sbg:content_hash": "afc2e83637bf47f9c38f360ad080c97d9c470817a1273c5313cacd818cb564c2b"
      },
      "label": "EBI vcf-validator",
      "sbg:x": 628.006591796875,
      "sbg:y": -399.0463562011719
    },
    {
      "id": "ebi_vcf_debugulator_cwl1",
      "in": [
        {
          "id": "input_vcf",
          "source": "vcftools_hwe_0_1_14_cwl1/output_file"
        },
        {
          "id": "errors_report_db_file",
          "source": "ebi_vcf_validator_cwl1/db_out"
        }
      ],
      "out": [
        {
          "id": "corrected_vcf_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.1",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/ebi-vcf-debugulator-cwl1/0",
        "baseCommand": [],
        "inputs": [
          {
            "id": "input_vcf",
            "type": "File",
            "label": "Input VCF file",
            "doc": "VCF file which debugulator should attempt to correct.",
            "sbg:fileTypes": "VCF, VCF.GZ"
          },
          {
            "sbg:altPrefix": "-l",
            "sbg:toolDefaultValue": "warning",
            "id": "validation_level",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "error",
                  "warning",
                  "stop"
                ],
                "name": "validation_level"
              }
            ],
            "inputBinding": {
              "prefix": "--level",
              "shellQuote": false,
              "position": 3
            },
            "label": "Validation level",
            "doc": "Validation level (error, warning, stop)."
          },
          {
            "sbg:altPrefix": "-e",
            "id": "errors_report_db_file",
            "type": "File",
            "inputBinding": {
              "prefix": "--errors",
              "shellQuote": false,
              "position": 3
            },
            "label": "Error reports file from EBI vcf-validator",
            "doc": "Path to the errors report produced by EBI vcf-validator from the input VCF.",
            "sbg:fileTypes": "DB"
          }
        ],
        "outputs": [
          {
            "id": "corrected_vcf_file",
            "doc": "Corrected VCF file.",
            "label": "Corrected VCF file",
            "type": "File?",
            "outputBinding": {
              "glob": "*_debugged.vcf",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_vcf)\n\n}"
            },
            "sbg:fileTypes": "VCF"
          }
        ],
        "doc": "**EBI vcf-debugulator** will attempt to correct inconsistencies in input VCF files detected using **EBI vcf-validator**. \n\n### Common Use Cases\n\n**EBI vcf-debugulator** can be used to attempt to bring misbehaving VCFs closer to the VCF specification, based on the error reports (DB format) generated using **EBI vcf-validator**.\n\n### Changes Introduced by Seven Bridges\n\n* `--output` flag is hardcoded in the tool wrapper.\n\n### Common Issues and Important Notes\n\n* Inputs **Input VCF file** and **Error reports file from EBI vcf-validator** are required inputs. \n* Input **Error reports file from EBI vcf-validator** should be a DB file (not a TXT error report or error summary).",
        "label": "EBI vcf-debugulator",
        "arguments": [
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    filename = inputs.input_vcf.path.split('/').pop()\n    ext = filename.split('.').pop()\n    if (ext == 'gz') {\n        return \"zcat \".concat(inputs.input_vcf.path, \" | /opt/vcf-validator-0.7/build/bin/vcf_debugulator\")\n    } {\n        return \"/opt/vcf-validator-0.7/build/bin/vcf_debugulator -i \".concat(inputs.input_vcf.path)\n    }\n}"
          },
          {
            "prefix": "-o",
            "shellQuote": false,
            "position": 102,
            "valueFrom": "${\n    filename = inputs.input_vcf.path.split('/').pop()\n    primename = filename.split('.vcf')[0]\n    return primename.concat('_debugged.vcf')\n}"
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
            "dockerPull": "images.sbgenomics.com/jrandjelovic/ebi-vcf-validator:v1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "/opt/vcf-validator-0.7/build/bin/vcf_debugulator -i /path/to/input_vcf.vcf --errors /path/to/errors_report_db_file.ext -o input_vcf_debugged.vcf",
        "sbg:toolkitVersion": "0.7",
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:license": "Apache License 2.0",
        "sbg:links": [
          {
            "label": "Source Code",
            "id": "https://github.com/EBIvariation/vcf-validator"
          }
        ],
        "sbg:toolkit": "EBI vcf-validator",
        "sbg:toolAuthor": "The European Bioinformatics Institute",
        "sbg:appVersion": [
          "v1.1"
        ],
        "sbg:id": "markoz/hgi/ebi-vcf-debugulator-cwl1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/ebi-vcf-debugulator",
        "sbg:modifiedOn": 1636982297,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982297,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:projectName": "HGI",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982297,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/ebi-vcf-debugulator"
          }
        ],
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a4a7c4d22fed3db15b0eda89fbf65cf5ec84b8471c0a3d41cc7abfe0579537f5b"
      },
      "label": "EBI vcf-debugulator",
      "sbg:x": 797.25830078125,
      "sbg:y": -481.4967041015625
    },
    {
      "id": "bcftools_index_cwl1",
      "in": [
        {
          "id": "input_file",
          "source": "tabix_bgzip_1_9_cwl1_0/output_file"
        },
        {
          "id": "output_vcf_with_index",
          "default": true
        }
      ],
      "out": [
        {
          "id": "output_file"
        },
        {
          "id": "index_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.1",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/bcftools-index-cwl1/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:category": "File Input",
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 43,
              "valueFrom": "${\n    fname = inputs.input_file.path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() != 'gz') {\n        fname = inputs.input_file.path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n        if (inputs.tbi_index) {\n            return fname + \".vcf.gz\"\n        } else {\n            return fname + \".vcf\"\n        }\n    } else {\n\n        return fname\n\n    }\n}"
            },
            "label": "Input file",
            "doc": "Input file.",
            "sbg:fileTypes": "VCF.GZ, VCF"
          },
          {
            "sbg:altPrefix": "-o",
            "sbg:category": "Indexing options",
            "id": "output_name",
            "type": "string?",
            "inputBinding": {
              "prefix": "--output-file",
              "shellQuote": false,
              "position": 7
            },
            "label": "Output file name",
            "doc": "Output file name. If not set, then the index will be created using the input file name plus a .csi or .tbi extension"
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
            "doc": "Number of threads. Number of output compression threads to use in addition to main thread. Only used when output type is CompressedBCF CompressedVCF."
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
            "sbg:altPrefix": "-c",
            "sbg:category": "Indexing options",
            "id": "csi_index",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--csi",
              "shellQuote": false,
              "position": 4
            },
            "label": "Generate CSI Index",
            "doc": "Generate CSI-format index for VCF/BCF files [default]."
          },
          {
            "sbg:altPrefix": "-f",
            "sbg:category": "Indexing options",
            "id": "force",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--force",
              "shellQuote": false,
              "position": 5
            },
            "label": "Overwrite index",
            "doc": "Overwrite index if it already exists"
          },
          {
            "sbg:toolDefaultValue": "14",
            "sbg:altPrefix": "-m",
            "sbg:category": "Indexing options",
            "id": "min_interval",
            "type": "int?",
            "inputBinding": {
              "prefix": "--min-shift",
              "shellQuote": false,
              "position": 6
            },
            "label": "Minimal interval",
            "doc": "Set minimal interval size for CSI indices to 2^INT."
          },
          {
            "sbg:altPrefix": "-t",
            "sbg:category": "Indexing options",
            "id": "tbi_index",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--tbi",
              "shellQuote": false,
              "position": 8
            },
            "label": "Generate TBI Index",
            "doc": "Generate TBI-format index for VCF files."
          },
          {
            "sbg:altPrefix": "-n",
            "sbg:category": "Stats options",
            "id": "nrecords",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--nrecords",
              "shellQuote": false,
              "position": 9
            },
            "label": "Number of records",
            "doc": "Print the number of records based on the CSI or TBI index files."
          },
          {
            "sbg:altPrefix": "-s",
            "sbg:category": "Stats options",
            "id": "stats",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--stats",
              "shellQuote": false,
              "position": 11
            },
            "label": "Print stats",
            "doc": "Print per contig stats based on the CSI or TBI index files. Output format is three tab-delimited columns listing the contig name, contig length (. if unknown) and number of records for the contig. Contigs with zero records are not printed."
          },
          {
            "sbg:toolDefaultValue": "False",
            "sbg:category": "General options",
            "id": "output_vcf_with_index",
            "type": "boolean?",
            "label": "Output VCF file with index",
            "doc": "Output VCF file with index."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "label": "Output file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    fname = inputs.input_file.path.replace(/^.*[\\\\\\/]/, '')\n    if (inputs.output_vcf_with_index == true) {\n        if (inputs.tbi_index) {\n            return \"*.vcf.gz\"\n        } else if (fname.split('.').pop().toLowerCase() == 'gz') {\n            return \"*.vcf.gz\"\n        } else {\n            return \"*.vcf\"\n        }\n    } else {\n        return \"\"\n    }\n\n}",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
            },
            "secondaryFiles": [
              {
                "pattern": ".tbi",
                "required": false
              },
              {
                "pattern": ".csi",
                "required": false
              }
            ],
            "sbg:fileTypes": "VCF, VCF.GZ"
          },
          {
            "id": "index_file",
            "doc": "Index file only",
            "label": "Index file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n\n    if (inputs.tbi_index) {\n\n        return \"*.tbi\"\n\n    } else {\n\n        return \"*.csi\"\n\n    }\n\n\n\n\n\n}",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
            },
            "sbg:fileTypes": "TBI, CSI"
          }
        ],
        "doc": "**BCFtools Index**: Creates index for bgzip-compressed VCF/BCF files for random access. CSI (coordinate-sorted index) is created by default. \n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n### Common Use Cases\n\nThe **CSI** format supports indexing of chromosomes up to 2^31 in length with **Generate CSI Index** (`--csi`). When loading an index file, BCFtools will try the **CSI** first and then the **TBI**.\n```\n$bcftools index --csi file.vcf.gz\n```\n**TBI** (tabix index) index files, which support chromosome lengths up to 2^29, can be created by using the **Generate TBI Index** (`--tbi`) option or using the tabix program packaged with **htslib**.\n```\n$bcftools index --tbi file.vcf.gz\n```\n\n### Changes Introduced by Seven Bridges\n\n* Added the **Output VCF file with index** option that allows you to get the index file along with the VCF on the output, to be used in tools that require VCFs with a corresponding secondary index file.\n\n### Common Issues and Important Notes\n\n * No common issues specific to the tool's execution on the Seven Bridges Platform have been detected.\n\n### Performance Benchmarking\n\nIt took 3 minutes to execute this tool on AWS c4.2xlarge instance using an input of 7 MB. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
        "label": "Bcftools Index",
        "arguments": [
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    fname = inputs.input_file.path.replace(/^.*[\\\\\\/]/, '')\n    if (fname.split('.').pop().toLowerCase() != 'gz') {\n        fname = inputs.input_file.path.replace(/^.*[\\\\\\/]/, '').replace(/\\.[^/.]+$/, \"\")\n        if (inputs.tbi_index) {\n            return \"bgzip -c -f \" + fname + \".vcf > \" + fname + \".vcf.gz &&\"\n        } else {\n            return \"\"\n        }\n    } else {\n\n        return \"\"\n\n    }\n}"
          },
          {
            "shellQuote": false,
            "position": 1,
            "valueFrom": "bcftools"
          },
          {
            "shellQuote": false,
            "position": 2,
            "valueFrom": "index"
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
        "abg:revisionNotes": "Changed output glob",
        "sbg:image_url": null,
        "sbg:license": "MIT License",
        "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:toolkit": "bcftools",
        "sbg:cmdPreview": "bgzip -c -f annotated_input_file.vcf > annotated_input_file.vcf.gz && bcftools index  annotated_input_file.vcf.gz",
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
        "sbg:appVersion": [
          "v1.1"
        ],
        "sbg:id": "markoz/hgi/bcftools-index-cwl1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/bcftools-index",
        "sbg:modifiedOn": 1636982200,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982200,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:projectName": "HGI",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982200,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/bcftools-index"
          }
        ],
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a103de0e186add7d66fcb4b39a3fcce5a0e29bfec552ef5e05694ccf22e7cf607"
      },
      "label": "Bcftools Index",
      "sbg:x": 1168.471435546875,
      "sbg:y": -332.7357177734375
    },
    {
      "id": "vcftools_sort_0_1_14_cwl2",
      "in": [
        {
          "id": "input_file",
          "source": "bcftools_concat_1_10_2/out_variants"
        },
        {
          "id": "parallel",
          "default": 7
        },
        {
          "id": "compressed",
          "default": true
        }
      ],
      "out": [
        {
          "id": "output_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.1",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/vcftools-sort-0-1-14-cwl1/0",
        "baseCommand": [
          "vcf-sort"
        ],
        "inputs": [
          {
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "separate": false,
              "shellQuote": false,
              "position": 1
            },
            "label": "Input file",
            "doc": "Input file (vcf or vcf.gz)",
            "sbg:fileTypes": "VCF, VCF.GZ"
          },
          {
            "id": "chromosomal_order",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "-c",
              "shellQuote": false,
              "position": 0
            },
            "label": "Chromosomal order",
            "doc": "Use natural ordering (1,2,10,MT,X) rather then the default (1,10,2,MT,X). This requires                                      new version of the unix \"sort\" command which supports the --version-sort option."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "1",
            "id": "parallel",
            "type": "int?",
            "inputBinding": {
              "prefix": "-p",
              "shellQuote": false,
              "position": 2
            },
            "label": "Parallel threads",
            "doc": "Change the number of sorts run concurrently to <int>."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "N/A",
            "id": "mem_mb",
            "type": "int?",
            "label": "Memory in MB",
            "doc": "Memory in MB for execution."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "FALSE",
            "id": "compressed",
            "type": "boolean?",
            "label": "Compressed output",
            "doc": "Check to make the output compressed (usually for further processing)."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "label": "Output file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    filepath = inputs.input_file.path\n\n    filename = filepath.split(\"/\").pop();\n\n    file_dot_sep = filename.split(\".\");\n    file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n\n    if (inputs.compressed) {\n        new_filename += \".gz\";\n    }\n\n    return new_filename;\n\n}",
              "outputEval": "${\n\n    for (var i = 0; i < self.length; i++) {\n        var out_metadata = {\n            '__inherit__': 'input_file'\n        };\n        self[i] = setMetadata(self[i], out_metadata)\n    };\n\n    return self\n\n}"
            },
            "sbg:fileTypes": "VCF, VCF.GZ"
          }
        ],
        "doc": "VCFtools sort sorts a VCF file.",
        "label": "VCFtools Sort",
        "arguments": [
          {
            "shellQuote": false,
            "position": 100,
            "valueFrom": "${\n\n\n    if (inputs.compressed) {\n\n        filepath = inputs.input_file.path\n\n        filename = filepath.split(\"/\").pop();\n\n        file_dot_sep = filename.split(\".\");\n        file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n        new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n        return '&& bgzip -c -f ' + new_filename + ' > ' + new_filename + '.gz'\n\n    }\n\n\n\n}"
          },
          {
            "shellQuote": false,
            "position": 50,
            "valueFrom": "${\n    filepath = inputs.input_file.path\n\n    filename = filepath.split(\"/\").pop();\n\n    file_dot_sep = filename.split(\".\");\n    file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    new_filename = filename.substr(0, filename.lastIndexOf(\".vcf\")) + \".sorted.vcf\";\n\n    return '> ' + new_filename;\n\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_mb) {\n\n        return inputs.mem_mb\n\n    } else {\n\n        return 1000\n\n    }\n}",
            "coresMin": "${\n    if (inputs.parallel) {\n        return inputs.parallel\n    } else {\n        return 1\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/ognjenm/vcftools:0.1.14"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "vcf-sort sample1.vcf  > sample1.sorted.vcf  && bgzip -c -f sample1.sorted.vcf > sample1.sorted.vcf.gz",
        "sbg:toolkitVersion": "0.1.14",
        "sbg:license": "GNU General Public License version 3.0 (GPLv3)",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "https://vcftools.github.io"
          },
          {
            "label": "Source code",
            "id": "https://github.com/vcftools/vcftools"
          },
          {
            "label": "Publications",
            "id": "http://bioinformatics.oxfordjournals.org/content/27/15/2156"
          }
        ],
        "sbg:toolkit": "VCFtools",
        "sbg:toolAuthor": "Adam Auton, Petr Danecek, Anthony Marcketta",
        "sbg:appVersion": [
          "v1.1"
        ],
        "sbg:id": "markoz/hgi/vcftools-sort-0-1-14-cwl1/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-sort-0-1-14",
        "sbg:modifiedOn": 1636982240,
        "sbg:modifiedBy": "marko_zecevic",
        "sbg:createdOn": 1636982240,
        "sbg:createdBy": "marko_zecevic",
        "sbg:project": "markoz/hgi",
        "sbg:projectName": "HGI",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "marko_zecevic"
        ],
        "sbg:latestRevision": 0,
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1636982240,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-sort-0-1-14"
          }
        ],
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a8be3a0977da128b3cf37fc288b3ed90823dc54be3bcc4992f7e5607de593c9a2"
      },
      "label": "VCFtools Sort",
      "sbg:x": 2200.52294921875,
      "sbg:y": -502.03973388671875
    },
    {
      "id": "vcftools_max_missing_cwl1",
      "in": [
        {
          "id": "input_file",
          "source": "bcftools_view/out_variants"
        },
        {
          "id": "max-missing",
          "source": "max-missing"
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
        "id": "markoz/hgi/vcftools-max-missing-cwl1/1",
        "baseCommand": [
          "vcftools"
        ],
        "inputs": [
          {
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 3,
              "valueFrom": "${\n    // sufix = \"_CNVs\";\n    // sufix_ext = \"txt\";\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    var file_dot_sep = filename.split(\".\");\n\n    var basename = file_dot_sep[0]\n    var file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    var prefix = \"\"\n\n    if (file_ext == \"vcf\") {\n        var prefix = \"--vcf \"\n    }\n    if (file_ext == \"bcf\") {\n        var prefix = \"--bcf --gatk \"\n    }\n    if (file_ext == \"gz\") {\n        var prefix = \"--gzvcf \"\n    }\n\n    // new_filename = basename + \".analyzed.hwe\";\n\n    return prefix + filepath;\n}"
            },
            "label": "Input file",
            "doc": "Input file (vcf, vcf.gz, bcf)",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF"
          },
          {
            "sbg:category": "Execution",
            "id": "max-missing",
            "type": "float",
            "inputBinding": {
              "prefix": "--max-missing",
              "shellQuote": false,
              "position": 2
            },
            "label": "Threshold variant call rate",
            "doc": "Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed)."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "doc": "Analyzed VCF file.",
            "label": "Analyzed output file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    if (filename.lastIndexOf(\".vcf.gz\") != -1) {\n        var basename = filename.substr(0, filename.lastIndexOf(\".vcf.gz\"))\n    } else {\n        var basename = filename.substr(0, filename.lastIndexOf(\".\"))\n    }\n\n    var new_filename = basename + \".filtered.recode.vcf\";\n\n    return new_filename;\n}",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
            },
            "sbg:fileTypes": "VCF"
          }
        ],
        "doc": "VCFtools max-missing filters a multi-sample VCF on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).",
        "label": "VCFtools max-missing",
        "arguments": [
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "--recode"
          },
          {
            "prefix": "--out",
            "shellQuote": false,
            "position": 101,
            "valueFrom": "${\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    if (filename.lastIndexOf(\".vcf.gz\") != -1) {\n        var basename = filename.substr(0, filename.lastIndexOf(\".vcf.gz\"))\n    } else {\n        var basename = filename.substr(0, filename.lastIndexOf(\".\"))\n    }\n\n    var new_filename = basename + \".filtered\";\n\n    return new_filename;\n}"
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
            "dockerPull": "images.sbgenomics.com/ognjenm/vcftools:0.1.14"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "vcftools --recode --hwe 0 --vcf sample.vcf --out sample.analyzed",
        "sbg:toolkitVersion": "0.1.14",
        "sbg:license": "GNU General Public License version 3.0 (GPLv3)",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "https://vcftools.github.io"
          },
          {
            "label": "Source code",
            "id": "https://github.com/vcftools/vcftools"
          },
          {
            "label": "Publications",
            "id": "http://bioinformatics.oxfordjournals.org/content/27/15/2156"
          }
        ],
        "sbg:toolkit": "VCFtools",
        "sbg:toolAuthor": "Adam Auton, Petr Danecek, Anthony Marcketta",
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1637012661,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-max-missing"
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637012837,
            "sbg:revisionNotes": ""
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/vcftools-max-missing-cwl1/1",
        "sbg:revision": 1,
        "sbg:revisionNotes": "",
        "sbg:modifiedOn": 1637012837,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1637012661,
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
        "sbg:content_hash": "a914eeb5d4f67723e132b441449ee7db56f9f21de5af937fbafcda4a47af5c2b0"
      },
      "label": "VCFtools max-missing",
      "sbg:x": 351.8326416015625,
      "sbg:y": -380
    },
    {
      "id": "vcftools_het",
      "in": [
        {
          "id": "input_file",
          "source": "bcftools_index_cwl1/output_file"
        }
      ],
      "out": [
        {
          "id": "het_file"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/vcftools-het/4",
        "baseCommand": [
          "vcftools"
        ],
        "inputs": [
          {
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 3,
              "valueFrom": "${\n    // sufix = \"_CNVs\";\n    // sufix_ext = \"txt\";\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    var file_dot_sep = filename.split(\".\");\n\n    var basename = file_dot_sep[0]\n    var file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    var prefix = \"\"\n\n    if (file_ext == \"vcf\") {\n        var prefix = \"--vcf \"\n    }\n    if (file_ext == \"bcf\") {\n        var prefix = \"--bcf --gatk \"\n    }\n    if (file_ext == \"gz\") {\n        var prefix = \"--gzvcf \"\n    }\n\n    // new_filename = basename + \".analyzed.hwe\";\n\n    return prefix + filepath;\n}"
            },
            "label": "Input file",
            "doc": "Input file (vcf, vcf.gz, bcf)",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF"
          }
        ],
        "outputs": [
          {
            "id": "het_file",
            "doc": "HET file",
            "label": "HET file",
            "type": "File?",
            "outputBinding": {
              "glob": "*.het"
            },
            "sbg:fileTypes": "HET"
          }
        ],
        "doc": "VCFtools het calculates a measure of heterozygosity on a per-individual basis. Specfically, the inbreeding coefficient, F, is estimated for each individual using a method of moments. The resulting file has the suffix \".het\".",
        "label": "VCFtools het",
        "arguments": [
          {
            "prefix": "--out",
            "shellQuote": false,
            "position": 101,
            "valueFrom": "${\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    if (filename.lastIndexOf(\".vcf.gz\") != -1) {\n        var basename = filename.substr(0, filename.lastIndexOf(\".vcf.gz\"))\n    } else {\n        var basename = filename.substr(0, filename.lastIndexOf(\".\"))\n    }\n\n    var new_filename = basename\n\n    return new_filename;\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 10,
            "valueFrom": "--het"
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
            "dockerPull": "images.sbgenomics.com/ognjenm/vcftools:0.1.14"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "vcftools --recode --hwe 0 --vcf sample.vcf --out sample.analyzed",
        "sbg:toolkitVersion": "0.1.14",
        "sbg:license": "GNU General Public License version 3.0 (GPLv3)",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "https://vcftools.github.io"
          },
          {
            "label": "Source code",
            "id": "https://github.com/vcftools/vcftools"
          },
          {
            "label": "Publications",
            "id": "http://bioinformatics.oxfordjournals.org/content/27/15/2156"
          }
        ],
        "sbg:toolkit": "VCFtools",
        "sbg:toolAuthor": "Adam Auton, Petr Danecek, Anthony Marcketta",
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637685185,
            "sbg:revisionNotes": "Copy of markoz/hgi/vcftools-max-missing-cwl1/1"
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637685811,
            "sbg:revisionNotes": ""
          },
          {
            "sbg:revision": 2,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637686665,
            "sbg:revisionNotes": ""
          },
          {
            "sbg:revision": 3,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637687590,
            "sbg:revisionNotes": ""
          },
          {
            "sbg:revision": 4,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1637688917,
            "sbg:revisionNotes": ""
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/vcftools-het/4",
        "sbg:revision": 4,
        "sbg:revisionNotes": "",
        "sbg:modifiedOn": 1637688917,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1637685185,
        "sbg:createdBy": "markoz",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "markoz"
        ],
        "sbg:latestRevision": 4,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a470ecf4d44a133ca5839ae6deab238dd9e10370caf6af12c6ebb0f6111f88cc2"
      },
      "label": "VCFtools het",
      "sbg:x": 1418.146240234375,
      "sbg:y": -201.28431701660156
    },
    {
      "id": "tabix_index_1_9_cwl1_0",
      "in": [
        {
          "id": "preset",
          "default": "vcf.gz"
        },
        {
          "id": "in_variants",
          "source": "vcftools_sort_0_1_14_cwl1/output_file"
        }
      ],
      "out": [
        {
          "id": "out_variants"
        },
        {
          "id": "out_index"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/tabix-index-1-9-cwl1-0/0",
        "baseCommand": [],
        "inputs": [
          {
            "sbg:toolDefaultValue": "gff.gz",
            "sbg:category": "Config Inputs",
            "id": "preset",
            "type": {
              "type": "enum",
              "symbols": [
                "gff.gz",
                "bed.gz",
                "sam.gz",
                "vcf.gz"
              ],
              "name": "preset"
            },
            "inputBinding": {
              "shellQuote": false,
              "position": 1,
              "valueFrom": "${\n  if (inputs.index_file==undefined && inputs.preset!=undefined && inputs.preset!='')\n  {\n    return ' -p ' + inputs.preset.substring(0, inputs.preset.length-3)\n  }\n  else return ''\n}"
            },
            "label": "Select input file format",
            "doc": "Select input file format."
          },
          {
            "sbg:category": "File Inputs",
            "id": "in_variants",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 99,
              "valueFrom": "${\n    inp = [].concat(inputs.in_variants)[0]\n  if (inputs.index_file==undefined)\n  \treturn inp.path.split('/').slice(-1)[0]\n  return ''\n}"
            },
            "label": "Input file",
            "doc": "Input file for tabix indexing.",
            "sbg:fileTypes": "GFF.GZ, BED.GZ, SAM.GZ, VCF.GZ, PSLTBL.GZ"
          },
          {
            "sbg:toolDefaultValue": "0",
            "sbg:category": "Config Inputs",
            "id": "skip_lines",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 5,
              "valueFrom": "${\n  if (inputs.index_file==undefined && inputs.skip_lines!=undefined)\n    return \" -S \" + inputs.skip_lines\n}"
            },
            "label": "Skip first N lines",
            "doc": "Skip first N lines in the data file."
          },
          {
            "sbg:toolDefaultValue": "5",
            "sbg:category": "Config Inputs",
            "id": "end",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 4,
              "valueFrom": "${\n  if (inputs.index_file!=undefined && inputs.end!=undefined)\n    return \" -e \" + inputs.end\n  return ''\n}"
            },
            "label": "Column number for region end",
            "doc": "Column number for region end (if no end, set INT to -b)."
          },
          {
            "sbg:toolDefaultValue": "1",
            "sbg:category": "Config Inputs",
            "id": "sequence",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 2,
              "valueFrom": "${\n  if (inputs.index_file!=undefined && inputs.sequence!=undefined)\n    return \" -s \" + inputs.sequence\n  return ''\n}"
            },
            "label": "Column number for sequence names",
            "doc": "Column number for sequence names (suppressed by -p)."
          },
          {
            "sbg:toolDefaultValue": "4",
            "sbg:category": "Config Inputs",
            "id": "begin",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 3,
              "valueFrom": "${\n  if (inputs.index_file!=undefined && inputs.begin!=undefined)\n    return \" -b \" + inputs.begin\n  return ''\n}"
            },
            "label": "Column number for region start",
            "doc": "Column number for region start."
          },
          {
            "sbg:toolDefaultValue": "#",
            "sbg:category": "Config Inputs",
            "id": "comment",
            "type": "string?",
            "inputBinding": {
              "shellQuote": false,
              "position": 6,
              "valueFrom": "${\n  if (inputs.index_file!=undefined && inputs.comment!=undefined && inputs.comment != '')\n    return \" -c \" + inputs.comment\n  return ''\n}"
            },
            "label": "Skip comment lines starting with character CHAR",
            "doc": "Skip comment lines starting with character CHAR."
          },
          {
            "sbg:toolDefaultValue": "Default is not defined in the tool.",
            "sbg:category": "Config Inputs",
            "id": "zero_based",
            "type": "boolean?",
            "inputBinding": {
              "shellQuote": false,
              "position": 7,
              "valueFrom": "${\n  if (inputs.index_file!=undefined && inputs.zero_based==true)\n    return \" -0 \"\n  return ''\n}"
            },
            "label": "Specify if the position in the data file is 0 based",
            "doc": "Specify if the position in the data file is 0 based."
          },
          {
            "sbg:toolDefaultValue": "2048",
            "sbg:category": "Platform Options",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Reserve N MB of RAM",
            "doc": "Reserve N MB of RAM for tool execution.."
          },
          {
            "sbg:category": "File inputs",
            "id": "index_file",
            "type": "File?",
            "label": "Index file",
            "doc": "Index file.",
            "sbg:fileTypes": "TBI, CSI"
          },
          {
            "sbg:toolDefaultValue": "False",
            "sbg:category": "Config Inputs",
            "id": "dont_output_data_file",
            "type": "boolean?",
            "label": "Don't output data file",
            "doc": "Don't output data file (only index file will be outputed)."
          },
          {
            "sbg:toolDefaultValue": "1",
            "sbg:category": "Platform Options",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "Cpus to be used on the platform",
            "doc": "Cpus to be used on the platform."
          },
          {
            "sbg:toolDefaultValue": "False",
            "sbg:category": "Config Inputs",
            "id": "csi",
            "type": "boolean?",
            "inputBinding": {
              "shellQuote": false,
              "position": 7,
              "valueFrom": "${\n    if (inputs.index_file==undefined && inputs.csi != undefined && inputs.csi == true)\n    {\n        return ' --csi '\n    }\n    return ''\n}"
            },
            "label": "Generate CSI index for VCF",
            "doc": "Generate CSI index for VCF (default is TBI)."
          },
          {
            "sbg:toolDefaultValue": "14",
            "sbg:category": "Config Inputs",
            "id": "min_shift",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 0,
              "valueFrom": "${\n    if(inputs.index_file==undefined && inputs.min_shift != undefined)\n    {\n        return ' -m ' + inputs.min_shift\n    }\n    return ''\n}"
            },
            "label": "Set minimal interval size for CSI indices to 2^INT",
            "doc": "Set minimal interval size for CSI indices to 2^INT."
          }
        ],
        "outputs": [
          {
            "id": "out_variants",
            "doc": "Tabix indexed file.",
            "label": "Tabix indexed file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n  if (inputs.dont_output_data_file==true)\n    return ''\n  return inputs.in_variants.path.split('/').slice(-1)[0]\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
            },
            "secondaryFiles": [
              ".tbi"
            ],
            "sbg:fileTypes": "GFF.GZ, BED.GZ, SAM.GZ, VCF.GZ, PSLTBL.GZ"
          },
          {
            "id": "out_index",
            "doc": "Tabix index file.",
            "label": "Tabix index",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    if (inputs.csi!=undefined && inputs.csi==true)\n    {\n        return '*.csi'\n    }\n    return '*.tbi'\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_variants))"
            },
            "sbg:fileTypes": "VCF.TBI"
          }
        ],
        "doc": "**Tabix Index** indexes a TAB-delimited genome position file\u00a0IN.TAB.BGZ\u00a0and creates an index file (IN.TAB.BGZ.TBI\u00a0or\u00a0IN.TAB.BGZ.CSI). The input data file must be position sorted and compressed by\u00a0bgzip\u00a0which has a\u00a0gzip-like interface[1].\n\nA list of all inputs and parameters with corresponding descriptions can be found at the bottom of this page.\n\n###Common Use Cases\n**Tabix Index** can be used to generate index or to pass-through data and data\u2019s index if the index file is provided as input.\n\nDepending on the usage of the tool, whether it is used as a stand-alone tool or as a part of a workflow, **Don't output data file** parameter is used as follows:\n 1. If the index is expected to be generated inside a workflow, **Don't output data file** should be set to *False* or left unset. \n 2. If the index is generated when the tool is run separately, **Don't output data file** should be set *True*. This way only index will be outputted and there will be no duplication of the files on the output. \n\n###Changes Introduced by Seven Bridges\nPlease note that in this tool, **Tabix** is wrapped only for indexing. Querying and other options are not supported.\n\n###Common Issues and Important Notes\nAs described by tool\u2019s official description, **Select input file format** ('-p', '-- preset STR') should not be used with any of the following options: **Column number for sequence names** ('-s', '-- sequence INT'), **Column number for region start** ('-b', '-- begin INT'), **Column number for region end** ('-e', '-- end INT'), **Skip comment lines starting with character CHAR** ('-c', '-- comment CHAR') and **Specify if the position in the data file is 0 based** ('-0', '-- zero-based'). As a consequence, wrapper will ignore '-s', '-b', '-e', '-c', '-0' options (will not be present in the command-line) if they are used together with '-p'.\n\n###Performance Benchmarking\n**Tabix Index** is not CPU/Memory intensive. The default c4.2 AWS instance can be used.\nCost can be significantly reduced by using spot instances. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.\n\n### References\n\n[1] [Tabix](http://www.htslib.org/doc/tabix.html)",
        "label": "Tabix Index CWL1.0",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n  if (inputs.index_file==undefined)\n    return \"/opt/htslib-1.9/tabix \"\n  else\n    return \"echo \\\"Passing inputs to outputs.\\\" \"\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n  if (inputs.index_file==undefined)\n    return \" -f \"\n  return ''\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n  if (inputs.mem_per_job)\n  {\n    return inputs.mem_per_job\n  }\n  else\n  {\n    return 2048\n  }\n}",
            "coresMin": "${\n  if (inputs.cpu_per_job)\n  {\n    return inputs.cpu_per_job\n  }\n  else\n  {\n    return 1\n  }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/nevenam/htslib-1-9:0"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              "$(inputs.in_variants)",
              "$(inputs.index_file)"
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:toolkitVersion": "1.9.0",
        "sbg:image_url": null,
        "sbg:links": [
          {
            "id": "http://www.htslib.org/",
            "label": "Homepage"
          },
          {
            "id": "https://github.com/samtools/htslib/tree/master",
            "label": "Sourcecode"
          },
          {
            "id": "http://www.htslib.org/doc/#manual-pages",
            "label": "Documentation"
          },
          {
            "id": "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/",
            "label": "Publication"
          },
          {
            "id": "http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.5.tar.bz2/download",
            "label": "Download"
          }
        ],
        "sbg:toolAuthor": "Heng Li -  Broad Institue",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639610254,
            "sbg:revisionNotes": "Copy of admin/sbg-public-data/tabix-index-1-9-cwl1-0/17"
          }
        ],
        "sbg:license": "The MIT/Expat License",
        "sbg:projectName": "HGI",
        "sbg:toolkit": "Tabix",
        "sbg:categories": [
          "Indexing",
          "CWL1.0",
          "Utilities",
          "VCF Processing"
        ],
        "sbg:appVersion": [
          "v1.0"
        ],
        "sbg:id": "markoz/hgi/tabix-index-1-9-cwl1-0/0",
        "sbg:revision": 0,
        "sbg:revisionNotes": "Copy of admin/sbg-public-data/tabix-index-1-9-cwl1-0/17",
        "sbg:modifiedOn": 1639610254,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1639610254,
        "sbg:createdBy": "markoz",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "markoz"
        ],
        "sbg:latestRevision": 0,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a9ce783911f24cf6f9397ee5da442eb54d1f715ca9a0201fd26b8ff1e79c1bff4",
        "sbg:copyOf": "admin/sbg-public-data/tabix-index-1-9-cwl1-0/17"
      },
      "label": "Tabix Index CWL1.0",
      "sbg:x": -265.5848693847656,
      "sbg:y": -479.6109924316406
    },
    {
      "id": "bcftools_view",
      "in": [
        {
          "id": "in_variants",
          "source": "bcftools_reheader_1_10_1/out_variants"
        },
        {
          "id": "include_expression",
          "source": "include_expression"
        },
        {
          "id": "samples_file",
          "source": "samples_file"
        },
        {
          "id": "keep_ids_file",
          "source": "keep_ids_file"
        }
      ],
      "out": [
        {
          "id": "out_variants"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/bcftools-view/1",
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
              {
                "pattern": "${ \n    if(self.basename.split('.').pop() == 'gz'){\n    if(self.nameroot.split('.').pop() == 'bcf'){\n        return self.nameroot + \".gz.csi\"}\n    else{\n        return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n    if(self.basename.split('.').pop() == 'bcf'){\n        return self.basename + \".csi\"\n    }\n    else{\n    return self.basename + \".tbi\"}\n}\n\n}",
                "required": true
              }
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
            "doc": "Restrict to comma-separated list of regions (e.g. chr|chr:pos|chr:from-to|chr:from-[,\u2026])."
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
          },
          {
            "id": "keep_ids_file",
            "type": "File?",
            "label": "Keep IDs file",
            "doc": "File containing IDs of variants to keep"
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
              "$(inputs.in_variants)",
              {
                "entry": "$(inputs.keep_ids_file)",
                "writable": false
              }
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
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639786213,
            "sbg:revisionNotes": null
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639786732,
            "sbg:revisionNotes": "import + keep_ids_file added"
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/bcftools-view/1",
        "sbg:revision": 1,
        "sbg:revisionNotes": "import + keep_ids_file added",
        "sbg:modifiedOn": 1639786732,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1639786213,
        "sbg:createdBy": "markoz",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "markoz"
        ],
        "sbg:latestRevision": 1,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a20dc63ef7ead22db76ce3db3fa2f73f65b82d0e2e4af4382396ff2af4cb6905d"
      },
      "label": "Bcftools View",
      "sbg:x": 140.64418029785156,
      "sbg:y": -467.8241882324219
    },
    {
      "id": "vcftools_hwe_0_1_14_cwl1",
      "in": [
        {
          "id": "not_chr",
          "source": [
            "not_chr"
          ]
        },
        {
          "id": "input_file",
          "source": "vcftools_max_missing_cwl1/output_file"
        },
        {
          "id": "chr",
          "source": [
            "chr"
          ]
        },
        {
          "id": "hwe",
          "source": "hwe"
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
        "id": "markoz/hgi/vcftools-hwe-0-1-14-cwl1/1",
        "baseCommand": [
          "vcftools"
        ],
        "inputs": [
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "None",
            "id": "not_chr",
            "type": "string[]?",
            "label": "Chromosomes to ommit",
            "doc": "Chromosomes to omit."
          },
          {
            "id": "input_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 3,
              "valueFrom": "${\n    // sufix = \"_CNVs\";\n    // sufix_ext = \"txt\";\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    var file_dot_sep = filename.split(\".\");\n\n    var basename = file_dot_sep[0]\n    var file_ext = file_dot_sep[file_dot_sep.length - 1];\n\n    var prefix = \"\"\n\n    if (file_ext == \"vcf\") {\n        prefix = \"--vcf \"\n    }\n    if (file_ext == \"bcf\") {\n        prefix = \"--bcf --gatk \"\n    }\n    if (file_ext == \"gz\") {\n        prefix = \"--gzvcf \"\n    }\n\n    // new_filename = basename + \".analyzed.hwe\";\n\n    return prefix + filepath;\n}"
            },
            "label": "Input file",
            "doc": "Input file (vcf, vcf.gz, bcf)",
            "sbg:fileTypes": "VCF, VCF.GZ, BCF"
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "All",
            "id": "chr",
            "type": "string[]?",
            "label": "Chromosomes to analyze",
            "doc": "Chromosomes to analyze."
          },
          {
            "sbg:category": "Execution",
            "sbg:toolDefaultValue": "0.01",
            "id": "hwe",
            "type": "float",
            "inputBinding": {
              "prefix": "--hwe",
              "shellQuote": false,
              "position": 2
            },
            "label": "Threshold p-value for HW equilibrium",
            "doc": "p-value for elimination via the Hardy-Weinberg Equilibrium filtering."
          }
        ],
        "outputs": [
          {
            "id": "output_file",
            "doc": "Analyzed VCF file.",
            "label": "Analyzed output file",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n\n    filepath = inputs.input_file.path\n    filename = filepath.split(\"/\").pop();\n\n    if (filename.lastIndexOf(\".vcf.gz\") != -1) {\n        basename = filename.substr(0, filename.lastIndexOf(\".vcf.gz\"))\n    } else {\n        basename = filename.substr(0, filename.lastIndexOf(\".\"))\n    }\n\n    new_filename = basename + \".analyzed.recode.vcf\";\n\n    return new_filename;\n}",
              "outputEval": "${\n    return inheritMetadata(self, inputs.input_file)\n\n}"
            },
            "sbg:fileTypes": "VCF"
          }
        ],
        "doc": "VCFtools hwe assesses sites for Hardy-Weinberg Equilibrium using an exact test, as defined by Wigginton, Cutler and Abecasis (2005). Sites with a p-value below the threshold defined by this option are taken to be out of HWE and therefore excluded.",
        "label": "VCFtools Hwe",
        "arguments": [
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "--recode"
          },
          {
            "prefix": "--out",
            "shellQuote": false,
            "position": 101,
            "valueFrom": "${\n\n    var filepath = inputs.input_file.path\n    var filename = filepath.split(\"/\").pop();\n\n    if (filename.lastIndexOf(\".vcf.gz\") != -1) {\n        var basename = filename.substr(0, filename.lastIndexOf(\".vcf.gz\"))\n    } else {\n        var basename = filename.substr(0, filename.lastIndexOf(\".\"))\n    }\n\n    var new_filename = basename + \".analyzed\";\n\n    return new_filename;\n}"
          },
          {
            "shellQuote": false,
            "position": 5,
            "valueFrom": "${\n    var out = \"\"\n    for (var i = 0; i < [].concat(inputs.not_chr).length; i++ ){\n        out += \" --not-chr \" + [].concat(inputs.not_chr)[i]\n    }    \n    return out\n}"
          },
          {
            "shellQuote": false,
            "position": 4,
            "valueFrom": "${\n    var out = \"\"\n    for (var i = 0; i < [].concat(inputs.chr).length; i++ ){\n        out += \" --chr \" + [].concat(inputs.chr)[i]\n    }    \n    return out\n}"
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
            "dockerPull": "images.sbgenomics.com/ognjenm/vcftools:0.1.14"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": []
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:categories": [
          "VCF-Processing"
        ],
        "sbg:image_url": null,
        "sbg:cmdPreview": "vcftools --recode --hwe 0 --vcf sample.vcf --out sample.analyzed",
        "sbg:toolkitVersion": "0.1.14",
        "sbg:license": "GNU General Public License version 3.0 (GPLv3)",
        "sbg:links": [
          {
            "label": "Homepage",
            "id": "https://vcftools.github.io"
          },
          {
            "label": "Source code",
            "id": "https://github.com/vcftools/vcftools"
          },
          {
            "label": "Publications",
            "id": "http://bioinformatics.oxfordjournals.org/content/27/15/2156"
          }
        ],
        "sbg:toolkit": "VCFtools",
        "sbg:toolAuthor": "Adam Auton, Petr Danecek, Anthony Marcketta",
        "sbg:projectName": "HGI",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "marko_zecevic",
            "sbg:modifiedOn": 1640108686,
            "sbg:revisionNotes": "Upgraded to v1.1 from markoz/hgi/vcftools-hwe-0-1-14"
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1642710071,
            "sbg:revisionNotes": ""
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/vcftools-hwe-0-1-14-cwl1/1",
        "sbg:revision": 1,
        "sbg:revisionNotes": "",
        "sbg:modifiedOn": 1642710071,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1640108686,
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
        "sbg:content_hash": "ab566725cfd3e09a42575bf668582688e4eb661c94795e9686d4951bc99f6c6a7",
        "sbg:workflowLanguage": "CWL"
      },
      "label": "VCFtools Hwe",
      "sbg:x": 501.0728454589844,
      "sbg:y": -496.1324768066406
    },
    {
      "id": "bcftools_view_1",
      "in": [
        {
          "id": "in_variants",
          "source": "vcftools_sort_0_1_14_cwl2/output_file"
        },
        {
          "id": "include_expression",
          "source": "include_expression_1"
        },
        {
          "id": "output_type",
          "default": "CompressedVCF"
        }
      ],
      "out": [
        {
          "id": "out_variants"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "markoz/hgi/bcftools-view/1",
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
              {
                "pattern": "${ \n    if(self.basename.split('.').pop() == 'gz'){\n    if(self.nameroot.split('.').pop() == 'bcf'){\n        return self.nameroot + \".gz.csi\"}\n    else{\n        return self.nameroot + \".gz.tbi\"\n    }\n}  else{\n    if(self.basename.split('.').pop() == 'bcf'){\n        return self.basename + \".csi\"\n    }\n    else{\n    return self.basename + \".tbi\"}\n}\n\n}",
                "required": true
              }
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
            "doc": "Restrict to comma-separated list of regions (e.g. chr|chr:pos|chr:from-to|chr:from-[,\u2026])."
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
          },
          {
            "id": "keep_ids_file",
            "type": "File?",
            "label": "Keep IDs file",
            "doc": "File containing IDs of variants to keep"
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
              "$(inputs.in_variants)",
              {
                "entry": "$(inputs.keep_ids_file)",
                "writable": false
              }
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
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639786213,
            "sbg:revisionNotes": null
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "markoz",
            "sbg:modifiedOn": 1639786732,
            "sbg:revisionNotes": "import + keep_ids_file added"
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "markoz/hgi/bcftools-view/1",
        "sbg:revision": 1,
        "sbg:revisionNotes": "import + keep_ids_file added",
        "sbg:modifiedOn": 1639786732,
        "sbg:modifiedBy": "markoz",
        "sbg:createdOn": 1639786213,
        "sbg:createdBy": "markoz",
        "sbg:project": "markoz/hgi",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "markoz"
        ],
        "sbg:latestRevision": 1,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a20dc63ef7ead22db76ce3db3fa2f73f65b82d0e2e4af4382396ff2af4cb6905d",
        "sbg:workflowLanguage": "CWL"
      },
      "label": "Bcftools View",
      "sbg:x": 2355.334716796875,
      "sbg:y": -332.3305358886719
    }
  ],
  "requirements": [
    {
      "class": "ScatterFeatureRequirement"
    },
    {
      "class": "MultipleInputFeatureRequirement"
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
      "sbg:modifiedOn": 1636922916,
      "sbg:revisionNotes": null
    },
    {
      "sbg:revision": 1,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1636923012,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 2,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1636923075,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 3,
      "sbg:modifiedBy": "marko_zecevic",
      "sbg:modifiedOn": 1636982063,
      "sbg:revisionNotes": "Workflow decomposed"
    },
    {
      "sbg:revision": 4,
      "sbg:modifiedBy": "marko_zecevic",
      "sbg:modifiedOn": 1636983353,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 5,
      "sbg:modifiedBy": "marko_zecevic",
      "sbg:modifiedOn": 1636983458,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 6,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637065880,
      "sbg:revisionNotes": "added VCFtools max-missing"
    },
    {
      "sbg:revision": 7,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637106821,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 8,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637106857,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 9,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637106901,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 10,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637180284,
      "sbg:revisionNotes": "removed subsetting samples (BCFview)"
    },
    {
      "sbg:revision": 11,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637660802,
      "sbg:revisionNotes": "output separate files"
    },
    {
      "sbg:revision": 12,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637686818,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 13,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1637687617,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 14,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1638870458,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 15,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1639436180,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 16,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1639610557,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 17,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1639610615,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 18,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1639788095,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 19,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1639952451,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 20,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1640108868,
      "sbg:revisionNotes": "added HWE filter"
    },
    {
      "sbg:revision": 21,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642709786,
      "sbg:revisionNotes": "HWE added"
    },
    {
      "sbg:revision": 22,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642710173,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 23,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642723445,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 24,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642723949,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 25,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642724497,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 26,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642882268,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 27,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1642890967,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 28,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1643546864,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 29,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1645112278,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 30,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1645984326,
      "sbg:revisionNotes": "added post imputation filtering (R2)"
    },
    {
      "sbg:revision": 31,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1645984356,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 32,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1646174463,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 33,
      "sbg:modifiedBy": "markoz",
      "sbg:modifiedOn": 1646442724,
      "sbg:revisionNotes": ""
    }
  ],
  "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/markoz/hgi/infinium-array-genotyping/33.png",
  "sbg:appVersion": [
    "v1.2",
    "v1.0",
    "v1.1"
  ],
  "id": "https://cgc-api.sbgenomics.com/v2/apps/markoz/hgi/infinium-array-genotyping/33/raw/",
  "sbg:id": "markoz/hgi/infinium-array-genotyping/33",
  "sbg:revision": 33,
  "sbg:revisionNotes": "",
  "sbg:modifiedOn": 1646442724,
  "sbg:modifiedBy": "markoz",
  "sbg:createdOn": 1636922916,
  "sbg:createdBy": "markoz",
  "sbg:project": "markoz/hgi",
  "sbg:sbgMaintained": false,
  "sbg:validationErrors": [],
  "sbg:contributors": [
    "marko_zecevic",
    "markoz"
  ],
  "sbg:latestRevision": 33,
  "sbg:publisher": "sbg",
  "sbg:content_hash": "a0f8da352ff32315ed1078fee19a3a427811f710322beba906bb027d8b6e05ea0",
  "sbg:workflowLanguage": "CWL"
}