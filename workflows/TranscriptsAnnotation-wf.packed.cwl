{
    "cwlVersion": "v1.0", 
    "$schemas": [
        "https://schema.org/docs/schema_org_rdfa.html", 
        "http://edamontology.org/EDAM_1.20.owl", 
        "http://edamontology.org/EDAM_1.16.owl"
    ], 
    "$graph": [
        {
            "inputs": [
                {
                    "doc": "Force tblastn to run on a single core and ignore the --cpu argument for\nthis step only. Useful if inconsistencies when using multiple threads are\nnoticed\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--blast_single_core"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/blastSingleCore", 
                    "label": "Force tblastn to run on a single core"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cpu"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#BUSCO-v3.cwl/cpu", 
                    "label": "Specify the number of threads/cores to use (default: 1)"
                }, 
                {
                    "doc": "Allowed formats: 0.001 or 1e-03 (default: 1e-03).\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--evalue"
                    }, 
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "id": "#BUSCO-v3.cwl/evalue", 
                    "label": "E-value cutoff for BLAST searches"
                }, 
                {
                    "doc": "Must be used when output files with the provided name already exist.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--force"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/force", 
                    "label": "Force rewriting of existing files/folders"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--help"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/help", 
                    "label": "Show this help message and exit"
                }, 
                {
                    "doc": "Specify location of the BUSCO lineage data to be used.\nVisit http://busco.ezlab.org/ for available lineages.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--lineage_path"
                    }, 
                    "type": "Directory", 
                    "id": "#BUSCO-v3.cwl/lineage", 
                    "label": "Location of the BUSCO lineage data to use (e.g. fungi_odb9)"
                }, 
                {
                    "doc": "Adds substantially to the run time!\nCan improve results for some non-model organisms.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--long"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/long", 
                    "label": "Turn on Augustus optimization mode for self-training (default: Off)"
                }, 
                {
                    "doc": "Specify which BUSCO analysis mode to run.\nThere are three valid modes:\n- geno or genome, for genome assemblies (DNA).\n- tran or transcriptome, for transcriptome assemblies (DNA).\n- prot or proteins, for annotated gene sets (protein).\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--mode"
                    }, 
                    "type": "#BUSCO-assessment_modes.yaml/assessment_modes", 
                    "id": "#BUSCO-v3.cwl/mode", 
                    "label": "Sets the assessment MODE: genome, proteins, transcriptome"
                }, 
                {
                    "doc": "Give your analysis run a recognisable short name.\nOutput folders and files will be labelled (prepended) with this name.\nWARNING: do not provide a path.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--out"
                    }, 
                    "type": "string", 
                    "id": "#BUSCO-v3.cwl/outputName", 
                    "label": "Name to use for the run and all temporary files (appended)"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--quiet"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/quiet", 
                    "label": "Disable the info logs, display only errors"
                }, 
                {
                    "doc": "NB: this limit is on scaffolds, chromosomes, or transcripts, not individual hit regions.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--limit"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#BUSCO-v3.cwl/regionLimit", 
                    "label": "How many candidate regions to consider (integer, default: 3)"
                }, 
                {
                    "doc": "NB: If all the required results files from previous steps are not all found then this will not be possible.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--restart"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/restart", 
                    "label": "Restart the BUSCO run from the last successfully-completed step"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--in"
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "Input sequence file in FASTA format (not compressed/zipped!).\nCan be an assembled genome (genome mode) or transcriptome (DNA,\ntranscriptome mode), or protein sequences from an annotated gene set\n(proteins mode).\nNB: select just one transcript/protein per gene for your input,\notherwise they will appear as 'Duplicated' matches.\n", 
                    "label": "Sequence file in FASTA format", 
                    "type": "File", 
                    "id": "#BUSCO-v3.cwl/sequenceFile"
                }, 
                {
                    "doc": "See Augustus documentation for available options.\nEach lineage has a default species (see below on assessment sets).\nSelecting a closely-related species usually produces better results.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--species"
                    }, 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#BUSCO-v3.cwl/species", 
                    "label": "Name of existing Augustus species gene finding parameters"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--tarzip"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/tarzip", 
                    "label": "Results folders with many files will be tarzipped"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--tmp"
                    }, 
                    "type": [
                        "null", 
                        "Directory"
                    ], 
                    "id": "#BUSCO-v3.cwl/tempPath", 
                    "label": "Where to store temporary files (default: ./tmp)"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--version"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#BUSCO-v3.cwl/version", 
                    "label": "Show this version information and exit"
                }
            ], 
            "requirements": [
                {
                    "coresMin": 1, 
                    "class": "ResourceRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "symbols": [
                                "#BUSCO-assessment_modes.yaml/assessment_modes/geno", 
                                "#BUSCO-assessment_modes.yaml/assessment_modes/prot", 
                                "#BUSCO-assessment_modes.yaml/assessment_modes/tran"
                            ], 
                            "type": "enum", 
                            "name": "#BUSCO-assessment_modes.yaml/assessment_modes", 
                            "id": "#BUSCO-assessment_modes.yaml"
                        }
                    ]
                }
            ], 
            "doc": "BUSCO v3 provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness, based on evolutionarily-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9. BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Benchmarking Universal Single-Copy Orthologs. These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.\nPlease visit http://busco.ezlab.org/ for full documentation.\nThe BUSCO assessment software distribution is available from the public GitLab project: https://gitlab.com/ezlab/busco where it can be downloaded or cloned using a git client (git clone https://gitlab.com/ezlab/busco.git). We encourage users to opt for the git client option in order to facilitate future updates.\nBUSCO is written for Python 3.x and Python 2.7+. It runs with the standard packages. We recommend using Python3 when available.\n", 
            "baseCommand": [
                "run_BUSCO.py"
            ], 
            "label": "BUSCO: assesses genome assembly and annotation completeness with single-copy orthologs", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "outputs": [
                {
                    "doc": "tBLASTn results, not created for assessment of proteins.\nFile: tblastn_XXXX.txt = tabular tBLASTn results\nFile: coordinates_XXXX.txt = locations of BUSCO matches (genome mode)\n", 
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/blast_output"
                    }, 
                    "type": "Directory", 
                    "id": "#BUSCO-v3.cwl/blastOutput"
                }, 
                {
                    "doc": "Contains the complete results in a tabular format with scores and lengths of BUSCO matches, and coordinates (for genome mode) or gene/protein IDs (for transcriptome or proteins mode).\n", 
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/full_table_*.tsv"
                    }, 
                    "type": "File", 
                    "id": "#BUSCO-v3.cwl/fullTable", 
                    "format": "https://www.iana.org/assignments/media-types/text/tab-separated-values"
                }, 
                {
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/hmmer_output"
                    }, 
                    "type": "Directory", 
                    "id": "#BUSCO-v3.cwl/hmmerOutput", 
                    "label": "Tabular format HMMER output of searches with BUSCO HMMs"
                }, 
                {
                    "format": "https://www.iana.org/assignments/media-types/text/tab-separated-values", 
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/missing_busco_list_*.tsv"
                    }, 
                    "type": "File", 
                    "id": "#BUSCO-v3.cwl/missingBUSCOs", 
                    "label": "Contains a list of missing BUSCOs"
                }, 
                {
                    "doc": "Contains a plain text summary of the results in BUSCO notation.\nAlso gives a brief breakdown of the metrics.\n", 
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/short_summary_*.txt"
                    }, 
                    "type": "File", 
                    "id": "#BUSCO-v3.cwl/shortSummary"
                }, 
                {
                    "outputBinding": {
                        "glob": "run_$(inputs.outputName)/translated_proteins"
                    }, 
                    "type": "Directory", 
                    "id": "#BUSCO-v3.cwl/translatedProteins", 
                    "label": "Transcript sequence translations, only created during transcriptome assessment"
                }
            ], 
            "id": "#BUSCO-v3.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "$namespaces": {
                "edam": "http://edamontology.org/", 
                "s": "http://schema.org/"
            }, 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "packages": [
                        {
                            "version": [
                                "3.0.2"
                            ], 
                            "specs": [
                                "https://identifiers.org/rrid/RRID:SCR_015008"
                            ], 
                            "package": "BUSCO"
                        }
                    ], 
                    "class": "SoftwareRequirement"
                }, 
                {
                    "dockerPull": "comics/busco:3.0.2", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "blastSingleCore"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "cpu"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "float", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "evalue"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "force"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "help"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "lineage"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "long"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#value": "tran", 
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#name": "mode"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#value": "TEST", 
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#name": "outputName"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "quiet"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "regionLimit"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "restart"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "sequenceFile"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "species"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "tarzip"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "tempPath"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "version"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--block-size"
                    }, 
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/blockSize", 
                    "label": "sequence block size in billions of letters (default=2.0)"
                }, 
                {
                    "doc": "Path to the DIAMOND database file.", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--db"
                    }, 
                    "type": "File", 
                    "id": "#Diamon.blastx-v0.9.21.cwl/databaseFile", 
                    "label": "DIAMOND database input file"
                }, 
                {
                    "doc": "0   = BLAST pairwise\n5   = BLAST XML\n6   = BLAST tabular\n100 = DIAMOND alignment archive (DAA)\n101 = SAM\n\nValue 6 may be followed by a space-separated list of these keywords", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--outfmt"
                    }, 
                    "type": [
                        "null", 
                        "#Diamond-output_formats.yaml/output_formats"
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/outputFormat", 
                    "label": "Format of the output file"
                }, 
                {
                    "doc": "Ignore translated sequences that do not contain an open reading frame of at least this length.\nBy default this feature is disabled for sequences of length below 30, set to 20 for sequences of length below 100, and set to 40 otherwise. Setting this option to 1 will disable this feature.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--min-orf"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/queryGeneticCode", 
                    "label": "Genetic code used for the translation of the query sequences"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--query"
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "Path to the query input file in FASTA or FASTQ format (may be gzip compressed). If this parameter is omitted, the input will be read from stdin\n", 
                    "label": "Query input file in FASTA", 
                    "type": "File", 
                    "id": "#Diamon.blastx-v0.9.21.cwl/queryInputFile"
                }, 
                {
                    "doc": "Set strand of query to align for translated searches. By default both strands are searched. Valid values are {both, plus, minus}", 
                    "inputBinding": {
                        "position": -3, 
                        "prefix": "--strand"
                    }, 
                    "type": [
                        "null", 
                        "#Diamond-strand_values.yaml/strand"
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/strand", 
                    "label": "Set strand of query to align for translated searches"
                }, 
                {
                    "doc": "Comma-separated list of NCBI taxonomic IDs to filter the database by. Any taxonomic rank can be used, and only reference sequences matching one of the specified taxon ids will be searched against. Using this option requires setting the --taxonmap and --taxonnodes parameters for makedb.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--taxonlist"
                    }, 
                    "type": [
                        "null", 
                        {
                            "items": "int", 
                            "type": "array"
                        }
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/taxonList", 
                    "label": "Protein accession to taxon identifier NCBI mapping file"
                }, 
                {
                    "doc": "Number of CPU threads. By default, the program will auto-detect and use all available virtual cores on the machine.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--threads"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#Diamon.blastx-v0.9.21.cwl/threads", 
                    "label": "Number of CPU threads"
                }
            ], 
            "requirements": [
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "symbols": [
                                "#Diamond-strand_values.yaml/strand/both", 
                                "#Diamond-strand_values.yaml/strand/plus", 
                                "#Diamond-strand_values.yaml/strand/minus"
                            ], 
                            "type": "enum", 
                            "name": "#Diamond-strand_values.yaml/strand", 
                            "id": "#Diamond-strand_values.yaml"
                        }, 
                        {
                            "symbols": [
                                "#Diamond-output_formats.yaml/output_formats/0", 
                                "#Diamond-output_formats.yaml/output_formats/5", 
                                "#Diamond-output_formats.yaml/output_formats/6", 
                                "#Diamond-output_formats.yaml/output_formats/100", 
                                "#Diamond-output_formats.yaml/output_formats/101", 
                                "#Diamond-output_formats.yaml/output_formats/102"
                            ], 
                            "type": "enum", 
                            "name": "#Diamond-output_formats.yaml/output_formats", 
                            "id": "#Diamond-output_formats.yaml"
                        }
                    ]
                }, 
                {
                    "ramMin": 1000, 
                    "class": "ResourceRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "DIAMOND is a sequence aligner for protein and translated DNA searches,\ndesigned for high performance analysis of big sequence data.\n\nThe key features are:\n      + Pairwise alignment of proteins and translated DNA at 500x-20,000x speed of BLAST.\n      + Frameshift alignments for long read analysis.\n      + Low resource requirements and suitable for running on standard desktops or laptops.\n      + Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.\n\nPlease visit https://github.com/bbuchfink/diamond for full documentation.\n\nReleases can be downloaded from https://github.com/bbuchfink/diamond/releases\n", 
            "baseCommand": [
                "diamond", 
                "blastx"
            ], 
            "label": "Diamond blastx: Aligns DNA query sequences against a protein reference database", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "arguments": [
                {
                    "position": 0, 
                    "prefix": "--out", 
                    "valueFrom": "$(inputs.queryInputFile.basename).diamond_matches"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.queryInputFile.basename).diamond_matches"
                    }, 
                    "type": "File", 
                    "id": "#Diamon.blastx-v0.9.21.cwl/matches", 
                    "format": "http://edamontology.org/format_2333"
                }
            ], 
            "id": "#Diamon.blastx-v0.9.21.cwl", 
            "http://schema.org/author": "Maxim Scheremetjews", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "buchfink/diamond:version0.9.21", 
                    "class": "DockerRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-T"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#phmmer-v3.2.cwl/bitscoreThreshold", 
                    "label": "report sequences >= this bit score threshold in output"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cpu"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#phmmer-v3.2.cwl/cpu", 
                    "label": "Number of parallel CPU workers to use for multithreads"
                }, 
                {
                    "inputBinding": {
                        "position": 1
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "Search one or more query protein sequences against a protein sequence database.\n", 
                    "label": "Query sequence(s) file", 
                    "type": "File", 
                    "id": "#phmmer-v3.2.cwl/seqFile"
                }, 
                {
                    "label": "Target database of sequences", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#phmmer-v3.2.cwl/seqdb", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "requirements": [
                {
                    "coresMin": 2, 
                    "ramMin": 1024, 
                    "class": "ResourceRequirement", 
                    "coresMax": 4
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "The phmmer and jackhmmer programs search a single protein sequence against a protein sequence database, akin to BLASTP and PSIBLAST, respectively. (Internally, they just produce a profile HMM from the query sequence, then run HMM searches.) Please visit https://github.com/EddyRivasLab/hmmer for full documentation.\nReleases can be downloaded from https://github.com/EddyRivasLab/hmmer/releases\n", 
            "baseCommand": [
                "phmmer"
            ], 
            "label": "HMMER: Search a single protein sequence against a protein sequence database. (BLASTP-like)", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "arguments": [
                {
                    "position": 0, 
                    "prefix": "--tblout", 
                    "valueFrom": "$(inputs.seqFile.basename).phmmer_matches.tblout"
                }, 
                {
                    "position": 0, 
                    "prefix": "-o", 
                    "valueFrom": "$(inputs.seqFile.basename).phmmer_matches.out"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.seqFile.basename).phmmer_matches.tblout"
                    }, 
                    "type": "File", 
                    "id": "#phmmer-v3.2.cwl/matches"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.seqFile.basename).phmmer_matches.out"
                    }, 
                    "type": "File", 
                    "id": "#phmmer-v3.2.cwl/programOutput"
                }
            ], 
            "id": "#phmmer-v3.2.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "quay.io/biocontainers/hmmer:3.2--hfc679d8_3", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "bitscoreThreshold"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "cpu"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "seqFile"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "seqdb"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": "File", 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/covariance_model_database"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cpu"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/cpu", 
                    "label": "Number of parallel CPU workers to use for multithreads"
                }, 
                {
                    "default": false, 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cut_ga"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/cut_ga", 
                    "label": "use CM's GA gathering cutoffs as reporting thresholds"
                }, 
                {
                    "doc": "This can greatly reduce the output volume.", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--noali"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/omit_alignment_section", 
                    "label": "Omit the alignment section from the main output."
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--hmmonly"
                    }, 
                    "default": false, 
                    "doc": "Only filter stages F1 through F3 will be executed, using strict P-value\nthresholds (0.02 for F1, 0.001 for F2 and 0.00001 for F3). Additionally\na bias composition filter is used after the F1 stage (with P=0.02\nsurvival threshold). Any hit that survives all stages and has an HMM\nE-value or bit score above the reporting threshold will be output.\n", 
                    "label": "Only use the filter profile HMM for searches, do not use the CM", 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/only_hmm"
                }, 
                {
                    "streamable": true, 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "type": "File", 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/query_sequences", 
                    "format": "http://edamontology.org/format_1929"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-Z"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/search_space_size", 
                    "label": "search space size in *Mb* to <x> for E-value calculations"
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "Infernal (\"INFERence of RNA ALignment\") is for searching DNA sequence databases for RNA structure and sequence similarities. It is an implementation of a special case of profile stochastic context-free grammars called covariance models (CMs). A CM is like a sequence profile, but it scores a combination of sequence consensus and RNA secondary structure consensus, so in many cases, it is more capable of identifying RNA homologs that conserve their secondary structure more than their primary sequence.\nPlease visit http://eddylab.org/infernal/ for full documentation.\nVersion 1.1.2 can be downloaded from http://eddylab.org/infernal/infernal-1.1.2.tar.gz\n", 
            "baseCommand": [
                "cmsearch"
            ], 
            "label": "Infernal: Search sequence(s) against a covariance model database", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "arguments": [
                {
                    "position": 0, 
                    "prefix": "--tblout", 
                    "valueFrom": "$(inputs.query_sequences.basename).cmsearch_matches.tbl"
                }, 
                {
                    "position": 0, 
                    "prefix": "-o", 
                    "valueFrom": "$(inputs.query_sequences.basename).cmsearch.out"
                }
            ], 
            "outputs": [
                {
                    "doc": "http://eddylab.org/infernal/Userguide.pdf#page=60", 
                    "outputBinding": {
                        "glob": "$(inputs.query_sequences.basename).cmsearch_matches.tbl"
                    }, 
                    "type": "File", 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/matches", 
                    "label": "target hits table, format 2"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.query_sequences.basename).cmsearch.out"
                    }, 
                    "type": "File", 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/programOutput", 
                    "label": "direct output to file, not stdout"
                }
            ], 
            "id": "#infernal-cmsearch-v1.1.2.cwl", 
            "http://schema.org/author": "Michael Crusoe, Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "packages": [
                        {
                            "version": [
                                "1.1.2"
                            ], 
                            "specs": [
                                "https://identifiers.org/rrid/RRID:SCR_011809"
                            ], 
                            "package": "infernal"
                        }
                    ], 
                    "class": "SoftwareRequirement"
                }, 
                {
                    "dockerPull": "quay.io/biocontainers/infernal:1.1.2--h470a237_1", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "covariance_model_database"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "cpu"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "cut_ga"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "omit_alignment_section"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "only_hmm"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "query_sequences"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#default_value": 1000, 
                            "http://galaxyproject.org/cwl#name": "search_space_size"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 8, 
                        "prefix": "--input"
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "Optional, path to fasta file that should be loaded on Master startup. Alternatively, in CONVERT mode, the InterProScan 5 XML file to convert.", 
                    "label": "Input file path", 
                    "type": "File", 
                    "id": "#InterProScan-v5.cwl/inputFile"
                }, 
                {
                    "doc": "Optional, comma separated list of analyses. If this option is not set, ALL analyses will be run.", 
                    "inputBinding": {
                        "position": 9, 
                        "prefix": "--applications", 
                        "itemSeparator": ","
                    }, 
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-apps.yaml/apps", 
                            "type": "array"
                        }
                    ], 
                    "id": "#InterProScan-v5.cwl/applications", 
                    "label": "Analysis"
                }, 
                {
                    "inputBinding": {
                        "position": 10, 
                        "prefix": "--formats", 
                        "itemSeparator": ","
                    }, 
                    "default": "TSV", 
                    "doc": "Optional, case-insensitive, comma separated list of output formats. Supported formats are TSV, XML, JSON, GFF3, HTML and SVG. Default for protein sequences are TSV, XML and GFF3, or for nucleotide sequences GFF3 and XML.", 
                    "label": "output format", 
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-protein_formats.yaml/protein_formats", 
                            "type": "array"
                        }
                    ], 
                    "id": "#InterProScan-v5.cwl/outputFormat"
                }, 
                {
                    "type": "Directory", 
                    "id": "#InterProScan-v5.cwl/databases"
                }, 
                {
                    "doc": "Optional, excludes sites from the XML, JSON output.", 
                    "inputBinding": {
                        "position": 11, 
                        "prefix": "--disable-residue-annot"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#InterProScan-v5.cwl/disableResidueAnnotation", 
                    "label": "Disables residue annotation"
                }, 
                {
                    "doc": "Optional, the type of the input sequences (dna/rna (n) or protein (p)). The default sequence type is protein.", 
                    "inputBinding": {
                        "position": 12, 
                        "prefix": "--seqtype"
                    }, 
                    "type": [
                        "null", 
                        {
                            "symbols": [
                                "#InterProScan-v5.cwl/seqtype/seqtype/p", 
                                "#InterProScan-v5.cwl/seqtype/seqtype/n"
                            ], 
                            "type": "enum", 
                            "name": "#InterProScan-v5.cwl/seqtype/seqtype"
                        }
                    ], 
                    "id": "#InterProScan-v5.cwl/seqtype", 
                    "label": "Sequence type"
                }
            ], 
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "symbols": [
                                "#InterProScan-apps.yaml/apps/TIGRFAM", 
                                "#InterProScan-apps.yaml/apps/SFLD", 
                                "#InterProScan-apps.yaml/apps/SUPERFAMILY", 
                                "#InterProScan-apps.yaml/apps/Gene3D", 
                                "#InterProScan-apps.yaml/apps/Hamap", 
                                "#InterProScan-apps.yaml/apps/Coils", 
                                "#InterProScan-apps.yaml/apps/ProSiteProfiles", 
                                "#InterProScan-apps.yaml/apps/SMART", 
                                "#InterProScan-apps.yaml/apps/CDD", 
                                "#InterProScan-apps.yaml/apps/PRINTS", 
                                "#InterProScan-apps.yaml/apps/PIRSF", 
                                "#InterProScan-apps.yaml/apps/ProSitePatterns", 
                                "#InterProScan-apps.yaml/apps/PfamA", 
                                "#InterProScan-apps.yaml/apps/ProDom", 
                                "#InterProScan-apps.yaml/apps/MobiDBLite", 
                                "#InterProScan-apps.yaml/apps/SignalP_GRAM_POSITIVE", 
                                "#InterProScan-apps.yaml/apps/SignalP_GRAM_NEGATIVE", 
                                "#InterProScan-apps.yaml/apps/SignalP_EUK", 
                                "#InterProScan-apps.yaml/apps/Phobius", 
                                "#InterProScan-apps.yaml/apps/TMHMM"
                            ], 
                            "type": "enum", 
                            "name": "#InterProScan-apps.yaml/apps", 
                            "id": "#InterProScan-apps.yaml"
                        }, 
                        {
                            "symbols": [
                                "#InterProScan-protein_formats.yaml/protein_formats/TSV", 
                                "#InterProScan-protein_formats.yaml/protein_formats/XML", 
                                "#InterProScan-protein_formats.yaml/protein_formats/GFF3", 
                                "#InterProScan-protein_formats.yaml/protein_formats/JSON", 
                                "#InterProScan-protein_formats.yaml/protein_formats/HTML", 
                                "#InterProScan-protein_formats.yaml/protein_formats/SVG"
                            ], 
                            "type": "enum", 
                            "name": "#InterProScan-protein_formats.yaml/protein_formats", 
                            "id": "#InterProScan-protein_formats.yaml"
                        }
                    ]
                }, 
                {
                    "coresMin": 3, 
                    "ramMin": 8192, 
                    "class": "ResourceRequirement"
                }, 
                {
                    "dockerPull": "biocontainers/interproscan:v5.30-69.0_cv1", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "InterProScan is the software package that allows sequences (protein and nucleic) to be scanned against InterPro's signatures. Signatures are predictive models, provided by several different databases, that make up the InterPro consortium.\n\nThis tool description is using a Docker container tagged as version v5.30-69.0.\n\nDocumentation on how to run InterProScan 5 can be found here: https://github.com/ebi-pf-team/interproscan/wiki/HowToRun", 
            "baseCommand": [], 
            "id": "#InterProScan-v5.cwl", 
            "arguments": [
                {
                    "shellQuote": false, 
                    "position": 0, 
                    "valueFrom": "cp -r /opt/interproscan $(runtime.outdir)/interproscan;"
                }, 
                {
                    "shellQuote": false, 
                    "position": 1, 
                    "valueFrom": "cp -rs $(inputs.databases.path)/data/* $(runtime.outdir)/interproscan/data;"
                }, 
                {
                    "shellQuote": false, 
                    "position": 2, 
                    "valueFrom": "$(runtime.outdir)/interproscan/interproscan.sh"
                }, 
                {
                    "position": 3, 
                    "valueFrom": "--disable-precalc"
                }, 
                {
                    "position": 4, 
                    "valueFrom": "--goterms"
                }, 
                {
                    "position": 5, 
                    "valueFrom": "--pathways"
                }, 
                {
                    "position": 6, 
                    "prefix": "--tempdir", 
                    "valueFrom": "$(runtime.tmpdir)"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.inputFile.nameroot).fasta.tsv"
                    }, 
                    "type": "File", 
                    "id": "#InterProScan-v5.cwl/i5Annotations"
                }
            ], 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/author": "Michael Crusoe, Aleksandra Ola Tarkowska, Maxim Scheremetjew", 
            "label": "InterProScan: protein sequence classifier", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "applications"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "proteinFile"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--gene_trans_map"
                    }, 
                    "format": "http://edamontology.org/format_3475", 
                    "doc": "gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return>)", 
                    "label": "gene-to-transcript mapping", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/geneToTranscriptMap"
                }, 
                {
                    "doc": "genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-G"
                    }, 
                    "type": [
                        "null", 
                        "#TransDecoder-v5-genetic_codes.yaml/genetic_codes"
                    ], 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/geneticCode", 
                    "label": "genetic code"
                }, 
                {
                    "doc": "minimum protein length (default: 100)", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-m"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/minimumProteinLength", 
                    "label": "minimum protein length"
                }, 
                {
                    "doc": "strand-specific (only analyzes top strand)", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-S"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/strandSpecific", 
                    "label": "strand-specific"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-t"
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "FASTA formatted sequence file containing your transcripts.", 
                    "label": "transcripts.fasta", 
                    "type": "File", 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/transcriptsFile"
                }
            ], 
            "requirements": [
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "symbols": [
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/universal", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Euplotes", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Tetrahymena", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Candida", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Acetabularia", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Canonical", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Vertebrates", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Arthropods", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Echinoderms", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Molluscs", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Ascidians", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Nematodes", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Platyhelminths", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Yeasts", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Euascomycetes", 
                                "#TransDecoder-v5-genetic_codes.yaml/genetic_codes/Mitochondrial-Protozoans"
                            ], 
                            "type": "enum", 
                            "name": "#TransDecoder-v5-genetic_codes.yaml/genetic_codes", 
                            "id": "#TransDecoder-v5-genetic_codes.yaml"
                        }
                    ]
                }, 
                {
                    "ramMin": 1024, 
                    "class": "ResourceRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments to the genome using Tophat and Cufflinks.\nTransDecoder identifies likely coding sequences based on the following criteria:\n      + a minimum length open reading frame (ORF) is found in a transcript sequence\n      + a log-likelihood score similar to what is computed by the GeneID software is > 0.\n      + the above coding score is greatest when the ORF is scored in the 1st reading frame\n      as compared to scores in the other 2 forward reading frames.\n      + if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF,\n      the longer one is reported. However, a single transcript can report multiple ORFs\n      (allowing for operons, chimeras, etc).\n      + a PSSM is built/trained/used to refine the start codon prediction.\n      + optional the putative peptide has a match to a Pfam domain above the noise cutoff score.\n\nPlease visit https://github.com/TransDecoder/TransDecoder/wiki for full documentation.\n\nReleases can be downloaded from https://github.com/TransDecoder/TransDecoder/releases", 
            "baseCommand": [
                "TransDecoder.LongOrfs"
            ], 
            "label": "TransDecoder.LongOrfs: Perl script, which extracts the long open reading frames", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.transcriptsFile.basename).transdecoder_dir"
                    }, 
                    "type": "Directory", 
                    "id": "#TransDecoder.LongOrfs-v5.cwl/workingDir"
                }
            ], 
            "id": "#TransDecoder.LongOrfs-v5.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "greatfireball/ime_transdecoder:5.0.2", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "geneToTranscriptMap"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "geneticCode"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "minimumProteinLength"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "strandSpecific"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "transcriptsFile"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-G"
                    }, 
                    "type": [
                        "null", 
                        "#TransDecoder-v5-genetic_codes.yaml/genetic_codes"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/geneticCode", 
                    "label": "genetic code"
                }, 
                {
                    "type": "Directory", 
                    "id": "#TransDecoder.Predict-v5.cwl/longOpenReadingFrames"
                }, 
                {
                    "doc": "Start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default.", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--no_refine_starts"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/noRefineStarts", 
                    "label": "No refine starts"
                }, 
                {
                    "doc": "blastp output in '-outfmt 6' format.\nAny ORF with a blast match will be retained in the final output.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--retain_blastp_hits"
                    }, 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/retainBlastpHits", 
                    "label": "Retain Blastp hits"
                }, 
                {
                    "doc": "Under 'strict' mode, retain all ORFs found that are equal or longer than these many nucleotides\neven if no other evidence marks it as coding (default: 1000000) so essentially turned off by default.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--retain_long_orfs_length"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/retainLongOrfsLength", 
                    "label": "Retain long ORFs length"
                }, 
                {
                    "doc": "'dynamic' (default) or 'strict'. In dynamic mode, sets range according to 1%FDR in random sequence of same GC content.", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--retain_long_orfs_mode"
                    }, 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/retainLongOrfsMode", 
                    "label": "Retain long ORFs mode"
                }, 
                {
                    "doc": "Domain table output file from running hmmscan to search Pfam (see transdecoder.github.io for info).\nAny ORF with a pfam domain hit will be retained in the final output.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--retain_pfam_hits"
                    }, 
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/retainPfamHits", 
                    "label": "Retain Pfam hits"
                }, 
                {
                    "doc": "Retain only the single best ORF per transcript (prioritized by homology then ORF length)", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--single_best_only"
                    }, 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/singleBestOnly", 
                    "label": "Single best only"
                }, 
                {
                    "doc": "If no --train, top longest ORFs to train Markov Model (hexamer stats) (default: 500)\nNote, 10x this value are first selected for removing redundancies,\nand then this -T value of longest ORFs are selected from the non-redundant set.\n", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-T"
                    }, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#TransDecoder.Predict-v5.cwl/train", 
                    "label": "minimum protein length"
                }, 
                {
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-t"
                    }, 
                    "format": "http://edamontology.org/format_1929", 
                    "doc": "FASTA formatted sequence file containing your transcripts.", 
                    "label": "transcripts.fasta", 
                    "type": "File", 
                    "id": "#TransDecoder.Predict-v5.cwl/transcriptsFile"
                }
            ], 
            "requirements": [
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "$import": "#TransDecoder-v5-genetic_codes.yaml"
                        }
                    ]
                }, 
                {
                    "ramMin": 1024, 
                    "class": "ResourceRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.transcriptsFile)", 
                        {
                            "writable": true, 
                            "entry": "$(inputs.longOpenReadingFrames)", 
                            "entryname": "$(inputs.transcriptsFile.basename).transdecoder_dir"
                        }
                    ]
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments to the genome using Tophat and Cufflinks.\nTransDecoder identifies likely coding sequences based on the following criteria:\n      + a minimum length open reading frame (ORF) is found in a transcript sequence\n      + a log-likelihood score similar to what is computed by the GeneID software is > 0.\n      + the above coding score is greatest when the ORF is scored in the 1st reading frame\n      as compared to scores in the other 2 forward reading frames.\n      + if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF,\n      the longer one is reported. However, a single transcript can report multiple ORFs \n      (allowing for operons, chimeras, etc).\n      + a PSSM is built/trained/used to refine the start codon prediction.\n      + optional the putative peptide has a match to a Pfam domain above the noise cutoff score.\n      \nPlease visit https://github.com/TransDecoder/TransDecoder/wiki for full documentation.\nReleases can be downloaded from https://github.com/TransDecoder/TransDecoder/releases\n", 
            "baseCommand": [
                "TransDecoder.Predict"
            ], 
            "label": "TransDecoder.Predict: Perl script, which predicts the likely coding regions", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.transcriptsFile.basename).transdecoder.bed"
                    }, 
                    "type": "File", 
                    "id": "#TransDecoder.Predict-v5.cwl/bed_output", 
                    "format": "http://edamontology.org/format_3003"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.transcriptsFile.basename).transdecoder.cds"
                    }, 
                    "type": "File", 
                    "id": "#TransDecoder.Predict-v5.cwl/coding_regions", 
                    "format": "http://edamontology.org/format_1929"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.transcriptsFile.basename).transdecoder.gff3"
                    }, 
                    "type": "File", 
                    "id": "#TransDecoder.Predict-v5.cwl/gff3_output", 
                    "format": "http://edamontology.org/format_1975"
                }, 
                {
                    "outputBinding": {
                        "glob": "$(inputs.transcriptsFile.basename).transdecoder.pep"
                    }, 
                    "type": "File", 
                    "id": "#TransDecoder.Predict-v5.cwl/peptide_sequences", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "id": "#TransDecoder.Predict-v5.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "greatfireball/ime_transdecoder:5.0.2", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "geneticCode"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "longOpenReadingFrames"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "noRefineStarts"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "retainBlastpHits"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "retainLongOrfsLength"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "retainLongOrfsMode"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "text", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "retainPfamHits"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "singleBestOnly"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "train"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "transcriptsFile"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "Not all models provided need to be a member of a clan", 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--clanin"
                    }, 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "id": "#cmsearch-deoverlap-v0.02.cwl/clan_information", 
                    "label": "clan information on the models provided"
                }, 
                {
                    "inputBinding": {
                        "position": 1, 
                        "valueFrom": "$(self.basename)"
                    }, 
                    "type": "File", 
                    "id": "#cmsearch-deoverlap-v0.02.cwl/cmsearch_matches"
                }
            ], 
            "requirements": [
                {
                    "envDef": [
                        {
                            "envName": "LC_ALL", 
                            "envValue": "C"
                        }
                    ], 
                    "class": "EnvVarRequirement"
                }, 
                {
                    "ramMin": 100, 
                    "class": "ResourceRequirement", 
                    "coresMax": 1
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.cmsearch_matches)"
                    ]
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "https://github.com/nawrockie/cmsearch_tblout_deoverlap/blob/master/00README.txt", 
            "baseCommand": [
                "cmsearch-deoverlap.pl"
            ], 
            "label": "Cmsearch-deoverlap: Remove lower scoring overlaps from cmsearch --tblout files.", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "outputs": [
                {
                    "doc": "http://eddylab.org/infernal/Userguide.pdf#page=60", 
                    "outputBinding": {
                        "glob": "*.deoverlapped"
                    }, 
                    "type": "File", 
                    "id": "#cmsearch-deoverlap-v0.02.cwl/deoverlapped_matches", 
                    "label": "target hits table, format 2"
                }
            ], 
            "id": "#cmsearch-deoverlap-v0.02.cwl", 
            "http://schema.org/author": "Michael Crusoe, Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "packages": [
                        {
                            "version": [
                                "0.02"
                            ], 
                            "specs": [
                                "https://github.com/nawrockie/cmsearch_tblout_deoverlap"
                            ], 
                            "package": "cmsearch_tblout_deoverlap"
                        }
                    ], 
                    "class": "SoftwareRequirement"
                }, 
                {
                    "dockerPull": "biocrusoe/cmsearch-deoverlap", 
                    "class": "DockerRequirement"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#optional": true, 
                            "http://galaxyproject.org/cwl#name": "clan_information"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#name": "cmsearch_matches"
                        }
                    ]
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "streamable": true, 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#concatenate.cwl/files"
                }, 
                {
                    "type": "string", 
                    "id": "#concatenate.cwl/outputFileName"
                }
            ], 
            "requirements": [
                {
                    "ramMin": 100, 
                    "class": "ResourceRequirement", 
                    "coresMax": 1
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "stdout": "$(inputs.outputFileName)", 
            "doc": "The cat (short for \u201cconcatenate\u201c) command is one of the most frequently used command in Linux/Unix like operating systems. cat command allows us to create single or multiple files, view contain of file, concatenate files and redirect output in terminal or files.\n", 
            "baseCommand": [
                "cat"
            ], 
            "label": "Redirecting Multiple Files Contain in a Single File", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.outputFileName)", 
                        "outputEval": "${ self[0].format = inputs.files[0].format;\n   return self; }\n"
                    }, 
                    "type": "File", 
                    "id": "#concatenate.cwl/result"
                }
            ], 
            "id": "#concatenate.cwl", 
            "http://schema.org/author": "Michael Crusoe, Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "alpine:3.7", 
                    "class": "DockerRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--replace", 
                        "valueFrom": "$(self.find):$(self.replace)"
                    }, 
                    "type": [
                        "null", 
                        "#esl-reformat-replace.yaml/replace"
                    ], 
                    "id": "#esl-reformat.cwl/replace"
                }, 
                {
                    "inputBinding": {
                        "position": 3
                    }, 
                    "type": "File", 
                    "id": "#esl-reformat.cwl/sequences", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "requirements": [
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "id": "#esl-reformat-replace.yaml", 
                            "fields": [
                                {
                                    "doc": "Must be equal length with \"replace\". Each character from the input file\nwill be replaced by its counterpart (at the same position) from \"replace\"\n", 
                                    "type": "string", 
                                    "name": "#esl-reformat-replace.yaml/replace/find", 
                                    "label": "find"
                                }, 
                                {
                                    "doc": "must be equal length with \"find\"", 
                                    "type": "string", 
                                    "name": "#esl-reformat-replace.yaml/replace/replace", 
                                    "label": "replace"
                                }
                            ], 
                            "type": "record", 
                            "name": "#esl-reformat-replace.yaml/replace", 
                            "label": "sequence token replacement"
                        }
                    ]
                }, 
                {
                    "ramMin": 100, 
                    "class": "ResourceRequirement", 
                    "coresMax": 1
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "stdout": "$(inputs.sequences.basename).reformatted_seqs", 
            "doc": "Normalizes input sequences to FASTA with fixed number of sequence characters\nper line using esl-reformat from https://github.com/EddyRivasLab/easel\n", 
            "baseCommand": [
                "esl-reformat"
            ], 
            "label": "Normalizes input sequences to FASTA using esl-reformat", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "arguments": [
                {
                    "position": 2, 
                    "valueFrom": "fasta"
                }
            ], 
            "outputs": [
                {
                    "type": "stdout", 
                    "id": "#esl-reformat.cwl/reformatted_sequences", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "id": "#esl-reformat.cwl", 
            "http://schema.org/author": "Michael Crusoe, Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2017", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "quay.io/biocontainers/hmmer:3.2--hfc679d8_3", 
                    "class": "DockerRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "default": 10, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#fasta_chunker.cwl/chunk_size"
                }, 
                {
                    "type": "File", 
                    "id": "#fasta_chunker.cwl/seqs", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "requirements": [
                {
                    "ramMin": 100, 
                    "class": "ResourceRequirement", 
                    "coresMax": 1
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "based upon code by developers from EMBL-EBI", 
            "baseCommand": [
                "python3"
            ], 
            "label": "split FASTA by number of records", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "arguments": [
                {
                    "position": 0, 
                    "prefix": "-c", 
                    "valueFrom": "from Bio import SeqIO\nimport os\ncurrentSequences = []\nos.mkdir(\"$(runtime.outdir)/chunks\")\nfor record in SeqIO.parse(\"$(inputs.seqs.path)\", \"fasta\"):\n    currentSequences.append(record)\n    if len(currentSequences) == $(inputs.chunk_size):\n        fileName = currentSequences[0].id + \"_\" + currentSequences[-1].id + \".fasta\"\n        for char in [ \"/\", \" \", \":\" ]:\n            fileName = fileName.replace(char, \"_\")\n        SeqIO.write(currentSequences, \"$(runtime.outdir)/chunks/\"+fileName, \"fasta\")\n        currentSequences = []\n\n# write any remaining sequences\nif len(currentSequences) > 0:\n    fileName = currentSequences[0].id + \"_\" + currentSequences[-1].id + \".fasta\"\n    for char in [ \"/\", \" \", \":\" ]:\n        fileName = fileName.replace(char, \"_\")\n    SeqIO.write(currentSequences, \"$(runtime.outdir)/chunks/\"+fileName, \"fasta\")\n"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "chunks/*_*.fasta"
                    }, 
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#fasta_chunker.cwl/chunks", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "id": "#fasta_chunker.cwl", 
            "http://schema.org/author": "Michael Crusoe, Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "CommandLineTool", 
            "hints": [
                {
                    "dockerPull": "biopython/biopython:latest", 
                    "class": "DockerRequirement"
                }, 
                {
                    "packages": [
                        {
                            "version": [
                                "1.65", 
                                "1.66", 
                                "1.69", 
                                "1.72"
                            ], 
                            "specs": [
                                "https://identifiers.org/rrid/RRID:SCR_007173"
                            ], 
                            "package": "biopython"
                        }
                    ], 
                    "class": "SoftwareRequirement"
                }
            ]
        }, 
        {
            "inputs": [
                {
                    "doc": "Optional, path to fasta file that should be loaded on Master startup. Alternatively, in CONVERT mode, the InterProScan 5 XML file to convert.", 
                    "label": "Input file path", 
                    "type": "File", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/inputFile", 
                    "format": "http://edamontology.org/format_1929"
                }, 
                {
                    "doc": "Optional, comma separated list of analyses. If this option is not set, ALL analyses will be run.", 
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-apps.yaml/apps", 
                            "type": "array"
                        }
                    ], 
                    "id": "#InterProScan-v5-chunked-wf.cwl/applications", 
                    "label": "Analysis"
                }, 
                {
                    "doc": "Optional, case-insensitive, comma separated list of output formats. Supported formats are TSV, XML, JSON, GFF3, HTML and SVG. Default for protein sequences are TSV, XML and GFF3, or for nucleotide sequences GFF3 and XML.", 
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-protein_formats.yaml/protein_formats", 
                            "type": "array"
                        }
                    ], 
                    "id": "#InterProScan-v5-chunked-wf.cwl/outputFormat", 
                    "label": "output format"
                }, 
                {
                    "type": "Directory", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/databases"
                }, 
                {
                    "default": 10000, 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#InterProScan-v5-chunked-wf.cwl/chunk_size"
                }, 
                {
                    "doc": "Optional, excludes sites from the XML, JSON output.", 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#InterProScan-v5-chunked-wf.cwl/disableResidueAnnotation", 
                    "label": "Disables residue annotation"
                }, 
                {
                    "doc": "Optional, the type of the input sequences (dna/rna (n) or protein (p)). The default sequence type is protein.", 
                    "type": [
                        "null", 
                        {
                            "symbols": [
                                "#InterProScan-v5-chunked-wf.cwl/seqtype/seqtype/p", 
                                "#InterProScan-v5-chunked-wf.cwl/seqtype/seqtype/n"
                            ], 
                            "type": "enum", 
                            "name": "#InterProScan-v5-chunked-wf.cwl/seqtype/seqtype"
                        }
                    ], 
                    "id": "#InterProScan-v5-chunked-wf.cwl/seqtype", 
                    "label": "Sequence type"
                }, 
                {
                    "default": "full_i5_annotations", 
                    "type": "string", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/catOutputFileName"
                }
            ], 
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }, 
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "$import": "#InterProScan-apps.yaml"
                        }, 
                        {
                            "$import": "#InterProScan-protein_formats.yaml"
                        }
                    ]
                }
            ], 
            "outputs": [
                {
                    "outputSource": "#InterProScan-v5-chunked-wf.cwl/combine_interproscan_results/result", 
                    "type": "File", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/i5Annotations"
                }
            ], 
            "label": "Runs InterProScan on batches of sequences to retrieve functional annotations.", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "steps": [
                {
                    "out": [
                        "#InterProScan-v5-chunked-wf.cwl/combine_interproscan_results/result"
                    ], 
                    "run": "#concatenate.cwl", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/combine_interproscan_results", 
                    "in": [
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/i5Annotations", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/combine_interproscan_results/files"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/catOutputFileName", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/combine_interproscan_results/outputFileName"
                        }
                    ]
                }, 
                {
                    "run": "#InterProScan-v5.cwl", 
                    "scatter": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/inputFile", 
                    "in": [
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/applications", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/applications"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/databases", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/databases"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/disableResidueAnnotation", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/disableResidueAnnotation"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/split_seqs/chunks", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/inputFile"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/outputFormat", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan/outputFormat"
                        }
                    ], 
                    "label": "Run InterProScan on chunked sequence files", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/run_interproscan", 
                    "out": [
                        "#InterProScan-v5-chunked-wf.cwl/run_interproscan/i5Annotations"
                    ]
                }, 
                {
                    "out": [
                        "#InterProScan-v5-chunked-wf.cwl/split_seqs/chunks"
                    ], 
                    "run": "#fasta_chunker.cwl", 
                    "id": "#InterProScan-v5-chunked-wf.cwl/split_seqs", 
                    "in": [
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/chunk_size", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/split_seqs/chunk_size"
                        }, 
                        {
                            "source": "#InterProScan-v5-chunked-wf.cwl/inputFile", 
                            "id": "#InterProScan-v5-chunked-wf.cwl/split_seqs/seqs"
                        }
                    ]
                }
            ], 
            "id": "#InterProScan-v5-chunked-wf.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "Workflow"
        }, 
        {
            "inputs": [
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/singleBestOnly"
                }, 
                {
                    "type": "File", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/transcriptsFile", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputSource": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/bed_output", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/bed_output"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/coding_regions", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/coding_regions"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/gff3_output", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/gff3_output"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/peptide_sequences", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/peptide_sequences"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "out": [
                        "#TransDecoder-v5-wf-2steps.cwl/extract_long_orfs/workingDir"
                    ], 
                    "in": [
                        {
                            "source": "#TransDecoder-v5-wf-2steps.cwl/transcriptsFile", 
                            "id": "#TransDecoder-v5-wf-2steps.cwl/extract_long_orfs/transcriptsFile"
                        }
                    ], 
                    "run": "#TransDecoder.LongOrfs-v5.cwl", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/extract_long_orfs", 
                    "label": "Extracts the long open reading frames"
                }, 
                {
                    "out": [
                        "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/bed_output", 
                        "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/coding_regions", 
                        "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/gff3_output", 
                        "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/peptide_sequences"
                    ], 
                    "in": [
                        {
                            "source": "#TransDecoder-v5-wf-2steps.cwl/extract_long_orfs/workingDir", 
                            "id": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/longOpenReadingFrames"
                        }, 
                        {
                            "source": "#TransDecoder-v5-wf-2steps.cwl/singleBestOnly", 
                            "id": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/singleBestOnly"
                        }, 
                        {
                            "source": "#TransDecoder-v5-wf-2steps.cwl/transcriptsFile", 
                            "id": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions/transcriptsFile"
                        }
                    ], 
                    "run": "#TransDecoder.Predict-v5.cwl", 
                    "id": "#TransDecoder-v5-wf-2steps.cwl/predict_coding_regions", 
                    "label": "Predicts the likely coding regions"
                }
            ], 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "label": "TransDecoder 2 step workflow, running TransDecoder.LongOrfs (step 1) followed by TransDecoder.Predict (step2)", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "id": "#TransDecoder-v5-wf-2steps.cwl"
        }, 
        {
            "inputs": [
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "id": "#main/blockSize"
                }, 
                {
                    "type": "Directory", 
                    "id": "#main/buscoLineage"
                }, 
                {
                    "type": "#BUSCO-assessment_modes.yaml/assessment_modes", 
                    "id": "#main/buscoMode"
                }, 
                {
                    "type": "string", 
                    "id": "#main/buscoOutputName"
                }, 
                {
                    "type": "File", 
                    "id": "#main/clanInfoFile"
                }, 
                {
                    "type": "int", 
                    "id": "#main/cmsearchCores"
                }, 
                {
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#main/covariance_models"
                }, 
                {
                    "type": "File", 
                    "id": "#main/diamondSeqdb"
                }, 
                {
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-apps.yaml/apps", 
                            "type": "array"
                        }
                    ], 
                    "id": "#main/i5Applications"
                }, 
                {
                    "type": "Directory", 
                    "id": "#main/i5Databases"
                }, 
                {
                    "type": [
                        "null", 
                        {
                            "items": "#InterProScan-protein_formats.yaml/protein_formats", 
                            "type": "array"
                        }
                    ], 
                    "id": "#main/i5OutputFormat"
                }, 
                {
                    "type": "File", 
                    "id": "#main/phmmerSeqdb", 
                    "format": "http://edamontology.org/format_1929"
                }, 
                {
                    "type": [
                        "null", 
                        "#esl-reformat-replace.yaml/replace"
                    ], 
                    "id": "#main/replace"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#main/singleBestOnly"
                }, 
                {
                    "type": "File", 
                    "id": "#main/transcriptsFile", 
                    "format": "http://edamontology.org/format_1929"
                }
            ], 
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                }, 
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "$import": "#esl-reformat-replace.yaml"
                        }, 
                        {
                            "$import": "#BUSCO-assessment_modes.yaml"
                        }, 
                        {
                            "$import": "#InterProScan-apps.yaml"
                        }, 
                        {
                            "$import": "#InterProScan-protein_formats.yaml"
                        }
                    ]
                }
            ], 
            "outputs": [
                {
                    "outputSource": "#main/identify_coding_regions/bed_output", 
                    "type": "File", 
                    "id": "#main/bed_output"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/blastOutput", 
                    "type": "Directory", 
                    "id": "#main/busco_blast_output"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/fullTable", 
                    "type": "File", 
                    "id": "#main/busco_full_table"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/hmmerOutput", 
                    "type": "Directory", 
                    "id": "#main/busco_hmmer_output"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/missingBUSCOs", 
                    "type": "File", 
                    "id": "#main/busco_missing_buscos"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/shortSummary", 
                    "type": "File", 
                    "id": "#main/busco_short_summary"
                }, 
                {
                    "outputSource": "#main/run_transcriptome_assessment/translatedProteins", 
                    "type": "Directory", 
                    "id": "#main/busco_translated_proteins"
                }, 
                {
                    "outputSource": "#main/identify_coding_regions/coding_regions", 
                    "type": "File", 
                    "id": "#main/coding_regions"
                }, 
                {
                    "outputSource": "#main/identify_nc_rna/deoverlapped_matches", 
                    "type": "File", 
                    "id": "#main/deoverlapped_matches"
                }, 
                {
                    "outputSource": "#main/calculate_diamond_matches/matches", 
                    "type": "File", 
                    "id": "#main/diamond_matches"
                }, 
                {
                    "outputSource": "#main/identify_coding_regions/gff3_output", 
                    "type": "File", 
                    "id": "#main/gff3_output"
                }, 
                {
                    "outputSource": "#main/functional_analysis/i5Annotations", 
                    "type": "File", 
                    "id": "#main/i5Annotations"
                }, 
                {
                    "outputSource": "#main/identify_coding_regions/peptide_sequences", 
                    "type": "File", 
                    "id": "#main/peptide_sequences"
                }, 
                {
                    "outputSource": "#main/calculate_phmmer_matches/matches", 
                    "type": "File", 
                    "id": "#main/phmmer_matches"
                }, 
                {
                    "outputSource": "#main/remove_asterisks_and_reformat/reformatted_sequences", 
                    "type": "File", 
                    "id": "#main/reformatted_sequences"
                }
            ], 
            "label": "Transcripts annotation workflow", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "steps": [
                {
                    "out": [
                        "#main/calculate_diamond_matches/matches"
                    ], 
                    "in": [
                        {
                            "source": "#main/blockSize", 
                            "id": "#main/calculate_diamond_matches/blockSize"
                        }, 
                        {
                            "source": "#main/diamondSeqdb", 
                            "id": "#main/calculate_diamond_matches/databaseFile"
                        }, 
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/calculate_diamond_matches/queryInputFile"
                        }
                    ], 
                    "run": "#Diamon.blastx-v0.9.21.cwl", 
                    "id": "#main/calculate_diamond_matches", 
                    "label": "Calculates Diamond matches"
                }, 
                {
                    "out": [
                        "#main/calculate_phmmer_matches/matches", 
                        "#main/calculate_phmmer_matches/programOutput"
                    ], 
                    "in": [
                        {
                            "source": "#main/identify_coding_regions/peptide_sequences", 
                            "id": "#main/calculate_phmmer_matches/seqFile"
                        }, 
                        {
                            "source": "#main/phmmerSeqdb", 
                            "id": "#main/calculate_phmmer_matches/seqdb"
                        }
                    ], 
                    "run": "#phmmer-v3.2.cwl", 
                    "id": "#main/calculate_phmmer_matches", 
                    "label": "Calculates phmmer matches"
                }, 
                {
                    "doc": "Matches are generated against predicted CDS, using a sub set of databases\nfrom InterPro.\n", 
                    "out": [
                        "#main/functional_analysis/i5Annotations"
                    ], 
                    "run": "#InterProScan-v5-chunked-wf.cwl", 
                    "id": "#main/functional_analysis", 
                    "in": [
                        {
                            "source": "#main/i5Applications", 
                            "id": "#main/functional_analysis/applications"
                        }, 
                        {
                            "source": "#main/i5Databases", 
                            "id": "#main/functional_analysis/databases"
                        }, 
                        {
                            "source": "#main/remove_asterisks_and_reformat/reformatted_sequences", 
                            "id": "#main/functional_analysis/inputFile"
                        }, 
                        {
                            "source": "#main/i5OutputFormat", 
                            "id": "#main/functional_analysis/outputFormat"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#main/identify_coding_regions/peptide_sequences", 
                        "#main/identify_coding_regions/coding_regions", 
                        "#main/identify_coding_regions/gff3_output", 
                        "#main/identify_coding_regions/bed_output"
                    ], 
                    "in": [
                        {
                            "source": "#main/singleBestOnly", 
                            "id": "#main/identify_coding_regions/singleBestOnly"
                        }, 
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/identify_coding_regions/transcriptsFile"
                        }
                    ], 
                    "run": "#TransDecoder-v5-wf-2steps.cwl", 
                    "id": "#main/identify_coding_regions", 
                    "label": "Identifies candidate coding regions within transcript sequences"
                }, 
                {
                    "out": [
                        "#main/identify_nc_rna/deoverlapped_matches"
                    ], 
                    "in": [
                        {
                            "source": "#main/clanInfoFile", 
                            "id": "#main/identify_nc_rna/clan_info"
                        }, 
                        {
                            "source": "#main/cmsearchCores", 
                            "id": "#main/identify_nc_rna/cores"
                        }, 
                        {
                            "source": "#main/covariance_models", 
                            "id": "#main/identify_nc_rna/covariance_models"
                        }, 
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/identify_nc_rna/query_sequences"
                        }
                    ], 
                    "run": "#cmsearch-multimodel-wf.cwl", 
                    "id": "#main/identify_nc_rna", 
                    "label": "Identifies non-coding RNAs using Rfams covariance models"
                }, 
                {
                    "out": [
                        "#main/remove_asterisks_and_reformat/reformatted_sequences"
                    ], 
                    "in": [
                        {
                            "source": "#main/replace", 
                            "id": "#main/remove_asterisks_and_reformat/replace"
                        }, 
                        {
                            "source": "#main/identify_coding_regions/peptide_sequences", 
                            "id": "#main/remove_asterisks_and_reformat/sequences"
                        }
                    ], 
                    "run": "#esl-reformat.cwl", 
                    "id": "#main/remove_asterisks_and_reformat", 
                    "label": "Removes asterisks characters from given peptide sequences"
                }, 
                {
                    "out": [
                        "#main/run_transcriptome_assessment/shortSummary", 
                        "#main/run_transcriptome_assessment/fullTable", 
                        "#main/run_transcriptome_assessment/missingBUSCOs", 
                        "#main/run_transcriptome_assessment/hmmerOutput", 
                        "#main/run_transcriptome_assessment/translatedProteins", 
                        "#main/run_transcriptome_assessment/blastOutput"
                    ], 
                    "in": [
                        {
                            "source": "#main/buscoLineage", 
                            "id": "#main/run_transcriptome_assessment/lineage"
                        }, 
                        {
                            "source": "#main/buscoMode", 
                            "id": "#main/run_transcriptome_assessment/mode"
                        }, 
                        {
                            "source": "#main/buscoOutputName", 
                            "id": "#main/run_transcriptome_assessment/outputName"
                        }, 
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/run_transcriptome_assessment/sequenceFile"
                        }
                    ], 
                    "run": "#BUSCO-v3.cwl", 
                    "id": "#main/run_transcriptome_assessment", 
                    "label": "Performs transcriptome assessment using BUSCO"
                }
            ], 
            "id": "#main", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "class": "Workflow"
        }, 
        {
            "inputs": [
                {
                    "default": "full_cmsearch_output", 
                    "type": "string", 
                    "id": "#cmsearch-multimodel-wf.cwl/catOutputFileName"
                }, 
                {
                    "type": "File", 
                    "id": "#cmsearch-multimodel-wf.cwl/clan_info"
                }, 
                {
                    "type": "int", 
                    "id": "#cmsearch-multimodel-wf.cwl/cores"
                }, 
                {
                    "type": {
                        "items": "File", 
                        "type": "array"
                    }, 
                    "id": "#cmsearch-multimodel-wf.cwl/covariance_models"
                }, 
                {
                    "type": "File", 
                    "id": "#cmsearch-multimodel-wf.cwl/query_sequences"
                }
            ], 
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputSource": "#cmsearch-multimodel-wf.cwl/remove_overlaps/deoverlapped_matches", 
                    "id": "#cmsearch-multimodel-wf.cwl/deoverlapped_matches"
                }
            ], 
            "label": "Identifies non-coding RNAs using Rfams covariance models", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "steps": [
                {
                    "run": "#infernal-cmsearch-v1.1.2.cwl", 
                    "scatter": "#cmsearch-multimodel-wf.cwl/cmsearch/covariance_model_database", 
                    "in": [
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/covariance_models", 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/covariance_model_database"
                        }, 
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/cores", 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/cpu"
                        }, 
                        {
                            "default": true, 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/omit_alignment_section"
                        }, 
                        {
                            "default": true, 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/only_hmm"
                        }, 
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/query_sequences", 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/query_sequences"
                        }, 
                        {
                            "default": 1000, 
                            "id": "#cmsearch-multimodel-wf.cwl/cmsearch/search_space_size"
                        }
                    ], 
                    "label": "Search sequence(s) against a covariance model database", 
                    "id": "#cmsearch-multimodel-wf.cwl/cmsearch", 
                    "out": [
                        "#cmsearch-multimodel-wf.cwl/cmsearch/matches", 
                        "#cmsearch-multimodel-wf.cwl/cmsearch/programOutput"
                    ]
                }, 
                {
                    "out": [
                        "#cmsearch-multimodel-wf.cwl/concatenate_matches/result"
                    ], 
                    "run": "#concatenate.cwl", 
                    "id": "#cmsearch-multimodel-wf.cwl/concatenate_matches", 
                    "in": [
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/cmsearch/matches", 
                            "id": "#cmsearch-multimodel-wf.cwl/concatenate_matches/files"
                        }, 
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/catOutputFileName", 
                            "id": "#cmsearch-multimodel-wf.cwl/concatenate_matches/outputFileName"
                        }
                    ]
                }, 
                {
                    "out": [
                        "#cmsearch-multimodel-wf.cwl/remove_overlaps/deoverlapped_matches"
                    ], 
                    "in": [
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/clan_info", 
                            "id": "#cmsearch-multimodel-wf.cwl/remove_overlaps/clan_information"
                        }, 
                        {
                            "source": "#cmsearch-multimodel-wf.cwl/concatenate_matches/result", 
                            "id": "#cmsearch-multimodel-wf.cwl/remove_overlaps/cmsearch_matches"
                        }
                    ], 
                    "run": "#cmsearch-deoverlap-v0.02.cwl", 
                    "id": "#cmsearch-multimodel-wf.cwl/remove_overlaps", 
                    "label": "Remove lower scoring overlaps from cmsearch --tblout files."
                }
            ], 
            "id": "#cmsearch-multimodel-wf.cwl", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "class": "Workflow"
        }
    ]
}