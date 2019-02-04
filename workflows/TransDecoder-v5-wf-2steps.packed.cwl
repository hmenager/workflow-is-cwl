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
            "$namespaces": {
                "edam": "http://edamontology.org/", 
                "s": "http://schema.org/"
            }, 
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
            "outputs": [
                {
                    "type": "File", 
                    "outputSource": "#main/predict_coding_regions/bed_output", 
                    "id": "#main/bed_output"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#main/predict_coding_regions/coding_regions", 
                    "id": "#main/coding_regions"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#main/predict_coding_regions/gff3_output", 
                    "id": "#main/gff3_output"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#main/predict_coding_regions/peptide_sequences", 
                    "id": "#main/peptide_sequences"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "out": [
                        "#main/extract_long_orfs/workingDir"
                    ], 
                    "in": [
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/extract_long_orfs/transcriptsFile"
                        }
                    ], 
                    "run": "#TransDecoder.LongOrfs-v5.cwl", 
                    "id": "#main/extract_long_orfs", 
                    "label": "Extracts the long open reading frames"
                }, 
                {
                    "out": [
                        "#main/predict_coding_regions/bed_output", 
                        "#main/predict_coding_regions/coding_regions", 
                        "#main/predict_coding_regions/gff3_output", 
                        "#main/predict_coding_regions/peptide_sequences"
                    ], 
                    "in": [
                        {
                            "source": "#main/extract_long_orfs/workingDir", 
                            "id": "#main/predict_coding_regions/longOpenReadingFrames"
                        }, 
                        {
                            "source": "#main/singleBestOnly", 
                            "id": "#main/predict_coding_regions/singleBestOnly"
                        }, 
                        {
                            "source": "#main/transcriptsFile", 
                            "id": "#main/predict_coding_regions/transcriptsFile"
                        }
                    ], 
                    "run": "#TransDecoder.Predict-v5.cwl", 
                    "id": "#main/predict_coding_regions", 
                    "label": "Predicts the likely coding regions"
                }
            ], 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/author": "Maxim Scheremetjew", 
            "label": "TransDecoder 2 step workflow, running TransDecoder.LongOrfs (step 1) followed by TransDecoder.Predict (step2)", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute, 2018", 
            "id": "#main"
        }
    ]
}