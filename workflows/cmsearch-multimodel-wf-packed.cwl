{
    "cwlVersion": "v1.0", 
    "$schemas": [
        "https://schema.org/docs/schema_org_rdfa.html", 
        "http://edamontology.org/EDAM_1.16.owl"
    ], 
    "$graph": [
        {
            "class": "CommandLineTool", 
            "baseCommand": [
                "cmsearch"
            ], 
            "inputs": [
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/covariance_model_database", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1
                    }
                }, 
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/cpu", 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cpu"
                    }, 
                    "label": "Number of parallel CPU workers to use for multithreads"
                }, 
                {
                    "default": false, 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/cut_ga", 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--cut_ga"
                    }, 
                    "label": "use CM's GA gathering cutoffs as reporting thresholds"
                }, 
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/omit_alignment_section", 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--noali"
                    }, 
                    "label": "Omit the alignment section from the main output.", 
                    "doc": "This can greatly reduce the output volume."
                }, 
                {
                    "default": false, 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/only_hmm", 
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--hmmonly"
                    }, 
                    "label": "Only use the filter profile HMM for searches, do not use the CM", 
                    "doc": "Only filter stages F1 through F3 will be executed, using strict P-value\nthresholds (0.02 for F1, 0.001 for F2 and 0.00001 for F3). Additionally\na bias composition filter is used after the F1 stage (with P=0.02\nsurvival threshold). Any hit that survives all stages and has an HMM\nE-value or bit score above the reporting threshold will be output.\n"
                }, 
                {
                    "format": "http://edamontology.org/format_1929", 
                    "id": "#infernal-cmsearch-v1.1.2.cwl/query_sequences", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "streamable": true
                }, 
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/search_space_size", 
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "-Z"
                    }, 
                    "label": "search space size in *Mb* to <x> for E-value calculations"
                }
            ], 
            "outputs": [
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/matches", 
                    "doc": "http://eddylab.org/infernal/Userguide.pdf#page=60", 
                    "label": "target hits table, format 2", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(inputs.query_sequences.basename).cmsearch_matches.tbl"
                    }
                }, 
                {
                    "id": "#infernal-cmsearch-v1.1.2.cwl/programOutput", 
                    "label": "direct output to file, not stdout", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(inputs.query_sequences.basename).cmsearch.out"
                    }
                }
            ], 
            "doc": "Infernal (\"INFERence of RNA ALignment\") is for searching DNA sequence databases for RNA structure and sequence similarities. It is an implementation of a special case of profile stochastic context-free grammars called covariance models (CMs). A CM is like a sequence profile, but it scores a combination of sequence consensus and RNA secondary structure consensus, so in many cases, it is more capable of identifying RNA homologs that conserve their secondary structure more than their primary sequence.\nPlease visit http://eddylab.org/infernal/ for full documentation.\nVersion 1.1.2 can be downloaded from http://eddylab.org/infernal/infernal-1.1.2.tar.gz\n", 
            "label": "Infernal: Search sequence(s) against a covariance model database", 
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
            "hints": [
                {
                    "class": "SoftwareRequirement", 
                    "packages": [
                        {
                            "specs": [
                                "https://identifiers.org/rrid/RRID:SCR_011809"
                            ], 
                            "version": [
                                "1.1.2"
                            ], 
                            "package": "infernal"
                        }
                    ]
                }, 
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biocontainers/infernal:v1.1.2-1-deb_cv1"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#name": "covariance_model_database", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#format": "txt"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "cpu", 
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#optional": true
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "cut_ga", 
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "omit_alignment_section", 
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "only_hmm", 
                            "http://galaxyproject.org/cwl#type": "boolean", 
                            "http://galaxyproject.org/cwl#optional": true
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "query_sequences", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#format": "txt"
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "search_space_size", 
                            "http://galaxyproject.org/cwl#type": "integer", 
                            "http://galaxyproject.org/cwl#default_value": 1000, 
                            "http://galaxyproject.org/cwl#optional": true
                        }
                    ]
                }
            ], 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "id": "#infernal-cmsearch-v1.1.2.cwl", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "https://www.sevenbridges.comwrapperAuthor": "Maxim Scheremetjew", 
            "$namespaces": {
                "edam": "http://edamontology.org/", 
                "s": "http://schema.org/", 
                "sbg": "https://www.sevenbridges.com"
            }
        }, 
        {
            "class": "CommandLineTool", 
            "baseCommand": [
                "cmsearch-deoverlap.pl"
            ], 
            "inputs": [
                {
                    "id": "#cmsearch-deoverlap-v0.02.cwl/clan_information", 
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "inputBinding": {
                        "position": 0, 
                        "prefix": "--clanin"
                    }, 
                    "label": "clan information on the models provided", 
                    "doc": "Not all models provided need to be a member of a clan"
                }, 
                {
                    "id": "#cmsearch-deoverlap-v0.02.cwl/cmsearch_matches", 
                    "type": "File", 
                    "inputBinding": {
                        "position": 1, 
                        "valueFrom": "$(self.basename)"
                    }
                }
            ], 
            "outputs": [
                {
                    "id": "#cmsearch-deoverlap-v0.02.cwl/deoverlapped_matches", 
                    "doc": "http://eddylab.org/infernal/Userguide.pdf#page=60", 
                    "label": "target hits table, format 2", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.deoverlapped"
                    }
                }
            ], 
            "doc": "https://github.com/nawrockie/cmsearch_tblout_deoverlap/blob/master/00README.txt", 
            "label": "Cmsearch-deoverlap: Remove lower scoring overlaps from cmsearch --tblout files.", 
            "requirements": [
                {
                    "class": "EnvVarRequirement", 
                    "envDef": [
                        {
                            "envValue": "C", 
                            "envName": "LC_ALL"
                        }
                    ]
                }, 
                {
                    "class": "ResourceRequirement", 
                    "ramMin": 100, 
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
            "hints": [
                {
                    "class": "SoftwareRequirement", 
                    "packages": [
                        {
                            "specs": [
                                "https://github.com/nawrockie/cmsearch_tblout_deoverlap"
                            ], 
                            "version": [
                                "0.02"
                            ], 
                            "package": "cmsearch_tblout_deoverlap"
                        }
                    ]
                }, 
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biocrusoe/cmsearch-deoverlap"
                }, 
                {
                    "class": "http://galaxyproject.org/cwl#interface", 
                    "http://galaxyproject.org/cwl#inputs": [
                        {
                            "http://galaxyproject.org/cwl#name": "clan_information", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#format": "txt", 
                            "http://galaxyproject.org/cwl#optional": true
                        }, 
                        {
                            "http://galaxyproject.org/cwl#name": "cmsearch_matches", 
                            "http://galaxyproject.org/cwl#type": "data", 
                            "http://galaxyproject.org/cwl#format": "txt"
                        }
                    ]
                }
            ], 
            "id": "#cmsearch-deoverlap-v0.02.cwl", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "https://www.sevenbridges.comwrapperAuthor": "Maxim Scheremetjew"
        }, 
        {
            "class": "CommandLineTool", 
            "baseCommand": [
                "cat"
            ], 
            "inputs": [
                {
                    "id": "#concatenate.cwl/files", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }, 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "streamable": true
                }
            ], 
            "outputs": [
                {
                    "id": "#concatenate.cwl/result", 
                    "type": "File", 
                    "outputBinding": {
                        "glob": "result", 
                        "outputEval": "${ self[0].format = inputs.files[0].format;\n   return self; }\n"
                    }
                }
            ], 
            "requirements": [
                {
                    "class": "ResourceRequirement", 
                    "ramMin": 100, 
                    "coresMax": 1
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "alpine:3.7"
                }
            ], 
            "stdout": "result", 
            "id": "#concatenate.cwl", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0"
        }, 
        {
            "class": "Workflow", 
            "inputs": [
                {
                    "id": "#main/clan_info", 
                    "type": "File"
                }, 
                {
                    "id": "#main/cores", 
                    "type": "int"
                }, 
                {
                    "id": "#main/covariance_models", 
                    "type": {
                        "type": "array", 
                        "items": "File"
                    }
                }, 
                {
                    "id": "#main/query_sequences", 
                    "type": "File"
                }
            ], 
            "outputs": [
                {
                    "id": "#main/deoverlapped_matches", 
                    "outputSource": [
                        "#main/remove_overlaps/deoverlapped_matches"
                    ], 
                    "type": "File"
                }
            ], 
            "steps": [
                {
                    "id": "#main/cmsearch", 
                    "in": [
                        {
                            "id": "#main/cmsearch/covariance_model_database", 
                            "source": "#main/covariance_models"
                        }, 
                        {
                            "id": "#main/cmsearch/cpu", 
                            "source": "#main/cores"
                        }, 
                        {
                            "id": "#main/cmsearch/omit_alignment_section", 
                            "default": true
                        }, 
                        {
                            "id": "#main/cmsearch/only_hmm", 
                            "default": true
                        }, 
                        {
                            "id": "#main/cmsearch/query_sequences", 
                            "source": "#main/query_sequences"
                        }, 
                        {
                            "id": "#main/cmsearch/search_space_size", 
                            "default": 1000
                        }
                    ], 
                    "out": [
                        {
                            "id": "#main/cmsearch/matches"
                        }, 
                        {
                            "id": "#main/cmsearch/programOutput"
                        }
                    ], 
                    "run": "#infernal-cmsearch-v1.1.2.cwl", 
                    "label": "Search sequence(s) against a covariance model database", 
                    "scatter": [
                        "#main/cmsearch/covariance_model_database"
                    ]
                }, 
                {
                    "id": "#main/concatenate_matches", 
                    "in": [
                        {
                            "id": "#main/concatenate_matches/files", 
                            "source": [
                                "#main/cmsearch/matches"
                            ]
                        }
                    ], 
                    "out": [
                        {
                            "id": "#main/concatenate_matches/result"
                        }
                    ], 
                    "run": "#concatenate.cwl"
                }, 
                {
                    "id": "#main/remove_overlaps", 
                    "in": [
                        {
                            "id": "#main/remove_overlaps/clan_information", 
                            "source": "#main/clan_info"
                        }, 
                        {
                            "id": "#main/remove_overlaps/cmsearch_matches", 
                            "source": "#main/concatenate_matches/result"
                        }
                    ], 
                    "out": [
                        {
                            "id": "#main/remove_overlaps/deoverlapped_matches"
                        }
                    ], 
                    "run": "#cmsearch-deoverlap-v0.02.cwl", 
                    "label": "Remove lower scoring overlaps from cmsearch --tblout files."
                }
            ], 
            "requirements": [
                {
                    "class": "ScatterFeatureRequirement"
                }
            ], 
            "id": "#main", 
            "http://schema.org/copyrightHolder": "EMBL - European Bioinformatics Institute", 
            "http://schema.org/license": "https://www.apache.org/licenses/LICENSE-2.0", 
            "https://www.sevenbridges.comwrapperAuthor": "Maxim Scheremetjew"
        }
    ]
}