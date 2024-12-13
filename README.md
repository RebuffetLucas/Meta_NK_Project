# Meta-NK consortium: A standardized terminology that provides a basis for subsequent studies of the role of human Natural killer (NK) cells in health and disease.

## Article information

**Title:** High-Dimensional Single-Cell Analysis of Human Natural Killer Cell Heterogeneity

**Authors:**
Lucas Rebuffet,1, Janine E. Melsen,2,3, Bertrand Escalière,1, Daniela Basurto-Lozada,4,5, Avinash Bhandoola,6, Niklas K. Björkström,7, Yenan T. Bryceson,8, Roberta Castriconi9,10, Frank Cichocki,11, Marco Colonna,12, Daniel M. Davis,13, Andreas Diefenbach,14,15, Yi Ding,6, Muzlifah Haniffa,4,5,16, Amir Horowitz,17,18,19, Lewis L. Lanier,20, Karl-Johan Malmberg,7,21,22, Jeffrey S. Miller,11, Lorenzo Moretta,23, Emilie Narni-Mancinelli,1, Luke A.J. O’Neill,24, Chiara Romagnani,25,26,27, Dylan G. Ryan,28, Simona Sivori,9,29, Dan Sun,7, Constance Vagne,30, Eric Vivier,1,30,31,32*

<small>
1 Aix Marseille Université, CNRS, INSERM, Centre d'Immunologie de Marseille-Luminy, Marseille, France.
2 Leiden University Medical Center, Willem-Alexander Children’s Hospital, Laboratory for Pediatric Immunology, Leiden, The Netherlands.
3 Leiden University Medical Center, Department of Immunology, Leiden, The Netherlands.
4 Wellcome Sanger Institute, Wellcome Genome Campus, Cambridge CB10 1SA, UK.
5 Biosciences Institute, Newcastle University, NE24HH, UK.
6 T Cell Biology and Development Unit, Laboratory of Genome Integrity, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA.
7 Center for Infectious Medicine, Department of Medicine Huddinge, Karolinska Institutet, Karolinska University Hospital, Stockholm, Sweden.
8 Center for Hematology and Regenerative Medicine, Department of Medicine Huddinge, Karolinska Institutet, Stockholm. Sweden.
9 Department of Experimental Medicine (DIMES), University of Genoa, 16132 Genoa, Italy.
10 Laboratory of Clinical and Experimental Immunology, IRCCS Istituto Giannina Gaslini, 16147 Genova, Italy.
11 Department of Medicine, University of Minnesota, Minneapolis, MN, USA.
12 Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, MO 63110, USA.
13 Department of Life Sciences, Imperial College London, Sir Alexander Fleming Building, South Kensington, London, SW7 2AZ, UK.
14 Laboratory of Innate Immunity, Institute of Microbiology, Infectious Diseases and Immunology (I-MIDI), Campus Benjamin Franklin, Charité - Universitätsmedizin Berlin, Corporate Member of Freie Universität Berlin and Humboldt-Universität zu Berlin, Berlin, Germany.
15 Mucosal and Developmental Immunology, Deutsches Rheuma-Forschungszentrum (DRFZ), an Institute of the Leibniz Association, Berlin, Germany.
16 Department of Dermatology and NIHR Biomedical Research Centre, Newcastle Hospitals NHS Foundation Trust, Newcastle upon Tyne NE1 4LP, UK.
17 Department of Immunology & Immunotherapy, The Marc and Jennifer Lipschultz Precision Immunology Institute, Icahn School of Medicine at Mount Sinai, New York, NY, USA.
18 Department of Oncological Sciences, The Tisch Cancer Institute, Icahn School of Medicine at Mount Sinai, New York, NY, USA.19 Department of Oncological Sciences, Icahn School of Medicine at Mount Sinai, New York, NY, USA.
20 Department of Microbiology and Immunology and the Parker Institute for Cancer Immunotherapy, University of California San Francisco, CA U.S.A.
21 Precision Immunotherapy Alliance, The University of Oslo, Norway.
22 The Institute for Cancer Research, Oslo University Hospital, Norway.
23 Tumor Immunology Unit, Bambino Gesù Children’s Hospital, IRCCS, Rome, Italy.
24 School of Biochemistry and Immunology, Trinity Biomedical Sciences Institute, Trinity College Dublin, Dublin, Ireland.
25 Institute of Medical Immunology, Charité Universitätsmedizin Berlin, Corporate Member of Freie Universität Berlin and Humboldt Universität zu Berlin, Berlin, Germany.
26 Innate Immunity, Deutsches Rheuma-Forschungszentrum Berlin (DRFZ), ein Leibniz Institut, Berlin, Germany.
27 Berlin University Alliance, Berlin, Germany.
28 MRC Mitochondrial Biology Unit, University of Cambridge, UK.
29 IRCCS Ospedale Policlinico San Martino, Genova, Italy
30 Innate Pharma Research Laboratories, Innate Pharma, Marseille, France.
31 APHM, Hôpital de la Timone, Marseille-Immunopôle, Marseille, France. 
32 Paris-Saclay Cancer Cluster, Le Kremlin-Bicêtre, France.
* Corresponding author: vivier@ciml.univ-mrs.fr
</small>


**Summary:**
Natural killer (NK) cells are innate lymphoid cells (ILCs) that serve as a first line of defense against pathogens and tumors. Historically, their classification hinged on a limited array of surface protein markers. However, recent advancements in single-cell technologies, such as single-cell RNA sequencing (scRNA-seq) and cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq), have unveiled a more complex and nuanced understanding of NK cells. This has led to variations in nomenclature and inconsistencies across scientific literature. In our study, we utilized these technologies to dissect the heterogeneity of NK cells. We identified three prominent NK cell subsets in healthy human blood: NK1, NK2 and NK3, further differentiated into six distinct subgroups. All subsets were present in all donors, irrespective of their human cytomegalovirus (HCMV) status. Our findings delineate the unique molecular characteristics, key transcription factors, principal biological functions, significant metabolic traits, and cytokine responses of each subgroup. These data also suggest two separate ontogenetic origins for NK cells, leading to divergent transcriptional trajectories. We also analyzed the distribution of NK cell subsets in the lung, tonsils and intraepithelial lymphocytes isolated from healthy individuals and in 22 tumor types. This standardized terminology aims at fostering clarity and consistency in future research, thereby improving cross-study comparisons.


## Goal of the github
This github project contains the instructions and material to reproduce the analyses reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. 


## Repository Structure: `03_Script`

The `03_Script` directory organizes code into modules, each addressing specific analyses or explorations: 

```plaintext
03_Script/
├── 01_GlobalHeterogeneity               # Emergence of a new NK classification from meta-analysis
├── 02_RegulatoryNetworkAnalysis         # SCENIC+ analysis to identify regulatory elements
├── 03_GO_KEGG                           # Functional profiling via gene ontology and KEGG enrichment
├── 04_RNAvelocity                       # RNA velocity analysis for differentiation trajectories
├── 05_Monocle                           # Monocle3 analysis for trajectory inference
├── 10_VerifOnV3Romagnani                # Verification of classification on external datasets
├── 11_Comparing_Tang_Paper              # Comparison with tumor datasets from Tang's paper
├── 12_Using_CITEseq_Seurat              # Analysis of CITE-seq data in Seurat
├── 13_Destiny                           # Diffusion map analysis for trajectory identification
├── 14_SCTonOurData                      # Single-cell transformation preprocessing
├── 17_Metabo_Investigation              # Metabolic profiling of NK cell subsets
├── 19_SetUp_Data_Portal                 # Setting up the data portal for public access
├── 21_FACS_Data_Additional_analysis     # Flow cytometry additional analyses
├── 22_Comparing_Colonna_Paper           # Verification of NK subsets in healthy tissues
├── 23_MarkerCount_Analysis              # Marker count analysis for subsets
├── 24_Check_LabelTransfer_Reliability   # Evaluating label transfer reliability
├── 25_Check_Signatures_Discrim_NK       # Assessing discriminative power of signatures
└── 26_Cytosig_Analysis                  # Analysis of NK response to cytokines

```

# Meta-NK Consortium: A Standardized Terminology for Studying Human Natural Killer (NK) Cells in Health and Disease

## Key Results

1. **A New Classification for NK Cells** 
   - Based on meta-analysis of single-cell data, this study introduces a novel classification of NK cells (*03_Script/01_GlobalHeterogeneity*), identifying three primary subsets (*NK1, NK2, NK3*) and six distinct subgroups with unique molecular and functional profiles. 

2. **Verification Across Datasets and Tissues** 
   - The proposed classification was validated across external datasets (*03_Script/10_VerifOnV3Romagnani*), confirming its robustness. 
   - The subsets were further identified in healthy tissues (*03_Script/22_Comparing_Colonna_Paper*) and tumors (*03_Script/11_Comparing_Tang_Paper*). 

3. **Gene Regulation Networks (GRNs)** 
   - Regulatory network analysis using SCENIC (*03_Script/02_RegulatoryNetworkAnalysis*) highlighted key transcription factors driving the functional specialization of NK subsets. 

4. **Functional and Metabolic Profiling** 
   - Functional annotation using gene ontology and KEGG enrichment (*03_Script/03_GO_KEGG*) characterized the unique biological roles of NK subsets. 
   - Metabolic investigation (*03_Script/17_Metabo_Investigation*) revealed significant metabolic traits associated with each subset. 
   - Cytokine response profiling (*03_Script/26_Cytosig_Analysis*) demonstrated the differential responsiveness of subsets to various cytokines. 

5. **Differentiation Trajectories** 
   - A three-step approach was employed to trace NK cell differentiation trajectories: 
     - RNA velocity analysis (*03_Script/04_RNAvelocity*). 
     - Pseudotime analysis with Monocle3 (*03_Script/05_Monocle*). 
     - Diffusion map analysis (*03_Script/13_Destiny*). 

---

## Data Accessibility

All data from this study are publicly accessible through the [Human Cell Atlas Portal](https://collections.cellatlas.io/meta-nk). 
These datasets are part of the Human Cell Atlas effort, enabling researchers to download or interactively explore the data directly on the portal. 






