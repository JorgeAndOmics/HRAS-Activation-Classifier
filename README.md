
# HRAS Activation Classifier


HRAS Activation Classifier is a complete pipeline for the generation of a Random Forest model, capable of discriminating a given protein structure between its active or inactive status, based on homology with similar structures available in databases. It's part of the FMP of Omics Data Analysis in the University of VIC.

## Table of Contents

-   Installation
-   Usage
-   Usage Example
-   Contributing
-   License

## Installation

To use this pipeline, you will need to have Python 3.8.

Once you have Python installed, you have to install the WCN Standalone Library from Martin Floor repository [[Martin-Floor/WCN (github.com)](https://github.com/Martin-Floor/WCN)]. Please follow the instalation instructions in the repository, but do not create a custom environment just yet.

After this, you will have to acquire a license for the Pyrosetta library [[PyRosetta - Licensing PyRosetta](https://www.pyrosetta.org/home/licensing-pyrosetta)]. Install as advised after this.

Then, use the following command to create a custom environment. Use of Conda is advised for environment management:

    conda env create -f hras_classifier.yaml
Activate the environment with:

    conda activate hras_classifier
After this, we download the repository to run :

git clone https://github.com/JorgeAndOmics/HRAS-Activation-Classifier
cd HRAS-Activation-Classifier
cd Pipeline

## Usage

The HRAS Activation Classifier pipeline takes the following command line arguments:

`--uniprot_code`: The Uniprot identifier of the desired protein
`--identity_threshold `: The identity selection threshold for the retrieved structures 
`--query_coverage_threshold`: The query coverage selection threshold for the retrieved structures
`--active_ligands`: A comma-separated list of ligands that unequivocally identify a protein as active
`--inactive_ligands`: A comma-separated list of ligands that unequivocally identify a protein as inactive
`--gap_open_penalty`: Gap open penalty for the pairwise alignment
`--gap_extend_penalty`: Gap extend penalty for the pairwise alignment
`--seed`: The seed for the pseudo-random processes


To run the script, open a terminal or command prompt, navigate to the downloaded directory  `HRAS-Activation-Classifier`  file, and enter the following command:

    nextflow protein_activity_classifier.nf \  
      --uniprot_code <UNIPROT_CODE> \  
      --identity_threshold <IDENTITY_THRESHOLD> \  
      --query_coverage_threshold <QUERY_COVERAGE_THRESHOLD> \  
      --active_ligands <ACTIVE_LIGANDS> \  
      --inactive_ligands <INACTIVE_LIGANDS> \  
      --gap_open_penalty <GAP_OPEN_PENALTY> \  
      --gap_extend_penalty <GAP_EXTEND_PENALTY> \  
      --seed <SEED>  

Replace  each of the options  with the appropriate values for your analysis.

## Usage Example

The pipeline was created for the study of *Human HRAS* (Uniprot code P01112), but it is usable for any kind of enzyme. Suppose you want to generate a model for the protein *Citrate synthase*, Uniprot code P09948:

    nextflow protein_activity_classifier.nf \  
          --uniprot_code = P09948 \  
          --identity_threshold = 300 \  
          --query_coverage_threshold = 50 \  
          --active_ligands = acetyl-CoA,oxaloacetate \  
          --inactive_ligands = CoA,H> \  
          --gap_open_penalty = -1 \  
          --gap_extend_penalty = -1 \  
          --seed = 123  


The script will commence to download the required files and extract the data for them. The average time for a protein of ~200 amino acids long in a single 6 core-processor with hyperthreading was ~1h to retrieve and preprocess the data and ~2 days for training the model.


## Contributing

If you would like to contribute to the HRAS Activation Classifier script, please fork the repository, make your changes, and submit a pull request.

## License

This project is licensed under the terms of the MIT license. See the  `LICENSE`  file for more information.

## References
All structures used for training the model are available at the end of the preprocessing phase in the generated metadata file.
- Floor, M., WCN. GitHub. Retrieved 6 September 2023, from https://github.com/Martin-Floor/WCN 
- Olsson, M. H. M., Søndergaard, C. R., Rostkowski, M., & Jensen, J. H. (2011). PROPKA3: Consistent 
Treatment of Internal and Surface Residues in Empirical p K a Predictions. Journal of Chemical Theory 
and Computation, 7(2), 525–537. https://doi.org/10.1021/ct100578z .
- Chaudhury, S., Lyskov, S., & Gray, J. J. (2010). PyRosetta: A script-based interface for implementing 
molecular modeling algorithms using Rosetta. Bioinformatics, 26(5), 689–691. 
https://doi.org/10.1093/bioinformatics/btq007
- All structures used for training the model are available at the end of the preprocessing phase in the generated metadata file.
