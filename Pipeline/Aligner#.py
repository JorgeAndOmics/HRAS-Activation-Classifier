import propka.molecular_container
from Bio.Align import PairwiseAligner
from Bio import Seq, SeqIO, SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBList, PDBParser, PDBIO, PPBuilder, Select, Superimposer
from Bio.Data import IUPACData
import WCN
from propka.lib import loadOptions
import propka.molecular_container
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol, fa_elec
from itertools import combinations_with_replacement
from collections import defaultdict
import pandas as pd
import numpy as np
import urllib.request
import requests
import warnings
import pickle
import argparse
import os






class Datahub:
    '''
        Datahub class contains all the data used by the Aligner class.
        '''

    @staticmethod
    def directory_generator():
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir','full_pdb_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'wcn_full_pdb_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'wcn_full_pdb_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'chain_fasta_files')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'chain_fasta_files'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes'))
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles')):
            os.mkdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles'))

    @staticmethod
    def parsers_and_builders():
        pdblist = PDBList()
        pdbparser = PDBParser()
        polybuilder = PPBuilder()
        pdbio = PDBIO()
        return pdblist, pdbparser, polybuilder, pdbio

    @staticmethod
    def warning_silencer():
        warnings.filterwarnings('ignore', message='.*discontinuous.*')
        warnings.filterwarnings('ignore', message='.*Used element.*')

    canonical_seq = ''
    dict_str_chain_id = {} # Dictionary of structure:chain_id
    dict_chain_residues_short = {} # Sequences only filtered by identity percentage
    dict_chain_residues_long = {} # Sequences filtered by query coverage
    dict_chain_residues_coverage = {} # Coverage percent of the selected sequences post-query coverage filtering
    dict_chain_identity = {} # Dictionary of structure:identity percentage
    dict_chain_residues_postcover = {} # Sequences filtered by query coverage for alignment dataframe
    dict_chain_residues_postcover_aligned = {} # Sequences filtered by query coverage for alignment dataframe with gaps
    dict_metadata = defaultdict(dict) # Dictionary dictionaries of lists of metadata
    dict_activity_status = {} # Dictionary of activity status
    dict_residues_ppb = {} # Dictionary of structures post-trimming for pdb creation
    dict_residues_ppb_3L = {} # Dictionary of structures post-trimming in 3-letter code
    dict_wcn_alpha_metrics = {} # Dictionary of WCN metrics for alpha carbon atoms
    dict_wcn_alpha_metrics_aligned = {} # Dictionary of WCN metrics for alpha carbon atoms with alignment gaps
    dict_wcn_full_metrics = {} # Dictionary of WCN metrics for all atoms
    dict_wcn_b_factors = {} # Dictionary of b-factors for alpha carbons for substitution
    dict_propka_metrics = {} # Dictionary of propka metrics
    dict_propka_metrics_aligned = {} # Dictionary of propka metrics with alignment gaps
    dict_pyrosetta_total_energy = defaultdict(list) # Dictionary of total energy values
    dict_pyrosetta_total_energy_aligned = {}  # Dictionary of total energy values
    dict_pyrosetta_vanderwaals1 = defaultdict(list) # Dictionary of vanderwaals1 energy values
    dict_pyrosetta_vanderwaals1_aligned = {}  # Dictionary of vanderwaals1 energy values
    dict_pyrosetta_vanderwaals2 = defaultdict(list) # Dictionary of vanderwaals2 energy values
    dict_pyrosetta_vanderwaals2_aligned = {}  # Dictionary of vanderwaals2 energy values
    dict_pyrosetta_solvation = defaultdict(list) # Dictionary of solvation energy values
    dict_pyrosetta_solvation_aligned = {}  # Dictionary of solvation energy values
    dict_pyrosetta_electrostatics = defaultdict(list) # Dictionary of electrostatics energy values
    dict_pyrosetta_electrostatics_aligned = {}  # Dictionary of electrostatics energy values
    dict_rmsd_values = {} # Dictionary of RMSD values
    dict_pdb_coord_alpha = {} # Dictionary of aligned PDB coordinates of c-alpha atoms
    uniprot_accession_code = 'P01112'
    list_active_ligands = ['GTP', 'GNP']
    list_inactive_ligands = ['GDP', 'GSP']
    query_coverage_threshold = 53
    identity_threshold = 30
    gap_open_penalty = -0.2
    gap_extend_penalty = -0.2
    dataframe_alignment = ''
    dataframe_alignment_substitutions = ''
    dataframe_metadata = ''
    dataframe_wcn = ''
    dataframe_propka = ''





class ChainSelector(Select):
    '''
        ChainSelector class is used to select a chain from a PDB file.
        Esto es cremísima. Se usa junto con PDBIO.io.save() para guardar un archivo PDB con una sola cadena. Amazing!
    '''
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if chain.get_id() == self.chain_id else 0

class ResidueSelector(Select):
    def __init__(self, residue_id):
        self.residue_id = residue_id

    def accept_residue(self, residue):
        return residue.get_id()[1] == self.residue_id


def download_sequence(uniprot_code):
    '''
        Downloads the target sequence from Uniprot by making an HTTP request.

        Returns:
            sequence sequence.seq (Sequence Object): Target sequence.
        '''
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_code}.fasta'  # Constructs URL for Aligner.uniprot_code
    response = urllib.request.urlopen(url)  # Sends an HTTP request to the URL
    sequence: object = response.read().decode('utf-8')  # Reads the response and decodes it to utf-8
    sequence = sequence.splitlines()  # Splits the sequence into lines
    sequence = ''.join(line for line in sequence if not line.startswith('>'))  # Constructs the sequence without header
    Datahub.canonical_seq = Seq.Seq(sequence) # Stores the canonical sequence as Seq object in Datahub class
    return Seq.Seq(sequence)



def blast_sequence(sequence):
    '''
        Performs a BLAST search on the target sequence.

        Returns:
            blast_record (BLAST Record): BLAST record.
        '''
    print('Performing BLAST search...')

    result_handle = NCBIWWW.qblast(
        'blastp',
        'pdb',
        sequence,
        expect=0.01,
        hitlist_size=10000,
        word_size=3
    )

    return NCBIXML.read(result_handle)  # Parses the BLAST result



def switch_blaster():
    try:
        if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'blast_seq.pickle')):
            print('BLAST Pickle file not found. Starting new query...')
            # Save the BLAST result to a pickle file
            with open(os.path.join('./temp_dir/pickles', 'blast_seq.pickle'), "wb") as f:
                query_sequence = download_sequence(Datahub.uniprot_accession_code)
                print(query_sequence)
                blast_sequence_result = blast_sequence(query_sequence)
                print('Dumping BLAST record to pickle file...')
                pickle.dump(blast_sequence_result, f)
                print('BLAST record saved to pickle file!')
                return blast_sequence_result
        elif os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'blast_seq.pickle')):
            print('BLAST Pickle file found! Reading...')
        # Load the BLAST result from a pickle file
            def read_pickle(file_name):
                '''
                Reads a pickle file and returns the object.
                :param file_name:
                :return:
                '''
                download_sequence(Datahub.uniprot_accession_code)
                with open(file_name, "rb") as f:
                    return pickle.load(f)
            return read_pickle(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'blast_seq.pickle'))
    except Exception as e:
        print(e)
        print('Error with BLAST record!')



def get_chain_ids(blast_record):
    '''
        Gets the chain IDs from the BLAST record.

        Returns:
            chain_ids (List): List of chain IDs.
        '''
    for alignment in blast_record.alignments:
        for hsps in alignment.hsps:
            real_perc_ident = (hsps.identities / hsps.align_length) * 100
            Datahub.dict_chain_identity[alignment.hit_id.split('|')[1]] = real_perc_ident
            if real_perc_ident >= Datahub.identity_threshold:
                Datahub.dict_str_chain_id[alignment.hit_id.split('|')[1]] = alignment.hit_id.split('|')[2] # First element is the PDB ID, second element is the chain ID







def get_pdb_files():
    '''
    Downloads the PDB files from the PDB database.
    :return:
    '''
    pdblist, pdbparser, polybuilder, pdbio = Datahub.parsers_and_builders()
    for key, value in Datahub.dict_str_chain_id.items():
        try:
            if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files', f'pdb{key.lower()}.ent')):
                pdblist.retrieve_pdb_file(key, pdir=os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files'), file_format='pdb')
                print(f'Saved {key}.ent at temporary directory')
                structure = pdbparser.get_structure(key, os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files', f'pdb{key.lower()}.ent'))
                pdbio.set_structure(structure)
                pdbio.save(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'), ChainSelector(value))
                print(f'Saved c_{key}.pdb from chain {value} at temporary directory')
        except FileNotFoundError as e:
            print(e)
            continue





def chain_residue_parser():
    '''
    Parses the PDB files to get the chain IDs and the residues.
    :return:
    '''
    pdblist, pdbparser, polybuilder, pdbio = Datahub.parsers_and_builders()
    for filename in os.listdir('./temp_dir/full_pdb_files'):
        num_chains = 0
        print(f'Parsing {filename} for chain sequence...')
        if filename.endswith('.pdb'):
            structure = pdbparser.get_structure(filename, os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', filename))
            for polypep in polybuilder.build_peptides(structure):
                sequence = polypep.get_sequence()
                Datahub.dict_chain_residues_short[f'{str(filename.split("_")[1].split(".")[0])}'] = Seq.Seq(str(sequence))
            for model in structure:
                num_chains += len(model)
        if num_chains > 1:  # If there are more than one chain, delete the sequence from the dictionary
            print(f'{filename} has {num_chains} chains')
            del Datahub.dict_chain_residues_short[f'{str(filename.split("_")[1].split(".")[0])}']
    for key, value in Datahub.dict_chain_residues_short.items():
        print(f'Parsing {key} for query coverage...')
        try:
            if (Datahub.query_coverage_threshold / 100) <= (len(value) / len(Datahub.canonical_seq)) <= 1: # If the query coverage is greater than the threshold, add the sequence to the dictionary
                Datahub.dict_chain_residues_long[key] = value
                Datahub.dict_chain_residues_coverage[key] = (len(value) / len(Datahub.canonical_seq)) * 100
        except ZeroDivisionError:
            continue
    print(f'{len(Datahub.dict_chain_residues_long.keys())} chains selected with {Datahub.query_coverage_threshold}% query coverage')


def metadata_generator():
    '''
    Generates the metadata for the dataframe.
    :return:
    '''
    print(Datahub.dict_chain_residues_long)
    def fetch_ligand_info(pdb_id):
        """Fetch ligand information from RCSB PDB API for a given PDB ID."""
        entry_info_url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
        try:
            response = requests.get(entry_info_url)
            if response.status_code == 200:
                entry_data = response.json()
                nonpolymer_entities = [entity for entity in
                                       entry_data['rcsb_entry_container_identifiers']['non_polymer_entity_ids']]

                ligands = []
                ligand_abbreviations = []
                for entity_id in nonpolymer_entities:
                    ligand_info_url = f'https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}'
                    ligand_response = requests.get(ligand_info_url)
                    ligand_data = ligand_response.json()

                    # Extract ligand name based on the observed structure
                    ligand_name = ligand_data.get('rcsb_nonpolymer_entity', {}).get('pdbx_description',
                                                                                    ligand_data.get(
                                                                                        'pdbx_entity_nonpoly', {}).get(
                                                                                        'name', 'Unknown Ligand'))

                    # Extract ligand abbreviation based on the observed structure
                    ligand_abbr = ligand_data.get('pdbx_entity_nonpoly', {}).get('comp_id', 'Unknown Ligand')

                    ligands.append(ligand_name)
                    ligand_abbreviations.append(ligand_abbr)

                return [", ".join(ligand_abbreviations), ", ".join(ligands)]

            else:
                return f"HTTP Error {response.status_code} when fetching entry data"

        except Exception as e:
            return f'Error: {str(e)}'


    def chain_fasta_generator():
        '''
        Generates a FASTA file for each chain.
        :return:
        '''
        # Loop over all files in the current directory
        for file_name in os.listdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files')):
            structure_name = file_name.split('_')[1].split('.')[0]
            # Check if this file is in the dictionary
            if structure_name in Datahub.dict_chain_residues_long.keys():
                # Parse the structure
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure(structure_name, os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', file_name))

                # Iterate over the residues in the structure and extract the sequence
                print(f'Parsing {structure_name} for FASTA file generation...')
                sequence = ''
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname().title() in IUPACData.protein_letters_3to1_extended.keys():  # residue_map is a dictionary that maps PDB residue names to single-letter codes
                                sequence += IUPACData.protein_letters_3to1_extended[residue.get_resname().title()]

                # Create a SeqRecord and write it to a FASTA file
                print(f'Writing {structure_name} to FASTA file...')
                seq_record = SeqRecord.SeqRecord(Seq.Seq(sequence), id=structure_name, description='')
                print(seq_record)
                SeqIO.write(seq_record, os.path.join(os.path.dirname(__file__), 'temp_dir', 'chain_fasta_files', f'{structure_name}.fasta'),
                            'fasta')
                print(f'{structure_name} exported to FASTA file!')

                # Append the sequence to the output file
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'chain_fasta_files', 'batch_file.fasta'), 'w') as f:
                    SeqIO.write(seq_record, f, 'fasta')
                print(f'{structure_name} sequence appended to the batch output file!')

    def check_activity(pdb_file):
        '''
            Reads the PDB file and checks whether it is active or inactive based on the keywords inserted as input.

            Args:
                pdb_file (str): Path to the PDB file.

            Returns:
                None
        '''

        def check_status(pdb_file):
            with open(pdb_file, 'r') as file:  # Opens the PDB file within context manager in reading mode for ease of
                # access and automatic closure
                pdb_text = file.read()  #
                # if Datahub.list_activity_ligands in pdb_text:  # Checks whether any of the keywords searched are in the PDB
                if any(keyword in Datahub.dict_metadata[structure_name]["RCSB Ligand Codes"] for keyword in
                       Datahub.list_active_ligands):
                    if any(keyword in Datahub.dict_metadata[structure_name]["RCSB Ligand Codes"] for keyword in
                       Datahub.list_inactive_ligands):
                        return 'unknown'
                    return 'active'
                elif any(keyword in Datahub.dict_metadata[structure_name]["RCSB Ligand Codes"] for keyword in
                       Datahub.list_inactive_ligands):
                    return 'inactive'
                else:
                    return 'unknown'

        print(f'Retrieving activity status from {structure_name}..')
        if structure_name in Datahub.dict_chain_residues_long.keys():
            status = check_status(pdb_file)
            return status

    # Define PDB parser
    parser = PDBParser(QUIET=True)

    # Get a list of all PDB files in the directory
    pdb_files = [f for f in os.listdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files')) if f.endswith('.pdb')]

    # Process each PDB file
    for pdb_file in pdb_files:
        try:
            structure_name = pdb_file.split('_')[1].split('.')[0]

            if structure_name in Datahub.dict_chain_residues_long.keys():
                # Load the structure
                print(f'Parsing {structure_name} for metadata...')
                structure = parser.get_structure(pdb_file[:-4], os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', pdb_file))


                # Counting number of chains, atoms, residues, hetero residues, waters, and ligands
                num_chains = 0
                num_atoms = 0
                num_residues = 0
                num_hetero_residues = 0
                num_ligands = 0
                amino_acid_count = defaultdict(int)
                ligands = []



                for model in structure:
                    for chain in model:
                        for residue in chain:
                            num_residues += 1
                            if residue.id[0].startswith('H_'):
                                num_hetero_residues += 1
                                if residue.resname != 'HOH':
                                    ligands.append(residue.resname)
                                    num_ligands += 1
                            else:
                                amino_acid_count[residue.resname] += 1
                            for atom in residue:
                                num_atoms += 1


                # Get metadata from RCSB
                try:
                    try:
                        print(f'Requesting {structure_name} metadata from RCSB...')
                        response = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{structure_name}')
                        response_species = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{structure_name}/1")
                        if response.status_code == 200:
                            data = response.json()
                            doi = data['rcsb_primary_citation'].get('pdbx_database_id_doi', None)
                            Datahub.dict_metadata[structure_name] = {}
                            Datahub.dict_metadata[structure_name]["Title"] = data['struct'].get('title', None)
                            Datahub.dict_metadata[structure_name]["Structure Details"] = data['struct'].get('pdbx_descriptor', None)
                        if response_species.status_code == 200:
                            data_species = response_species.json()
                            organism_info = data_species.get('rcsb_entity_source_organism', [{}])
                            if organism_info:
                                Datahub.dict_metadata[structure_name]["Source Organism"] = organism_info[0].get('ncbi_scientific_name')
                                Datahub.dict_metadata[structure_name]["Taxonomy ID"] = organism_info[0].get(
                                    'ncbi_taxonomy_id')
                            else:
                                Datahub.dict_metadata[structure_name]["Source Organism"] = 'No organism information found'
                                Datahub.dict_metadata[structure_name]["Taxonomy ID"] = 'No taxonomy ID found'
                    except Exception as e:
                        print(e)
                        continue

                    try:
                        # Check if a DOI is available before making the request
                        if doi is not None:
                            print(f'Requesting abstract for DOI: {doi} from Europe PMC...')
                            response_abstract = requests.get(
                                f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}&resultType=core&format=json')
                            results_abstract = response_abstract.json().get('resultList', {}).get('result', [])
                            if not results_abstract:
                                print(f"No results found for DOI: {doi}")
                                Datahub.dict_metadata[structure_name]["Abstract"] = 'No abstract found'
                            else:
                                Datahub.dict_metadata[structure_name]["Abstract"] = results_abstract[0].get('abstractText',
                                                                                                    None)
                        else:
                            print(f"No DOI found for {structure_name}")
                            Datahub.dict_metadata[structure_name]["Abstract"] = 'No DOI found'
                    except Exception as e:
                        print(e)
                        continue
                except Exception as e:
                    print(e)
                    pass


                # Extract the features

                structure_ent = parser.get_structure(pdb_file[:-4], os.path.join(os.path.dirname(__file__), 'temp_dir', 'ent_files', f'pdb{structure_name.lower()}.ent'))
                print(f'Parsing {structure_name} original ENT file for content...')
                Datahub.dict_metadata[structure_name]["Method"] = structure_ent.header["structure_method"]
                Datahub.dict_metadata[structure_name]["Resolution"] = structure_ent.header["resolution"]
                Datahub.dict_metadata[structure_name]["Original Number of Models"] = len(structure_ent)

                for model in structure_ent:
                    num_chains += len(model)



                print(f'Parsing {structure_name} PDB chain for file content...')
                Datahub.dict_metadata[structure_name]["Original Number of Chains"] = num_chains
                Datahub.dict_metadata[structure_name]["Selected Chain"] = Datahub.dict_str_chain_id[structure_name]
                Datahub.dict_metadata[structure_name]["Canonical Sequence Coverage"] = \
                f'{Datahub.dict_chain_residues_coverage[structure_name]}%'
                Datahub.dict_metadata[structure_name]["Identity Percentage"] = f'{Datahub.dict_chain_identity[structure_name]}%'
                Datahub.dict_metadata[structure_name]["Number of Atoms"] = num_atoms
                Datahub.dict_metadata[structure_name]["Number of Residues"] = num_residues
                Datahub.dict_metadata[structure_name]["Number of Hetero Residues"] = num_hetero_residues
                Datahub.dict_metadata[structure_name]["Number of Ligands in Chain"] = num_ligands
                Datahub.dict_metadata[structure_name]["Ligands in Chain"] = list(set(ligands)) if ligands else 'No ligands found'
                Datahub.dict_metadata[structure_name]["RCSB Ligand Codes"] = fetch_ligand_info(structure_name)[0]
                Datahub.dict_metadata[structure_name]["RCSB Ligands"] = fetch_ligand_info(structure_name)[1]
                print(f'RCSB Ligand Codes: {Datahub.dict_metadata[structure_name]["RCSB Ligand Codes"]}')
                print(f'RCSB Ligands: {Datahub.dict_metadata[structure_name]["RCSB Ligands"]}')
                Datahub.dict_metadata[structure_name]["Read Activity Status"] = check_activity(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', pdb_file))
                Datahub.dict_metadata[structure_name]["Number of Waters"] = amino_acid_count["HOH"]

                for amino_acid, count in amino_acid_count.items():
                    if amino_acid != "HOH":
                        Datahub.dict_metadata[structure_name][f"Number of {amino_acid}"] = count


        except Exception as e:
            print(e)
            continue


    # Generate the FASTA file for the batch
    chain_fasta_generator()


    # Convert the data into a dataframe
    df = pd.DataFrame(Datahub.dict_metadata).T
    Datahub.dataframe_metadata = df.sort_index()

    # Trim dataframe and selected post-query coverage chains from unknown activation status
    unknown_indices = df[df['Read Activity Status'] == 'unknown'].index.to_list()
    Datahub.dataframe_metadata = Datahub.dataframe_metadata.drop(unknown_indices)
    Datahub.dict_chain_residues_long = {key: value for key, value in Datahub.dict_chain_residues_long.items() if key in Datahub.dataframe_metadata.index.to_list()}
    print(f'{len(unknown_indices)} chains with unknown activation status removed from analysis')
    # print(f'{len(Datahub.dict_chain_residues_long)} chains remaining in analysis')

    Datahub.dataframe_metadata.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'metadata_dataframe.csv'))
    print('Metadata exported to dataframe!')


def sequence_adjuster():
    metadata = pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'metadata_dataframe.csv'), index_col=0)
    Datahub.dict_chain_residues_long = {key: value for key, value in Datahub.dict_chain_residues_long.items() if
                                            key in metadata.index.to_list()}
    print(f'{len(Datahub.dict_chain_residues_long)} chains remaining in analysis')

def dataframe_generator():
    columns = [f'pos{i + 1}:{Datahub.canonical_seq[i]}' for i in range(len(Datahub.canonical_seq))] # Create the columns for the dataframe
    Datahub.dataframe_alignment = pd.DataFrame(columns=columns) # Initialize the dataframe
    for key, value in Datahub.dict_chain_residues_long.items(): # Iterate over the sequences pre-alignment
        try:
            aligner = PairwiseAligner()  # Initialize the pairwise aligner
            aligner.mode = 'global'  # Set the alignment mode to local
            aligner.open_gap_score = Datahub.gap_open_penalty  # Increase the gap open penalty
            aligner.extend_gap_score = Datahub.gap_extend_penalty  # Increase the gap extend penalty
            #### Alignment ####
            alignments = aligner.align(Datahub.canonical_seq,
                                       value)  # Align canonical sequence with the sequence of the chain
            print(key)

            index_list = [  # Get the index of the target alignment where there are no gaps
                i
                for i, aa in enumerate(alignments[0][0])
                if aa != '-'
                   and i <= (len(Datahub.dict_chain_residues_long[key]) - 1)
            ]
            index_gap = [
                i
                for i, aa in enumerate(alignments[0][1])
                if aa == '-'
                   and i <= (len(Datahub.dict_chain_residues_long[key]) - 1)

            ]  # Get the index of the query alignment where there are gaps

            def insert_nas(list, indexes, value):
                for i in sorted(indexes):
                    list.insert(i, value)

            Datahub.dict_chain_residues_postcover = Datahub.dict_chain_residues_long.copy()
            Datahub.dict_chain_residues_postcover[key] = list(Datahub.dict_chain_residues_postcover[key])
            insert_nas(Datahub.dict_chain_residues_postcover[key], index_gap,
                       np.nan)  # Insert NaN's where there are gaps in the alignment query

            print('non aligned:')
            print(f'Length: {len(Datahub.dict_chain_residues_postcover[key])}: {Datahub.dict_chain_residues_postcover[key]}')
            Datahub.dict_chain_residues_postcover_aligned[key] = [Datahub.dict_chain_residues_postcover[key][i] for i in
                                                           index_list]  # Trim the aminoacids to the alignment
            print('aligned:')
            print(f'Length: {len(Datahub.dict_chain_residues_postcover_aligned[key])}: {Datahub.dict_chain_residues_postcover_aligned[key]}')
            Datahub.dict_chain_residues_postcover_aligned[key] = Datahub.dict_chain_residues_postcover_aligned[key][
                                                          : len(Datahub.canonical_seq)]
            print(len(Datahub.dict_chain_residues_postcover_aligned[key]))

            print('Alignment Query:')
            print(f'length: {len(alignments[0][1])}, {alignments[0][1]}')
            Datahub.dict_chain_residues_postcover_aligned[key].extend([np.nan] * (
                        len(Datahub.canonical_seq) - len(Datahub.dict_chain_residues_postcover_aligned[key])))
            print(f'Length post-extend: {len(Datahub.dict_chain_residues_postcover_aligned[key])}: {Datahub.dict_chain_residues_postcover_aligned[key]}')# Extend the alignment to the length of the canonical sequence with gaps
            Datahub.dataframe_alignment.loc[key] = Datahub.dict_chain_residues_postcover_aligned[key]  # Add the alignment to the dataframe in order
            Datahub.dataframe_alignment = Datahub.dataframe_alignment.sort_index()
        except Exception as e:
            print(e)
            continue

    print(Datahub.dataframe_alignment)
    Datahub.dataframe_alignment.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'alignment_dataframe.csv'))

    print('Alignment Exported to Dataframe')
            # print(key)
            # print(alignments[-1])


#################¿Al final?#####################
def dataframe_trimmer():
    '''
    Removes non-perfect alignments between columns
    :return:
    '''
    columns = Datahub.dataframe.columns.to_list()
    # Initialize pointers
    left = 0
    right = len(columns) - 1
    # Iterate over the columns towards the center
    while left < right:
        if (Datahub.dataframe[columns[left]] == '-').any():
            left += 1
        elif (Datahub.dataframe[columns[right]] == '-').any():
            right -= 1
        else:
            break
    Datahub.dataframe = Datahub.dataframe.iloc[:, left:right + 1]
    #### Poner aquí los otros dataframes para recortar igual





def dataframe_parser():
    '''
    Parses dataframe and extracts trimmed sequences to dictionary
    :return:
    '''
    # Temporarily convert dataframe columns to string sequence for extraction via lambda function
    dict_letters = IUPACData.protein_letters_1to3_extended
    df = Datahub.dataframe.apply(lambda row: ''.join(row.values.astype(str)), axis=1) # Convert dataframe to string row by row in the columns axis
    # Extract sequences to dictionary
    ppb_dict = df.to_dict()
    # Convert sequences to list of letters for pdb extraction
    # for key, value in ppb_dict.items():
    #     ppb_dict[key] = list(value)
    #     ppb_dict[key] = [dict_letters[letter].upper() for letter in value]


#########################¿Al final?#########################


def wcn_function():

    def wcn_alfa_parser():
        '''
        Extract the C-alpha metrics from the PDB files.
        '''
        aligner = PairwiseAligner()  # Initialize aligner
        aligner.mode = 'global'  # Set alignment mode
        aligner.open_gap_score = Datahub.gap_open_penalty  # Increase the gap open penalty
        aligner.extend_gap_score = Datahub.gap_extend_penalty  # Increase the gap extend penalty
        for key, value in Datahub.dict_chain_residues_long.items():  # Iterate over the dictionary post-alignment
            instance = WCN.WCNObject(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'))  # Initialize WCN object
            Datahub.dict_wcn_alpha_metrics[key] = instance.calculateCAlpha(normalized=False,
                                                                           reciprocal=False)  # Calculate C-alpha metrics
            #### Alignment ####
            alignments = aligner.align(Datahub.canonical_seq,
                                       value)  # Align canonical sequence with the sequence of the chain
            print(key)
            print(alignments[0])
            index_list = [  # Get the index of the target alignment where there are no gaps
                i
                for i, aa in enumerate(alignments[0][0])
                if aa != '-'
                   and i <= (len(Datahub.dict_wcn_alpha_metrics[key]) - 1)
            ]
            index_gap = [
                i
                for i, aa in enumerate(alignments[0][1])
                if aa == '-'
                   and i <= (len(Datahub.dict_wcn_alpha_metrics[key]) - 1)

            ]  # Get the index of the query alignment where there are gaps

            def insert_nas(list, indexes, value):
                for i in sorted(indexes):
                    list.insert(i, value)

            Datahub.dict_wcn_alpha_metrics[key] = list(Datahub.dict_wcn_alpha_metrics[key])
            insert_nas(Datahub.dict_wcn_alpha_metrics[key], index_gap,
                       np.nan)  # Insert NA's where there are gaps in the alignment query
            print('non aligned:')
            print(Datahub.dict_wcn_alpha_metrics[key])
            Datahub.dict_wcn_alpha_metrics_aligned[key] = [Datahub.dict_wcn_alpha_metrics[key][i] for i in
                                                           index_list]  # Trim the C-alpha metrics to the alignment
            print('aligned:')
            print(Datahub.dict_wcn_alpha_metrics_aligned[key])
            Datahub.dict_wcn_alpha_metrics_aligned[key] = Datahub.dict_wcn_alpha_metrics_aligned[key][
                                                          : len(Datahub.canonical_seq)]

    wcn_alfa_parser()
    def wcn_dataframe_generator():
        '''
        Generates a dataframe with the C-alpha metrics of the chains
        '''
        columns = [f'pos{i + 1}:{Datahub.canonical_seq[i]}' for i in
                   range(len(Datahub.canonical_seq))]  # Generate column names
        Datahub.dataframe_wcn = pd.DataFrame(columns=columns)  # Initialize dataframe
        for key, value in Datahub.dict_wcn_alpha_metrics_aligned.items():  # Iterate over the dictionary of calpha metrics
            value.extend([np.nan] * (len(Datahub.canonical_seq) - len(
                value)))  # Extend the list of calpha metrics to the length of the canonical sequence with gaps
            Datahub.dataframe_wcn.loc[key] = value  # Add the list of calpha metrics to the dataframe in order
        Datahub.dataframe_wcn = Datahub.dataframe_wcn.sort_index()
        print(Datahub.dataframe_wcn)
        Datahub.dataframe_wcn.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'wcn_dataframe.csv'))

    wcn_dataframe_generator()

    def wcn_full_recorder():
        '''
        Generates and stores all-atom WCN profiles for PDB files
        '''
        try:
            if not os.path.exists(os.path.join('./temp_dir/pickles', 'wcn_all_atom.pickle')):
                for file in os.listdir('./temp_dir/full_pdb_files'):
                    if file.endswith('.pdb'):
                        key = file.split('.')[0].split('_')[1]
                        instance = WCN.WCNObject(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', file))
                        Datahub.dict_wcn_full_metrics[key] = instance.calculateAllAtom(normalized=False)
                        print(f'All-atom WCN profile for {key} calculated!')
                # Save the all-atom WCN result to a pickle file
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'wcn_all_atom.pickle'), "wb") as f:
                    print('Dumping all-atom WCN profile to pickle file...')
                    pickle.dump(Datahub.dict_wcn_full_metrics, f)
                    print('All-atom WCN profile saved to pickle file!')
            elif os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'wcn_all_atom.pickle')):
                # Load the all-atom WCN profile from a pickle file
                def read_pickle(file_name):
                    '''
                    Reads a pickle file and returns the object.
                    :param file_name:
                    :return:
                    '''
                    with open(file_name, "rb") as f:
                        return pickle.load(f)

                Datahub.dict_wcn_full_metrics = read_pickle(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'wcn_all_atom.pickle'))
        except Exception as e:
            print(e)
            print('Error with all-atom WCN profile record!')

    def wcn_to_pdb():
        '''
        Generates a PDB file with the C-alpha WCN metrics of the chains instead of b-factors
        '''
        pdblist, pdbparser, polybuilder, pdbio = Datahub.parsers_and_builders()
        dict_iupac = IUPACData.protein_letters_1to3_extended
        dict_iupac = {key: value.upper() for key, value in dict_iupac.items()}
        for key, value in Datahub.dict_wcn_alpha_metrics.items():
            if not os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'wcn_full_pdb_files', f'w_{key}.pdb')):
                structure = pdbparser.get_structure(key, os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'))

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[1] <= len(Datahub.dict_wcn_alpha_metrics[key]) and not np.isnan(
                                    Datahub.dict_wcn_alpha_metrics[key][residue.id[1] - 1]):
                                wcn_value = Datahub.dict_wcn_alpha_metrics[key][residue.id[1] - 1]
                                for atom in residue:
                                    atom.set_bfactor(wcn_value)

                pdbio.set_structure(structure)
                pdbio.save(os.path.join(os.path.dirname(__file__), 'temp_dir', 'wcn_full_pdb_files', f'w_{key}.pdb'))
                print(f'Saved w_{key}.pdb at temporary directory')
    # wcn_full_recorder()
    # wcn_to_pdb()


######################### WCN ##############################


######################### WCN ##############################


######################### PROPKA #########################

def propka_function():
    def propka_file_generator():
        '''
        Parses the pdb files for their propka profile and stores the pKa values in a dictionary
        '''

        # Generate the propka files
        for key, value in Datahub.dict_chain_residues_long.items():  # Iterate over the dictionary post-alignment
            try:
                if os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files', f'p_{key}.pka')):
                    print(f'p_{key}.pka file already exists at temporary directory')
                    continue
                else:
                    options = [os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb')]
                    args = loadOptions(options)
                    parameters = read_parameter_file(args.parameters, Parameters())
                    molecule = propka.molecular_container.MolecularContainer(parameters, args)
                    molecule = read_molecule_file(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'), molecule)
                    molecule.calculate_pka()
                    molecule.write_pka(
                        os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files', f'p_{key}.pka'))  # It performs local alignments on its own
                    print(f'p_{key}.pka file saved at temporary directory')
            except Exception as e:
                continue

    propka_file_generator()

    def propka_parser():
        '''
        Parses the propka files and stores the pKa values in a dictionary
        :return:
        '''
        # Extract the pKa values from the propka files
        tuple_values = tuple(value.upper() for value in
                             IUPACData.protein_letters_1to3_extended.values())  # Generate a tuple of the full set of amino acids for pKa extraction from files
        tuple_values = tuple('   ' + value for value in tuple_values)
        for key, value in Datahub.dict_chain_residues_long.items():
            Datahub.dict_propka_metrics[key] = [np.nan] * len(value)
            if os.path.exists(os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files', f'p_{key}.pka')):
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'propk_files', f'p_{key}.pka'), 'r') as p:
                    text = p.readlines()
                    for line in text:
                        if line.startswith(tuple_values):
                            print('Problem')
                            print(line.split())
                            if len(line.split()[0]) > 3:
                                process_line = [line.split()[0][0:3], line.split()[0][3:], line.split()[1:]]
                            else:
                                process_line = line.split()
                            print('process_line')
                            print(process_line)
                            if int(process_line[1]) + 1 <= len(Datahub.dict_propka_metrics[key]) - 1:
                                Datahub.dict_propka_metrics[key][int(process_line[1]) - 1] = float(process_line[3])

        # Convert propka metrics to alignment
        aligner = PairwiseAligner()  # Initialize aligner
        aligner.mode = 'local'  # Set alignment mode
        for key, value in Datahub.dict_chain_residues_long.items():  # Iterate over the dictionary post-alignment
            try:
                #### Alignment ####
                alignments = aligner.align(Datahub.canonical_seq,
                                           value)  # Align canonical sequence with the sequence of the chain
                alignment = alignments[0][0]  # Get the first alignment
                # print(key)
                # print(alignments[0])

                Datahub.dict_propka_metrics_aligned[key] = Datahub.dict_propka_metrics[
                    key]  # For keeping naming consistency
                # print(key)
                # print(Datahub.dict_propka_metrics_aligned[key])
                Datahub.dict_propka_metrics_aligned[key] = Datahub.dict_propka_metrics_aligned[key][
                                                           :len(Datahub.canonical_seq)]
            except Exception as e:
                continue

    propka_parser()

    def propka_dataframe_generator():
        '''
        Generates a dataframe with the propka metrics of the chains
        :return:
        '''
        # Generate a dataframe with the propka metrics of the chains
        columns = [f'pos{i + 1}:{Datahub.canonical_seq[i]}' for i in
                   range(len(Datahub.canonical_seq))]  # Generate column names
        Datahub.dataframe_propka = pd.DataFrame(columns=columns)  # Initialize dataframe
        for key, value in Datahub.dict_propka_metrics_aligned.items():  # Iterate over the dictionary of propka metrics
            value.extend([np.nan] * (len(Datahub.canonical_seq) - len(
                value)))  # Extend the list of propka metrics to the length of the canonical sequence with gaps
            Datahub.dataframe_propka.loc[key] = value  # Add the list of propka metrics to the dataframe in order
        Datahub.dataframe_propka = Datahub.dataframe_propka.sort_index()
        print(Datahub.dataframe_propka)
        Datahub.dataframe_propka.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'propka_dataframe.csv'))

    propka_dataframe_generator()





###################### PROPKA #######################


###################### PYROSETTA #######################

def pyrosetta_function():
    def pyrosetta_energy_parser():
        '''
        Extract the Pyrosetta metrics from the PDB files.
        '''
        aligner = PairwiseAligner()  # Initialize aligner
        aligner.mode = 'global'  # Set alignment mode
        aligner.open_gap_score = Datahub.gap_open_penalty  # Increase the gap open penalty
        aligner.extend_gap_score = Datahub.gap_extend_penalty  # Increase the gap extend penalty
        for key, value in Datahub.dict_chain_residues_long.items():  # Iterate over the dictionary post-alignment

            # List of substrings to check for
            unrecognized_residues = ["HETATM"]

            def pyrosetta_pdb_generator():
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'), 'r') as f:
                    lines = f.readlines()
                print(f'Saving PyRosetta {key} PDB file without unrecognized residues')
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files', f'r_{key}.pdb'), 'w') as f:
                    for line in lines:
                        # Only write the line if none of the substrings are in it
                        if not any(residue in line for residue in unrecognized_residues):
                            f.write(line)

            def pyrosetta_pdb_corrector():
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files', f'r_{key}.pdb'), 'r') as f:
                    lines = f.readlines()
                print('Saving PyRosetta PDB file without Updated Unrecognized Residues')
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files', f'r_{key}.pdb'), 'w') as f:
                    for line in lines:
                        # Only write the line if none of the substrings are in it
                        if not any(residue in line for residue in unrecognized_residues):
                            f.write(line)

            pyrosetta_pdb_generator()

            while True:
                try:
                    # Initialize PyRosetta
                    init()

                    # Load the structure
                    pose = pose_from_pdb(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pyrosetta_pdb_files', f'r_{key}.pdb'))

                    # Create a score function
                    scorefxn = create_score_function('ref2015')

                    # Score the pose
                    scorefxn(pose)

                    # Get the per-residue energies
                    energies = pose.energies()

                    # Iterate over the residues
                    for i in range(1, pose.total_residue() + 1):
                        residue = pose.residue(i)
                        # Get the energy for the alpha carbon
                        energy = energies.residue_total_energy(i)
                        vanderwaals1 = energies.residue_total_energies(i)[fa_atr]
                        vanderwalls2 = energies.residue_total_energies(i)[fa_rep]
                        solvation = energies.residue_total_energies(i)[fa_sol]
                        electrostatics = energies.residue_total_energies(i)[fa_elec]
                        Datahub.dict_pyrosetta_total_energy[key].append(energy)
                        Datahub.dict_pyrosetta_vanderwaals1[key].append(vanderwaals1)
                        Datahub.dict_pyrosetta_vanderwaals2[key].append(vanderwalls2)
                        Datahub.dict_pyrosetta_solvation[key].append(solvation)
                        Datahub.dict_pyrosetta_electrostatics[key].append(electrostatics)
                        # print(f"Total Energy -> Residue {i}: {energy}")
                        # print(f"Vanderwaals1 -> Residue {i}: {vanderwaals1}")
                        # print(f"Vanderwaals2 -> Residue {i}: {vanderwalls2}")
                        # print(f"Solvation -> Residue {i}: {solvation}")
                        # print(f"Electrostatics -> Residue {i}: {electrostatics}")

                    #### Alignment ####
                    alignments = aligner.align(Datahub.canonical_seq,
                                               value)  # Align canonical sequence with the sequence of the chain
                    print(key)
                    print(alignments[0])
                    index_list = [  # Get the index of the target alignment where there are no gaps
                        i
                        for i, aa in enumerate(alignments[0][0])
                        if aa != '-'
                           and i <= (len(Datahub.dict_chain_residues_long[key]) - 1)
                    ]
                    index_gap = [
                        i
                        for i, aa in enumerate(alignments[0][1])
                        if aa == '-'
                           and i <= (len(Datahub.dict_chain_residues_long[key]) - 1)

                    ]  # Get the index of the query alignment where there are gaps

                    def insert_nas(list, indexes, value):
                        for i in sorted(indexes):
                            list.insert(i, value)

                    def total_energy():
                        insert_nas(Datahub.dict_pyrosetta_total_energy[key], index_gap,
                                   np.nan)  # Insert NA's where there are gaps in the alignment query
                        print('Total energy non aligned:')
                        print(Datahub.dict_pyrosetta_total_energy[key])
                        Datahub.dict_pyrosetta_total_energy_aligned[key] = [Datahub.dict_pyrosetta_total_energy[key][i]
                                                                            for i in
                                                                            index_list]  # Trim the C-alpha metrics to the alignment
                        print('Total energy aligned:')
                        print(Datahub.dict_pyrosetta_total_energy_aligned[key])
                        Datahub.dict_pyrosetta_total_energy_aligned[key] = Datahub.dict_pyrosetta_total_energy_aligned[
                                                                               key][: len(Datahub.canonical_seq)]

                    def vanderwaals1():
                        insert_nas(Datahub.dict_pyrosetta_vanderwaals1[key], index_gap, np.nan)
                        print('Van der Waals 1 non aligned:')
                        print(Datahub.dict_pyrosetta_vanderwaals1[key])
                        Datahub.dict_pyrosetta_vanderwaals1_aligned[key] = [Datahub.dict_pyrosetta_vanderwaals1[key][i]
                                                                            for i in
                                                                            index_list]  # Trim the C-alpha metrics to the alignment
                        print('Van der Waals 1 aligned:')
                        print(Datahub.dict_pyrosetta_vanderwaals1_aligned[key])
                        Datahub.dict_pyrosetta_vanderwaals1_aligned[key] = Datahub.dict_pyrosetta_vanderwaals1_aligned[
                                                                               key][
                                                                           : len(Datahub.canonical_seq)]

                    def vanderwaals2():
                        insert_nas(Datahub.dict_pyrosetta_vanderwaals2[key], index_gap, np.nan)
                        print('Van der Waals 2 non aligned:')
                        print(Datahub.dict_pyrosetta_vanderwaals2[key])
                        Datahub.dict_pyrosetta_vanderwaals2_aligned[key] = [Datahub.dict_pyrosetta_vanderwaals2[key][i]
                                                                            for i in
                                                                            index_list]
                        print('Van der Waals 2 aligned:')
                        print(Datahub.dict_pyrosetta_vanderwaals2_aligned[key])
                        Datahub.dict_pyrosetta_vanderwaals2_aligned[key] = Datahub.dict_pyrosetta_vanderwaals2_aligned[
                                                                               key][
                                                                           : len(Datahub.canonical_seq)]

                    def solvation():
                        insert_nas(Datahub.dict_pyrosetta_solvation[key], index_gap, np.nan)
                        print('Solvation non aligned:')
                        print(Datahub.dict_pyrosetta_solvation[key])
                        Datahub.dict_pyrosetta_solvation_aligned[key] = [Datahub.dict_pyrosetta_solvation[key][i] for i
                                                                         in
                                                                         index_list]
                        print('Solvation aligned:')
                        print(Datahub.dict_pyrosetta_solvation_aligned[key])
                        Datahub.dict_pyrosetta_solvation_aligned[key] = Datahub.dict_pyrosetta_solvation_aligned[key][
                                                                        : len(Datahub.canonical_seq)]

                    def electrostatics():
                        insert_nas(Datahub.dict_pyrosetta_electrostatics[key], index_gap, np.nan)
                        print('Electrostatics non aligned:')
                        print(Datahub.dict_pyrosetta_electrostatics[key])
                        Datahub.dict_pyrosetta_electrostatics_aligned[key] = [
                            Datahub.dict_pyrosetta_electrostatics[key][i] for i in
                            index_list]
                        print('Electrostatics aligned:')
                        print(Datahub.dict_pyrosetta_electrostatics_aligned[key])
                        Datahub.dict_pyrosetta_electrostatics_aligned[key] = \
                        Datahub.dict_pyrosetta_electrostatics_aligned[key][
                        : len(Datahub.canonical_seq)]

                    total_energy()
                    vanderwaals1()
                    vanderwaals2()
                    solvation()
                    electrostatics()

                    break

                except Exception as e:
                    e = str(e)
                    print(f'Error message in {key}:')
                    print(e)
                    if "Unrecognized residue" in e:
                        residue = e.split(':')[-1].strip()
                        unrecognized_residues.append(residue)
                        pyrosetta_pdb_corrector()
                        print(unrecognized_residues)
                    else:
                        continue

    pyrosetta_energy_parser()

    def pyrosetta_dataframe_generator():
        '''
        Generates a dataframe with the pyrosetta metrics of the chains
        '''
        columns = [f'pos{i + 1}:{Datahub.canonical_seq[i]}' for i in
                   range(len(Datahub.canonical_seq))]  # Generate column names

        def total_energy():
            Datahub.dataframe_pyrosetta_total_energy = pd.DataFrame(columns=columns)  # Initialize dataframe
            for key, value in Datahub.dict_pyrosetta_total_energy_aligned.items():  # Iterate over the dictionary of calpha metrics
                value.extend([np.nan] * (len(Datahub.canonical_seq) - len(
                    value)))  # Extend the list of calpha metrics to the length of the canonical sequence with gaps
                Datahub.dataframe_pyrosetta_total_energy.loc[
                    key] = value  # Add the list of calpha metrics to the dataframe in order
            Datahub.dataframe_pyrosetta_total_energy = Datahub.dataframe_pyrosetta_total_energy.sort_index()
            print(Datahub.dataframe_pyrosetta_total_energy)
            Datahub.dataframe_pyrosetta_total_energy.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'pyrosetta_total_energy_dataframe.csv'))

        def vanderwaals1():
            Datahub.dataframe_pyrosetta_vanderwaals1 = pd.DataFrame(columns=columns)
            for key, value in Datahub.dict_pyrosetta_vanderwaals1_aligned.items():
                value.extend([np.nan] * (len(Datahub.canonical_seq) - len(value)))
                Datahub.dataframe_pyrosetta_vanderwaals1.loc[key] = value
            Datahub.dataframe_pyrosetta_vanderwaals1 = Datahub.dataframe_pyrosetta_vanderwaals1.sort_index()
            print(Datahub.dataframe_pyrosetta_vanderwaals1)
            Datahub.dataframe_pyrosetta_vanderwaals1.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'pyrosetta_vanderwaals1_dataframe.csv'))

        def vanderwaals2():
            Datahub.dataframe_pyrosetta_vanderwaals2 = pd.DataFrame(columns=columns)
            for key, value in Datahub.dict_pyrosetta_vanderwaals2_aligned.items():
                value.extend([np.nan] * (len(Datahub.canonical_seq) - len(value)))
                Datahub.dataframe_pyrosetta_vanderwaals2.loc[key] = value
            Datahub.dataframe_pyrosetta_vanderwaals2 = Datahub.dataframe_pyrosetta_vanderwaals2.sort_index()
            print(Datahub.dataframe_pyrosetta_vanderwaals2)
            Datahub.dataframe_pyrosetta_vanderwaals2.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'pyrosetta_vanderwaals2_dataframe.csv'))

        def solvation():
            Datahub.dataframe_pyrosetta_solvation = pd.DataFrame(columns=columns)
            for key, value in Datahub.dict_pyrosetta_solvation_aligned.items():
                value.extend([np.nan] * (len(Datahub.canonical_seq) - len(value)))
                Datahub.dataframe_pyrosetta_solvation.loc[key] = value
            Datahub.dataframe_pyrosetta_solvation = Datahub.dataframe_pyrosetta_solvation.sort_index()
            print(Datahub.dataframe_pyrosetta_solvation)
            Datahub.dataframe_pyrosetta_solvation.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'pyrosetta_solvation_dataframe.csv'))

        def electrosatics():
            Datahub.dataframe_pyrosetta_electrostatics = pd.DataFrame(columns=columns)
            for key, value in Datahub.dict_pyrosetta_electrostatics_aligned.items():
                value.extend([np.nan] * (len(Datahub.canonical_seq) - len(value)))
                Datahub.dataframe_pyrosetta_electrostatics.loc[key] = value
            Datahub.dataframe_pyrosetta_electrostatics = Datahub.dataframe_pyrosetta_electrostatics.sort_index()
            print(Datahub.dataframe_pyrosetta_electrostatics)
            Datahub.dataframe_pyrosetta_electrostatics.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'pyrosetta_electrostatics_dataframe.csv'))

        total_energy()
        vanderwaals1()
        vanderwaals2()
        solvation()
        electrosatics()

    pyrosetta_dataframe_generator()





###################### PYROSETTA #######################

###################### RMSD #######################

def calculate_rmsd():
    '''
    Calculates the RMSD values between all possible pairs of pdb files
    :return:
    '''

    # Function to get atom coordinates from a PDB file
    def get_atoms(pdb_file, positions):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        atom_list = []

        for model in structure:
            for chain in model:
                for i, residue in enumerate(chain):
                    if i + 1 in positions:  # +1 since positions is 1-indexed
                        for atom in residue:
                            if atom.get_name() == 'CA':  # only consider alpha-carbon atoms
                                atom_list.append(atom)

        return atom_list

    pdb_list_long = [pdb for pdb in os.listdir(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files'))
                     if pdb.split('_')[1].split('.')[0] in Datahub.dict_chain_residues_long.keys()] # Get the list of pdb files that are in the dictionary of query-selected pdb files
    pdb_pairs = list(combinations_with_replacement(pdb_list_long, 2))  # Generate all possible pairs of pdb files
    print(pdb_pairs)

    super = Superimposer()
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = Datahub.gap_open_penalty  # Increase the gap open penalty
    aligner.extend_gap_score = Datahub.gap_extend_penalty  # Increase the gap extend penalty


    for pair in pdb_pairs:
        target_seq = pair[0].split('_')[1].split('.')[0]
        query_seq = pair[1].split('_')[1].split('.')[0]
        alignment = aligner.align(Datahub.dict_chain_residues_long[target_seq], Datahub.dict_chain_residues_long[query_seq])
        print(pair)
        print(alignment[0])
        # print(alignment[0][0])
        # print(alignment[0][1])
        alignment_data = {target_seq: alignment[0][0], query_seq: alignment[0][1]}
        print(f'alignment dictionary: {alignment_data}')
        # Initialize the mappings of residues to positions


        def get_coincidences_no_gaps(target, query):
            '''
            Calculates the positions of the residues in the alignment without gaps, so they take into account
            their real position in their sequence rather than in the alignment.
            Otherwise, gaps distort their positions.
            :param target:
            :param query:
            :return:
            '''
            position_target = []
            position_query = []
            pos_target_counter = 0
            pos_query_counter = 0

            for i in range(max(len(target), len(query))):
                if target[i] != '-':
                    pos_target_counter += 1
                if query[i] != '-':
                    pos_query_counter += 1
                if target[i] == query[i] and target[i] != '-':
                    position_target.append(pos_target_counter)
                    position_query.append(pos_query_counter)

            return position_target, position_query

        residue_positions_target_seq, residue_positions_query_seq = get_coincidences_no_gaps(alignment_data[target_seq],
                                                                                             alignment_data[query_seq])


        # print(residue_positions_target_seq)
        # print(residue_positions_query_seq)

        # print(residue_positions_target_seq)
        # print(residue_positions_query_seq)
        # print(mutual_positions)
        # print(alignment[0])
        # print(alignment[0][0])
        # print(alignment[0][1])

        # Get atom coordinates for both proteins
        atoms_target_seq = get_atoms(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', pair[0]), residue_positions_target_seq)
        atoms_query_seq = get_atoms(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', pair[1]), residue_positions_query_seq)


        # if len(atoms_target_seq) < len(atoms_query_seq):         # There are some pdbs with small deficiencies and discontinuities. This removes mostly the last alfa carbon for complete alignment
        #     atoms_query_seq = atoms_query_seq[:len(atoms_target_seq)]
        # elif len(atoms_query_seq) < len(atoms_target_seq):
        #     atoms_target_seq = atoms_target_seq[:len(atoms_query_seq)]


        print(f'residues target: {atoms_target_seq}')
        print(f'total len: {len(atoms_target_seq)}')
        # ca_list_target_seq = [atom for atom in atoms_target_seq if atom.get_name() == 'CA']
        # print(f'ca list len: {len(ca_list_target_seq)}')
        print("Residue positions for target_seq: ", residue_positions_target_seq)
        print("Len of residue positions for target_seq: ", len(residue_positions_target_seq))
        print(f'residues query: {atoms_query_seq}')
        print(f'total len: {len(atoms_query_seq)}')
        # ca_list_query_seq = [atom for atom in atoms_query_seq if atom.get_name() == 'CA']
        # print(f'ca list len: {len(ca_list_query_seq)}')
        print("Residue positions for query_seq: ", residue_positions_query_seq)
        print("Len of residue positions for query_seq: ", len(residue_positions_query_seq))



        try:
            if pair[0] == pair[1]:
                Datahub.dict_rmsd_values[pair] = 0
                print(f'rmsd: {0}')
            else:
                super.set_atoms(atoms_target_seq, atoms_query_seq)
                rmsd = super.rms
                Datahub.dict_rmsd_values[pair] = rmsd
                print(f'rmsd: {rmsd}')
        except Exception as e:
            Datahub.dict_rmsd_values[pair] = np.nan
            print(e)
            continue

        with open(os.path.join('./temp_dir/pickles', 'rmsd_dictionary.pickle'), "wb") as f:
            pickle.dump(Datahub.dict_rmsd_values, f)
            print('Dumping RMSD dictionary record to pickle file...')
            print('RMSD dictionary record saved to pickle file!')



def calculate_pdb_coordinates():
    '''
    Calculates the coordinates of the alpha-carbon atoms in the PDB files.
    :return:
    '''

    # Function to get atom coordinates from a PDB file
    def get_atoms(pdb_file, positions):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        atom_list = []

        for model in structure:
            for chain in model:
                for i, residue in enumerate(chain):
                    if i + 1 in positions:  # +1 since positions is 1-indexed
                        for atom in residue:
                            if atom.get_name() == 'CA':  # only consider alpha-carbon atoms
                                atom_list.append(atom.get_coord())

        return atom_list

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = Datahub.gap_open_penalty  # Increase the gap open penalty
    aligner.extend_gap_score = Datahub.gap_extend_penalty  # Increase the gap extend penalty


    for key, value in Datahub.dict_chain_residues_long.items():
        target_seq = Datahub.canonical_seq
        query_seq = Datahub.dict_chain_residues_long[key]
        alignment = aligner.align(target_seq, query_seq)
        print(key)
        print(alignment[0])
        # print(alignment[0][0])
        # print(alignment[0][1])
        alignment_data = {'Canonical': alignment[0][0], 'Query': alignment[0][1]}
        print(f'alignment dictionary: {alignment_data}')
        # Initialize the mappings of residues to positions

        def get_coincidences_no_gaps(target, query):
            '''
            Calculates the positions of the residues in the alignment without gaps, so they take into account
            their real position in their sequence rather than in the alignment.
            Otherwise, gaps distort their positions.
            :param target:
            :param query:
            :return:
            '''
            position_target = []
            position_query = []
            pos_target_counter = 0
            pos_query_counter = 0

            for i in range(max(len(target), len(query))):
                if target[i] != '-':
                    pos_target_counter += 1
                if query[i] != '-':
                    pos_query_counter += 1
                if target[i] == query[i] and target[i] != '-':
                    position_target.append(pos_target_counter)
                    position_query.append(pos_query_counter)

            return position_target, position_query

        residue_positions_target_seq, residue_positions_query_seq = get_coincidences_no_gaps(alignment_data['Canonical'],
                                                                                             alignment_data['Query'])


        # print(residue_positions_target_seq)
        # print(residue_positions_query_seq)

        # print(residue_positions_target_seq)
        # print(residue_positions_query_seq)
        # print(mutual_positions)
        # print(alignment[0])
        # print(alignment[0][0])
        # print(alignment[0][1])

        # Get atom coordinates for both proteins
        # atoms_target_seq = get_atoms(f'./temp_dir/wcn_full_pdb_files/{key}', residue_positions_target_seq)
        atoms_query_seq = get_atoms(os.path.join(os.path.dirname(__file__), 'temp_dir', 'full_pdb_files', f'c_{key}.pdb'), residue_positions_query_seq)
        Datahub.dict_pdb_coord_alpha[key] = atoms_query_seq


        with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'calpha_coord_dictionary.pickle'), "wb") as f:
            pickle.dump(Datahub.dict_pdb_coord_alpha, f)
            print('Dumping C-alpha coordinates dictionary record to pickle file...')
            print('C-alpha coordinates dictionary record saved to pickle file!')



        # if len(atoms_target_seq) < len(atoms_query_seq):         # There are some pdbs with small deficiencies and discontinuities. This removes mostly the last alfa carbon for complete alignment
        #     atoms_query_seq = atoms_query_seq[:len(atoms_target_seq)]
        # elif len(atoms_query_seq) < len(atoms_target_seq):
        #     atoms_target_seq = atoms_target_seq[:len(atoms_query_seq)]


        # print(f'residues target: {atoms_target_seq}')
        # print(f'total len: {len(atoms_target_seq)}')
        # ca_list_target_seq = [atom for atom in atoms_target_seq if atom.get_name() == 'CA']
        # print(f'ca list len: {len(ca_list_target_seq)}')
        # print("Residue positions for target_seq: ", residue_positions_target_seq)
        # print("Len of residue positions for target_seq: ", len(residue_positions_target_seq))
        print(f'residues query: {atoms_query_seq}')
        print(f'total len: {len(atoms_query_seq)}')
        # ca_list_query_seq = [atom for atom in atoms_query_seq if atom.get_name() == 'CA']
        # print(f'ca list len: {len(ca_list_query_seq)}')
        print("Residue positions for query_seq: ", residue_positions_query_seq)
        print("Len of residue positions for query_seq: ", len(residue_positions_query_seq))



###################### RMSD #######################

if __name__ == '__main__':
    # DATA AND PREPROCESSING
    Datahub.directory_generator()
    Datahub.warning_silencer()
    get_chain_ids(switch_blaster())
    get_pdb_files()
    chain_residue_parser()
    metadata_generator()
    sequence_adjuster()
    dataframe_generator()
    # dataframe_trimmer()
    # dataframe_parser()

    # FEATURES
    wcn_function()
    propka_function()
    pyrosetta_function()

    # # METRICS AND COORDINATES
    calculate_rmsd()
    calculate_pdb_coordinates()

    # SCRIPT CALL
    parser = argparse.ArgumentParser(description="Script for Datahub configuration")

    # Required arguments
    parser.add_argument("--uniprot_code", type=str, required=True, help="UniProt accession code")
    parser.add_argument("--active_ligands", type=str, required=True, help="Comma-separated list of active ligands")
    parser.add_argument("--inactive_ligands", type=str, required=True, help="Comma-separated list of inactive ligands")

    # Optional arguments with default values
    parser.add_argument("--query_coverage_threshold", type=int, default=50, help="Query coverage threshold (default: 50)")
    parser.add_argument("--identity_threshold", type=int, default=30, help="Identity threshold (default: 30)")
    parser.add_argument("--gap_open_penalty", type=float, default=-0.2, help="Gap open penalty (default: -0.2)")
    parser.add_argument("--gap_extend_penalty", type=float, default=-0.2, help="Gap extend penalty (default: -0.2)")

    args = parser.parse_args()

    Datahub.uniprot_accession_code = args.uniprot_code
    Datahub.list_active_ligands = args.active_ligands.split(',')
    Datahub.list_inactive_ligands = args.inactive_ligands.split(',')
    Datahub.query_coverage_threshold = args.query_coverage_threshold
    Datahub.identity_threshold = args.identity_threshold
    Datahub.gap_open_penalty = args.gap_open_penalty
    Datahub.gap_extend_penalty = args.gap_extend_penalty

    print(f'Finished Analysing {len(Datahub.dict_chain_residues_long.keys())} Structures for Selected Features')
