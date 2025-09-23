"""
Data Import and Processing Module

This module handles all data import, processing, and conversion operations for sequence data.
It provides comprehensive functionality for:

- Downloading sequences from GenBank using accession numbers
- Converting between different file formats (GenBank, FASTA, JSON, TSV)
- Processing sequence features and annotations
- Managing input sequences from various sources
- Extracting and organizing genomic features

The module supports flexible input handling, allowing users to provide sequences as:
- GenBank accession numbers (requires email for NCBI compliance)
- Local GenBank (.gb) files
- Local FASTA (.fasta) files

Key Features:
- Automatic format detection and conversion
- Feature extraction from GenBank annotations
- JSON serialization of sequence metadata
- Comprehensive error handling for network operations
- NCBI Entrez integration with proper email handling

Dependencies:
    - BioPython (Entrez, SeqIO)
    - pandas
    - json
    - NCBI Entrez utilities

Author: [Author Name]
Date: [Date]
Version: 1.0
"""

from email.policy import default
import shutil
from Bio import Entrez, SeqIO, SearchIO
import os
import sys
import pandas as pd
import json
from Bio.Blast.Applications import NcbiblastnCommandline
from numpy.core.records import record


def download_sequence(accession, email, format, mode="text"):
    """
    Download sequences from GenBank using accession numbers.
    
    This function uses NCBI's Entrez system to download sequence data in various formats.
    It handles both text and XML modes and supports different return types.
    
    Args:
        accession (str): GenBank accession number for the sequence
        email (str): Valid email address (required by NCBI)
        format (str): Return format ('fasta', 'gb', etc.)
        mode (str, optional): Return mode ('text' or 'xml'). Defaults to 'text'
        
    Returns:
        Bio.SeqRecord or dict: Sequence record object or parsed XML data
        
    Raises:
        Exception: If download fails due to network issues, invalid accession, etc.
        
    Note:
        Email address is required by NCBI's usage policies. Invalid or missing
        email addresses may result in blocked access.
    """
    Entrez.email = email

    try:
        # Fetch sequence data from NCBI nucleotide database
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype=format,
            retmode=mode
        )

        # Parse response based on mode
        if mode == 'xml':
            record = Entrez.read(handle)
        else:
            record = SeqIO.read(handle, 'genbank')

        handle.close()
        return record

    except Exception as e:
        print(f"### Error downloading sequence {accession}: {e}")
        raise


def get_json_data(accession, email, gb_file):
    """
    Download GenBank data and save in JSON format.
    
    This function downloads sequence data from GenBank in XML format and converts
    it to JSON for easier processing. The JSON filename is derived from the
    GenBank filename for consistency.
    
    Args:
        accession (str): GenBank accession number
        email (str): Valid email address for NCBI
        gb_file (str): GenBank filename (used to generate JSON filename)
        
    Returns:
        str: Path to the created JSON file
        
    Raises:
        Exception: If download or file writing fails
        
    Note:
        The JSON file contains the complete GenBank record structure
        including all annotations and features.
    """
    # Generate JSON filename from GenBank filename
    if gb_file.endswith('.gb'):
        json_file = gb_file.replace('.gb', '.json')
    else:
        json_file = gb_file + '.json'
    
    # Download sequence data in XML format
    record = download_sequence(accession, email, "gb", "xml")
    
    # Convert to JSON and save
    json_data = json.dumps(record[0], indent=2)
    with open(json_file, "w") as f:
        f.write(json_data)

    return json_file


def gb_to_fasta(gb_file):
    """
    Convert GenBank file to FASTA format.
    
    This function reads a GenBank file and extracts the sequence data,
    saving it in FASTA format. Multiple sequences in the GenBank file
    are all converted.
    
    Args:
        gb_file (str): Path to input GenBank file
        
    Returns:
        str: Path to the created FASTA file
        
    Raises:
        FileNotFoundError: If GenBank file doesn't exist
        IOError: If output file cannot be written
        
    Note:
        Output filename is derived by replacing .gb extension with .fasta
    """
    # Generate FASTA filename
    if gb_file.endswith('.gb'):
        fasta_file = gb_file.replace('.gb', '.fasta')
    else:
        fasta_file = gb_file + '.fasta'
    
    # Convert GenBank to FASTA
    with open(fasta_file, "w") as output_handle:
        for record in SeqIO.parse(gb_file, "genbank"):
            SeqIO.write(record, output_handle, "fasta")
    
    print(f'{gb_file} converted to {fasta_file}')
    return fasta_file


def dict_from_row(row):
    """
    Convert a pandas DataFrame row to a dictionary.
    
    This utility function converts a pandas Series (DataFrame row) into
    a dictionary with column names as keys.
    
    Args:
        row (pandas.Series): DataFrame row to convert
        
    Returns:
        dict: Dictionary representation of the row
    """
    dat = {}
    cols = row.index
    for col in cols:
        dat[col] = row[col]
    return dat


def get_features_for_feature(data, row, typ):
    """
    Extract sub-features for a given genomic feature.
    
    This function finds all features of a specific type that belong to a parent
    feature (e.g., exons belonging to a gene). It uses pattern matching on
    feature names to establish relationships.
    
    Args:
        data (pandas.DataFrame): DataFrame containing all features
        row (pandas.Series): Parent feature row
        typ (str): Type of sub-features to find ('exon', 'intron', etc.)
        
    Returns:
        dict: Dictionary of sub-features indexed by their DataFrame index
        
    Note:
        Uses pattern matching with '[' character to find related features
    """
    # Find features of specified type that contain the parent feature name
    f = data[(data['type'] == typ) & (data['name'].str.contains(row['name'] + "[", regex=False))]
    subdat = {}
    
    if f.shape[0] > 0:
        nr = 0
        for index, feature_row in f.iterrows():
            nr += 1
            feature_dict = dict_from_row(feature_row)
            subdat[index] = feature_dict

    return subdat


def find_next_gene(df, i):
    """
    Find the index of the next gene feature in a DataFrame.
    
    This function searches forward from a given position to find the next
    feature of type 'gene'. Used for grouping features that belong to
    the same gene.
    
    Args:
        df (pandas.DataFrame): DataFrame containing genomic features
        i (int): Starting index to search from
        
    Returns:
        int: Index of the next gene feature, or length of DataFrame if none found
    """
    next_gene_index = len(df)
    for j in range(i + 1, len(df)):
        if df.iloc[j]['type'] == 'gene':
            next_gene_index = j
            break
    return next_gene_index


def make_json_data(data: pd.DataFrame, gb_file, seq_description):
    """
    Create structured JSON data from genomic features.
    
    This function processes a DataFrame of genomic features and creates a hierarchical
    JSON structure that groups features by genes and organizes sub-features (CDS, tRNA, rRNA)
    under their parent genes. It also handles misc_features separately.
    
    Args:
        data (pandas.DataFrame): DataFrame containing genomic features with columns:
                                type, name, start, end, etc.
        gb_file (str): Path to GenBank file (used to generate JSON filename)
        seq_description (dict): Dictionary containing sequence metadata
        
    Returns:
        str: Path to the created JSON file
        
    Raises:
        IOError: If JSON file cannot be written
        
    Note:
        The resulting JSON has a hierarchical structure:
        - description: sequence metadata
        - features: organized by gene with nested CDS, tRNA, rRNA features
    """
    # Generate JSON filename
    if gb_file.endswith('.gb'):
        json_file = gb_file.replace('.gb', '.json')
    else:
        json_file = gb_file + '.json'
    
    json_data = {}
    json_all_data = {}
    nr = 0
    
    # Add sequence description
    json_all_data['description'] = {
        'seq_description': seq_description
    }
    if gb_file.endswith('.gb'):
        json_file = gb_file.replace('.gb', '.json')
    else:
        json_file = gb_file + '.json'
    genes = data[data['type'] == 'gene']
    for i, row in genes.iterrows():
        nr += 1
        gene = dict_from_row(row)
        next_gene_index = find_next_gene(data, i)

        # Fix for CDS
        cdss = data[(data['type'] == 'CDS') & (data.index > i) & (data.index <= next_gene_index)]
        if cdss.shape[0] > 0:
            cds_d = {}
            cd_nr = 0
            for index2, row2 in cdss.iterrows():
                cd_nr += 1
                cds = dict_from_row(row2)
                # Add exons and introns to CDS
                cds['exons'] = get_features_for_feature(data, row, 'exon')
                cds['introns'] = get_features_for_feature(data, row, 'introns')
                cds_d[cd_nr] = cds
            gene['CDS'] = cds_d

        # Fix for tRNA
        trnas = data[(data['type'] == 'tRNA') & (data.index > i) & (data.index <= next_gene_index)]
        if trnas.shape[0] > 0:
            for index2, row2 in trnas.iterrows():
                trna = dict_from_row(row2)
                trna['tRNA'] = get_features_for_feature(data, row, 'tRNA')
                gene['tRNA'] = trna

        # Fix for rRNA
        rrnas = data[(data['type'] == 'rRNA') & (data.index > i) & (data.index <= next_gene_index)]
        if rrnas.shape[0] > 0:
            for index2, row2 in rrnas.iterrows():
                rrna = dict_from_row(row2)
                rrna['rRNA'] = get_features_for_feature(data, row, 'tRNA')
                gene['rRNA'] = rrna
                
        json_data[nr] = gene
    
    # Process misc_feature entries separately
    misc_features = data[data['type'] == 'misc_feature']
    for i, row in misc_features.iterrows():
        nr += 1
        misc_feature = dict_from_row(row)
        json_data[nr] = misc_feature

    # Save complete data structure to JSON
    json_all_data['features'] = json_data
    with open(json_file, 'w') as res_file:
        json.dump(json_all_data, res_file, indent=4)
    
    return json_file


def download_files_for_organism(acc, email, output_dir):
    """
    Download and process all files for a given organism.
    
    This function downloads a sequence from GenBank, extracts organism information,
    and creates all necessary files (GenBank, FASTA, JSON, TSV) with consistent naming.
    
    Args:
        acc (str): GenBank accession number
        email (str): Valid email address for NCBI
        output_dir (str): Directory to save output files
        
    Returns:
        tuple: A tuple containing:
            - gb_file (str): Path to GenBank file
            - fasta_file (str): Path to FASTA file
            - json_file (str): Path to JSON file
            - seq_description (dict): Sequence metadata
            
    Raises:
        Exception: If download or file processing fails
        
    Note:
        Files are named using organism name and accession number for clarity
    """
    # Download sequence record
    record = download_sequence(acc, email, 'gb')

    # Extract organism name for file naming
    organism = record.annotations.get('organism', 'Organism name not found')
    print(f'Getting data for organism: {organism}')
    
    if organism == 'Organism name not found':
        output_file_name = f'{acc}'
    else:
        output_file_name = f'{organism.replace(" ", "_")}-{acc}'
    # Save in sequence in gb format
    gb_file = os.path.join(f'{output_dir}',f'{output_file_name}.gb')
    SeqIO.write(record, gb_file, "genbank")
    print(f'\tSaved {gb_file} file.')
    
    # Convert to FASTA
    fasta_file = gb_to_fasta(gb_file)
    
    # Create JSON file
    json_file = get_json_data(acc, email, gb_file)
    
    # Extract features to TSV and create structured JSON
    seq_description, json_file = export_features_from_gb_to_tsv(gb_file)
    
    return gb_file, fasta_file, json_file, seq_description


def download_all_sequences(query_acc, subject_acc, email):
    """
    Download both query and subject sequences.
    
    Convenience function to download both query and subject sequences
    using their accession numbers.
    
    Args:
        query_acc (str): Query sequence accession number
        subject_acc (str): Subject sequence accession number
        email (str): Valid email address for NCBI
        
    Returns:
        tuple: Paths to query and subject GenBank and FASTA files
    """
    query_gb_file, query_fasta_file, json_file, seq_description = download_files_for_organism(query_acc, email)
    subject_gb_file, subject_fasta_file, json_file, seq_description = download_files_for_organism(subject_acc, email)
    return query_gb_file, query_fasta_file, subject_gb_file, subject_fasta_file


def copy_file(output_dir, in_file):
    """
    Copy a file to the output directory.
    
    Args:
        output_dir (str): Destination directory
        in_file (str): Source file path
        
    Returns:
        str: Path to the copied file
        
    Raises:
        IOError: If file cannot be copied
    """
    copy_file_path = os.path.join(output_dir, in_file)
    shutil.copy(in_file, copy_file_path)
    print(f'File {in_file} copied to {copy_file_path}.')
    return copy_file_path


def use_gb_file(output_dir, gb_file):
    """
    Process an existing GenBank file.
    
    This function copies a GenBank file to the output directory and creates
    all derived files (FASTA, JSON, TSV).
    
    Args:
        output_dir (str): Output directory
        gb_file (str): Path to GenBank file
        
    Returns:
        tuple: Paths to copied GenBank file, FASTA file, sequence description, and JSON file
    """
    copy_gb_file = copy_file(output_dir, gb_file)
    fasta_file = gb_to_fasta(copy_gb_file)
    seq_description, json_file = export_features_from_gb_to_tsv(copy_gb_file)
    return copy_gb_file, fasta_file, seq_description, json_file


def process_input_sequences(type_seq, acc_acc, acc_gb_file, acc_fasta_file, output_dir, email):
    """
    Process input sequences from various sources.
    
    This is the main function for handling different types of sequence input.
    It can process sequences from:
    - GenBank accession numbers (downloads from NCBI)
    - Local GenBank files
    - Local FASTA files
    
    Args:
        type_seq (str): Type of sequence ('query' or 'subject')
        acc_acc (str): GenBank accession number (if provided)
        acc_gb_file (str): Path to GenBank file (if provided)
        acc_fasta_file (str): Path to FASTA file (if provided)
        output_dir (str): Output directory for processed files
        email (str): Email address for NCBI (required for accession numbers)
        
    Returns:
        tuple: A tuple containing:
            - gb_file (str): Path to GenBank file
            - fasta_file (str): Path to FASTA file
            - json_file (str): Path to JSON file
            
    Raises:
        SystemExit: If no valid input is provided or email is missing when required
        
    Note:
        Exactly one input method must be provided per sequence
    """
    print(f'\n----------------------------\nProcessing {type_seq}')
    
    if acc_acc == '':
        print(f'No accession number for {type_seq} sequence provided. Checking for files.')

        if acc_gb_file != '':
            if os.path.isfile(acc_gb_file):
                print(f'{type_seq} gb file {acc_gb_file} for {type_seq} sequence is provided. I will try to use it.')
                acc_gb_file, acc_fasta_file, seq_description, json_file = use_gb_file(output_dir, acc_gb_file)
            else:
                print(f'File {acc_gb_file} cannot be found.')
                sys.exit(1)
        else:
            print(f'No gb file for {type_seq} sequence is provided.')
            if acc_fasta_file != '':
                if os.path.isfile(acc_fasta_file):
                    print(f'{type_seq} fasta file {acc_fasta_file} for {type_seq} sequence is provided. I will try to use it, ' +
                          f'but no sequence features for {type_seq} will be shown.')
                    acc_fasta_file_copy = copy_file(output_dir, acc_fasta_file)

                    """
                    with open(acc_fasta_file_copy, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                description = line.replace('>','').rstrip()
                                """
                    seq_descr = {}
                    for record in SeqIO.parse(acc_fasta_file_copy, "fasta"):
                        print("Description of sequence:", record.description)
                        print("Length of sequence:", len(record.seq))
                        seq_descr['seq_len'] = len(record.seq)
                        seq_descr['description'] = record.description
                        seq_descr['name'] = record.description
                    json_all_data = {}


                    json_all_data['description'] = {
                        'seq_description': seq_descr
                    }
                    json_all_data['features'] = {}
                    json_file = os.path.join(output_dir, f'{type_seq}.json')
                    with open(json_file, 'w') as res_file:
                        json.dump(json_all_data, res_file, indent=4)
                else:
                    print(f'File {acc_fasta_file} cannot be found.')
                    sys.exit(1)
            else:
                print(f'ERROR: No accession number, gb file, or fasta file was provided for {type_seq} sequence. At least one of them is required.')
                sys.exit(1)
    else:
        check_email(email, type_seq)
        print(f'Accession number {acc_acc} for {type_seq} sequence was provided. I will try to download data from GenBank.')
        acc_gb_file, acc_fasta_file, json_file, seq_description = download_files_for_organism(acc_acc, email, output_dir)
        print(f'{type_seq} files: {acc_gb_file}, {acc_fasta_file}, {json_file}')

    return os.path.basename(acc_gb_file), os.path.basename(acc_fasta_file), json_file
def check_email(email, type_seq):
    """
    Validate email address for NCBI compliance.
    
    Args:
        email (str): Email address to validate
        type_seq (str): Type of sequence being processed (for error messages)
        
    Raises:
        SystemExit: If email is missing or invalid
        
    Note:
        NCBI requires a valid email address for all Entrez queries
    """
    if email == '':
        print(f'ERROR: Email address is required when using accession numbers.')
        print(f'Please provide your email address with the -e option.')
        print(f'This is required by NCBI for {type_seq} sequence download.')
        sys.exit(1)


def export_features_from_gb_to_tsv(gb_file):
    """
    Extract genomic features from GenBank file and export to TSV and JSON formats.
    
    This function parses a GenBank file, extracts all genomic features and annotations,
    and creates both a TSV file for tabular analysis and a structured JSON file for
    visualization purposes.
    
    Args:
        gb_file (str): Path to input GenBank file
        
    Returns:
        tuple: A tuple containing:
            - tsv_file (str): Path to created TSV file
            - json_file (str): Path to created JSON file
            
    Raises:
        FileNotFoundError: If GenBank file doesn't exist
        IOError: If output files cannot be written
        
    Note:
        Processes features: exon, intron, gene, CDS, tRNA, rRNA, misc_feature
        Creates hierarchical JSON structure for visualization
    """
    # Initialize DataFrame for feature data
    data = pd.DataFrame(columns=[
        'nr',
        'type',
        'name',
        'start',
        'end',
        'note',
        'qualifier',
        'product',
        'qualifier_note',
        'plastid-derived',
    ])
    
    seq_description = {}
    
    # Parse GenBank file
    for record in SeqIO.parse(gb_file, "genbank"):
        nr = 0
        
        # Extract sequence metadata
        seq_description = {
            'seq_len': len(record.seq),
            'name': record.name,
            'id': record.id,
            'description': record.description,
            'source': record.annotations['source'],
            'organism': record.annotations['organism'],
            'taxonomy': record.annotations['taxonomy'],
        }

        # Process each feature in the record
        for feature in record.features:
            nr += 1

            # Initialize feature record
            rec = {
                'nr': nr,
                'type': '',
                'name': '',
                'start': '',
                'end': '',
                'note': '',
                'qualifier': '',
                'product': '',
                'qualifier_note': '',
                'plastid-derived': '',
                'strand': ''
            }
            
            # Process relevant feature types
            if feature.type in ['exon', 'intron', 'gene', 'CDS', 'tRNA', 'rRNA', 'misc_feature']:
                rec['type'] = feature.type
                rec['start'] = int(feature.location.start)
                rec['end'] = int(feature.location.end)
                rec['strand'] = feature.location.strand
                
                if 'note' in feature.qualifiers:
                    rec['note'] = feature.qualifiers['note'][0]

                # Check for plastid-derived misc_features
                is_plastid = (feature.type == 'misc_feature' and
                              'note' in feature.qualifiers and
                              'plastid-derived' in feature.qualifiers['note'][0])
                rec['plastid-derived'] = is_plastid
                
                # Extract feature name from various qualifier types
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0]
                    if 'number' in feature.qualifiers:
                        name += f"[{feature.qualifiers['number'][0]}]"
                elif 'product' in feature.qualifiers:
                    rec['product'] = feature.qualifiers['product'][0]
                    name = feature.qualifiers['product'][0]
                elif 'note' in feature.qualifiers:
                    name = feature.qualifiers['note'][0]  # First 20 characters from note
                    rec['note'] = feature.qualifiers['note'][0][:20]
                else:
                    name = feature.type
                    
                rec['name'] = name
                new_df = pd.DataFrame([rec])
                data = pd.concat([data, new_df])
    
    # Save features to TSV file
    tsv_file = gb_file.replace(".gb", ".tsv")
    data.reset_index(drop=True, inplace=True)
    data.to_csv(tsv_file, sep='\t', index=False)
    print(f'Features from {gb_file} exported to {tsv_file}')
    
    # Create structured JSON file
    json_file = make_json_data(data, gb_file, seq_description)

    return tsv_file, json_file