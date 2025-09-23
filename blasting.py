"""
BLAST Operations Module

This module handles all BLAST-related operations including database creation, 
BLAST execution, and result processing. It provides functions to run local BLAST
searches and generate output files in multiple formats.

The module supports:
- Creating BLAST databases from FASTA files
- Running blastn searches with customizable parameters
- Generating output in XML, TSV, and text formats
- Extracting aligned sequences to separate FASTA files
- Automatic cleanup of temporary database files

Dependencies:
    - NCBI BLAST+ tools (makeblastdb, blastn)
    - BioPython
    - pandas

Author: [Author Name]
Date: [Date]
Version: 1.0
"""

import sys
import glob
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import fileinput
import pandas as pd


def make_blast_db(subject_file):
    """
    Create a BLAST database from a FASTA file.
    
    This function uses the makeblastdb command-line tool to create a nucleotide
    BLAST database from the provided subject sequence file. The database is
    created with sequence ID parsing enabled.
    
    Args:
        subject_file (str): Path to the FASTA file containing subject sequences
        
    Returns:
        str: Name of the created database ('blast_db')
        
    Raises:
        OSError: If makeblastdb command fails or is not found
        
    Note:
        Requires NCBI BLAST+ tools to be installed and in PATH
    """
    db_name = 'blast_db'
    # Create nucleotide database with sequence ID parsing
    os.system(f"makeblastdb -in {subject_file} -dbtype nucl -parse_seqids -out {db_name}")
    print(f"BLAST database {db_name} for {subject_file} created")
    return db_name


def run_local_blast(query_file, subject_file, e_value=1e-10):
    """
    Execute local BLAST search and generate output files in multiple formats.
    
    This function performs a complete BLAST workflow:
    1. Creates a BLAST database from the subject file
    2. Runs blastn search in three different output formats
    3. Processes TSV output to add headers
    4. Extracts aligned sequences to FASTA files
    5. Cleans up temporary database files
    
    Args:
        query_file (str): Path to FASTA file containing query sequences
        subject_file (str): Path to FASTA file containing subject sequences  
        e_value (float, optional): E-value threshold for BLAST. Defaults to 1e-10
        
    Returns:
        tuple: A tuple containing:
            - xml_file (str): Path to XML format BLAST results
            - tsv_file (str): Path to TSV format BLAST results
            
    Raises:
        OSError: If BLAST commands fail or files cannot be created
        FileNotFoundError: If input files don't exist
        
    Note:
        Creates temporary database files that are automatically cleaned up
    """
    # Create BLAST database for subject sequences
    blast_db = make_blast_db(subject_file)

    print(f"Executing blastn using {blast_db} database and saving results...")

    # Run BLAST in default text format
    txt_file = "blast_results.tsv"
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=blast_db,
        evalue=e_value,
        out=txt_file
    )
    stdout, stderr = blastn_cline()

    # Run BLAST in XML format for detailed results
    xml_file = "blast_results.xml"
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=blast_db,
        evalue=e_value,
        outfmt=5,  # XML format
        out=xml_file
    )
    stdout, stderr = blastn_cline()

    # Run BLAST in custom TSV format with specific columns
    tsv_file = "blast_results.tsv"
    tsv_format = '6 length qlen slen qcovs qcovhsp pident mismatch gapopen evalue bitscore ' +\
                 'gaps qseqid sseqid qstart qend sstart send qseq sseq'
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=blast_db,
        evalue=e_value,
        outfmt=tsv_format,
        out=tsv_file
    )
    stdout, stderr = blastn_cline()
    
    # Add column headers to TSV file for easier processing
    headers = 'length\tqlen\tslen\tqcovs\tqcovhsp\tpident\tmismatch\tgapopen\tevalue\tbitscore'+\
               '\tgaps\tqseqid\tsseqid\tqstart\tqend\tsstart\tsend\tqseq\tsseq'
    with fileinput.input(tsv_file, inplace=True) as file:
        is_first = True
        for line in file:
            if is_first:
                print(headers)  # Add header as first line
                is_first = False
            print(line.rstrip())

    print(f'BLAST results saved in files: {txt_file}, {xml_file}, {tsv_file}')
    
    # Extract aligned sequences to separate FASTA files
    save_to_fasta(query_file, subject_file, tsv_file)
    
    # Clean up temporary database files
    files_to_remove = glob.glob("blast_db*")
    for file in files_to_remove:
        os.remove(file)

    return xml_file, tsv_file


def create_fasta_file(output_fasta, data, stype):
    """
    Create a FASTA file from aligned sequence data.
    
    This function extracts either query or subject sequences from BLAST results
    and writes them to a FASTA file with descriptive headers including coordinates.
    
    Args:
        output_fasta (str): Path for the output FASTA file
        data (pandas.DataFrame): DataFrame containing BLAST results with sequence data
        stype (str): Type of sequences to extract ('query' or 'subject')
        
    Raises:
        SystemExit: If stype is not 'query' or 'subject'
        IOError: If output file cannot be written
        
    Note:
        Sequence headers include sequence ID and alignment coordinates
    """
    # Determine column names based on sequence type
    if stype == 'query':
        name_col = 'qseq_name'
        seq_col = 'qseq'
    elif stype == 'subject':
        name_col = 'sseq_name'
        seq_col = 'sseq'
    else:
        print(f'ERROR: Invalid sequence type "{stype}". Only "query" or "subject" allowed!')
        sys.exit(1)
    
    # Write sequences to FASTA file
    with open(output_fasta, 'w') as outfile:
        for index, row in data.iterrows():
            outfile.write(f'{row[name_col].strip()}\n')  # Header line
            outfile.write(f'{row[seq_col].strip()}\n')   # Sequence line
    
    print(f'FASTA file for aligned {stype} sequence(s) created: {output_fasta}')


def save_to_fasta(query_file, subject_file, tsv_file):
    """
    Extract aligned sequences from BLAST results and save to FASTA files.
    
    This function reads BLAST results in TSV format and creates separate FASTA files
    for query and subject sequences that participated in alignments. Each sequence
    header includes the original sequence ID and alignment coordinates.
    
    Args:
        query_file (str): Original query filename (used in output filename)
        subject_file (str): Original subject filename (used in output filename)
        tsv_file (str): Path to TSV file containing BLAST results
        
    Raises:
        FileNotFoundError: If TSV file doesn't exist
        pandas.errors.EmptyDataError: If TSV file is empty
        
    Note:
        Creates two output files:
        - aligned_query-{query_file}: Contains query sequences from alignments
        - aligned_subject-{subject_file}: Contains subject sequences from alignments
    """
    # Read BLAST results from TSV file
    data = pd.read_table(tsv_file)
    
    # Create descriptive sequence names with coordinates
    # Subject sequences: >sequence_id_start-end
    data['sseq_name'] = '>' + data['sseqid'].astype(str) + '_' + \
                       data['sstart'].astype(str) + '-' + data['send'].astype(str)
    
    # Query sequences: >sequence_id_start-end  
    data['qseq_name'] = '>' + data['qseqid'].astype(str) + '_' + \
                       data['qstart'].astype(str) + '-' + data['qend'].astype(str)
    
    # Create separate FASTA files for query and subject sequences
    create_fasta_file(f'aligned_query-{query_file}', data, 'query')
    create_fasta_file(f'aligned_subject-{subject_file}', data, 'subject')