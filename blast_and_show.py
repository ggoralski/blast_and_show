"""
BLAST and Show - Sequence Alignment Visualization Tool

This module provides the main entry point for performing BLAST searches between query and subject
sequences and generating visualizations of the results. It supports multiple input formats including
GenBank files, FASTA files, and GenBank accession numbers.

The program performs the following steps:
1. Parse command line arguments for input sequences and parameters
2. Process input sequences (download from GenBank if needed)
3. Run local BLAST search
4. Generate visualizations in multiple formats (SVG, PDF, PNG)

Author: [Author Name]
Date: [Date]
Version: 1.0
"""

from argparse import ArgumentParser
import dataimport as di
import blasting as bl
import draw as dr
import os
import shutil


def parse_arguments():
    """
    Parse command line arguments for the BLAST and visualization pipeline.
    
    This function sets up the argument parser with all necessary options for input sequences,
    output settings, and BLAST parameters. It supports flexible input methods including
    GenBank accession numbers, GenBank files, and FASTA files.
    
    Returns:
        tuple: A tuple containing all parsed arguments in the following order:
            - query_acc (str): Query accession number
            - subject_acc (str): Subject accession number  
            - query_fasta_file (str): Path to query FASTA file
            - query_gb_file (str): Path to query GenBank file
            - subject_fasta_file (str): Path to subject FASTA file
            - subject_gb_file (str): Path to subject GenBank file
            - title (str): Title for graphics
            - min_len (int): Minimum alignment length to display
            - email (str): Email address for NCBI queries
            - output_dir (str): Output directory path
            - e_value (float): E-value threshold for BLAST
    
    Note:
        Either accession numbers OR files must be provided for both query and subject.
        Email is required when using accession numbers.
    """
    parser = ArgumentParser(
        usage='Usage:\npython blast_and_show.py [options]',
        description='Performs BLAST search between query and subject sequences and generates visualizations.\n'
                   'Supports GenBank files, FASTA files, and GenBank accession numbers as input.\n'
                   'Outputs include BLAST results and sequence alignment visualizations in multiple formats.'
    )

    # Query sequence options
    parser.add_argument("-q", "--query", dest="query", default="",
                        help="Query accession number in GenBank. Required if query file is not provided.")
    parser.add_argument("-a", "--query_gb_file", dest="query_gb_file", default="",
                        help="GenBank format (gb) file for query sequence. Required if query accession is not provided.")
    parser.add_argument("-b", "--query_fasta_file", dest="query_fasta_file", default="",
                        help="FASTA file with query sequence(s). Can be used with query accession or gb file.")

    # Subject sequence options
    parser.add_argument("-s", "--subject", dest="subject", default="",
                        help="Subject accession number in GenBank. Required if subject file is not provided.")
    parser.add_argument("-c", "--subject_gb_file", dest="subject_gb_file", default="",
                        help="GenBank format (gb) file for subject sequence. Required if subject accession is not provided.")
    parser.add_argument("-d", "--subject_fasta_file", dest="subject_fasta_file", default="",
                        help="FASTA file with subject sequence. Can be used with subject accession or gb file.")

    # Output and visualization options
    parser.add_argument("-t", "--title", dest="title", default="Results of blastn",
                        help="Title of graphics, also used (without spaces) as the name of graphics files.")
    parser.add_argument("-o", "--output_dir", dest="output_dir", default='output',
                        help="Output directory. Default: output")
    
    # BLAST parameters
    parser.add_argument("-v", "--e_value", dest="e_value", default=1e-10,
                        help="E-value threshold for BLAST. Default: 1e-10")
    parser.add_argument("-m", "--minlen", dest="minlen", default=100,
                        help="Minimum length (bp) of alignments to draw. Default: 100")
    
    # NCBI requirements
    parser.add_argument("-e", "--email", dest="email", default='',
                        help="Your email address (required by NCBI when using accession numbers).")

    args = parser.parse_args()
    
    # Extract all arguments
    query_acc = args.query
    subject_acc = args.subject
    query_fasta_file = args.query_fasta_file
    query_gb_file = args.query_gb_file
    subject_fasta_file = args.subject_fasta_file
    subject_gb_file = args.subject_gb_file
    title = args.title
    min_len = int(args.minlen)
    email = args.email
    output_dir = args.output_dir
    e_value = args.e_value

    return query_acc, subject_acc, query_fasta_file, query_gb_file, subject_fasta_file, subject_gb_file, title, min_len, email, output_dir, e_value


def main():
    """
    Main function that orchestrates the entire BLAST and visualization pipeline.
    
    This function:
    1. Parses command line arguments
    2. Creates output directory if needed
    3. Processes input sequences (downloads from GenBank if necessary)
    4. Copies settings file to output directory
    5. Runs BLAST search
    6. Generates visualizations
    
    The function handles all the coordination between different modules and ensures
    proper file management and directory structure.
    
    Raises:
        SystemExit: If required arguments are missing or invalid
        OSError: If output directory cannot be created
        Various exceptions from called modules for file I/O, network, or BLAST errors
    """
    # Parse command line arguments
    query_acc, subject_acc, query_fasta_file, query_gb_file, \
    subject_fasta_file, subject_gb_file, title, min_len, email, \
    output_dir, e_value = parse_arguments()
    
    # Display parsed options for user verification
    print(f'''Options:
        {query_acc = }
        {subject_acc = }
        {query_fasta_file = }
        {query_gb_file = }
        {subject_fasta_file = }
        {subject_gb_file = }
        {title = }
        {min_len = }
        {email = }
        {output_dir = }
        {e_value = }
    ''')
    
    # Create output directory if it doesn't exist
    if os.path.exists(output_dir):
        print(f'Output directory {output_dir} exists. Using existing directory.')
    else:
        print(f'Creating output directory: {output_dir}')
        os.makedirs(output_dir)

    # Process input sequences - download from GenBank if needed, convert formats
    print("Processing query sequence...")
    query_gb_file, query_fasta_file, query_json_file = di.process_input_sequences(
        'query', query_acc, query_gb_file, query_fasta_file, output_dir, email
    )
    
    print("Processing subject sequence...")
    subject_gb_file, subject_fasta_file, subject_json_file = di.process_input_sequences(
        'subject', subject_acc, subject_gb_file, subject_fasta_file, output_dir, email
    )
    
    # Copy settings file to output directory for visualization configuration
    shutil.copy("settings.json", output_dir)
    
    # Change to output directory for BLAST operations
    main_dir = os.getcwd()
    os.chdir(output_dir)
    
    # Run BLAST search
    print("Running BLAST search...")
    xml_file, tsv_file = bl.run_local_blast(query_fasta_file, subject_fasta_file, e_value)
    
    # Generate visualizations
    print("Generating visualizations...")
    dr.run(os.path.basename(query_json_file), os.path.basename(subject_json_file), tsv_file)
    
    print("Pipeline completed successfully!")


if __name__ == "__main__":
    main()