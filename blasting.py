import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import fileinput
import pandas as pd


def _run_cmd(cmd):
    """Helper: uruchamia komendę jako lista argumentów i zwraca (stdout, stderr).
    Rzuca subprocess.CalledProcessError, jeśli proces zwróci kod != 0."""
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout, result.stderr


def make_blast_db(subject_file):
    """Creating database for subject sequence"""
    db_name = 'blast_db'
    _run_cmd([
        "makeblastdb",
        "-in", subject_file,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", db_name,
    ])
    print(f"BLAST database {db_name} for {subject_file} created")
    return db_name


def _run_blastn(query_file, db, evalue, out, outfmt=None):
    """Wywołanie programu blastn z linii poleceń (zamiennik NcbiblastnCommandline)."""
    cmd = [
        "blastn",
        "-query", query_file,
        "-db", db,
        "-evalue", str(evalue),
        "-out", out,
    ]
    if outfmt is not None:
        cmd.extend(["-outfmt", str(outfmt)])
    return _run_cmd(cmd)


def run_local_blast(query_file, subject_file, e_value=1e-10):
    """
    Executing BLAST locally and saving output files in XML and TSV formats

    Parameters:
    query_file - fasta query file
    subject_file - fasta subject file
    e_value - E-value (default 1e-10)
    """
    # Creating database for subject sequence
    blast_db = make_blast_db(subject_file)

    # Executing blastn and saving results in xml file
    print(f"Executing blastn using {blast_db} database and saving results in xml file...")

    # Text format (domyślny pairwise)
    txt_file = "blast_results.tsv"
    _run_blastn(query_file, blast_db, e_value, txt_file)

    # XML format
    xml_file = "blast_results.xml"
    _run_blastn(query_file, blast_db, e_value, xml_file, outfmt=5)

    # TSV format
    tsv_file = "blast_results.tsv"
    tsv_format = '6 length qlen slen qcovs qcovhsp pident mismatch gapopen evalue bitscore ' + \
                 'gaps qseqid sseqid qstart qend sstart send qseq sseq'
    _run_blastn(query_file, blast_db, e_value, tsv_file, outfmt=tsv_format)

    headers = 'length\tqlen\tslen\tqcovs\tqcovhsp\tpident\tmismatch\tgapopen\tevalue\tbitscore' + \
              '\tgaps\tqseqid\tsseqid\tqstart\tqend\tsstart\tsend\tqseq\tsseq'
    with fileinput.input(tsv_file, inplace=True) as file:
        is_first = True
        for line in file:
            if is_first:
                print(headers)
                is_first = False
            print(line.rstrip())

    print(f'BLAST results saved in files: {txt_file}, {xml_file}, {tsv_file}')
    save_to_fasta(query_file, subject_file, tsv_file)
    # Deleting database files
    files_to_remove = glob.glob("blast_db*")
    for file in files_to_remove:
        os.remove(file)

    return xml_file, tsv_file


def create_fasta_file(output_fasta, data, stype):
    if stype == 'query':
        name = 'qseq_name'
        seq = 'qseq'
    elif stype == 'subject':
        name = 'sseq_name'
        seq = 'sseq'
    else:
        print(f'ERROR: bad value for type of sequence (only "query" or "subject" allowed)!')
        sys.exit(1)
    with open(output_fasta, 'w') as outfile:
        for index, row in data.iterrows():
            outfile.write(f'{row[name].strip()}\n')
            outfile.write(f'{row[seq].strip()}\n')
    print(f'Fasta file for aligned {stype} sequence(s) created: {output_fasta}.')


def save_to_fasta(query_file, subject_file, tsv_file):
    data = pd.read_table(tsv_file)
    data['sseq_name'] = '>' + data['sseqid'].astype(str) + '_' + data['sstart'].astype(str) + '-' + data['send'].astype(str)
    data['qseq_name'] = '>' + data['qseqid'].astype(str) + '_' + data['qstart'].astype(str) + '-' + data['qend'].astype(str)
    create_fasta_file(f'aligned_query-{query_file}', data, 'query')
    create_fasta_file(f'aligned_subject-{subject_file}', data, 'subject')
