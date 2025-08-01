from email.policy import default
import shutil
from Bio import Entrez, SeqIO, SearchIO
import os
import sys
#import math
import pandas as pd
import json
from Bio.Blast.Applications import NcbiblastnCommandline
from numpy.core.records import record


def download_sequence(accession, email, format, mode="text"):
    """
    Downloads sequences from GenBank in fasta or gb format
    """
    Entrez.email = email

    try:
        # Getting the sequence
        handle = Entrez.efetch(db="nucleotide",
                             id=accession,
                             rettype=format,
                             retmode=mode)

        if mode=='xml':
            record = Entrez.read(handle)
        else:
            record = SeqIO.read(handle, 'genbank')

        handle.close()
        return record

    except Exception as e:
        print(f"### Error: {e}")

def get_json_data(accession, email, gb_file):
    """ Gets data form GenBank and saves in json format.
    Name of gb file is needed only to create json file name to be consistent with other file names."""
    if gb_file.endswith('.gb'):
        json_file =  gb_file.replace('.gb', '.json')
    else:
        json_file = gb_file+'.json'
    record = download_sequence(accession, email, "gb", "xml")
    # Convert to JSON
    #print(type(record))
    json_data = json.dumps(record[0], indent=2)
    # Save to file
    with open(json_file, "w") as f:
        f.write(json_data)

    return json_file

def gb_to_fasta(gb_file):

    if gb_file.endswith('.gb'):
        fasta_file =  gb_file.replace('.gb', '.fasta')
    else:
        fasta_file = gb_file+'.fasta'
    #fasta_file = os.path.join(output_dir, fasta_file_name)
    with open(fasta_file, "w") as output_handle:
        for record in SeqIO.parse(gb_file, "genbank"):
            SeqIO.write(record, output_handle, "fasta")
    print(f'{gb_file} converted to {fasta_file}')
    return fasta_file

def dict_from_row(row):
    dat = {}
    cols = row.index
    for col in cols:
        dat[col] = row[col]
    return dat

def get_features_for_feature(data, row, typ):
    f = data[(data['type'] == typ) & (data['name'].str.contains(row['name'] + "[", regex=False))]
    subdat = {}
    if f.shape[0] > 0:
        nr = 0
        for index, row in f.iterrows():
            nr += 1
            # print(f"\t\t exon: {row3['name']}")
            exon = dict_from_row(row)
            subdat[index] = exon

    return subdat
def find_next_gene(df, i):
    next_gene_index = len(df)
    for j in range(i + 1, len(df)):
        if df.iloc[j]['type'] == 'gene':
            next_gene_index = j
            break
    return next_gene_index


def make_json_data(data: pd.DataFrame, gb_file, seq_description):
    json_data = {}
    json_all_data = {}
    nr = 0
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

    misc_features = data[data['type'] == 'misc_feature']
    for i, row in misc_features.iterrows():
        nr += 1
        misc_feature = dict_from_row(row)
        json_data[nr] = misc_feature

    json_all_data['features'] = json_data
    with open(json_file, 'w') as res_file:
        json.dump(json_all_data, res_file, indent=4)
    return json_file
"""
def make_json_data(data : pd.DataFrame, gb_file, seq_description):
    json_data = {}
    json_all_data = {}
    nr = 0
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
        cdss = data[i+1:next_gene_index+1][(data['type'] == 'CDS')]
        if cdss.shape[0] > 0:
            #print(f"Gene: {row['name']} {i}:")
            cds_d = {}
            cd_nr = 0
            for index2, row2 in cdss.iterrows():
                cd_nr += 1
                #print(f"\t CDS: {row2['name']} {index2}")
                cds = dict_from_row(row2)
                cds['exons'] = get_features_for_feature(data, row, 'exon')
                cds['introns'] = get_features_for_feature(data, row, 'introns')
                cds_d[cd_nr] = cds
            gene['CDS'] = cds_d
        trnas = data[i+1:next_gene_index+1][(data['type'] == 'tRNA')]
        if trnas.shape[0] > 0:
            for index2, row2 in trnas.iterrows():
                #print(f"\t tRNA: {row2['name']}")
                trna = dict_from_row(row2)
                trna['tRNA'] = get_features_for_feature(data, row, 'tRNA')
                gene['tRNA'] = trna
        rrnas = data[i + 1:next_gene_index + 1][(data['type'] == 'rRNA')]
        if rrnas.shape[0] > 0:
            for index2, row2 in rrnas.iterrows():
                # print(f"\t rRNA: {row2['name']}")
                rrna = dict_from_row(row2)
                rrna['rRNA'] = get_features_for_feature(data, row, 'tRNA')
                gene['rRNA'] = rrna
        json_data[nr] = gene
    misc_features = data[data['type'] == 'misc_feature']
    for i, row in misc_features.iterrows():
        nr += 1
        misc_feature = dict_from_row(row)
        json_data[nr] = misc_feature

    json_all_data['features'] = json_data
    with open(json_file, 'w') as res_file:
        json.dump(json_all_data, res_file, indent=4)
    return json_file
"""

def download_files_for_organism(acc, email, output_dir):
    # Getting gb file and Organism name
    record = download_sequence(acc, email, 'gb')

    # Get organism name
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
    # Converting to FASTA file
    fasta_file = gb_to_fasta(gb_file)
    json_file = get_json_data(acc, email, gb_file)
    seq_description, json_file = export_features_from_gb_to_tsv(gb_file)
    return gb_file, fasta_file, json_file, seq_description

def download_all_sequences(query_acc, subject_acc, email):
    query_gb_file, query_fasta_file, json_file, seq_description  = download_files_for_organism(query_acc, email)
    subject_gb_file, subject_fasta_file, json_file, seq_description = download_files_for_organism(subject_acc, email)
    return query_gb_file, query_fasta_file, subject_gb_file, subject_fasta_file

def copy_file(output_dir, in_file):
    copy_file = os.path.join(output_dir, in_file)
    shutil.copy(in_file, copy_file)
    print(f'File {in_file} copied to {copy_file}.')
    return copy_file

def use_gb_file(output_dir, gb_file):
    copy_gb_file = copy_file(output_dir, gb_file)
    fasta_file = gb_to_fasta(copy_gb_file)
    seq_description, json_file = export_features_from_gb_to_tsv(copy_gb_file)
    return copy_gb_file, fasta_file, seq_description, json_file

def process_input_sequences(type_seq, acc_acc, acc_gb_file, acc_fasta_file, output_dir, email):
    print(f'\n----------------------------\nProcessing {type_seq}')
    if acc_acc == '':
        print(f'No acc. number for {type_seq} sequence is provided. Checking if files are provided.')
        if acc_gb_file != '':
            if os.path.isfile(acc_gb_file):
                print(f'{type_seq} gb file {acc_gb_file} for {type_seq} sequence is provided, I will try to use it!')
                acc_gb_file, acc_fasta_file, seq_description, json_file = use_gb_file(output_dir, acc_gb_file)
            else:
                print(f'File {acc_gb_file} cannot be found.')
                sys.exit(1)
        else:
            print(f'No gb file for {type_seq} sequence is provided.')
            if acc_fasta_file != '':
                if os.path.isfile(acc_fasta_file):
                    print(f'{type_seq} fasta file {acc_fasta_file} for {type_seq} sequence is provided. I will try to use it,' +
                          f' but no sequence features for {type_seq} will be shown!')
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
                print(f'ERROR: No acc. number or gb or fasta file was provided for {type_seq} sequence. I need at least one of them!')
                sys.exit(1)
    else:
        check_email(email, type_seq)
        print(f'Acc. number {acc_acc} for {type_seq} sequence is given. I will try to download data from GenBank.')
        acc_gb_file, acc_fasta_file, json_file, seq_description = download_files_for_organism(acc_acc, email, output_dir)
        print(f'{type_seq} files: {acc_gb_file}, {acc_fasta_file}, {json_file}')

    return os.path.basename(acc_gb_file), os.path.basename(acc_fasta_file), json_file
def check_email(email, type_seq):
    if email == '':
        print(f'ERROR: You did not provided gb or fasta file provided for {type_seq} sequence.'+
              ' I need your e-mail address (-e argument) to fetch data from GenBank.'+
              ' It is required by NCBI.')
        sys.exit(1)

def export_features_from_gb_to_tsv(gb_file):
    """Extracts features from gb file and saves them to tsv file."""
    features = []
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
    for record in SeqIO.parse(gb_file, "genbank"):
        nr = 0
        #print(f'%%%%%% RECORD {dir(record)}')
        #print(record.annotations)
        seq_description = {
            'seq_len': len(record.seq),
            'name': record.name,
            'id' : record.id,
            'description': record.description,
            'source': record.annotations['source'],
            'organism': record.annotations['organism'],
            'taxonomy': record.annotations['taxonomy'],
        }


        for feature in record.features:
            nr += 1

            # New record
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
            }
            if feature.type in ['exon', 'intron', 'gene', 'CDS', 'tRNA', 'rRNA', 'misc_feature']:
                rec['type'] = feature.type
                rec['start'] = int(feature.location.start)
                rec['end'] = int(feature.location.end)
                rec['strand'] = feature.location.strand
                if 'note' in feature.qualifiers:
                    rec['note'] = feature.qualifiers['note'][0]

                # Check if plastid-derived misc_feature
                is_plastid = (feature.type == 'misc_feature' and
                              'note' in feature.qualifiers and
                              'plastid-derived' in feature.qualifiers['note'][0])
                rec['plastid-derived'] = is_plastid
                # Getting the name of the feature
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0]
                    if 'number' in feature.qualifiers:
                        name += f"[{feature.qualifiers['number'][0]}]"

                elif 'product' in feature.qualifiers:
                    rec['product'] = feature.qualifiers['product'][0]
                    name = feature.qualifiers['product'][0]
                elif 'note' in feature.qualifiers:
                    name = feature.qualifiers['note'][0]  # Pierwsze 20 znaków z note
                    rec['note'] = feature.qualifiers['note'][0][:20]
                else:
                    name = feature.type
                rec['name'] = name
                new_df = pd.DataFrame([rec])
                data = pd.concat([data, new_df])
    tsv_file = gb_file.replace(".gb", ".tsv")
    data.reset_index(drop=True, inplace=True)
    data.to_csv(tsv_file, sep='\t', index=False)
    print(f'Features from {gb_file} exported to {tsv_file}')
    json_file = make_json_data(data, gb_file, seq_description)

    return tsv_file, json_file

