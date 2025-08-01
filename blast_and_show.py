from argparse import ArgumentParser
import dataimport as di
import blasting as bl
import draw as dr
import os
import shutil

def parse_arguments():
    parser = ArgumentParser(usage='Usage:\npython visualise_blast_results.py -x blast_results.xml -g gb-file.gb -t "Title for graphics" -m 200',
                        description='Prepares graphics for blastn results, with alignments on target file and described features of target sequence.\n'+
                                    'It needs:\n'+
                                    '  - results of blastn in xml file'+
                                    '  - target sequence in gb format with described features.\n'+
                                    'The blastn should be previously executed as e.g.:\n'+
                                    '\tmakeblastdb -in dna_target.fasta -dbtype nucl -out database_name\n'+
                                    '\tblastn -query query.fasta -db database_name -out blast_results.xml -outfmt 5'+
                                    'The target sequence in gb format may be fetched from GenBank e.g using command line efetch tool:'+
                                    '\tefetch -db nuccore -id acc_number -format gb > target.gb'
                            )

    parser.add_argument("-q", "--query", dest="query", default="",
                        help="Query accession number in GenBank. Required if (and only if) FASTA file for the query "+
                        "(-a argument) is not given.")
    parser.add_argument("-s", "--subject", dest="subject", default="",
                        help="Subject accession number in GenBank.  Required if (and only if) FASTA and gb files for "+
                             "the subject (-b and c- arguments) are not given.")
    parser.add_argument("-a", "--query_gb_file", dest="query_gb_file", default="",
                        help="GenBank format (gb) file for query sequence. Required if (and only if) accession number "
                        "for the subject (-q argument) is not given.")
    parser.add_argument("-b", "--query_fasta_file", dest="query_fasta_file", default="",
                        help="FASTA file with query sequence(s). Required if (and only if) accession number for the query "+
                         "(-q argument) is not given.")

    parser.add_argument("-c", "--subject_gb_file", dest="subject_gb_file", default="",
                        help="GenBank format (gb) file for subject sequence. Required if (and only if) accession number "
                        "for the subject (-s argument) is not given.")
    parser.add_argument("-d", "--subject_fasta_file", dest="subject_fasta_file", default="",
                        help="FASTA file with subject sequence. Required if (and only if) accession number for the subject "+
                        "(-s argument) is not given.")

    parser.add_argument("-t", "--title", dest="title", default="Results of blastn",
                    help="Title of graphics, also used (without spaces) as the name of graphics files.")
    parser.add_argument("-m", "--minlen", dest="minlen", default=100,
                        help="Minimal length (bp) of alignments to draw (default 100).")
    parser.add_argument("-e", "--email", dest="email", default='',
                        help="Your true e-mail address (required by NCBI!).")
    parser.add_argument("-o", "--output_dir", dest="output_dir", default='output',
                        help="Output directory. Default: output")
    parser.add_argument("-v", "--e_value", dest="e_value", default=1e-10,
                        help="e_value for blast. Default: 1e-10")


    args = parser.parse_args()
    query_acc = args.query
    subject_acc = args.subject
    query_fasta_file = args.query_fasta_file
    query_gb_file = args.query_gb_file
    subject_fasta_file = args.subject_fasta_file
    subject_gb_file = args.subject_gb_file
    title = args.title
    min_len = int(args.minlen)
    email = args.email
    output_dir =  args.output_dir
    e_value = args.e_value

    return query_acc, subject_acc, query_fasta_file, query_gb_file, subject_fasta_file, subject_gb_file, title, min_len, email, output_dir, e_value

def main():
    query_acc, subject_acc, query_fasta_file, query_gb_file, \
    subject_fasta_file, subject_gb_file, title, min_len, email, \
    output_dir, e_value = parse_arguments()
    #query_gb_file, query_fasta_file, subject_gb_file, subject_fasta_file = download_all_sequences(query_acc, subject_acc, email)
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
    if os.path.exists(output_dir):
        print(f'Output directory {output_dir} exists. I will use it.')
    else:
        print(f'Creating output directory: {output_dir}')
        os.makedirs(output_dir)

    query_gb_file, query_fasta_file, query_json_file = di.process_input_sequences('query',query_acc, query_gb_file, query_fasta_file,
                                                              output_dir, email)
    subject_gb_file, subject_fasta_file, subject_json_file = di.process_input_sequences('subject', subject_acc, subject_gb_file,
                                                                  subject_fasta_file, output_dir, email)
    shutil.copy("settings.json", output_dir)
    main_dir = os.getcwd()
    os.chdir(output_dir)
    xml_file, tsv_file = bl.run_local_blast(query_fasta_file, subject_fasta_file, e_value)
    #os.chdir(main_dir)
    dr.run(os.path.basename(query_json_file), os.path.basename(subject_json_file), tsv_file)
if __name__ == "__main__":
  main()
