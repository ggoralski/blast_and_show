from argparse import ArgumentParser, RawDescriptionHelpFormatter
import dataimport as di
import blasting as bl
import draw as dr
import os
import sys
import shutil

def parse_arguments():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        usage='python blast_and_show.py [-q QUERY_ACC | -Q QUERY_GB | -a QUERY_FASTA]\n'
              '                          [-s SUBJECT_ACC | -S SUBJECT_GB | -b SUBJECT_FASTA]\n'
              '                          [-o OUTPUT_DIR] [-e EMAIL] [-c CONFIG]\n'
              '                          [-t TITLE] [-m MINLEN] [-v E_VALUE]',
        description=
            'Runs a local blastn comparison between two DNA sequences (query and subject) '
            'and produces a graphical visualisation of the alignments overlaid on both '
            'sequences with their annotated features (genes, CDS, tRNA, rRNA, misc_features).\n\n'
            'For each sequence you must provide ONE of:\n'
            '  - GenBank accession number (downloaded automatically; requires -e EMAIL)\n'
            '  - local GenBank file (.gb) — provides full feature annotation\n'
            '  - local FASTA file — no features will be drawn for that sequence\n\n'
            'Output (in OUTPUT_DIR) includes:\n'
            '  - blast_results.tsv / .xml — BLAST results in tabular and XML formats\n'
            '  - out_draw.svg / .pdf / .png — visualisations of the alignments\n'
            '  - aligned_query-*.fasta, aligned_subject-*.fasta — aligned regions\n'
            '  - copy of the configuration file used (settings.json)\n'
    )

    parser.add_argument("-q", "--query", dest="query", default="",
                        help="Query accession number in GenBank. Use INSTEAD of -Q/-a. "
                             "Requires -e EMAIL.")
    parser.add_argument("-Q", "--query_gb_file", dest="query_gb_file", default="",
                        help="Local GenBank (.gb) file for the query sequence. "
                             "Use INSTEAD of -q/-a. Provides feature annotation.")
    parser.add_argument("-a", "--query_fasta_file", dest="query_fasta_file", default="",
                        help="Local FASTA file for the query sequence. "
                             "Use INSTEAD of -q/-Q. No features will be drawn.")

    parser.add_argument("-s", "--subject", dest="subject", default="",
                        help="Subject accession number in GenBank. Use INSTEAD of -S/-b. "
                             "Requires -e EMAIL.")
    parser.add_argument("-S", "--subject_gb_file", dest="subject_gb_file", default="",
                        help="Local GenBank (.gb) file for the subject sequence. "
                             "Use INSTEAD of -s/-b. Provides feature annotation.")
    parser.add_argument("-b", "--subject_fasta_file", dest="subject_fasta_file", default="",
                        help="Local FASTA file for the subject sequence. "
                             "Use INSTEAD of -s/-S. No features will be drawn.")

    parser.add_argument("-e", "--email", dest="email", default='',
                        help="Your real e-mail address. Required by NCBI when downloading "
                             "sequences by accession number (-q / -s).")
    parser.add_argument("-o", "--output_dir", dest="output_dir", default='output',
                        help="Output directory. Will be created if it does not exist. "
                             "Default: output")
    parser.add_argument("-c", "--config", dest="config_file", default="settings.json",
                        help="Path to JSON configuration file with drawing settings "
                             "(geometry, fonts, colours, feature filters). "
                             "Default: settings.json. Two ready-made profiles are provided: "
                             "settings_short.json (1-15 kb sequences) and "
                             "settings_long.json (200-800 kb sequences).")

    parser.add_argument("-t", "--title", dest="title", default="Results of blastn",
                        help='Title displayed on the graphic. Default: "Results of blastn".')
    parser.add_argument("-m", "--minlen", dest="minlen", default=100,
                        help="Minimum alignment length (bp) to draw. Shorter alignments "
                             "will be filtered out. Default: 100.")
    parser.add_argument("-v", "--e_value", dest="e_value", default=1e-10,
                        help="E-value threshold for BLAST. Default: 1e-10.")

    args = parser.parse_args()

    return (args.query, args.subject,
            args.query_fasta_file, args.query_gb_file,
            args.subject_fasta_file, args.subject_gb_file,
            args.title, int(args.minlen), args.email,
            args.output_dir, args.e_value, args.config_file)


def main():
    (query_acc, subject_acc,
     query_fasta_file, query_gb_file,
     subject_fasta_file, subject_gb_file,
     title, min_len, email,
     output_dir, e_value, config_file) = parse_arguments()

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
        {config_file = }
    ''')

    # Sprawdzenie pliku konfiguracyjnego
    if not os.path.isfile(config_file):
        print(f'ERROR: Configuration file "{config_file}" not found.')
        sys.exit(1)

    if os.path.exists(output_dir):
        print(f'Output directory {output_dir} exists. I will use it.')
    else:
        print(f'Creating output directory: {output_dir}')
        os.makedirs(output_dir)

    query_gb_file, query_fasta_file, query_json_file = di.process_input_sequences(
        'query', query_acc, query_gb_file, query_fasta_file, output_dir, email)
    subject_gb_file, subject_fasta_file, subject_json_file = di.process_input_sequences(
        'subject', subject_acc, subject_gb_file, subject_fasta_file, output_dir, email)

    # Kopiowanie pliku konfiguracyjnego do output_dir jako settings.json
    # (draw.py szuka pliku o tej nazwie w bieżącym katalogu)
    target_config = os.path.join(output_dir, "settings.json")
    if os.path.abspath(config_file) != os.path.abspath(target_config):
        shutil.copy(config_file, target_config)
        print(f'Configuration "{config_file}" copied to "{target_config}".')
    else:
        print(f'Using configuration in place: {target_config}')

    main_dir = os.getcwd()
    os.chdir(output_dir)
    xml_file, tsv_file = bl.run_local_blast(query_fasta_file, subject_fasta_file, e_value)
    #os.chdir(main_dir)
    dr.run(os.path.basename(query_json_file), os.path.basename(subject_json_file), tsv_file)

if __name__ == "__main__":
    main()
