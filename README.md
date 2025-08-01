

# Examples of use cases:

- `python blast_and_show.py -a Laburnum_anagyroides-OZ176120.gb -c Cuscuta_pedicellata-NC_052872.gb -o results` - gb files are provided both for query and subject sequences
- `python blast_and_show.py -a Laburnum_anagyroides-OZ176120.gb -s NC_052872 -o results -e your@email.address` - gb file for query is provided and accession number for subject sequence are provided
- `python blast_and_show.py -q OZ176120 -s NC_052872 -o results_gb -e ggoralski@gmail.com` - accession numbers for query and subject sequence are provided
- `python blast_and_show.py -q OZ176120 -s BK059237 -o results_gb -e ggoralski@gmail.com` - mitochondrial genomes accession numbers for query and subject sequence are provided
Cuscuta_epilinum-BK059237.gb
- `python blast_and_show.py -a Laburnum_anagyroides-OZ176120.gb -c Cuscuta_epilinum-BK059237.gb -o results` - gb files are provided both for query and subject sequences

- `python blast_and_show.py -a Cuscuta_pedicellata-NC_052872.gb -c Cuscuta_epilinum-BK059237.gb -o results` - gb files are provided both for query and subject sequences
- `python blast_and_show.py -b 223-LaF4R4-223-LaF1R1.fasta -c Cuscuta_epilinum-BK059237.gb -o results` - fasta and gb files are provided for query and subject sequences
- `python blast_and_show.py -b Cuscuta_epithymum_mt_328700-340700.fasta -s NC_070194 -o results -v 100 -e ggoralski@gmail.com` - fasta file and acc. number are provided 

