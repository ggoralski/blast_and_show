# BLAST and Show (v. 0.1.0)

## Overview

This program performs a BLAST search between query and subject sequences, then visualizes the results in multiple formats (SVG, PDF, PNG). See the `EXAMPLE_results` folder for example. The program accepts various input types including FASTA files, GenBank files, and GenBank accession numbers.
**Warning**: It is created for plastid/mitochondrial genomes or shorter sequences. It was not tested e.g., for eukaryotic chromosomes.
## Features

- **Flexible Input**: Accepts query and subject sequences from GenBank files, FASTA files, or accession numbers
- **BLAST Execution**: Performs local BLAST searches with customizable parameters
- **Visualization**: Generates detailed visualizations showing sequence alignments and features
- **Multiple Output Formats**: Creates SVG, PDF, and PNG images
- **Customizable Settings**: Configuration through `settings.json` file

## Requirements

### Python Dependencies

#### Option 1: Using pip
```bash
pip install biopython pandas matplotlib pycairo
```

#### Option 2: Using conda (recommended for Cairo dependencies)
```bash
# Create a new conda environment (optional but recommended)
conda create -n blast-and-show python=3.12
conda activate blast-and-show

# Install dependencies
conda install -c conda-forge biopython pandas matplotlib pycairo

# Alternative: Install all at once
conda install -c conda-forge -c bioconda biopython pandas matplotlib pycairo
```

#### Option 3: Using conda-forge and bioconda channels
```bash
# Add channels (one-time setup)
conda config --add channels conda-forge
conda config --add channels bioconda

# Install dependencies
conda install biopython pandas matplotlib pycairo
```

### System Requirements
- Python 3.6 or higher (Python 3.12+ recommended)
- BLAST+ command line tools (makeblastdb, blastn)
- Cairo graphics library (automatically handled by conda)

### Installing BLAST+

#### Using conda (recommended)
```bash
conda install -c bioconda blast
```

#### Using system package managers
- **Ubuntu/Debian**: `sudo apt-get install ncbi-blast+`
- **macOS**: `brew install blast`
- **Windows**: Download from [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### Complete conda setup (all dependencies)
```bash
# Create environment and install everything at once
conda create -n blast-and-show -c conda-forge -c bioconda python=3.12 biopython pandas matplotlib pycairo blast

# Activate environment
conda activate blast-and-show
```

## Usage

### Basic Syntax
```bash
python blast_and_show.py [options]
```

### Command Line Options

- `-q, --query`: Query accession number in GenBank
- `-s, --subject`: Subject accession number in GenBank  
- `-a, --query_gb_file`: GenBank format file for query sequence
- `-b, --query_fasta_file`: FASTA file with query sequence(s)
- `-c, --subject_gb_file`: GenBank format file for subject sequence
- `-d, --subject_fasta_file`: FASTA file with subject sequence
- `-o, --output_dir`: Output directory (default: 'output')
- `-e, --email`: Your email address (required by NCBI for downloads)
- `-v, --e_value`: E-value threshold for BLAST (default: 1e-10)
- `-m, --minlen`: Minimum length of alignments to display (default: 100 bp)
- `-t, --title`: Title for graphics (default: "Results of blastn")

### Input Requirements

You must provide **one** of the following for both query and subject:
- GenBank accession number (requires email)
- GenBank file (.gb)
- FASTA file (.fasta)

## Examples

### 1. Using GenBank Files
```bash
python blast_and_show.py -a query.gb -c subject.gb -o results
```
Uses GenBank files for both query and subject sequences.

### 2. Mixed Input Types
```bash
python blast_and_show.py -a query.gb -s subject_acc.nr. -o results -e your@email.address
```
Uses a GenBank file for query and downloads subject sequence using accession number.

### 3. Using Accession Numbers
```bash
python blast_and_show.py -q query_acc.nr. -s subject_acc.nr. -o results -e your@email.address
```
Downloads both sequences using GenBank accession numbers.

### 4. FASTA and GenBank Combination
```bash
python blast_and_show.py -b query_sequence.fasta -c subject_genome.gb -o results
```
Uses a FASTA file for query and GenBank file for subject.

### 5. Custom Parameters
```bash
python blast_and_show.py -q query_acc.nr. -s subject_acc.nr. -o results -v 1e-5 -m 50 -e your@email.address
```
Uses custom E-value (1e-5) and minimum alignment length (50 bp).

## Configuration

The `settings.json` file contains visualization parameters including:
- Image dimensions and margins
- Font settings
- Color schemes for different sequence features
- Alignment visualization parameters

**Note**: Your first visualization may not meet expectations, so adjustments to `settings.json` may be necessary.

## Output Files

The program generates:
- **Visualization files**: SVG, PDF, and PNG formats
- **BLAST results**: XML and TSV formats
- **Sequence files**: FASTA files for aligned sequences
- **Feature data**: JSON files with sequence annotations

## Email Requirement

When using GenBank accession numbers, you must provide a valid email address as required by NCBI's usage policies.


The tool generates several output files in the specified output directory:

### Sequence Files
- `*.gb`: GenBank format files (downloaded or copied)
- `*.fasta`: FASTA format files for BLAST input
- `*.json`: Structured JSON files with sequence features
- `*.tsv`: Tab-separated files with extracted features

### BLAST Results
- `blast_results.xml`: BLAST results in XML format
- `blast_results.tsv`: BLAST results in tabular format with custom columns
- `aligned_query-*.fasta`: FASTA file with aligned query sequences
- `aligned_subject-*.fasta`: FASTA file with aligned subject sequences

### Visualizations
- `out_draw.svg`: Scalable vector graphics visualization
- `out_draw.pdf`: PDF format visualization
- `out_draw.png`: PNG raster image visualization

### Configuration
- `settings.json`: Copy of configuration file used for the analysis

## Configuration - settings.json

The `settings.json` file controls all aspects of the visualization appearance and behavior. It contains two main sections:

### Settings Section

#### Drawing Dimensions and Layout
- `ratio` (50): Scale ratio for converting base pairs to pixels (1:50 means 1 pixel per 50 bp)
- `left_margin` (100): Left margin in pixels
- `right_margin` (300): Right margin in pixels  
- `top_margin` (100): Top margin in pixels
- `bottom_margin` (20): Bottom margin in pixels
- `main_height` (600): Total canvas height in pixels
- `query_y_position` (300): Y-coordinate for query sequence line
- `subject_y_position` (500): Y-coordinate for subject sequence line

#### Tick Marks and Coordinates
- `tick` (1000): Interval between tick marks in base pairs
- `tick_line_width` (0.5): Width of tick mark lines
- `tick_length` (20): Length of tick marks in pixels
- `tick_font` ("Sans-Serif"): Font family for coordinate labels

#### Title Formatting
- `title_font` ("Sans-Serif"): Font family for main title
- `title_font_size` (12): Font size for main title

#### Feature Visualization
- `feature_height` (20): Height of feature rectangles in pixels
- `feature_margin` (65): Vertical spacing between sequence line and features
- `features_alpha` (50): Transparency level for features (0-100)
- `features_rendered`: Array of feature types to display
  - `["CDS", "exon", "intron", "tRNA", "rRNA", "misc_feature"]`
- `special_features`: Array of special feature types with custom positioning
  - `["plastid-derived"]`
- `feature_data`: Array of data fields to include in feature labels
  - `["name", "type", "start", "end"]`

#### Alignment Visualization
- `alignment_height` (10): Height of alignment rectangles in pixels
- `alignment_alpha` (0.9): Transparency for alignment rectangle fill
- `alignment_line_alpha` (0.3): Transparency for alignment rectangle outlines
- `alignment_line_width` (1.5): Width of alignment rectangle outline
- `alignment_font` ("Sans-Serif"): Font for alignment labels
- `alignment_font_size` (10): Font size for alignment labels
- `alignments_alpha` (50): Overall transparency for alignments
- `query_alignment_data`: Data fields shown for query alignments
  - `["pident", "length", "gaps", "qstart", "qend"]`
- `subject_alignment_data`: Data fields shown for subject alignments
  - `["pident", "length", "gaps", "sstart", "send"]`

#### Color Scale
- `min_scale` (70): Minimum identity percentage for color scale
- `scale_move_left` (100): Horizontal offset for scale position
- `scale_move_down` (50): Vertical offset for scale position
- `scale_step_height` (10): Height of each scale step
- `scale_step_width` (10): Width of each scale step
- `scale_font` ("Sans-Serif"): Font for scale labels
- `scale_font_size` (8): Font size for scale labels
- `scale_font_color` ([0.1, 0.1, 0.1]): RGB color for scale text

#### Background and Spacing
- `background` ([1, 1, 1, 1]): Background color as RGBA values (white)
- `space_description` (100): Spacing for description text

### Feature Colors Section

Defines colors for different genomic features using RGB values (0-1 range):

- `alpha` (0.5): Default transparency level for all features
- `CDS` ([0, 0.8, 0]): Bright green for coding sequences
- `exon` ([0, 0.8, 0]): Bright green for exons
- `intron` ([0.6, 0.8, 0.6]): Light green for introns
- `tRNA` ([0, 0, 0.9]): Blue for transfer RNA genes
- `rRNA` ([0.4, 0.4, 0]): Olive for ribosomal RNA genes
- `misc_feature` ([0.8, 0.8, 0.8]): Light gray for miscellaneous features
- `plastid-derived` ([1, 0, 0]): Red for plastid-derived sequences

### Customizing Settings

You can modify `settings.json` to customize the visualization:

#### Changing Colors
```json
"feature_colors": {
    "CDS": [1, 0, 0],        // Red instead of green
    "tRNA": [0.5, 0, 0.5],   // Purple instead of blue
    "alpha": 0.8             // More opaque features
}
```

#### Adjusting Scale and Layout
```json
"settings": {
    "ratio": 25,             // Higher resolution (1:25)
    "feature_height": 30,    // Taller features
    "min_scale": 80,         // Show only high-identity alignments
    "tick": 500              // Fewer tick marks
}
```

#### Filtering Features
```json
"settings": {
    "features_rendered": ["CDS", "tRNA", "rRNA"]  // Show only these features
}
```


## Troubleshooting

- Ensure BLAST+ tools are installed and in your PATH
- Verify that all required Python packages are installed
- Check that input files exist and are in the correct format
- Provide a valid email address when using accession numbers
