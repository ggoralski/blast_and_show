# BLAST Visualiser

A Python tool that runs a local **blastn** comparison between two DNA sequences
(query and subject) and produces a graphical visualisation of the resulting
alignments, overlaid on both sequences together with their annotated features
(genes, CDSs, exons, introns, tRNAs, rRNAs, misc_features).

The program is designed for comparative analyses of related genomes, with
particular focus on:

- Detection of horizontal gene transfer (HGT) between parasitic plants and
  their hosts (e.g. *Cuscuta*, *Orobanche*, *Phelipanche*).
- Comparison of organellar genomes (plastid, mitochondrial).
- Highlighting plastid-derived regions in mitochondrial genomes (intracellular
  gene transfer).
- Visualisation of conserved and divergent regions between related species.

---

## Features

- Local **blastn** execution — no on-line BLAST needed.
- Visualisation of alignments on both sequences with connecting lines.
- Annotated features drawn as coloured rectangles above/below each sequence
  line.
- Colour of each alignment encodes percentage identity (`pident`) using a
  rainbow gradient with a built-in legend.
- Special highlighting of user-defined feature classes (default:
  `plastid-derived` misc_features, drawn farther from the sequence line and in
  a dedicated colour).
- Three output formats: **SVG** (scalable), **PDF** (printable), **PNG**
  (raster).
- Fully configurable through a JSON file (geometry, fonts, colours, feature
  filters).
- Sequences can be supplied in three different ways and can be mixed freely
  between query and subject.

---

## Requirements

- Python 3.8 or newer
- **NCBI BLAST+** (`makeblastdb`, `blastn`) installed and available in `PATH`
- Python packages: `biopython`, `pandas`, `pycairo`, `matplotlib`

### Installation on Fedora

```bash
sudo dnf install ncbi-blast+ python3-cairo
pip install biopython pandas matplotlib
```

### Installation on Debian / Ubuntu

```bash
sudo apt install ncbi-blast+ libcairo2-dev pkg-config python3-dev
pip install biopython pandas pycairo matplotlib
```

### Installation with conda / mamba (any OS)

```bash
mamba install -c bioconda blast
mamba install -c conda-forge biopython pandas pycairo matplotlib
```

---

## Files in this repository

| File | Role |
|---|---|
| `blast_and_show.py` | Main entry-point script (command-line interface). |
| `blasting.py` | Builds the BLAST database and runs `blastn`. |
| `dataimport.py` | Downloads sequences from GenBank, parses `.gb` files, exports features. |
| `draw.py` | Builds the graphic (SVG / PDF / PNG) using Cairo. |
| `settings.json` | Default drawing configuration. |
| `settings_short.json` | Configuration tuned for short sequences (1–15 kb). |
| `settings_long.json` | Configuration tuned for long sequences (200–800 kb). |

---

## Usage

```
python blast_and_show.py [arguments]
```

### Input specification

For each sequence (query and subject) you must provide **exactly one** of:

1. **GenBank accession number** (`-q` / `-s`) — the sequence is downloaded
   automatically; an e-mail address (`-e`) is required by NCBI.
2. **Local GenBank file** (`.gb`) (`-a` / `-c`) — gives full feature
   annotation in the drawing.
3. **Local FASTA file** (`-b` / `-d`) — features will not be drawn for that
   sequence.

### Command-line arguments

| Short | Long | Description | Default |
|:---:|---|---|---|
| `-q` | `--query` | Query accession number in GenBank. | |
| `-Q` | `--query_gb_file` | Local GenBank file for query (gives features). | |
| `-a` | `--query_fasta_file` | Local FASTA file for query (no features). | |
| `-s` | `--subject` | Subject accession number in GenBank. | |
| `-S` | `--subject_gb_file` | Local GenBank file for subject (gives features). | |
| `-b` | `--subject_fasta_file` | Local FASTA file for subject (no features). | |
| `-e` | `--email` | Real e-mail address. Required by NCBI when using `-q`/`-s`. | |
| `-o` | `--output_dir` | Output directory (created if absent). | `output` |
| `-c` | `--config` | JSON configuration file with drawing settings. | `settings.json` |
| `-t` | `--title` | Title displayed on the graphic. | `"Results of blastn"` |
| `-m` | `--minlen` | Minimum alignment length (bp) to draw. | `100` |
| `-v` | `--e_value` | E-value threshold for BLAST. | `1e-10` |

> **Note on flag conventions:** the symmetry between query and subject is
> built into the short flags — lower-case for query (`-q`, `-Q`, `-a`),
> the matching switch in the same column for subject (`-s`, `-S`, `-b`).
> Accession numbers and GenBank files use the same letter
> (`-q`/`-Q`, `-s`/`-S`), FASTA files use neighbouring letters (`-a`/`-b`).

Run `python blast_and_show.py --help` to see the full inline help.

---

## Examples

```bash
# Both sequences provided as local GenBank files (host genome vs. parasite plastid)
python blast_and_show.py \
    -Q Laburnum_anagyroides-OZ176120.gb \
    -S Cuscuta_pedicellata-NC_052872.gb \
    -o results

# Query from a local GB file, subject downloaded by accession number
python blast_and_show.py \
    -Q Laburnum_anagyroides-OZ176120.gb \
    -s NC_052872 \
    -o results \
    -e your@email.address

# Both sequences downloaded by accession number
python blast_and_show.py \
    -q OZ176120 -s NC_052872 \
    -o results_gb \
    -e your@email.address

# Comparison of two mitochondrial genomes (Cuscuta) by accession
python blast_and_show.py \
    -q OZ176120 -s BK059237 \
    -o results_gb \
    -e your@email.address

# Two Cuscuta plastid genomes from local GB files
python blast_and_show.py \
    -Q Cuscuta_pedicellata-NC_052872.gb \
    -S Cuscuta_epilinum-BK059237.gb \
    -o results

# Custom FASTA query, GB subject
python blast_and_show.py \
    -a 223-LaF4R4-223-LaF1R1.fasta \
    -S Cuscuta_epilinum-BK059237.gb \
    -o results

# A mitochondrial sub-region in FASTA against a GenBank accession,
# raising the e-value threshold to 100 to keep weaker hits
python blast_and_show.py \
    -a Cuscuta_epithymum_mt_328700-340700.fasta \
    -s NC_070194 \
    -o results \
    -v 100 \
    -e your@email.address

# Using the long-sequence configuration profile
python blast_and_show.py \
    -q NC_052872 -s BK059237 \
    -o results \
    -e your@email.address \
    -c settings_long.json

# Using the short-sequence configuration profile and a custom title
python blast_and_show.py \
    -a region_5kb.fasta -b region_5kb_other.fasta \
    -o results \
    -c settings_short.json \
    -t "Pytheas region: Cuscuta vs. Orobanche"
```

---

## Output files

All output is written to `OUTPUT_DIR` (default: `output/`):

| File pattern | Content |
|---|---|
| `*.gb`, `*.fasta`, `*.json` | Sequence files for query and subject (downloaded or copied), and feature data in JSON. |
| `*.tsv` | Tabular export of features extracted from each GenBank file. |
| `blast_results.tsv` | BLAST results in tabular format (`-outfmt 6`) with a header row. |
| `blast_results.xml` | BLAST results in XML format (`-outfmt 5`) for further parsing. |
| `aligned_query-*.fasta` | Aligned regions of the query sequence in FASTA format. |
| `aligned_subject-*.fasta` | Aligned regions of the subject sequence in FASTA format. |
| `out_draw.svg`, `out_draw.pdf`, `out_draw.png` | The visualisation in three formats. |
| `settings.json` | Copy of the configuration used for this run. |

Temporary BLAST database files (`blast_db*`) are created and removed
automatically.

---

## Configuration

The drawing layout, fonts, colours and which feature types to display are all
controlled through a JSON configuration file. Three configuration files are
provided out of the box:

- **`settings.json`** — default profile, suitable for sequences in the order
  of tens of kb.
- **`settings_short.json`** — tuned for short sequences (1–15 kb), with bigger
  fonts, taller features and a denser tick spacing in base pairs.
- **`settings_long.json`** — tuned for long sequences (200–800 kb), such as
  plant mitochondrial genomes; smaller features, sparser ticks, larger image
  height and reduced label clutter.

Pass a configuration file with `-c FILE`. Default is `settings.json`.

### File structure

The configuration file has two top-level sections:

- **`settings`** — geometry, fonts, margins, what data to display.
- **`feature_colors`** — colours for each type of feature (CDS, intron, tRNA, etc.).

### General conventions

- **Colours** — lists of three numbers `[R, G, B]` in the range `0.0–1.0`
  (Cairo format, not `0–255`). The `background` parameter may have a fourth
  element — the alpha channel (transparency), e.g. `[1, 1, 1, 1]` =
  opaque white.
- **Transparency (alpha)** — values from `0.0` (fully transparent) to `1.0`
  (fully opaque).
- **Length units** — all geometric dimensions (margins, heights, lengths)
  are given in **Cairo pixels**. The actual width of the drawing in pixels
  depends additionally on the `ratio` parameter.
- **Fonts** — `"Sans-Serif"` means the system's default sans-serif font.
  You can also use a specific name like `"DejaVu Sans"` or `"Arial"`. If the
  named font is not available, Cairo substitutes a similar one.

### Parameters — `settings` section

**Scale and length conversion**

| Parameter | Default | Meaning |
|---|---|---|
| `ratio` | `50` | **Base pairs per pixel.** Higher = shorter drawing. For plant mitochondrial genomes (~500 kb) typically 50–200; for short sequences (~10 kb) go down to 5–10. |
| `tick` | `100` | **Tick spacing on the axis in base pairs.** Typical values: 1 000, 10 000, 100 000, depending on sequence length. |
| `tick_length` | `20` | Length of tick marks (in pixels), perpendicular to the sequence line. |
| `tick_line_width` | `0.5` | Tick line thickness. |
| `tick_font` | `"Sans-Serif"` | Font for tick labels. |
| `tick_font_size` | `2` | Tick label font size. **Note:** very small by default — tick labels are intentionally discrete. Increase if illegible. |

**Margins and overall dimensions**

The drawing dimensions are computed as:

```
draw_width  = left_margin + (sequence_length / ratio) + right_margin
draw_height = top_margin  + main_height               + bottom_margin
```

| Parameter | Default | Meaning |
|---|---|---|
| `left_margin` | `100` | Left margin — space for `"query"`/`"subject"` labels. |
| `right_margin` | `300` | Right margin — also holds the colour gradient legend. |
| `top_margin` | `100` | Top margin — above the title and first sequence line. |
| `bottom_margin` | `20` | Bottom margin. Increase if features below the subject line are cut off. |
| `main_height` | `600` | Height of the main drawing area. Determines room for both sequence lines and all features. |

**Query and subject line positions**

| Parameter | Default | Meaning |
|---|---|---|
| `query_y_position` | `300` | Y coordinate (from top) of the **query** sequence line. Features drawn above the line. |
| `subject_y_position` | `500` | Y coordinate of the **subject** line. Features drawn below the line. |

> The gap between the two lines (`subject_y_position − query_y_position`)
> determines the space for the connecting alignment lines. Too small → lines
> overlap and become unreadable.

**Features (genes, CDS, exons, tRNA, rRNA, misc_feature)**

| Parameter | Default | Meaning |
|---|---|---|
| `feature_height` | `20` | Height of the rectangle for a single feature. |
| `feature_margin` | `65` | Vertical distance from the sequence line to the feature track. Features listed in `special_features` are drawn at three times this distance (`feature_margin × 3`) for emphasis. |
| `features_rendered` | `["CDS", "exon", "intron", "tRNA", "rRNA", "misc_feature"]` | **List of feature types to draw.** Remove a type to hide it. Examples: `["CDS", "tRNA", "rRNA"]` hides exons, introns, misc_features; `["CDS"]` shows only coding sequences; `["misc_feature"]` is useful for viewing only plastid-derived regions. |
| `special_features` | `["plastid-derived"]` | List of features treated as "special" — drawn farther from the sequence line and coloured with a dedicated key in `feature_colors`. |
| `feature_data` | `["name", "type", "start", "end"]` | **Fields shown in the feature label** above each rectangle. The order matters. The `end` field is prefixed with a dash (`start-end`); others with a space. |

**Alignments (BLAST hits)**

| Parameter | Default | Meaning |
|---|---|---|
| `alignment_height` | `10` | Rectangle height for a single alignment. |
| `alignment_alpha` | `0.9` | Alignment rectangle opacity (`0.0`–`1.0`). |
| `alignment_line_alpha` | `0.3` | Opacity of the line connecting query and subject alignment. Lower = more discrete. |
| `alignment_line_width` | `0.5` | Width of the connecting line. |
| `alignment_font` | `"Sans-Serif"` | Font for alignment labels. |
| `alignment_font_size` | `5` | Font size for alignment labels. |
| `query_alignment_data` | `["pident", "length", "gaps", "qstart", "qend"]` | **Fields shown in the query-side alignment label.** `pident` is automatically prefixed with `%`; `qend`/`send` with a dash. |
| `subject_alignment_data` | `["pident", "length", "gaps", "sstart", "send"]` | Same, for the subject side. |

> **Available fields** in `query_alignment_data` / `subject_alignment_data`
> match BLAST tabular columns: `length`, `qlen`, `slen`, `qcovs`, `qcovhsp`,
> `pident`, `mismatch`, `gapopen`, `evalue`, `bitscore`, `gaps`, `qseqid`,
> `sseqid`, `qstart`, `qend`, `sstart`, `send`.

**Colour scale legend (gradient of percent identity)**

Alignments are coloured by `pident` (percent identity) using a rainbow
gradient from `min_scale` to 100 %.

| Parameter | Default | Meaning |
|---|---|---|
| `min_scale` | `70` | **Lower bound of the colour scale** (percent pident). Alignments with `pident < min_scale` have no assigned colour — must be ≤ the minimum pident in your BLAST result, otherwise the script raises `KeyError`. |
| `scale_move_left` | `100` | Horizontal position of the legend — distance from the **right edge**. |
| `scale_move_down` | `50` | Vertical position of the legend — distance from the top edge. |
| `scale_step_height` | `10` | Height of one rectangle in the legend. |
| `scale_step_width` | `10` | Width of one rectangle in the legend. |
| `scale_font` | `"Sans-Serif"` | Font for legend numbers. |
| `scale_font_size` | `8` | Font size for legend numbers. |
| `scale_font_color` | `[0.1, 0.1, 0.1]` | Colour of the legend numbers (almost black). |

**Title and background**

| Parameter | Default | Meaning |
|---|---|---|
| `title_font` | `"Sans-Serif"` | Title font. |
| `title_font_size` | `12` | Title font size. |
| `background` | `[1, 1, 1, 1]` | Background colour in RGBA `[R, G, B, alpha]`. Default white. For a transparent PNG use `[1, 1, 1, 0]`. |

### Parameters — `feature_colors` section

Feature colours in RGB format (`[R, G, B]`, each `0.0–1.0`).

| Key | Default | Meaning |
|---|---|---|
| `alpha` | `0.5` | **Common transparency** for all features. Applied to both the rectangle fill and its label. |
| `CDS` | `[0, 0.8, 0]` | Green — coding sequences. |
| `exon` | `[0, 0.8, 0]` | Green — exons (same as CDS). |
| `intron` | `[0.6, 0.8, 0.6]` | Pale green — introns. |
| `tRNA` | `[0, 0, 0.9]` | Blue — tRNA genes. |
| `rRNA` | `[0.4, 0.4, 0]` | Olive — rRNA genes. |
| `misc_feature` | `[0.8, 0.8, 0.8]` | Light grey — other annotations. |
| `plastid-derived` | `[1, 0, 0]` | Red — **plastid-derived misc_features** (HGT). Drawn farther from the sequence line. |

> **Adding new feature types:** to colour a new feature type (e.g. `repeat_region`),
> add a key under `feature_colors` and add the type to `features_rendered`.
> You also need to extend the list of accepted feature types in the
> `export_features_from_gb_to_tsv()` function in `dataimport.py`.

### Common editing scenarios

**Drawing is too wide**
Increase `ratio` (e.g. `50` → `100` — 2× horizontal compression). Then scale
`tick` up accordingly to avoid too many tick labels.

**Drawing is too short / lacks detail**
Decrease `ratio` (e.g. `50` → `20`). Decrease `tick` so a denser tick mark
spacing is still readable.

**Query and subject lines overlap with features**
Increase `feature_margin` (e.g. `65` → `100`), or move the lines further apart
by adjusting `query_y_position` and `subject_y_position`.

**Features below the subject line are cut off**
Increase `bottom_margin` (e.g. `20` → `100`), or increase `main_height` if you
want to keep more room for features above as well.

**Alignment labels unreadable**
Increase `alignment_font_size` (e.g. `5` → `8`), or shorten labels by removing
fields from `query_alignment_data` / `subject_alignment_data`.

**Too many alignments → visual clutter**
Filter results at load time (`-m` / `--minlen` argument in `blast_and_show.py`).
Lower `alignment_line_alpha` (e.g. to `0.1`) so connecting lines become more discrete.

**Gradient legend extends outside the image**
Increase `right_margin` (e.g. `300` → `400`), or decrease `scale_move_left`
(e.g. `100` → `60`) so the legend sits closer to the edge.

**I want to show only selected feature types**
Edit `features_rendered`. For example, to show only plastid-derived regions
without the background of all genes: `"features_rendered": ["misc_feature"]`.

**Error `KeyError: <number>` while drawing alignments**
At least one alignment has `pident < min_scale`. Lower `min_scale` (e.g.
`70` → `60`).

### Editing the file

The file is standard JSON — edit it in any text editor. Validate the syntax
after editing:

```bash
python -m json.tool settings.json
```

The command pretty-prints the file and reports any syntax errors (missing
comma, unclosed bracket, etc.).

---

## Troubleshooting

**`sh: makeblastdb: command not found`**
NCBI BLAST+ is not installed or not on `PATH`. Install with
`sudo dnf install ncbi-blast+` (Fedora) or
`sudo apt install ncbi-blast+` (Debian/Ubuntu).

**`ModuleNotFoundError: No module named 'cairo'`**
PyCairo is missing. Install with `sudo dnf install python3-cairo` (Fedora)
or `pip install pycairo` after installing the Cairo development headers
(`libcairo2-dev` on Debian/Ubuntu, `cairo-devel` on Fedora).

**`KeyError: <number>` while drawing alignments**
At least one alignment has `pident` below `min_scale` (default 70). Lower
`min_scale` in the configuration file.

**Features below the subject line are cut off**
Increase `bottom_margin` (or `main_height`) in the configuration file.

---

## How to cite

If you use this software in your research or it contributes to a publication,
we would appreciate it if you cited the following paper:

> Denysenko-Bennett, M., Kwolek, D., Góralski, G., Szklarczyk, M.,
> Piwowarczyk, R., Stefanović, S., Schneider, A. C. & Joachimiak, A. J.
> Horizontal gene transfer of the Pytheas sequence from *Cuscuta* to
> *Orobanche* via a host-mediated pathway. *Scientific Reports* **16**,
> 2056 (2026). https://doi.org/10.1038/s41598-025-31853-x

BibTeX:

```bibtex
@article{DenysenkoBennett2026Pytheas,
  author  = {Denysenko-Bennett, Magdalena and Kwolek, Dagmara and
             G{\'o}ralski, Grzegorz and Szklarczyk, Marek and
             Piwowarczyk, Renata and Stefanovi{\'c}, Sa{\v s}a and
             Schneider, Adam C. and Joachimiak, Andrzej J.},
  title   = {Horizontal gene transfer of the {Pytheas} sequence from
             \textit{Cuscuta} to \textit{Orobanche} via a host-mediated pathway},
  journal = {Scientific Reports},
  volume  = {16},
  number  = {1},
  pages   = {2056},
  year    = {2026},
  doi     = {10.1038/s41598-025-31853-x}
}
```
