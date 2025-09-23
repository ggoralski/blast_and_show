"""
Visualization and Drawing Module

This module handles the creation of sequence alignment visualizations using Cairo graphics.
It generates publication-quality images showing:

- Sequence alignments between query and subject sequences
- Genomic features (genes, CDS, tRNA, rRNA, exons, introns)
- Color-coded alignment identity scales
- Coordinate systems and annotations

The module supports multiple output formats:
- SVG (Scalable Vector Graphics)
- PDF (Portable Document Format)
- PNG (Portable Network Graphics)

Key Features:
- Automatic scaling and positioning of sequences
- Color-coded visualization of alignment identity
- Feature-based annotation display
- Customizable appearance through settings.json
- Support for complex genomic structures

Visualization Components:
- Sequence lines with coordinate ticks
- Feature rectangles with type-specific colors
- Alignment blocks connecting query and subject
- Identity scale legend
- Descriptive labels and annotations

Dependencies:
    - Cairo graphics library (pycairo)
    - matplotlib (for color generation)
    - pandas (for data processing)
    - json (for settings and data)

Author: [Author Name]
Date: [Date]
Version: 1.0
"""

import json
import sys
from copy import copy
import os
import pandas as pd
import math
import cairo
import matplotlib.pyplot as plt


def generate_gradient(x):
    """
    Generate a color gradient dictionary for alignment identity visualization.

    This function creates a rainbow color gradient from a minimum percentage
    identity value to 100%. The colors are used to represent different levels
    of sequence identity in alignment visualizations.

    Args:
        x (int): Minimum percentage identity value (0-100)

    Returns:
        dict: Dictionary mapping percentage values to RGB color tuples
              Keys: integers from x to 100
              Values: RGB tuples with values 0.0-1.0

    Raises:
        ValueError: If x is not between 0 and 100

    Example:
        >>> gradient = generate_gradient(80)
        >>> gradient[90]  # Returns RGB tuple for 90% identity
        (0.2, 0.8, 0.6)

    Note:
        Uses matplotlib's rainbow colormap for smooth color transitions
    """
    # Validate input range
    x = int(x)
    if x < 0 or x > 100:
        raise ValueError("Minimum value must be between 0 and 100.")

    # Create range of percentage values
    keys = list(range(x, 101))
    num_colors = len(keys)

    # Generate rainbow color scale
    cmap = plt.cm.get_cmap("rainbow", num_colors)
    colors = [tuple(round(c, 2) for c in cmap(i)[:3]) for i in range(num_colors)]

    # Create percentage-to-color mapping
    color_dict = {key: color for key, color in zip(keys, colors)}
    return color_dict


def read_data(input_file):
    """
    Read and parse JSON data from a file.

    This function loads sequence and feature data from JSON files created
    during the data processing pipeline.

    Args:
        input_file (str): Path to JSON file containing sequence data

    Returns:
        dict: Parsed JSON data containing sequence information and features

    Raises:
        FileNotFoundError: If input file doesn't exist
        json.JSONDecodeError: If file contains invalid JSON

    Note:
        Expected JSON structure includes 'description' and 'features' sections
    """
    print(f'Reading JSON data from input file: {input_file}')
    with open(input_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    return data


def prepare_draw_settings(query_file, subject_file, alignments_file, settings):
    """
    Prepare drawing settings and load all necessary data for visualization.

    This function loads sequence data, alignment results, and calculates drawing
    parameters including scaling, positioning, and color gradients. It determines
    the optimal layout for visualizing sequences of different lengths.

    Args:
        query_file (str): Path to query sequence JSON file
        subject_file (str): Path to subject sequence JSON file
        alignments_file (str): Path to BLAST results TSV file
        settings (dict): Visualization settings from settings.json

    Returns:
        tuple: A tuple containing:
            - draw_settings (dict): Calculated drawing parameters including:
                * scaling ratios
                * positioning offsets
                * image dimensions
                * color gradient
                * title and layout info
            - query_data (dict): Query sequence data and features
            - subject_data (dict): Subject sequence data and features
            - alignments_data (pandas.DataFrame): BLAST alignment results

    Raises:
        FileNotFoundError: If any input file doesn't exist
        KeyError: If required data fields are missing

    Note:
        Automatically adjusts layout based on sequence length differences
        and calculates appropriate scaling factors for visualization
    """
    # Load sequence and alignment data
    query_data = read_data(query_file)
    subject_data = read_data(subject_file)
    alignments_data = pd.read_csv(alignments_file, sep='\t')

    # Extract sequence lengths
    query_length, subject_length = 0, 0
    if query_data is not None:
        query_length = query_data["description"]["seq_description"]["seq_len"]
        print(f'Query length: {query_length}')
    if subject_data is not None:
        subject_length = subject_data["description"]["seq_description"]["seq_len"]
        print(f'Subject length: {subject_length}')

    # Calculate positioning adjustments for different sequence lengths
    left_move = round(abs(subject_length - query_length) / 2)

    # Determine which sequence is shorter and calculate positioning
    if query_length < subject_length:
        shorter = 'query'
        left_move_query = left_move
        left_move_subject = 0
    else:
        shorter = 'subject'
        left_move_query = 0
        left_move_subject = left_move

    # Calculate alignment identity range for color scaling
    # Calculate positioning adjustments for different sequence lengths
    left_move = round(abs(subject_length - query_length) / 2)

    # Determine which sequence is shorter and calculate positioning offsets
    # This ensures both sequences are visually centered when lengths differ
    if query_length < subject_length:
        shorter = 'query'
        left_move_query = left_move  # Move query right to center it
        left_move_subject = 0        # Subject stays at original position
    else:
        shorter = 'subject'
        left_move_query = 0          # Query stays at original position
        left_move_subject = left_move # Move subject right to center it

    # Calculate alignment identity range for color scaling
    min_pident = round(min(alignments_data['pident']))
    max_pident = round(max(alignments_data['pident']))
    print(f'Alignment identity range: {min_pident}% - {max_pident}%')

    # Compile all drawing settings
    draw_settings = {
        'main_length': math.ceil((max([query_length, subject_length]))),
        'query_length': math.ceil(query_length),
        'subject_length': math.ceil(subject_length),
        'gradient': generate_gradient(settings["min_scale"]),
        'min_pident': min_pident,
        'max_pident': max_pident,
        'query_data': query_data,
        'subject_data': subject_data,
        'left_move_query': round(left_move_query/settings['ratio']),
        'left_move_subject': round(left_move_subject/settings['ratio']),
        'shorter': shorter,
        'query_description': query_data["description"]["seq_description"]["description"],
        'subject_description': subject_data["description"]["seq_description"]["description"],
        'title': f'Query: {query_data["description"]["seq_description"]["description"]} '+
                 f'{query_data["description"]["seq_description"]["name"]} - ' +
                 f'Subject: {subject_data["description"]["seq_description"]["description"]} '+
                 f'{subject_data["description"]["seq_description"]["name"]}'
    }

    # Calculate final image dimensions
    draw_settings['draw_width'] = (settings['left_margin'] +
                                  math.ceil(draw_settings['main_length'] / settings['ratio']) +
                                  settings['right_margin'])
    draw_settings['draw_height'] = settings['top_margin'] + settings['main_height'] + settings['right_margin']

    return draw_settings, query_data, subject_data, alignments_data


def draw_sequence_line(context, settings, draw_settings, feature_colors, seq_data, y_position, seq_type):
    """
    Draw a sequence line with coordinate ticks and labels.

    This function draws the main sequence line representing either query or subject
    sequence, including coordinate ticks, labels, and calls the feature drawing function.
    It handles positioning adjustments for sequences of different lengths.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for different feature types
        seq_data (dict): Sequence data including length and features
        y_position (int): Vertical position for drawing the sequence line
        seq_type (str): Type of sequence ('query' or 'subject')

    Returns:
        cairo.Context: Updated drawing context

    Note:
        Automatically adjusts horizontal positioning based on sequence length differences.
        Draws coordinate ticks at regular intervals with rotated labels.
    """
    print(f'Drawing {seq_type} sequence line...')
    original_left_margin = copy(settings['left_margin'])

    # Apply positioning adjustments for shorter sequences
    if draw_settings['shorter'] == 'query' and seq_type == 'query':
        settings['left_margin'] += draw_settings['left_move_query']
    elif draw_settings['shorter'] == 'subject' and seq_type == 'subject':
        settings['left_margin'] += draw_settings['left_move_subject']

    # Draw sequence type label
    context.move_to(settings['left_margin']-45, y_position)
    context.select_font_face('Sans-serif')
    context.set_font_size(12)
    context.set_source_rgb(0, 0, 0)
    context.show_text(seq_type)

    # Draw main sequence line
    context.select_font_face(settings['tick_font'])
    context.set_font_size(settings['tick_font_size'])
    context.move_to(settings['left_margin'], y_position)
    length = math.ceil(seq_data['description']['seq_description']['seq_len']/settings['ratio'])
    context.line_to(length+settings['left_margin'], y_position)
    context.set_source_rgb(0.2, 0.2, 0.2)
    context.set_line_width(2)
    context.stroke()

    # Draw coordinate ticks and labels
    for x in range(0, length, math.ceil(settings['tick']/settings['ratio'])):
        x_pos = x + settings['left_margin']

        # Draw tick mark
        context.move_to(x_pos, y_position)
        context.line_to(x_pos, y_position+settings['tick_length'])
        context.set_source_rgb(0.2, 0.2, 0.2)
        context.set_line_width(settings['tick_line_width'])
        context.stroke()

        # Draw rotated coordinate label
        context.move_to(x_pos, y_position+settings['tick_length']+5)
        context.set_source_rgb(0.7, 0, 0)
        context.save()
        context.rotate(math.pi / 2)  # Rotate 90 degrees
        context.show_text(str(x*settings['ratio']))
        context.restore()

    # Draw genomic features on this sequence line
    draw_features(context, settings, draw_settings, feature_colors, seq_data, y_position, seq_type)

    # Restore left margin and return updated context
    settings['left_margin'] = original_left_margin
    return context


def draw_feature(context, settings, draw_settings, feature_colors, feature, y_position, i, seq_type, special=None):
    """
    Draw a single genomic feature as a colored rectangle with label.

    This function draws individual genomic features (genes, CDS, tRNA, etc.) as
    colored rectangles positioned relative to the sequence line. Features are
    drawn above or below the sequence line depending on sequence type, with
    alternating vertical positions to avoid overlap.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for different feature types
        feature (dict): Feature data including start, end, type, name
        y_position (int): Vertical position of the sequence line
        i (int): Feature index for alternating positioning
        seq_type (str): Type of sequence ('query' or 'subject')
        special (str, optional): Special feature type for custom coloring

    Returns:
        cairo.Context: Updated drawing context

    Note:
        - Query features are drawn above the sequence line
        - Subject features are drawn below the sequence line
        - Special features (like plastid-derived) get enhanced positioning
        - Feature labels are rotated and positioned at feature centers
    """
    # Calculate vertical offset for alternating feature positions
    i = i * settings['feature_height']/2

    # Set positioning and rotation based on sequence type
    if seq_type == 'subject':
        feature_margin = settings['feature_margin']
        rotation = math.pi / 2  # 90 degrees clockwise
        text_move_y = settings['feature_height'] + 5
    elif seq_type == 'query':
        feature_margin = -(settings['feature_margin'] + settings['feature_height'])
        rotation = (3*math.pi) / 2  # 270 degrees clockwise
        text_move_y = -5
        i = -i  # Reverse alternating direction for query
    else:
        print(f'ERROR: Invalid sequence type: {seq_type}')
        sys.exit(1)

    # Apply special positioning for highlighted features
    if special in settings['special_features']:
        feature_margin *= 3

    # Calculate feature dimensions and position
    length = math.ceil((abs(feature['end'] - feature['start']) / settings['ratio']))
    coordinates = (
        settings['left_margin'] + math.ceil(feature['start']/settings['ratio']),
        y_position + feature_margin + i,
        length,
        settings['feature_height']
    )

    # Draw feature rectangle
    context.rectangle(*coordinates)

    # Determine feature color (special features override default colors)
    if special is not None:
        if feature_colors.get(special) is not None:
            color = feature_colors.get(special)
        else:
            color = feature_colors[feature['type']]
    else:
        color = feature_colors[feature['type']]

    context.set_source_rgba(*color, feature_colors['alpha'])
    context.fill()

    # Draw feature label
    context.move_to(settings['left_margin'] + math.ceil(feature['start']/settings['ratio']) + length/2,
                    y_position + feature_margin + text_move_y + i)
    context.set_source_rgba(*color, feature_colors['alpha'])
    context.save()
    context.rotate(rotation)

    # Build feature label from configured data fields
    name = ''
    for d in settings['feature_data']:
        if d == 'end':
            name += "-" + str(feature[d])
        else:
            name += " " + str(feature[d])

    context.show_text(name)
    context.restore()

    return context


def draw_features(context, settings, draw_settings, feature_colors, seq_data, y_position, seq_type):
    """
    Draw all genomic features for a sequence.

    This function iterates through all features in a sequence and draws them
    using the draw_feature function. It handles the hierarchical structure of
    genomic features, drawing genes and their sub-features (CDS, tRNA, rRNA)
    as well as misc_features.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for different feature types
        seq_data (dict): Complete sequence data including features
        y_position (int): Vertical position of the sequence line
        seq_type (str): Type of sequence ('query' or 'subject')

    Returns:
        cairo.Context: Updated drawing context

    Note:
        Processes features in hierarchical order:
        1. Genes with their CDS, tRNA, rRNA sub-features
        2. Standalone misc_features
        3. Special handling for plastid-derived features
    """
    features = seq_data['features']
    n = 1

    # Process each feature in the sequence
    for key, feature in features.items():
        i = n % 2  # Alternate positioning for visual clarity
        n += 1

        # Handle gene features and their sub-features
        if feature['type'] == 'gene':
            if feature.get('CDS') is not None:
                cdss = feature.get('CDS')
                for key2, cds in cdss.items():
                    # Draw exons if present
                    if len(cds['exons']) > 0:
                        exons = cds['exons']
                        for key3, exon in exons.items():
                            draw_feature(context, settings, draw_settings, feature_colors,
                                       exon, y_position, i, seq_type)
                    # Draw introns if present
                    elif len(cds['introns']) > 0:
                        introns = cds['introns']
                        for key4, intron in introns.items():
                            draw_feature(context, settings, draw_settings, feature_colors,
                                       intron, y_position, i, seq_type)
                    # Draw CDS directly if no sub-features
                    else:
                        draw_feature(context, settings, draw_settings, feature_colors,
                                   cds, y_position, i, seq_type)
            # Handle tRNA features
            elif feature.get('tRNA') is not None:
                tRNA = feature.get('tRNA')
                draw_feature(context, settings, draw_settings, feature_colors,
                           tRNA, y_position, i, seq_type)
            # Handle rRNA features
            elif feature.get('rRNA') is not None:
                rRNA = feature.get('rRNA')
                draw_feature(context, settings, draw_settings, feature_colors,
                           rRNA, y_position, i, seq_type)
            else:
                print(f'Warning: Gene feature without recognized sub-features: {feature}')

        # Handle misc_feature entries
        elif feature['type'] == 'misc_feature':
            if feature['plastid-derived']:
                # Special handling for plastid-derived features
                draw_feature(context, settings, draw_settings, feature_colors, feature,
                           y_position, i, seq_type, 'plastid-derived')
            else:
                draw_feature(context, settings, draw_settings, feature_colors, feature,
                           y_position, i, seq_type)

    return context


def draw_scale(context, settings, draw_settings):
    """
    Draw the color scale legend showing alignment identity percentages.

    This function creates a color-coded scale that shows the relationship
    between colors and percentage identity values in the alignments. The
    scale is positioned in the right margin of the image.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters including gradient

    Returns:
        cairo.Context: Updated drawing context

    Note:
        Scale shows percentage values from min_scale to 100% with corresponding colors
    """
    x_pos = draw_settings['draw_width'] - settings['scale_move_left']
    y_pos = settings['scale_move_down']
    gradient = draw_settings['gradient']

    context.select_font_face(settings['scale_font'])
    context.set_font_size(settings['scale_font_size'])

    j = 0
    # Draw color blocks with percentage labels
    for i in range(settings['min_scale'], 101):
        y = y_pos + j * settings['scale_step_height'] + 5 * j
        coordinates = (
            x_pos,
            y,
            settings['scale_step_width'],
            settings['scale_step_height']
        )
        color = gradient[i]

        # Draw colored rectangle
        context.rectangle(*coordinates)
        context.set_source_rgba(*color, settings['alignment_alpha'])
        context.fill()

        # Draw percentage label
        context.move_to(x_pos + settings['scale_step_width'] + 5, y + settings['scale_step_height'])
        context.set_source_rgb(*settings['scale_font_color'])
        context.show_text(str(i))
        j += 1

    return context


def draw_alignment(context, settings, draw_settings, row):
    """
    Draw a single alignment block connecting query and subject sequences.

    This function draws an individual alignment as colored rectangles on both
    query and subject sequences, connected by a line. The color represents
    the percentage identity of the alignment.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        row (pandas.Series): Single alignment data from BLAST results

    Returns:
        cairo.Context: Updated drawing context

    Note:
        Alignment color is determined by percentage identity (pident column)
        Rectangles are positioned based on start/end coordinates in both sequences
    """
    # Determine horizontal margins based on which sequence is shorter.
    # This ensures proper alignment when sequences have different lengths by
    # applying an offset (left_move_query/left_move_subject) to the shorter one.
    if draw_settings['shorter'] == 'query':
        query_margin = settings['left_margin'] + draw_settings['left_move_query']
        subject_margin = settings['left_margin']
    elif draw_settings['shorter'] == 'subject':
        subject_margin = settings['left_margin'] + draw_settings['left_move_subject']
        query_margin = settings['left_margin']
    else:
        # Both sequences are treated as the same length — no offset needed.
        query_margin = settings['left_margin']
        subject_margin = settings['left_margin']

    # Calculate vertical positions for alignment rectangles.
    # Subtract half the alignment height so the rectangle is centered on the
    # sequence line (which is positioned at query_y_position/subject_y_position).
    y_query = settings['query_y_position'] - settings['alignment_height'] / 2
    y_subject = settings['subject_y_position'] - settings['alignment_height'] / 2

    # Handle alignments reported on the reverse strand by ensuring start < end.
    # BLAST and similar tools can report coordinates in reverse order; normalizing
    # them simplifies subsequent length and position calculations.
    if row['qend'] - row['qstart'] < 0:
        qstart = row['qend']
        qend = row['qstart']
    else:
        qstart = row['qstart']
        qend = row['qend']

    if row['send'] - row['sstart'] < 0:
        sstart = row['send']
        send = row['sstart']
    else:
        sstart = row['sstart']
        send = row['send']

    # Convert sequence coordinates to pixel positions using the scaling ratio.
    # math.ceil is used to avoid zero-width rectangles for very short alignments.
    x_start_query = query_margin + math.ceil(qstart / settings['ratio'])
    x_start_subject = subject_margin + math.ceil(sstart / settings['ratio'])

    # Compute the width of the alignment rectangles (in pixels) by scaling the
    # nucleotide/amino-acid span by the ratio. Use abs to be robust to ordering.
    length_query = math.ceil(abs(qend - qstart) / settings['ratio'])
    length_subject = math.ceil(abs(send - sstart) / settings['ratio'])

    # Pack rectangle coordinates as (x, y, width, height) for cairo.
    query_coordinates = (
        x_start_query,
        y_query,
        length_query,
        settings['alignment_height']
    )
    subject_coordinates = (
        x_start_subject,
        y_subject,
        length_subject,
        settings['alignment_height']
    )

    # Choose a color from the precomputed gradient based on percent identity.
    # round() maps the floating pident to an integer index into the gradient.
    color = draw_settings['gradient'][round(row['pident'])]

    # Draw the query alignment rectangle with the configured transparency.
    context.rectangle(*query_coordinates)
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.fill()
    put_alignment_description(context, settings, row, "query", x_start_query, y_query, color, length_query)

    # Draw the subject alignment rectangle with the configured transparency.
    context.rectangle(*subject_coordinates)
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.fill()
    put_alignment_description(context, settings, row, "subject", x_start_subject, y_subject, color, length_subject)

    # Draw a connecting line between the bottom-left of the query rectangle and
    # the top-left of the subject rectangle (or vice-versa depending on y positions).
    # This visually links the two corresponding aligned regions.
    context.move_to(x_start_query, y_query + settings['alignment_height'])
    context.line_to(x_start_subject, y_subject)
    context.set_source_rgba(*color, settings['alignment_line_alpha'])
    context.set_line_width(settings['alignment_line_width'])
    context.stroke()
    return context


def put_alignment_description(context, settings, row, seqtype, x, y, color, length):
    """
    Add descriptive text labels to alignment blocks.

    This function adds text labels to alignment rectangles showing relevant
    information like coordinates, identity percentage, and other BLAST statistics.
    Labels are positioned and rotated appropriately for query vs subject sequences.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings including label configuration
        row (pandas.Series): Alignment data from BLAST results
        seqtype (str): Type of sequence ('query' or 'subject')
        x (int): Horizontal position for label
        y (int): Vertical position for label
        color (tuple): RGB color tuple for text
        length (int): Length of alignment block

    Returns:
        cairo.Context: Updated drawing context

    Raises:
        SystemExit: If seqtype is not 'query' or 'subject'

    Note:
        Text rotation and positioning differ between query and subject sequences
        Label content is configurable through settings
    """
    # Set rotation and positioning based on sequence type
    if seqtype == 'query':
        rotation = (3 * math.pi) / 2  # 270 degrees
        x = x + length / 2
        data = settings['query_alignment_data']
    elif seqtype == 'subject':
        rotation = math.pi / 2  # 90 degrees
        x = x + length / 2
        y = y + settings['alignment_height'] + 1
        data = settings['subject_alignment_data']
    else:
        print(f'ERROR: Invalid sequence type: {seqtype} - must be "query" or "subject"')
        sys.exit(1)

    # Build label text from configured data fields
    name = ''
    for i in range(len(data)):
        if isinstance(row[data[i]], float):
            dat = str(round(row[data[i]], 1))
        else:
            dat = str(row[data[i]])

        # Add percentage sign for identity values
        if data[i] == 'pident':
            dat = dat + '%'

        # Format label based on field position
        if i == 0:
            name = dat
        elif data[i] == 'qend' or data[i] == 'send':
            name = name + f'-{dat}'  # Range format for end coordinates
        else:
            name = name + f' {dat}'

    # Draw rotated text label
    context.move_to(x, y)
    context.select_font_face(settings['alignment_font'])
    context.set_font_size(settings['alignment_font_size'])
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.save()
    context.rotate(rotation)
    context.show_text(name)
    context.restore()

    return context


def draw_alignments(context, settings, draw_settings, feature_colors, alignments_data):
    """
    Draw all alignment blocks connecting query and subject sequences.

    This function iterates through all BLAST alignments and draws them as
    colored rectangles on both sequences, connected by lines. Each alignment
    is colored according to its percentage identity.

    Args:
        context (cairo.Context): Cairo drawing context
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for features
        alignments_data (pandas.DataFrame): Complete BLAST results data

    Returns:
        cairo.Context: Updated drawing context

    Note:
        Processes each alignment row and calls draw_alignment for individual blocks
    """
    for index, row in alignments_data.iterrows():
        context = draw_alignment(context, settings, draw_settings, row)

    return context


def draw_picture(surface, settings,  draw_settings, feature_colors,
                 query_data, subject_data, alignments_data):
    """
    Create the complete visualization by orchestrating all drawing components.

    This is the main drawing function that coordinates all visualization elements:
    - Background
    - Title
    - Query sequence line with features
    - Subject sequence line with features
    - Alignment blocks and connections
    - Color scale legend

    Args:
        surface (cairo.Surface): Cairo surface (SVG, PDF, or PNG)
        settings (dict): Visualization settings from settings.json
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for genomic features
        query_data (dict): Query sequence data and features
        subject_data (dict): Subject sequence data and features
        alignments_data (pandas.DataFrame): BLAST alignment results

    Returns:
        cairo.Surface: Completed surface with all visualization elements

    Note:
        This function sets up the Cairo context and calls all other drawing
        functions in the correct order to create the final visualization
    """
    context = cairo.Context(surface)
    background = settings['background']

    # Set background color
    context.set_source_rgba(*background)
    context.paint()

    # Draw title
    context.move_to(10, 20)
    context.select_font_face(settings['title_font'])
    context.set_font_size(settings['title_font_size'])
    context.set_source_rgb(0, 0, 0)
    context.show_text(draw_settings['title'])

    # Draw sequence lines with features
    context = draw_sequence_line(context, settings, draw_settings, feature_colors,
                                 query_data, settings['query_y_position'], 'query')
    context = draw_sequence_line(context, settings, draw_settings, feature_colors,
                                 subject_data, settings['subject_y_position'], 'subject')

    # Draw alignments connecting the sequences
    context = draw_alignments(context, settings,  draw_settings, feature_colors, alignments_data)

    # Draw color scale legend
    context = draw_scale(context, settings, draw_settings)

    return surface

def draw_svg(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.svg'):
    """
    Create SVG format visualization of sequence alignments.

    Args:
        settings (dict): Visualization settings
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for features
        query_data (dict): Query sequence data
        subject_data (dict): Subject sequence data
        alignments_data (pandas.DataFrame): BLAST alignment results
        file_out (str, optional): Output filename. Defaults to 'out_draw.svg'
    """
    surface = cairo.SVGSurface(file_out, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.finish()
    print(f'SVG visualization saved as: {file_out}')


def draw_pdf(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.pdf'):
    """
    Create PDF format visualization of sequence alignments.

    Args:
        settings (dict): Visualization settings
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for features
        query_data (dict): Query sequence data
        subject_data (dict): Subject sequence data
        alignments_data (pandas.DataFrame): BLAST alignment results
        file_out (str, optional): Output filename. Defaults to 'out_draw.pdf'
    """
    surface = cairo.PDFSurface(file_out, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.finish()
    print(f'PDF visualization saved as: {file_out}')


def draw_png(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.png'):
    """
    Create PNG format visualization of sequence alignments.

    Args:
        settings (dict): Visualization settings
        draw_settings (dict): Calculated drawing parameters
        feature_colors (dict): Color scheme for features
        query_data (dict): Query sequence data
        subject_data (dict): Subject sequence data
        alignments_data (pandas.DataFrame): BLAST alignment results
        file_out (str, optional): Output filename. Defaults to 'out_draw.png'
    """
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.write_to_png(file_out)
    print(f'PNG visualization saved as: {file_out}')


def read_settings(settings_file):
    """
    Read visualization settings from JSON configuration file.

    This function loads the settings.json file containing visualization parameters
    and color schemes for different feature types.

    Args:
        settings_file (str): Path to JSON settings file

    Returns:
        tuple: A tuple containing:
            - settings (dict): General visualization settings
            - feature_colors (dict): Color scheme for genomic features

    Raises:
        FileNotFoundError: If settings file doesn't exist
        json.JSONDecodeError: If settings file contains invalid JSON

    Note:
        Settings file should contain 'settings' and 'feature_colors' sections
    """
    with open(settings_file, 'r', encoding='utf-8') as file:
        all_settings = json.load(file)
    settings = all_settings['settings']
    feature_colors = all_settings['feature_colors']
    return settings, feature_colors


def run(query_file, subject_file, alignments_file):
    """
    Main function to generate sequence alignment visualizations.

    This function orchestrates the entire visualization pipeline:
    1. Loads settings from configuration file
    2. Prepares drawing parameters and loads data
    3. Generates visualizations in SVG, PDF, and PNG formats

    Args:
        query_file (str): Path to query sequence JSON file
        subject_file (str): Path to subject sequence JSON file
        alignments_file (str): Path to BLAST results TSV file

    Raises:
        FileNotFoundError: If any input file doesn't exist
        Various exceptions from drawing functions

    Note:
        Creates three output files: out_draw.svg, out_draw.pdf, out_draw.png
        All files are saved in the current working directory
    """
    print(f'Starting visualization pipeline...')
    print(f'Current directory: {os.getcwd()}')

    # Load visualization settings
    settings, feature_colors = read_settings('settings.json')

    # Prepare drawing parameters and load data
    draw_settings, query_data, subject_data, alignments_data = \
        prepare_draw_settings(query_file, subject_file, alignments_file, settings)

    # Generate visualizations in multiple formats
    print('Generating visualizations...')
    draw_svg(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    draw_pdf(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    draw_png(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)

    print('Visualization pipeline completed successfully!')



