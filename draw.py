import json
import sys
from copy import copy
import os
import pandas as pd
import math
import cairo
import matplotlib.pyplot as plt


def generate_gradient(x):
    # Sprawdzenie poprawności wartości x
    x = int(x)
    if x < 0 or x > 100:
        raise ValueError("Min. vale must be between 0 and 100.")

    # Tworzymy listę liczb od x do 100
    keys = list(range(x, 101))

    # Liczba kolorów do wygenerowania
    num_colors = len(keys)

    # Używamy skali barw "rainbow" z matplotlib
    cmap = plt.cm.get_cmap("rainbow", num_colors)

    # Generujemy kolory w formacie RGB
    colors = [tuple(round(c, 2) for c in cmap(i)[:3]) for i in range(num_colors)]

    # Tworzymy słownik
    color_dict = {key: color for key, color in zip(keys, colors)}
    #print(f'{color_dict = }')
    return color_dict

"""
def visualize_color_scale(color_dict):
    fig, ax = plt.subplots(figsize=(10, 1))
    colors = [color for color in color_dict.values()]
    ax.imshow([colors], extent=[min(color_dict.keys()), max(color_dict.keys()), 0, 1])
    ax.set_yticks([])
    ax.set_xticks(list(color_dict.keys())[::5])  # Co 5 kluczy
    plt.show()
"""


def read_data(input_file):
    print(f'Reading json data from input file: {input_file}')
    with open(input_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    #print(data)
    return data

def prepare_draw_settings(query_file, subject_file, alignments_file, settings):
    query_data = read_data(query_file)
    subject_data = read_data(subject_file)
    alignments_data = pd.read_csv(alignments_file, sep='\t')
    query_length, subject_length = 0, 0
    if query_data is not None:
        query_length = query_data["description"]["seq_description"]["seq_len"]
        print(f'Query length: {query_length}')
    if subject_data is not None:
        subject_length = subject_data["description"]["seq_description"]["seq_len"]
        print(f'Subject length: {subject_length}')

    left_move = round(abs(subject_length - query_length) / 2)

    if query_length < subject_length:
        shorter = 'query'
        left_move_query = left_move
        left_move_subject = 0
    else:
        shorter = 'subject'
        left_move_query = 0
        left_move_subject = left_move

    min_pident = round(min(alignments_data['pident']))
    max_pident = round(max(alignments_data['pident']))

    print(f'>>> {min_pident = } - {max_pident = }')

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
    draw_settings['draw_width'] = (settings['left_margin'] + math.ceil(draw_settings['main_length'] / settings['ratio']) +
                   settings['right_margin'])
    draw_settings['draw_height'] = settings['top_margin'] + settings['main_height'] + settings['right_margin']

    return draw_settings, query_data, subject_data, alignments_data

def draw_sequence_line(context, settings,  draw_settings, feature_colors, seq_data, y_position, seq_type):
    print('------')
    original_left_margin = copy(settings['left_margin'])

    if draw_settings['shorter'] == 'query' and seq_type == 'query':
        settings['left_margin'] += draw_settings['left_move_query']
        #print(
        #    f'{seq_type = } {original_left_margin = } {settings['left_margin'] = } {draw_settings['left_move_query'] = }')
    elif draw_settings['shorter'] == 'subject' and seq_type == 'subject':
        settings['left_margin'] += draw_settings['left_move_subject']
        #print(
        #    f'{seq_type = } {original_left_margin = } {settings['left_margin'] = } {draw_settings['left_move_subject'] = }')
    context.move_to(settings['left_margin']-45, y_position)
    context.select_font_face('Sans-serif')
    context.set_font_size(12)
    context.set_source_rgb(0, 0, 0)
    context.show_text(seq_type)
   # print(f'>>>> {draw_settings['shorter'] = }  {seq_type = } {settings['left_margin'] = }')
    context.select_font_face(settings['tick_font'])
    context.set_font_size(settings['tick_font_size'])
    context.move_to(settings['left_margin'], y_position)
    length = math.ceil(seq_data['description']['seq_description']['seq_len']/settings['ratio'])
    context.line_to(length+settings['left_margin'],
                    y_position)
    context.set_source_rgb(0.2, 0.2, 0.2)
    context.set_line_width(2)
    context.stroke()

    # Drawing ticks
    for x in range(0, length, math.ceil(settings['tick']/settings['ratio'])):
        x_pos = x + settings['left_margin']
        context.move_to(x_pos, y_position)
        context.line_to(x_pos,
                        y_position+settings['tick_length'])
        context.set_source_rgb(0.2, 0.2, 0.2)
        context.set_line_width(settings['tick_line_width'])
        context.stroke()
        context.move_to(x_pos, y_position+settings['tick_length']+5)
        context.set_source_rgb(0.7, 0, 0)
        context.save()
        context.rotate(math.pi / 2)
        context.show_text(str(x*settings['ratio']))
        context.restore()
    draw_features(context, settings, draw_settings, feature_colors, seq_data, y_position, seq_type)
    settings['left_margin'] = original_left_margin
    return context


def draw_feature(context, settings, draw_settings, feature_colors, feature, y_position, i, seq_type, special=None):

    i = i * settings['feature_height']/2
    if seq_type == 'subject':
        feature_margin = settings['feature_margin']
        rotation = math.pi / 2
        text_move_y = settings['feature_height'] + 5
    elif seq_type == 'query':
        feature_margin = -(settings['feature_margin'] + settings['feature_height'])
        rotation = (3*math.pi) / 2
        text_move_y = -5
        i = -i
    else:
        print(f'Wrong type of sequence: {seq_type}')
        sys.exit(1)
    if special in settings['special_features']:
        feature_margin *= 3
    #print(f">> {feature = }")
    #print(f"Range: {feature['end'] = } - {feature['start'] = }")
    length = math.ceil((abs(feature['end'] - feature['start']) / settings['ratio']))
    coordinates = (
        settings['left_margin'] + math.ceil(feature['start']/settings['ratio']),
        y_position + feature_margin + i,
        length,
        settings['feature_height']
    )
    context.rectangle(*coordinates)
    if special is not None:
        if feature_colors.get(special) is not None:
            color = feature_colors.get(special)
        else:
            color = feature_colors[feature['type']]
    else:
        color = feature_colors[feature['type']]
    context.set_source_rgba(*color, feature_colors['alpha'])
    context.fill()
    # Description:

    context.move_to(settings['left_margin'] + math.ceil(feature['start']/settings['ratio']) + length/2,
                    y_position + feature_margin + text_move_y + i)
    context.set_source_rgba(*color, feature_colors['alpha'])
    context.save()
    context.rotate(rotation)
    name = ''
    for d in settings['feature_data']:
        if d == 'end':
            name += "-" + str(feature[d])
        else:
            name += " " + str(feature[d])

    context.show_text(name)
    context.restore()

    return context

def draw_features(context, settings,  draw_settings, feature_colors, seq_data, y_position, seq_type):
    features = seq_data['features']
    #print(f'{features = }')
    n = 1
    for key, feature in features.items():
        i = n%2
        n += 1
        #print(f'{feature["type"] = } {feature["name"] = }')
        if feature['type'] == 'gene':
            if feature.get('CDS') is not None:
                cdss = feature.get('CDS')
                for key2, cds in cdss.items():
                    #print(f'\t{cds = }')

                    if len(cds['exons']) > 0:
                        exons = cds['exons']
                        #print(f'\t\t{len(exons)} {exons = }')
                        for key3, exon in exons.items():
                            #print(f'\t\t\t{key3} {exon = }')
                            draw_feature(context, settings, draw_settings, feature_colors, exon, y_position, i, seq_type)

                    elif len(cds['introns']) > 0:
                        introns = cds['introns']
                        #print(f'\t\t{len(introns)} {introns = }')
                        for key4, intron in introns.items():
                            draw_feature(context, settings, draw_settings, feature_colors, intron, y_position, i, seq_type)
                    else:
                        #print(f'Only cds: {cds = }')
                        draw_feature(context, settings, draw_settings, feature_colors, cds, y_position, i, seq_type)
            elif feature.get('tRNA') is not None:
                tRNA = feature.get('tRNA')
                #print(f'{tRNA = }')
                draw_feature(context, settings, draw_settings, feature_colors, tRNA, y_position, i, seq_type)
            elif feature.get('rRNA') is not None:
                rRNA = feature.get('rRNA')
                #print(f'{rRNA = }')
                draw_feature(context, settings, draw_settings, feature_colors, rRNA, y_position, i, seq_type)
            else:
                print(f'!!!! Gene: {feature = }')
                #draw_feature(context, settings, draw_settings, feature_colors, feature, y_position, i, seq_type)

        elif feature['type'] == 'misc_feature':
           # print(f'Misc derived: {feature}')
            if feature['plastid-derived']:
               # print('plastid derived!')
                draw_feature(context, settings, draw_settings, feature_colors, feature, y_position, i, seq_type,
                'plastid-derived')
            else:
                draw_feature(context, settings, draw_settings, feature_colors, feature, y_position, i, seq_type)
    return context

def draw_scale(context, settings, draw_settings):
    x_pos = draw_settings['draw_width'] - settings['scale_move_left']
    y_pos = settings['scale_move_down']
    gradient = draw_settings['gradient']
    context.select_font_face(settings['scale_font'])
    context.set_font_size(settings['scale_font_size'])
    j = 0
    for i in range(settings['min_scale'], 101):
        y = y_pos + j * settings['scale_step_height'] + 5 * j
        coordinates = (
            x_pos,
            y,
            settings['scale_step_width'],
            settings['scale_step_height']
        )
        color = gradient[i]
        context.rectangle(*coordinates)
        context.set_source_rgba(*color, settings['alignment_alpha'])
        context.fill()
        context.move_to(x_pos + settings['scale_step_width'] + 5, y+settings['scale_step_height'])
        context.set_source_rgb(*settings['scale_font_color'])
        context.show_text(str(i))
        j += 1


    return context



def draw_alignment(context, settings, draw_settings, row):

    if draw_settings['shorter'] == 'query':
        query_margin =  settings['left_margin'] + draw_settings['left_move_query']
        subject_margin = settings['left_margin']
    elif draw_settings['shorter'] == 'subject':
        subject_margin =  settings['left_margin'] + draw_settings['left_move_subject']
        query_margin = settings['left_margin']
    else:
        query_margin = settings['left_margin']
        subject_margin = settings['left_margin']

    y_query = settings['query_y_position'] - settings['alignment_height']/2
    y_subject = settings['subject_y_position'] - settings['alignment_height']/2
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
    x_start_query = query_margin + math.ceil(qstart / settings['ratio'])
    x_start_subject = subject_margin + math.ceil(sstart / settings['ratio'])


    length_query = math.ceil((abs(qend - qstart) / settings['ratio']))
    length_subject = math.ceil((abs(send - sstart) / settings['ratio']))
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
    color = draw_settings['gradient'][round(row['pident'])]

    context.rectangle(*query_coordinates)
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.fill()
    put_alignment_description(context, settings, row, "query", x_start_query, y_query, color, length_query)

    context.rectangle(*subject_coordinates)
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.fill()
    put_alignment_description(context, settings, row, "subject", x_start_subject, y_subject, color, length_subject)

    # Draw line between query and subject alignment
    context.move_to(x_start_query,  y_query + settings['alignment_height'])
    context.line_to(x_start_subject, y_subject,)
    context.set_source_rgba(*color, settings['alignment_line_alpha'])
    context.set_line_width(settings['alignment_line_width'])
    context.stroke()
    return context

def put_alignment_description(context, settings, row, seqtype, x, y, color, length):
    if seqtype == 'query':
        rotation = (3 * math.pi) / 2
        x=x+length/2
        data = settings['query_alignment_data']
    elif seqtype == 'subject':
        rotation = math.pi / 2
        x = x + length / 2
        y = y + settings['alignment_height']+1
        data = settings['subject_alignment_data']
    else:
        print(f'Wrong type of sequence: {seqtype} - must be "query" or "subject"')
        sys.exit(1)
    name = ''

    for i in range(len(data)):
        if isinstance(row[data[i]], float):
            dat = str(round(row[data[i]],1))
        else:
            dat = str(row[data[i]])
        if data[i] == 'pident':
            dat =  dat + '%'
        if i == 0:
            name = dat
        elif data[i] == 'qend' or data[i]  == 'send':
            name = name + f'-{dat}'
        else:
            name = name +  f' {dat}'

    context.move_to(x, y)
    context.select_font_face(settings['alignment_font'])
    context.set_font_size(settings['alignment_font_size'])
    context.set_source_rgba(*color, settings['alignment_alpha'])
    context.save()
    context.rotate(rotation)

    context.show_text(name)
    context.restore()
    return context



def draw_alignments(context, settings,  draw_settings, feature_colors, alignments_data):
    """ Draws alignments on lines representing query and subject sequences"""
    for index, row in alignments_data.iterrows():
        context = draw_alignment(context, settings, draw_settings, row)

    return context

def draw_picture(surface, settings,  draw_settings, feature_colors,
                 query_data, subject_data, alignments_data):
    context = cairo.Context(surface)
    background = settings['background']
    #print(f'{background = }')

    context.set_source_rgba(*background)
    context.paint()

    context.move_to(10, 20)
    context.select_font_face(settings['title_font'])
    context.set_font_size(settings['title_font_size'])
    context.set_source_rgb(0, 0, 0)
    context.show_text(draw_settings['title'])
    # Line representing query genome

    context = draw_sequence_line(context, settings, draw_settings, feature_colors,
                                 query_data, settings['query_y_position'], 'query')
    context = draw_sequence_line(context, settings, draw_settings, feature_colors,
                                 subject_data, settings['subject_y_position'], 'subject')
    context = draw_alignments(context, settings,  draw_settings, feature_colors, alignments_data)
    context = draw_scale(context, settings, draw_settings)
    return surface

def draw_svg(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.svg'):
    surface = cairo.SVGSurface(file_out, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.finish()

def draw_pdf(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.pdf'):
    surface = cairo.PDFSurface(file_out, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.finish()

def draw_png(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data,
             file_out='out_draw.png'):
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, draw_settings['draw_width'], draw_settings['draw_height'])
    surface = draw_picture(surface, settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    surface.write_to_png(file_out)

def read_settings(settings_file):
    """ Reads settings from .json file"""
    with open(settings_file, 'r', encoding='utf-8') as file:
        all_settings = json.load(file)
    settings = all_settings['settings']
    feature_colors = all_settings['feature_colors']
    return settings, feature_colors

def run(query_file, subject_file, alignments_file):
    print(f'Current dir: {os.getcwd()}')
    #os.chdir(output_dir)
    settings, feature_colors = read_settings('settings.json')
    draw_settings, query_data, subject_data, alignments_data =\
        prepare_draw_settings(query_file, subject_file, alignments_file, settings)

    #print(f'SETTINGS LOADED:\n\tsettings:{settings}\n\tfeature_colors: {feature_colors}\n\tdraw settings: {draw_settings}')
    draw_svg(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    draw_pdf(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)
    draw_png(settings, draw_settings, feature_colors, query_data, subject_data, alignments_data)

def main():
    run('results/Cuscuta_pedicellata-NC_052872.json',
                              'results/Cuscuta_epilinum-BK059237.json',
                              'results/blast_results.tsv')

if __name__ == "__main__":
  main()

