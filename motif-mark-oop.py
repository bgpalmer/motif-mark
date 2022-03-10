#!/usr/bin/env python

import argparse
import os 
from math import pi
import cairo
import re
import string
from itertools import chain
import random


def oneline_fasta(fasta_fp) -> str:
    '''
    Create a fasta with the sequence in one line
    '''
    seq = ''
    for line in fasta_fp:
        if len(seq) > 0:
            if line.startswith('>'): 
                yield seq
                seq = ''
            else: 
                line = line.strip()
        seq += line
    
    yield seq


class NucleicAcidSequence():
    def __init__(self, sequence):
        self._sequence = sequence
        self._length = len(self._sequence)


    def __len__(self):
        return self._length


    def __getitem__(self, key):
        return self._sequence[key]


    def __str__(self):
        return self._sequence


class Motif(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence)


class IUPAC_Matcher():
    def __init__(self):
        self.__rules = {
            'A': 'A',  # Adenine
            'C': 'C',  # Cytosine
            'G': 'G',  # Guanine
            'T': 'T',  # Thymine
            'U': 'U',  # Uracil
            'W': '[A|T]', # Weak
            'S': '[C|G]', # Strong
            'M': '[A|C]', # Amino
            'K': '[G|T]', # Ketone
            'R': '[A|G]', # Purine
            'Y': '[C|T]', # Pyrimidine
            'B': '[C|G|T]', # Not A
            'D': '[A|G|T]', # Not C
            'H': '[A|C|T]', # Not G
            'V': '[A|C|G]', # Not T
            'N': '[A|C|G|T]', # Any one base
        }

    def match(self, motif, sequence):
        return re.match(''.join(self.__rules[x.upper()] for x in motif), sequence.upper()) 


class Gene(NucleicAcidSequence):
    def __init__(self, name, sequence):
        super().__init__(sequence)
        self._name = name
        self._exons = re.split(f"[{string.ascii_lowercase}]+", self._sequence)
        self._introns = re.split(f"[{string.ascii_uppercase}]+", self._sequence)
        self._splices = list(chain.from_iterable(zip(self._exons, self._introns)))
        self._splices = list(filter(None, self._splices))
        self._motifs = []
        self._overlaps = None  


    def __str__(self):
        return self._name


    def __nkmers(self, motif):
        return len(self) - len(motif) + 1


    def find_motifs(self, motifs, searcher):
        self._motifs = [ list() for _ in range(len(motifs)) ]
        for idx, motif in enumerate(motifs):
            for start in range(self.__nkmers(motif)):             
                sequence = gene[start:len(motif) + start]
                match = searcher.match(motif, sequence)
                if match: self._motifs[idx].append((start, motif, sequence))

        return self._motifs


class MotifMarkFigure():
    def __init__(self, title, genes, motifs, padding = (0, 0), dims = (720, 540)):
        self._surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, *dims)
        self._title = title
        self._dims = dims
        self._offset = (padding, padding) # ((x pad, y pad), (width pad, height pad))
        self._xy_scale = (1., 1.)
        self._plots = []
        self._plot_area = None
        self._n_genes = len(genes)
        self._n_motifs = len(motifs)
        self._motif_names = [str(s) for s in motifs]

        self._colors = ()
        for _ in range(self._n_motifs):
            motif_color = ()
            for _ in range(3):
                motif_color += (random.uniform(0, 1),)
            self._colors += (motif_color,)


    def append_gene_plot(self, gene, motifs):

        # determine if we will need to scale x dim
        x = 1. if len(gene) <= self._dims[0] else self._dims[0] / len(gene)
        self._xy_scale = (min(self._xy_scale[0], x), 1.)  # (gene length closest to figure width, height unchanged)
        (x_pad, y_pad), (width_pad, height_pad) = self._offset
        rect = None
        if 0 == len(self._plots):
            # rect = (x_pad, y_pad, len(gene) - width_pad, y_pad + (self._dims[1] / (1 + self._n_genes)))
            rect = (x_pad, y_pad, len(gene), y_pad + (self._dims[1] / (1 + self._n_genes)))

        else:
            prev_area, _, _ = self._plots[-1]
            # rect = (x_pad, prev_area[3], len(gene) - width_pad, prev_area[3] + (self._dims[1] / (1 + self._n_genes)))
            rect = (x_pad, prev_area[3], len(gene), prev_area[3] + (self._dims[1] / (1 + self._n_genes)))

        self._plots.append((rect, gene, motifs))

    def draw(self):
        ctx = cairo.Context(self._surface)
        for rect, gene, motifs in self._plots:
            x, y = (rect[0]), (rect[1])
            # Draw text 

            ctx.set_source_rgb(*(.1, .1, .1))
            ctx.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(13)
            (_, _, _, text_height, _, _) = ctx.text_extents(f"{str(gene)}")
            ctx.move_to(x, y + text_height)
            ctx.show_text(f"{str(gene)}")
            ctx.stroke()
            ctx.save()

            # Draw gene

            rgb = (.3, .4, .5)
            line_width = 3.
            ctx.set_source_rgb(*rgb)
            ctx.set_line_width(line_width)
            ctx.scale(*self._xy_scale)
            # an area with coordinates of
            # (left, right, top, bottom) edges in absolute coordinates: 
            gene_left, gene_top = x, y + text_height
            (x, y), (x1, y1) = (gene_left, gene_top), (rect[2], rect[3])
            r = 10

            ctx.arc(x + r, y + r, r, pi, 3 * pi / 2)   # note: pi is where I pictured 0 would be?
            ctx.arc(x1 - r, y + r, r, 3 * pi / 2, 0)
            ctx.arc(x1 - r, y1 - r, r, 0, pi / 2)
            ctx.arc(x + r, y1 - r, r, pi / 2, pi)
            ctx.close_path() # adds line segment to back to beginning
            ctx.stroke()
            
            # Draw Exons & Introns

            last_x = x
            for splice in gene._splices:
                ctx.arc(last_x + r, y + r, r, pi, 3 * pi / 2)   # note: pi is where i pictured 0 would be?
                ctx.arc((last_x + len(splice)) - r, y + r, r, 3 * pi / 2, 0)
                ctx.arc((last_x + len(splice)) - r, y1 - r, r, 0, pi / 2)
                ctx.arc(last_x + r, y1 - r, r, pi / 2, pi)
                ctx.close_path()
                rgb = (.6, .2, .2) if splice.isupper() else (.2, .6, .2)
                ctx.set_source_rgb(*rgb)
                ctx.fill()
                ctx.stroke()
                last_x += len(splice)

            # Draw Motifs

            gene_height = y1 - y
            motif_height = gene_height / len(motifs)
            motif_y = (motif_height / 2.) + gene_top
            ctx.set_source_rgb(*(.1, .1, .1))
            ctx.set_line_width(line_width * 3.)
            color_idx = 0
            for motif in motifs:
                for start, m, _ in motif:
                    motif_x = rect[0] + start
                    ctx.move_to(motif_x, motif_y)
                    ctx.line_to(motif_x + len(m), motif_y) 
                    ctx.set_source_rgb(*self._colors[color_idx])
                    ctx.stroke()
                color_idx += 1
                motif_y += motif_height

            ctx.restore()

        # Legend

        color_idx = 0
        sq = 20
        x, y, x1, y1 = 0, self._dims[1] - sq, sq, self._dims[1]
        for splice, color in (('exon', (.6, .2, .2)), ('intron', (.2, .6, .2))):
            ctx.arc(x + r, y + r, r, pi, 3 * pi / 2)   # note: pi is where I pictured 0 would be?
            ctx.arc(x1 - r, y + r, r, 3 * pi / 2, 0)
            ctx.arc(x1 - r, y1 - r, r, 0, pi / 2)
            ctx.arc(x + r, y1 - r, r, pi / 2, pi)
            ctx.set_source_rgb(*color)
            ctx.fill()
            ctx.stroke()

            ctx.set_source_rgb(*(.1, .1, .1))
            ctx.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(13)
            (_, _, text_width, text_height, _, _) = ctx.text_extents(f"{splice}")

            ctx.move_to(x1, (y1 - ((y1 - (y1 - text_height)) / 2.)))
            ctx.show_text(f"{splice}")
            ctx.stroke()

            x += text_width + (2*r)
            x1 += text_width + (2*r) 

        for motif in self._motif_names:
            ctx.arc(x + r, y + r, r, pi, 3 * pi / 2)   # note: pi is where I pictured 0 would be?
            ctx.arc(x1 - r, y + r, r, 3 * pi / 2, 0)
            ctx.arc(x1 - r, y1 - r, r, 0, pi / 2)
            ctx.arc(x + r, y1 - r, r, pi / 2, pi)
            ctx.set_source_rgb(*self._colors[color_idx])
            ctx.fill()
            ctx.stroke()
            color_idx += 1

            ctx.set_source_rgb(*(.1, .1, .1))
            ctx.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(13)
            (_, _, text_width, text_height, _, _) = ctx.text_extents(f"{motif}")

            ctx.move_to(x1, (y1 - ((y1 - (y1 - text_height)) / 2.)))
            ctx.show_text(f"{motif}")
            ctx.stroke()

            x += text_width + (2*r)
            x1 += text_width + (2*r) 

    def __getitem__(self, idx):
        return self._plots[idx]


    def save(self, filename):
        self._surface.write_to_png(f"{filename}.png")


if __name__ == "__main__":
    random.seed(7)

    parser = argparse.ArgumentParser(description="This program finds motifs and plots them.")
    parser.add_argument("-f", "--fasta-file", type=argparse.FileType('r'), action='store', help="fasta file", dest="FASTA", required=True)
    parser.add_argument("-m", "--motifs-file", type=argparse.FileType('r'), action='store', help="motifs file", dest="MOTIFS", required=True)
  
    options = parser.parse_args()
    name, ext = os.path.splitext(options.FASTA.name)
    if ext != ".fa" and ext != ".fasta":
        raise ValueError("Argument for fasta parameter is not a fasta file (*.fa, *.fasta)")

    OUTPUT_PNG = f"{os.getcwd()}/{name}.png"

    ''' Read in the motifs '''

    motifs = [ Motif(motif.strip()) for motif in options.MOTIFS ]

    ''' Read in each gene '''

    genes = []
    for fasta in oneline_fasta(options.FASTA):
        header, sequence = fasta.split('\n')
        gene_name = header[1:header.find(' ')]
        genes.append(Gene(gene_name, sequence))

    ''' Make Figure '''

    figure = MotifMarkFigure('Main', genes, motifs)
    matcher = IUPAC_Matcher()

    matches = []
    for gene in genes:
        matches = gene.find_motifs(motifs, matcher)
        figure.append_gene_plot(gene, matches)

    figure.draw()
    figure.save('output')

