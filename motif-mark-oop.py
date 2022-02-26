#!/usr/bin/env python

import argparse
import os 
from math import pi
import cairo


def oneline_fasta(file) -> tuple:
    '''
    Create a fasta with the sequence in one line
    '''
    seq = ''
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                line.strip()
            seq += line


class Motif():
    def __init__(self, sequence):
        self._sequence = sequence
        self._start = -1
        self._length = len(self._sequence)


    def __len__(self):
        return self._length


class Gene():
    def __init__(self, name, sequence):
        self._name = name
        self._sequence = sequence
        self._length = len(self._sequence)


    def __len__(self):
        return self._length


    def find(motif):
        '''
        Use gibbs sampling to find the motifs in the sequence
        '''
        pass


class Canvas():
    def __init__(self, gene, layout = (720, 540)):
        self._surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, *shape)
        self._context = cairo.Context(_surface)
        self._offset = offset
        self._layout = layout


    def draw_rectangle_rounded_corner(dim = (720, 540), offset = 50, rgb = (1, 1, 1), width = 3, radius = 75):
        self._context.set_line_width(width)
        self._context.set_source_rgb(*rgb)

        h, w = dim
        # an area with coordinates of
        # (top, bottom, left, right) edges in absolute coordinates:
        a,b,c,d = (offset, w - offset, offset, h - offset)
        self._context.arc(a + radius, c + radius, radius, 2*(pi/2), 3*(pi/2))
        self._context.arc(b - radius, c + radius, radius, 3*(pi/2), 4*(pi/2))
        self._context.arc(b - radius, d - radius, radius, 0*(pi/2), 1*(pi/2)) 
        self._context.arc(a + radius, d - radius, radius, 1*(pi/2), 2*(pi/2))
        self._context.close_path()
        self._context.stroke()


    def save(filename):
        self._surface.write_to_png(f"{filename}.png")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="This program finds motifs and plots them.")
    parser.add_argument("-f", "--fasta-file", type=argparse.FileType('r'), action='store', help="fasta file", dest="FASTA", required=True)
    parser.add_argument("-m", "--motifs-file", type=argparse.FileType('r'), action='store', help="motifs file", dest="MOTIFS", required=True)
  
    options = parser.parse_args()
    name, ext = os.path.splitext(options.FASTA.name)
    if ext != ".fa" or ext != ".fasta":
        raise ValueError("Argument for fasta parameter is not a fasta file (*.fa, *.fasta)")

    OUTPUT_PNG = f"{os.getcwd()}/{name}.png"

    ''' Read in the motifs '''

    motifs = {}
    with open(options.MOTIFS) as motif_file: 
        for motif in motif_file:
            motif = motif.strip()
            motifs[motif] = Motif(motif)

    ''' Read in each gene '''

    with open(options.FASTA) as fasta_file:
        fasta = oneline_fasta(file, output):
        header, sequence = fasta.split('\n')
        gene_name = header[1:header.find(' ')]
        gene = Gene(gene_name, sequence)

