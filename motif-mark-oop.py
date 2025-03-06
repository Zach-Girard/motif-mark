#!/usr/bin/env python

import cairo
import argparse
import re
import random



# argparse
def parse_args():
    parser = argparse.ArgumentParser(description="Process a FASTA file and a TXT file.")
    parser.add_argument('--f', type=str, help="Path to the input FASTA file", required=True)
    parser.add_argument('--m', type=str, help="Path to the input motif TXT file", required=True)
    return parser.parse_args()

args = parse_args()





def get_colors(motifs): 
    '''This function will be used to create a random rgb color code for each motif sequence.'''
    colors = dict()
    random.seed(5)
    for motif in motifs:
        while True:
            # Generate random RGB values between 0 and 255 and convert to decimal
            r = random.randint(0, 255) / 255
            g = random.randint(0, 255) / 255
            b = random.randint(0, 255) / 255
            
            # Ensure the color is not pure white (1, 1, 1) or pure black (0, 0, 0)
            if (r, g, b) != (1, 1, 1) and (r, g, b) != (0, 0, 0):
                break
        
        # Add the RGB values to the dictionary
        colors[motif] = [r, g, b]
    return colors



list_motifs = []
# Read motif file and put in each motif sequence in a list
dict_of_motif_color = {}
with open(args.m, mode ="r") as motif_file:
    for line in motif_file:
        motif_seq = line.strip()
        list_motifs.append(motif_seq)


# Create the rgb value for each motif
dict_of_motif_color = get_colors(list_motifs)




# IUPAC is the International Union of Pure and Applied Chemistry. This table is the shorthand codes for ambiguous bases.
IUPAC_REGEX = {'A': '[Aa]', 'a': '[Aa]',
               'C': '[Cc]', 'c': '[Cc]',
               'T': '[Tt]', 't': '[Tt]', 
               'G': '[Gg]', 'g': '[Gg]',
               'U': '[UuTt]', 'u': '[UuTt]',
               'W': '[AaTt]', 'w': '[AaTt]',
               'S': '[CcGg]', 's': '[CcGg]',
               'M': '[AaCc]', 'm': '[AaCc]',
               'K': '[GgTt]', 'k': '[GgTt]',
               'R': '[AaGg]', 'r': '[AaGg]',
               'Y': '[CcTt]', 'y': '[CcTt]',
               'B': '[CcGgTt]', 'b': '[CcGgTt]',
               'D': '[AaGgTt]', 'd': '[AaGgTt]',
               'H': '[AaCcTt]', 'h': '[AaCcTt]',
               'V': '[AaCcGg]', 'v': '[AaCcGg]',
               'N': '[AaCcTtUuGg]', 'n': '[AaCcTtUuGg]'}



def convert_regex(motif_seq):
    '''This function creates the regex equivalent of each motif. It intakes a single motif sequence and
    returns the equivalent IUPAC regex expression.'''
    regex = ""
    for char in motif_seq:
        # Look up the corresponding IUPAC regex pattern
        regex += IUPAC_REGEX.get(char)
    return regex



# Inititalize a dictionary holding motif sequences as keys and regex expressions as values
motif_regex = {}

# Do the regex conversion for each motif and add to a dictionary
for motif_seq in list_motifs:
    motif_regex[motif_seq] = convert_regex(motif_seq)





def find_motifs(record, dict_of_motif_color, regex_dict):
    '''This function will locaate every motif in the fasta record. It will then place each match into a motif object,
    with its start position, stop position, and corresponding color.'''
    motif_matches = []  # List to hold all motif match objects

    # Iterate through each motif and its regex pattern
    for motif_seq, regex_pattern in regex_dict.items():
        color = dict_of_motif_color[motif_seq]
        
        # Find all matches of the regex pattern in the sequence
        for match in re.finditer(regex_pattern, record):
            start, stop = match.start(), match.end() - 1  # Adjust for zero-based indexing
            # Create a Motif object and add it to the list
            motif_matches.append(motif(start, stop, color))

    return motif_matches






# Define Classes

# GeneGroup class
# Attributes: gene, exon, motifs, rank
# Methods: Draw
class geneGroup:
    def __init__(self, gene, exon, motifs, rank):
        self.gene = gene
        self.exon = exon
        self.motifs = motifs
        self.rank = rank*100+50

    def draw(self, context):
        self.gene.draw(context, self.rank)
        self.exon.draw(context, self.rank)
        for m in self.motifs:
            m.draw(context, self.rank)


# Gene class
# Attributes: start, stop, rank
# Methods: Draw
class gene:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop

    def draw(self, context, rank):
        context.set_source_rgb(0, 0, 0)  # black line
        context.set_line_width(2)
        context.move_to(self.start+200, rank)  #(x,y)
        context.line_to(self.stop+200,rank)
        context.stroke()



# Exon class
# Attributes: start, stop, color
# Methods: Draw
class exon:
    def __init__(self, start, stop, name):
        self.start = start
        self.stop = stop
        self.name = name

    def draw(self, context, rank):
        # Draw exon
        context.set_source_rgba(0,0,0, 0.5) 
        context.rectangle(self.start+200,rank -25,self.stop-self.start,50)
        context.fill()
        # Add gene name text next to gene
        context.set_source_rgba(0, 0, 0, 1)
        context.set_font_size(25) 
        context.move_to(50, rank+10)
        context.show_text(f"{self.name}")



# Motif class
# Attributes: start, stop, color, rank
# Methods: Draw
class motif:
    def __init__(self,start, stop, color):
        self.start = start
        self.stop = stop
        self.color = color

    def draw(self, context, rank):
        context.set_source_rgb(self.color[0], self.color[1], self.color[2]) 
        context.rectangle(self.start+200,rank-15,self.stop-self.start,30)  #height set to 30
        context.fill()





####   NEED TO CONVERT TO ONELINE FASTA FIRST
def oneline_fasta(filename: str, oneline):
    with open(filename, "r") as fasta:
        with open(oneline, "w") as oneline:
            first_time = True
            for line in fasta:
                if first_time:
                    oneline.write(line)
                    first_time = False
                elif ">" in line:
                    oneline.write(f"\n{line}")
                else:
                    oneline.write(line.strip())


oneline_fasta(args.f, "oneline.fa")





rank = 0 # The gene number in the fasta record
list_of_geneGroups = []  # Holds all GeneGroup objects created

# read through fasta file and do everything
with open("oneline.fa", mode = "r") as fasta:
    for line in fasta:
        line = line.strip()
        # if line is empty, break
        if line == "":
            break

        # read header line to get gene name
        if line.startswith(">"):
            name = line.split('\t')[0] 
            name = name.split('>')[1] 
            name = name.split(' ')[0] # pull out gene name

        # read each fasta record
        else:
            record = line

            # Increase rank with each record being processed
            rank += 1

            # define the start and stop of the gene and create gene object
            stop = len(record)
            current_Gene = gene(0, stop)

            # Loop through the string to find the first and last uppercase characters
            for i, char in enumerate(record):
                if char.isupper():
                    exon_start = i
                    break

            for i in range(len(record) - 1, -1, -1):  # Start from the end of the string and go backwards
                if record[i].isupper():
                    exon_stop = i
                    break

            current_Exon = exon(exon_start, exon_stop, name)



            # find all motifs and create motif objects
            # put motif objects into a list
    
            current_motifs = find_motifs(record, dict_of_motif_color, motif_regex)
        

            # create GeneGroup from gene, exon, motifs, rank, name
            current_GeneGroup = geneGroup(current_Gene,current_Exon,current_motifs,rank)
            # add GeneGroup to list
            list_of_geneGroups.append(current_GeneGroup)




filename = args.f
filename = filename.split(".")[0]

# Define the pycairo surface to draw on.
# Adjust width based on the fact seqs < 1000, height is sclaed by number of GeneGroups and motif sequences
width, height = 1100, (len(list_of_geneGroups)*100) + 150 + (len(list_motifs) * 30)
#create the coordinates to display your graphic, desginate output
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,width, height)
#create the coordinates you will be drawing on 
context = cairo.Context(surface)
context.set_source_rgb(1, 1, 1)  # RGB values for white surface
context.rectangle(0, 0, width, height)
context.fill()

# draw plot title
context.set_source_rgba(0, 0, 0, 1)
context.set_font_size(40) 
context.move_to(375, 50)
context.show_text(f"{filename} Motif Plot") 

# draw all genes and names
for GeneGroup in list_of_geneGroups:
    GeneGroup.draw(context)


# Draw a legend
for i, m in enumerate(dict_of_motif_color):
    # motif labels
    context.set_source_rgba(0, 0, 0, 1)
    context.set_font_size(20) 
    context.move_to(100, i*30 + len(list_of_geneGroups)*100 + 150)
    context.show_text(f"{m}")
    # Colored line for each motif
    color = dict_of_motif_color[m]
    context.set_source_rgb(color[0],color[1],color[2])  
    context.set_line_width(10)
    context.move_to(20, i*30 + (len(list_of_geneGroups)*100) + 145)  #(x,y)
    context.line_to(80,i*30 + (len(list_of_geneGroups)*100) + 145)
    context.stroke()



surface.write_to_png(f"{filename}.png")

