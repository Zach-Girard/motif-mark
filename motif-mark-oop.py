#!/usr/bin/env python

import cairo
import argparse
import re



# argparse
def parse_args():
    parser = argparse.ArgumentParser(description="Process a FASTA file and a TXT file.")
    parser.add_argument('--f', type=str, help="Path to the input FASTA file", required=True)
    parser.add_argument('--m', type=str, help="Path to the input motif TXT file", required=True)
    return parser.parse_args()

args = parse_args()





def get_colors(motifs: dict) -> dict: 
    colors: dict = dict()
    r = 0; g = 0; b = 0
    for i, m in enumerate(motifs): 
        match i % 3:
            case 0: 
                r += 1
            case 1:
                g += 1
            case 2: 
                b += 1
        if r == g and r == b and r != 0: 
                r = 0; g = 0
        color_list = [r, g, b, 0]
        colors.update({m: color_list})
    return colors


list_motifs = []
# read motif and put in dict with a color
dict_of_motif_color = {}
with open(args.m, mode ="r") as motif_file:
    for line in motif_file:
        motif_seq = line.strip()
        list_motifs.append(motif_seq)

dict_of_motif_color = get_colors(list_motifs)





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


# creating the regex equivalent of each motif
def convert_regex(motif_seq):
    regex = ""
    for char in motif_seq:
        # Look up the corresponding IUPAC regex pattern
        regex += IUPAC_REGEX.get(char)
    return regex


motif_regex = {}

for motif_seq in list_motifs:
    motif_regex[motif_seq] = convert_regex(motif_seq)






def find_motifs(record, dict_of_motif_color, regex_dict):
    motif_matches = []  # List to hold all motif match objects

    # Iterate through each motif and its regex pattern
    for motif_seq, regex_pattern in regex_dict.items():
        color = dict_of_motif_color[motif_seq]
        
        # Find all matches of the regex pattern in the sequence
        for match in re.finditer(regex_pattern, record):
            start, stop = match.start(), match.end() - 1  # Adjust for zero-based indexing
            # Create a Motif object and add it to the list
            print(f"Matched motifs: {match.group()} at position {match.start()} to {match.end()}")
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
        context.set_source_rgb(0,0,0) 
        context.rectangle(self.start+200,rank -25,self.stop-self.start,50)
        context.fill()
        context.set_source_rgba(0, 0, 0, 1)
        context.set_font_size(25) 
        context.move_to(50, rank)
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





rank = 0
list_of_geneGroups = []
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
            name = name.split('>')[1] # pull out gene name
            name = name.split(' ')[0] 

        # read each fasta record
        else:
            record = line
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







# Define the pycairo surface to draw on.
# Adjust width based on the fact seqs < 1000, height is number of (genes*100)+300
width, height = 1200, len(list_of_geneGroups)*100+300
#create the coordinates to display your graphic, desginate output
surface = cairo.SVGSurface(f"{args.f}.svg",width, height)
#create the coordinates you will be drawing on 
context = cairo.Context(surface)
context.set_source_rgb(1, 1, 1)  # RGB values for white surface
context.rectangle(0, 0, width, height)
context.fill()

  # draw plot title
context.set_source_rgba(0, 0, 0, 1)
context.set_font_size(40) 
context.move_to(500, 50)
context.show_text("Motif Plot") 

# draw all genes and names
for GeneGroup in list_of_geneGroups:
    GeneGroup.draw(context)


# Draw a legend
for i, m in enumerate(dict_of_motif_color):
    context.set_source_rgba(0, 0, 0, 1)
    context.set_font_size(20) 
    context.move_to(100, i*30 + len(list_of_geneGroups)*100 + 200)
    context.show_text(f"{m}")

    color = dict_of_motif_color[m]
    context.set_source_rgb(color[0],color[1],color[2])  
    context.set_line_width(5)
    context.move_to(20, i*30 + len(list_of_geneGroups)*100 + 200)  #(x,y)
    context.line_to(80,i*30 + len(list_of_geneGroups)*100 + 200)
    context.stroke()






surface.write_to_png(f"{args.f}.png")

