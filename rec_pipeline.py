#!/usr/local/bin/python3
from Bio import SeqIO
from itertools import islice
import sys
from itertools import chain
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import csv
from collections import defaultdict
from pprint import pprint

child = sys.argv[1]
parents = sys.argv[2]
window_size = sys.argv[3]
step = sys.argv[4]
blastDB = sys.argv[5]
parent1_name = sys.argv[6]
parent2_name = sys.argv[7]


class RecAnalysis(object):
    """
    input:  rec_genome: single .fasta file with the recombinant (child) genome
            parent_genomes: multi.fasta file containing both parent genomes
            win_size: integer (desired window size Recomended: 1000)
            step_size: integer (desired window overlap size Recomended: 100)
            parent_blastDB: local blast database from both parents
            parent1: name of parent1 exactly as it is in multifasta file for blastdb without '>'
            parent1: name of parent2 exactly as it is in multifasta file for blastdb without '>'

    """

    def __init__(self, rec_genome, parent_genomes, win_size, step_size, parent_blastDB, parent1, parent2):
        self.parse_rec = self.parsed_sequence(rec_genome)
        self.parent_genomes = parent_genomes
        self.win_size = win_size
        self.step_size = step_size
        self.parent_blastDB = parent_blastDB
        self.parent1 = parent1
        self.parent2 = parent2

    @staticmethod
    def parsed_sequence(genome):
        with open(genome, "rU") as sequence:
            for record in SeqIO.parse(sequence, "fasta"):
                return {
                    'id': record.id,
                    'sequence': record.seq
                        }

    @staticmethod
    def parse_blast(blastFile):
        d = defaultdict(list)
        with open(blastFile, 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
            for row in spamreader:
                d[row[0]].append((row[1], float(row[2])))
            newDict = {}
            for k, v in d.items():
                newDict[k] = dict(v)
        return newDict

    def sliding_window(self):
        """Returns a mulit.fasta file with chunks of input fasta.  Input sequence
        must be iterable."""
        step = self.win_size - self.step_size
        id = self.parse_rec['id']
        outfile = open(id + "_window_" + str(self.win_size) + ".fasta", 'w')
        sequence = self.parse_rec['sequence']

        # Pre-compute number of chunks to emit
        numOfChunks = int(((len(sequence) - self.win_size) / step) + 1)

        # Do the work
        win_count = 0
        coord_count = 0
        for i in range(0, numOfChunks * step, step):
            win_count += 1
            coord_count += step
            print(">" + id + "_window_" + str(win_count) + "_coord_" + str(coord_count), file=outfile)
            print(sequence[i:i + self.win_size], file=outfile)

    def blast_recombinant(self):

        strainMap = {
            self.parent1: 1,
            self.parent2: 2
        }

        def parse_blast(blastFile):
            d = defaultdict(list)
            with open(blastFile, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
                for row in spamreader:
                    d[row[0]].append((row[1], float(row[2])))
                newDict = {}
                for k, v in d.items():
                    newDict[k] = dict(v)
            return newDict

        def compute_score(hitDict):

            highest = max(hitDict.values())
            highestList = [k for k,v in hitDict.items() if v == highest]
            if len(highestList) > 1:
                return 0
            else:
                return (strainMap[highestList[0]])

        id = self.parse_rec['id']
        outBlast = 'blasted_windows.tsv'
        windows_file = id + "_window_" + str(self.win_size) + ".fasta"
        blastCline = NcbiblastnCommandline(query=windows_file, db=self.parent_blastDB, evalue=0.001, outfmt=6, out=outBlast)
        print(blastCline)
        cl_command = str(blastCline).split("-")
        call = cl_command[0]
        cl_command = ["-" + x.strip() for x in cl_command[1:]]
        cl_command = [x.split(" ") for x in cl_command]
        cl_command = list(chain.from_iterable(cl_command))
        cl_call = [call.strip()]
        cl_call.extend(cl_command)
        subprocess.call(cl_call)
        hits = parse_blast(outBlast)
        outTable = open('recombinant_table.tsv', 'w')
        print('index', 'window', 'score', file=outTable)
        for k, v in hits.items():
            print(k, compute_score(v), file=outTable)

tester = RecAnalysis(rec_genome=child, parent_genomes=parents, win_size=1000, step_size=10, parent_blastDB=blastDB, parent1=parent1_name, parent2=parent2_name)
chunks = tester.sliding_window()
blast = tester.blast_recombinant()
