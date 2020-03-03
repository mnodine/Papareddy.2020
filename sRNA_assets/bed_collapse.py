import os
import re
import sys
import argparse
import fasta_utils as fu

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    dest='BED',
    help="Input BED file.",
    type=str,
    nargs='?'
)

args = parser.parse_args()

#####################
# ENVIRONMENT SETUP #
#####################

current_chrom = None
current_start = None
class ReadStack():
    '''Temporary stack of reads for compiling identical species'''
    def __init__(self,chrom=None,pos=None,strand=None,species={}):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand
        self.species = species
    
    def __repr__(self):
        return '<{}:{}>'.format(
            self.chrom,
            self.pos
        )
    
    def populate(self,strand,endpos,string,value):
        '''Add a read to the species dict'''
        if strand not in self.species:
            self.species[strand] = {}
        if endpos not in self.species[strand]:
            self.species[strand][endpos] = {}
        self.species[strand][endpos][string] = \
            self.species[strand][endpos].get(string, float(0)) + value
    
    def dump(self,counter_plus,counter_minus,digits = 5):
        '''Dump the species dict to stdout in a sorted BED format'''
        new_plus = counter_plus
        new_minus = counter_minus
        # Print all species as BED lines
        if '+' in self.species:
            entries = sorted([(k,v) for k,v in self.species['+'].items()])
            for pos,posdict in entries:
                species = sorted([(k,v) for k,v in posdict.items()])
                for string,score in species:
                    head,sequence,tail,length = string.split(':')
                    new_plus += 1
                    print(
                        '\t'.join(
                            [
                                str(i) for i in [
                                    self.chrom,
                                    self.pos,
                                    pos,
                                    'p.{}'.format(new_plus),
                                    round(score,digits),
                                    '+',
                                    head,
                                    sequence,
                                    tail,
                                    length
                                ]
                            ]
                        )
                    )
        
        if '-' in self.species:
            entries = sorted([(k,v) for k,v in self.species['-'].items()])
            for pos,posdict in entries:
                species = sorted([(k,v) for k,v in posdict.items()])
                for string,score in species:
                    head,sequence,tail,length = string.split(':')
                    new_minus += 1
                    print(
                        '\t'.join(
                            [
                                str(i) for i in [
                                    self.chrom,
                                    self.pos,
                                    pos,
                                    'm.{}'.format(new_minus),
                                    round(score,digits),
                                    '-',
                                    head,
                                    sequence,
                                    tail,
                                    length
                                ]
                            ]
                        )
                    )
        
        return new_plus,new_minus
    

####################
# PARSE INPUT FILE #
####################
if args.BED:
    if args.BED.split('.')[-1].lower() != 'bed':
        print("\nERROR: input file must be BED format.")
        parser.print_help()
        sys.exit(1)
    file = open(args.BED)
elif not sys.stdin.isatty():
    file = sys.stdin
else:
    print("\nERROR: requires BED file as input.")
    parser.print_help()
    sys.exit(1)

current_stack = None
dump_stack = False
counter_plus = 0
counter_minus = 0

# chrom start end name reads strand	head sequence tail length
for line in file:
    if line[0] == '#':
        continue
    chrom,start,end,name,reads,strand,head,sequence,tail,length = line.rstrip().split('\t')
    string = ':'.join([head,sequence,tail,length])
    
    if not current_stack:
        current_stack = ReadStack()
        current_stack.chrom = chrom
        current_stack.pos = start
        current_stack.species = {}
    
    if chrom == current_stack.chrom and start == current_stack.pos:
        # If the same position, add the current read to the stack
        current_stack.populate(
            strand,
            end,
            string,
            float(reads)
        )
    else:
        dump_stack = True
    
    if dump_stack:
        # Output all reads in the stack and start a new one
        dump_stack = False
        counter_plus,counter_minus = current_stack.dump(counter_plus,counter_minus)
        current_stack = ReadStack()
        current_stack.chrom = chrom
        current_stack.pos = start
        current_stack.species = {}
        current_stack.populate(
            strand,
            end,
            string,
            float(reads)
        )

file.close()