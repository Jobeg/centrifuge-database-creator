#!/usr/bin/env python3
"""
=========================================================================================================================
 centrifuge_database_creator: Create all files required to build a centrifuge database using a newick and fastas files.
 
 Required: 	- Python 3.5.3+
			- BioPython [1]

[1] http://biopython.org

 REQUIRED OPTIONS:
 Newick file. Phylogeny computed from input sequences.
 Fasta reference folder. All fastas files for reference.
 
 Version: 1.0.1 (Janv. 2018)
 Author: Johann Beghain
 Contact: johann.beghain@inserm.fr
=========================================================================================================================
"""
import argparse
import logging
import os
from Bio import SeqIO
from Bio import Phylo

###############################################################################
#############                CLASS                            #################
###############################################################################

class Taxon:
    def __init__(self, sequences):
        self.sequences = sequences

class Sequence:
    def __init__(self, length, seq):
        self.seq = seq
        self.length = length

###############################################################################
#############                FUNCTIONS                            #############
###############################################################################

# Argument selection and help
def ArgumentsParser():
  parser=argparse.ArgumentParser(description = '''centrifuge_database_creator: Create all required files for building a custom centrifuge database.''',
                                 epilog = '''Require a taxonomy in newick format.''')
  parser.add_argument('-n', '--name', type = str, required = True,
                      help = 'Name of the output database.')
  parser.add_argument('-w', '--newick', type = str, required = True,
                      help = 'Path to a valid .newick file. This is the taxonomy used for building the db. The name of the leaves are the names of the fasta files without extension.')
  parser.add_argument('-f', '--fasta', type = str, required = True,
                      help = 'Path to a folder containing all .fasta files. Extension must be .fasta')
  logger_args = parser.add_argument_group('Logger configuration')
  logger_args.add_argument( '--log_level',
                         type = str,
                         nargs = '?',
                         default = 'INFO',
                         help = 'log level',
                         choices = ['ERROR', 'error', 'WARNING', 'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                         required = False )
  logger_args.add_argument( '--log_file',
                         type = str,
                         nargs = '?',
                         help = 'log file (use the stderr by default)',
                         required = False )
  return(parser.parse_args())

def extract_sequence_from_fasta(fasta_folder_name, db_name):
    #Concatenate all fastas file from folder in a single fasta
    db_fasta = db_name+'.fasta'
    output_handle = open(db_fasta, 'w')
    #Dict of class Taxon
    taxons = dict()
    for fasta_file_name in os.listdir(fasta_folder_name):
        if fasta_file_name.endswith(".fasta"):
            taxon_id = os.path.splitext(fasta_file_name)[0]
            #Dict of class Sequence
            sequences = dict()
            for seq_record in SeqIO.parse(os.path.join(fasta_folder_name, fasta_file_name), "fasta"):
                seq_record.id = taxon_id+"_"+seq_record.id
                sequences[seq_record.id] = Sequence(len(seq_record), seq_record.seq)
                output_handle.write('>'+seq_record.id+"\n")
                seq = [ seq_record.seq[i:i+60] for i in range(0, len(seq_record), 60) ]
                for seq_block in seq:
                    output_handle.write(str(seq_block)+"\n")
            taxons[taxon_id] = Taxon(sequences)
    return(taxons)

def extract_taxonomy_from_newick(newick_file_name, taxons, db_name):
    tree = Phylo.read(newick_file_name, "newick")
    #Create nodes and names.dmp
    db_nodes = db_name+'_nodes.dmp'
    db_names = db_name+'_names.dmp'
    db_conversion = db_name+'_seqid2taxid.map'
    node_handle = open(db_nodes, 'w')
    name_handle = open(db_names, 'w')
    conversion_handle = open(db_conversion, 'w')
    #Create the root
    node_handle.write("1\t|\t1\t|\tno rank\t|\n")
    name_handle.write("1\t|\tall\t|\t|\tsynonym\t|\n1\t|\troot\t|\t|\tscientific name\t|\n")
    
    #Current taxid
    taxid = 1
    last_taxid = get_subtree(tree.clade, taxid, taxons, node_handle, name_handle, conversion_handle)
    print("Centrifuge database can be now created using the following command line:\n")
    print("centrifuge-build --conversion-table "+db_conversion+" --taxonomy-tree "+db_nodes+" --name-table "+db_names+" "+db_name+".fasta "+db_name+"\n")

def get_subtree(clades, taxid, taxons, node_handle, name_handle, conversion_handle):
    parent_taxid = taxid
    for clade in clades:
        taxid = taxid+1
        if clade.is_terminal():
            name_handle.write(str(taxid)+"\t|\t"+clade.name+"\t|\t"+clade.name.lower()+"\t|\tscientific name\t|\n")
            node_handle.write(str(taxid)+"\t|\t"+str(parent_taxid)+"\t|\tno rank\t|\n")
            #Print every sequence from taxon associated with taxid
            for sequence in taxons[clade.name].sequences:
                conversion_handle.write(sequence+"\t"+str(taxid)+"\n")
        else:
            node_name = "node"+str(taxid)
            if clade.name != None:
                node_name = clade.name
            #name the nodes by the composition of all leaves composing it
            node_name2 = node_name.lower()
            name_handle.write(str(taxid)+"\t|\t"+node_name+"\t|\t"+node_name2+"\t|\tscientific name\t|\n")
            node_handle.write(str(taxid)+"\t|\t"+str(parent_taxid)+"\t|\tno rank\t|\n")
            taxid = get_subtree(clade.clades, taxid, taxons, node_handle, name_handle, conversion_handle)
    return(taxid)

###############################################################################
#############                MAIN                                 #############
###############################################################################
if __name__ == "__main__":
    args=ArgumentsParser()
    # Logger config
    logging_std_format = '[%(levelname)s] %(message)s'
    logging_debug_format = '%(asctime)s [%(levelname)s] [%(threadName)s - %(name)s] %(message)s'
    log_level = args.log_level.upper()
    if ( log_level == 'DEBUG' ):
        logging_std_format = logging_debug_format
    logging_datefmt = '%d/%m/%Y - %H:%M:%S'
    if ( args.log_file != None ):
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             filename = args.log_file,
                             filemode = 'w',
                             level = log_level )
    else:
        logging.basicConfig( format = logging_std_format,
                             datefmt = logging_datefmt,
                             level = log_level )
    logger = logging.getLogger( 'main' )

    #Parse fasta file
    taxons = extract_sequence_from_fasta(args.fasta, args.name)
    #Parse newick file
    extract_taxonomy_from_newick(args.newick, taxons, args.name)
