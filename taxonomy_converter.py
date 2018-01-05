#!/usr/bin/env python3
"""
=========================================================================================================================
 taxonomy_converter: Create a newick file from nodes and names.dmp
 
 Required: 	- Python 3.5.3+
			- BioPython [1]

[1] http://biopython.org

 REQUIRED OPTIONS:
 nodes.dmp and names.dmp from NCBI taxonomy
 
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
    def __init__(self, name, alternative_name):
        self.name = name
        self.alternative_name = alternative_name

class Sequence:
    def __init__(self, length, seq):
        self.seq = seq
        self.length = length

###############################################################################
#############                FUNCTIONS                            #############
###############################################################################

# Argument selection and help
def ArgumentsParser():
  parser=argparse.ArgumentParser(description = '''taxonomy_converter: Create a newick file fro mnames.dmp and nodes.dmp.''',
                                 epilog = '''Secondary script.''')
  parser.add_argument('-n', '--name', type = str, required = False,
                      help = 'Name of the output database.')
  parser.add_argument('-a', '--names', type = str, required = True,
                      help = 'Path to a valid names.dmp ')
  parser.add_argument('-d', '--nodes', type = str, required = True,
                      help = 'Path to a valid nodes.dmp ')
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

def extract_taxonomy(names_file, nodes_file, db_name):
    tree = Phylo.BaseTree.Tree()
    #Record names
    names = dict()
    names_handle = open(names_file, 'r')
    for line in names_handle:
        fields = line.split("|")
        for i in range(len(fields)):
            fields[i] = "".join(fields[i].split())
        names[int(fields[0])] = Taxon(fields[1], fields[2])
    names_handle.close()
    #Construct the tree
    tree.clade.name = names[1].name
    tree.clade.branch_length = 0.0
    nodes_handle = open(nodes_file, 'r')
    for line in nodes_handle:
        fields = line.split("|")
        for i in range(len(fields)):
            fields[i] = "".join(fields[i].split())
        #Check if node doesn't already exists
        if find_in_tree(tree, names[int(fields[0])].name) != None:
            logger.warning('Duplicates taxid in nodes.dmp: '+fields[0]+' ! Skipping.')
            continue
        #Create the node
        new_node = Phylo.BaseTree.Clade(name=names[int(fields[0])].name, branch_length=1.0)
        #Find the parent object in tree
        bounding_node = find_in_tree(tree, names[int(fields[1])].name)
        bounding_node.clades.append(new_node)
    nodes_handle.close()
    #Output the tree in newick format
    Phylo.write(tree, db_name+'.nwk', "newick")
    

def find_in_tree(tree, name):
    all_corresponding_nodes = tree.find_clades(name)
    nb_of_matchs = 0
    matching_node = Phylo.BaseTree.Clade()
    for node in all_corresponding_nodes:
        nb_of_matchs = nb_of_matchs+1
        matching_node = node
    if nb_of_matchs > 1:
        logger.warning('Duplicates names in tree: '+name+'. Taking only the last one in list: '+matching_node.name)
    elif nb_of_matchs == 0:
        return(None)
    return(matching_node)

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
    output_name = 'taxonomy'
    if args.name:
        output_name = args.name
    #Parse newick file
    extract_taxonomy(args.names, args.nodes, output_name)
