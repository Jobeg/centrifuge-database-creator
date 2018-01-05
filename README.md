# centrifuge-database-creator
Create all files required to build a centrifuge database using a newick and fastas files.

The creation of a custom centrifuge database can be painful. Indeed, you have to add all the sequences you need, sometimes you will need to find a new taxid, sometimes you will remove some species from your db, etc. Working with both names.dmp, nodes.dmp and seqid2taxid.map files can be confusing. For this reason I created this script that works with a single newick file for the taxonomy.


# Required
* Python 3.5.3 or later
* Biopython (http://biopython.org)


# Usage

```
% centrifuge_database_creator.py
optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Name of the output database.
  -w NEWICK, --newick NEWICK
                        Path to a valid .newick file. This is the taxonomy
                        used for building the db. The name of the leaves are
                        the names of the fasta files without extension.
  -f FASTA, --fasta FASTA
                        Path to a folder containing all .fasta files.
                        Extension must be .fasta
```



# Licence
