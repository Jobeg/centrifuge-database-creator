# centrifuge-database-creator
Create all files required to build a centrifuge database using a newick and fastas files.

If you don't know what I am talking about you can read the centrifuge documentation on [creating a custom dabase with centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database).

/!\ This script is for **centrifuge-1.0.3-beta**. I don't know if in the future there will be another way to create a centrifuge database.

The creation of a custom centrifuge database can be painful. Indeed, you have to add all the sequences you need in at least 3  files, sometimes you will need to find a new taxid, sometimes you will remove some species from your db, etc. Working with both names.dmp, nodes.dmp and seqid2taxid.map files can be confusing. For this reason I created this script that works with a single newick file for the taxonomy.

## What is this script for?

The main program is centrifuge_database_creator. It takes the wanted taxonomy in a newick format, a bunch of fasta files (for your database) and creates for you the corresponding names.dmp, nodes.dmp and seqid2taxid.map required for centrifuge-build.

## Use-case: ok how can I use this?

### Exemple 1: adding new sequences from official taxonomy

The centrifuge-download gave you the official taxonomy from NCBI (`centrifuge-download taxonomy`). You have the names.dmp and nodes.dmp . You can choose to edit theses files manually or to usethe secondary script-taxonomy converter.py : it will make the opposite of the main script, taking a names.dmp and nodes.dmp and generate a newick format.

From your newick format you can edit using a tree editor (I recommand [treegraph2](http://treegraph.bioinfweb.info/) ) to add your custom nodes, ie your own sequences. Add in the tree the name of the fastas files without the .fasta extension.
ex: you want to add the organism myVirus-strainB722.fasta, you add a node myVirus-strainB722 in the tree under the "Virus" clade.

With your newick file and fastas files you can now use the centrifuge database creator script to generate the names.dmp, nodes.dmp and seqid2taxid.map .
And finally youy can build your database following the correct command line in [centrifuge documentation](https://ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database)

### Exemple 2: Fully custom taxonomy

You can build your taxonomy in newick format using any tool you want. For exemple, I did a phylogeny on the core-genome from a set of 100 bacterias.

With your newick file and fastas files you can now use the centrifuge database creator script to generate the names.dmp, nodes.dmp and seqid2taxid.map .
And finally youy can build your database following the correct command line in [centrifuge documentation](https://ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database)

# Required
* Python 3.5.3 or later
* Biopython (http://biopython.org)


# Usage

For the main script:

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

For the secondary script:

```
% taxonomy_converter.py
taxonomy_converter: Create a newick file from names.dmp and nodes.dmp.

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Name of the output database.
  -a NAMES, --names NAMES
                        Path to a valid names.dmp
  -d NODES, --nodes NODES
                        Path to a valid nodes.dmp
```



# Licence
