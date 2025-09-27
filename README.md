# Sortadate_IQTree2_helper

STILL TESTING SCRIPT!!!!

When using the [iqtree](https://iqtree.github.io/) 2.3.6 `-S` flag to generate locus/gene trees, this script will split the locus trees into a new directory and reroot them for down stream use for using [sortadate](https://github.com/FePhyFoFum/phyx). The section in the `*.iqtree` file describing basic statistics of the resulting trees indicates what locus/gene they came from (ID column).

Notably, the script assumes the partition file came from AMAS, where the locus trees are appended with a `p*_{locus_name}`

To check that you have called things correctly, you can include the path to the directory that includes the fasta files and the script will print out if you have named the correctly

In development, I removed the full file name thinking I would be slick, but I had to re-add the `_trimalauto` suffix to the file name for consistency, so that is why that flag exists.

This script will also prepare a rerooted tree file of each using [phyx](https://github.com/FePhyFoFum/phyx) so you can run sortadate. A python v3.10 conda environment with `bioconda::phyx` installed is required for the rerooting to work. 
An additional flag called --keep_only_outgroup, when called will *not* reroot any locus trees that do not include an outgroup.

see the --help flag for full descript/explination of all this script can do


Example script
```
python iqtree_sortadate_helper.py \
  --logfile locustrees.iqtree \
  --input_locus_trees locustrees.treefile \
  --output locus_trees \
  -suffix _trimalauto \
  --fasta_directory loci_keep/ \
  --reroot GCA_048127345.1 \
  --keep_only_outgroup
```
