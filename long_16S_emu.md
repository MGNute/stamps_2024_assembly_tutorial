## Full-length and error prone 16S rRNA microbiome analysis using Emu

### Within the terminal

A few necessities you should know for using the command line
```
# this will make a new directory named new_folder
mkdir new_folder

# this changes the current directory to that folder
cd new_folder 

# if you're uncertain where you are this will tell you 
pwd

# make a file named sup_guys.txt 
vim sup_guys.sh

# let's try to use vim

# what do we have in our folder
ls

# to delete a file
rm sup_guys.txt

# to get back one level (recent directory)
cd ..

# to remove a folder
rm -r new_folder
```
That's a good start, we'll work on our skills more as we go!

#### Making our conda env 
By making a **conda** environment of course!

![conda](https://github.com/MGNute/stamps_2024_assembly_tutorial/blob/main/img/conda.png)

```
# Now make the conda environment
conda create -n stamps_emu 

# Activate the conda environment
conda activate stamps_emu

# Couple of dependencies
pip install flatten-dict
conda install -y -c conda-forge libgcc-ng

# Then install the programs
conda install -y -c bioconda emu chopper seqkit nanoplot

# and test
emu -h
```
#### Downloading test data
```
# Makes a new directory
mkdir emu_data

# changes to that directory
cd emu_data

# Grabs the culture data
wget https://zenodo.org/records/12790714/files/subsamp_colonies.tar.gz

# Unpacks the file
tar -xf subsamp_colonies.tar.gz

# Grabs the mock mixed colony data
wget https://zenodo.org/records/12790714/files/subsamp_mock.tar.gz

# Unpacks the file
tar -xf subsamp_mock.tar.gz
```

#### Raw sequence data observation
Once our dependencies are installed and data downloaded we will use [NanoPlot](https://github.com/wdecoster/NanoPlot). This is a program within the NanoPack2 suite of tools.

```
# For individual files
NanoPlot -t 15 --fastq sample.fastq --maxlength 40000 --plots dot sample_name
```
For many files
```
# Make the nanoplot file
vim nanoplotter.sh

# To get into edit mode press i

# Paste everything below to make your new file

#! /bin/bash

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate stamps_emu

echo "Happy NanoPlotting"

DIR="$PWD"
for f in "${DIR}"/*.fastq; do filename="${f%*.fastq}";
    NanoPlot -t 15 --fastq $filename".fastq" --maxlength 40000 --plots dot -o $filename
done

# Then exit vim using escape then :wq to save
```
```
# To run this script 

chmod +x nanoplotter.sh

./nanoplotter.sh
```

This should allow you to view your files in the bottom right window of Rstudio if you would like. If you want to learn more about vim, check the [FAQ](https://vimhelp.org/vim_faq.txt.html) link.

At around which quality score was this data prefiltered at?


#### Read filtering and trimming
There are many options for this step but for long reads use [chopper](https://github.com/wdecoster/chopper) or [filtlong](https://github.com/rrwick/Filtlong) (if you don't know if you still have primers on your reads choose filtlong but most basecalling steps currently remove them)

![mg_coolguy](https://github.com/MGNute/stamps_2024_assembly_tutorial/blob/main/img/mg_coolguy.png)
**Mission Impossible Mike**

To use chopper on a single file
```
# trimming to only keep 16S sequences 1,350 < x < 1,650 bp with 15 threads
chopper --minlength 1350 --maxlength 1650 --threads 15 -i filename.fastq > filename_chop.fastq
```
For many files
```
# Make the chopper file
vim choppa.sh

# To get into edit mode press i

# Then paste everything below to make your new file

#! /bin/bash

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate stamps_emu

echo "get to the choppa"

DIR="$PWD"
for f in "${DIR}"/*.fastq; do filename="${f%*.fastq}";
    chopper --minlength 1350 --maxlength 1650 --threads 15 -i $filename".fastq" > $filename"_chop.fastq"
done

# Then exit vim by clicking escape then :wq to save
```
```
# To run this script 

chmod +x choppa.sh

./choppa.sh
```
Challenge:
Try repeating the NanoPlot step on these new files to see how the data was affected by that step.

#### Emu classification
Now we can take the outputs of the chopper step and use these sequences within [emu](https://github.com/treangenlab/emu)!

![emu](https://github.com/MGNute/stamps_2024_assembly_tutorial/blob/main/img/emu.png)

Now to install the Emu database 
```
# Install the Emu database

# Where file housed at (make sure you're still in the stamps_emu conda env)
pip install osfclient

# Make directory for the Emu database
mkdir emu_db

# Go to that directory
cd emu_db

# Get that path
pwd

# Set the Emu DB dir with the path just found
export EMU_DATABASE_DIR=<path_to_database>

# You should still be in that Emu DB dir
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar

# Then just extract the tar file containing our DB
tar -xvf emu.tar
```

For single files
```
emu abundance filename_chop.fastq --db $EMU_DATABASE_DIR --threads 15 --keep-counts --output-dir emu_output
```
For many files
```
# Make the chopper file
vim emu_4_u.sh

# To get into edit mode press i

# Then paste everything below to make your new file

#! /bin/bash

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate stamps_emu

echo "the coolest flightless bird around"

DIR="$PWD"
for f in "${DIR}"/*_chop.fastq; do filename="${f%*_chop.fastq}"; 
    emu abundance $f --db EMU_DATABASE_DIR --threads 15 --keep-counts --output-dir emu_output ;  
done

# This takes the output from emu abundance and then creates an OTU tables for all your samples at the genus level
emu combine-outputs --split-tables --counts $DIR/emu_output genus

# This does the same thing but at the species level
emu combine-outputs --split-tables --counts $DIR/emu_output species

# And exit vim by clicking escape then :wq to save
```
```
# To run this script 

chmod +x emu_4_u.sh

./emu_4_u.sh
```
Once we get all of our files processed with emu, we will get a filename_rel-abundance.tsv file for each sample. For ease of processing we can use the last two lines convert that to an OTU-type table for genus level and species level counts.


### Now moving to Rstudio
#### Visualization

Let's make a phyloseq object from our newly made OTU table in R. (MIKE EDIT AS YOU PLEASE)

```
# Metadata file
meta = read.delim("metadata_all_PACS.tsv", check.names = FALSE) |>
  select(metadata_columns_to_add) |>
  tibble::column_to_rownames("sample")

# Load and prepare taxonomy file
tax_file <- read.delim("emu-combined-taxonomy-species.tsv") %>%
  select(superkingdom, phylum, class, order, family, genus, species) %>%
  as.matrix()

# Ensure all elements in tax_file are characters
tax_file <- apply(tax_file, 2, as.character)

# Load in Emu results
emu <- read.delim(file = "emu-combined-abundance-species-counts.tsv", check.names = FALSE) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("species")

# Filter out columns ending with "-threshold-0.0001"
filt_cols <- grep("-threshold-0.0001$", names(emu), value = TRUE, ignore.case = TRUE)
filt_emu <- emu[, !(names(emu) %in% filt_cols)]

# Convert to matrix and ensure numeric type for otu_table
emu_otu_mat <- as.matrix(filt_emu)
emu_otu_mat <- apply(emu_otu_mat, 2, as.numeric)

# Create the OTU, sample data, and taxonomy table objects
emu_otu <- otu_table(emu_otu_mat, taxa_are_rows = TRUE)
emu_sam <- sample_data(final_meta)
emu_tax <- tax_table(as.matrix(tax_file))

# Combine into a phyloseq object
physeq <- phyloseq(emu_otu, emu_sam, emu_tax)
```
Props to anyone who gets this working

There is a bit more to do if you want to make trees from this data so we can use UniFrac distances but its on the way.

To give credit where credit is due I used this [cool tool](https://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20) for the ASCII art. 

And I just want to thank my dog for being the best lil guy
![coopy](https://github.com/MGNute/stamps_2024_assembly_tutorial/blob/main/img/coopy.png)
