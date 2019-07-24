# Supplementary Note 1: RDP DB, Greengenes and SILVA statistics

The goal here is to count the number of unique genera, species and strains in all three databases for named high-quality named sequences. For RDP we use the training set for Naïve Bayes classifier. As the initial trusted sets (NOT representative set, but a set used for initial multiple sequence alignment) for Greengenes and SILVA are not available, we use Greengenes sequences clustered at 99% identity and SILVA “high-quality full length Ref dataset”. Then, we select only sequences that have species or strain annotations. Here we provide Linux shell commands to compare these databases.

## RDP DB Training Set v16

For this part we will need to download both the training set and the full database (because the species / strain names are only present in the full database and that is the only piece of information used from it):

```sh
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz
gunzip *.gz
```
Download the original training set:

```sh
wget https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo16_rawtrainingdata.zip/download
mv download rdp_trainset16.zip
unzip rdp_trainset16.zip
cd RDPClassifier_16S_trainsetNo16_rawtrainingdata
```

Extract RDP ID for each training set sequence, then use that to filter the entire RDP FASTA file:

```sh
grep "^>" trainset16_022016.fa | sed -E 's/.*\|([^\s]+)\sRoot.*$/\1/' > trainset16_rdpids.txt
```

Get all strains:
```sh
grep "^>" current_Bacteria_unaligned.fa | cut -f1 | sed -E 's/^>([A-Z0-9]+) /\1\t/' | sed -E 's/;//g' > rdp_v16_all_strains_bacteria.txt

grep "^>" current_Archaea_unaligned.fa | cut -f1 | sed -E 's/^>([A-Z0-9]+) /\1\t/' | sed -E 's/;//g' > rdp_v16_all_strains_archaea.txt

cat rdp_v16_all_strains_bacteria.txt rdp_v16_all_strains_archaea.txt > rdp_v16_all_strains.txt
```

Now filter out strains with unknown names:

```sh
awk '$2 !~ "^[a-z]" && $2 !~ "^Bacterium" { print $0 }' rdp_v16_all_strains.txt > rdp_v16_all_strains_with_names.txt
awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1]
  }' rdp_v16_all_strains_with_names.txt \
  RDPClassifier_16S_trainsetNo16_rawtrainingdata/trainset16_rdpids.txt > \
  rdp_v16_train_strain_names.txt
```

Subtract 1 because the first line is blank to get number of unique strains

```sh
cut -f2 rdp_v16_train_strain_names.txt | sort | uniq > rdp_v16_trainset_unique_strains.txt
wc -l rdp_v16_trainset_unique_strains.txt
# Unique species:
tail -n +2 rdp_v16_trainset_unique_strains.txt | cut -d" " -f1-2 | sort | uniq | wc -l
# Unique genera:
tail -n +2 rdp_v16_trainset_unique_strains.txt | cut -d" " -f1 | sort | uniq | wc -l
# Unique: 12089 strains, 11675 species, 2465 genera
```

## Greengenes 13.8

Download and extract:

```sh
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -xvzf gg_13_8_otus.tar.gz
```

Generate names with full species and count from sequences clustered at 99% identity:

```sh
grep -v "s__$" gg_13_8_otus/taxonomy/99_otu_taxonomy.txt | cut -d" " -f6,7 | sed -E 's/([gs]__|\;)//g' | sort | uniq > gg_13_8_99otu_taxonomy_species.txt
sort gg_13_8_99otu_taxonomy_species.txt | uniq -c | wc -l
# 3114
```

Now count unique genera (6th field, separated by space)

```sh
grep -v "g__;" gg_13_8_otus/taxonomy/99_otu_taxonomy.txt | cut -d " " -f6 | sed -E 's/([gs]__|\;)//g' | grep "^[A-Z][a-z]" | sort | uniq -c | wc -l
# 1889
```

There is no strain information here. Uniques: NA strains, 3114 species, 1889 genera

## SILVA v132

Generating statistics based on the SILVA_132_SSURef_tax_silva.fasta. First extract all entries that have species identification:

```sh
wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_tax_silva.fasta.gz
gunzip SILVA_132_SSURef_tax_silva.fasta.gz
mkdir silva
mv SILVA_132_SSURef_tax_silva.fasta.gz silva

awk '$0 ~ "^>" { 
    split($0, a, ";"); 
    split(a[1], king, " ");
    if (a[7] !~ "[Uu][nidentified|ncultured]" && a[7] !~ " [Bb]acterium" && 
        a[7] ~ " " && a[7] !~ "^[a-z]" && (king[2] == "Bacteria" || king[2] == "Archaea") ) { 
            sub("Candidatus[ ]?", "", a[7]); 
            sub("\\[", "", a[7]); 
            sub("\\]", "", a[7]); 
            print a[7] 
    }
}' SILVA_132_SSURef_tax_silva.fasta | sed -E s/\'//g | sort | uniq > SILVA_132_SSURef_tax_silva_species.txt
```

Also generate a single-line FASTA file with the same species:

```sh
awk ' { 
    if ($0 ~ "^>") {
        split($0, a, ";"); 
        if (a[7] !~ "[Uu][nidentified|ncultured]" && a[7] !~ " [Bb]acterium" && 
            a[7] !~ "Bacterium" && a[7] !~ "^[Bb]acterium" && a[7] ~ " " && 
            a[7] !~ "^[a-z]") { 
                sub("Candidatus[ ]?", "", a[7]); 
                sub("\\[", "", a[7]); 
                sub("\\]", "", a[7]);
                numrec += 1
                lastgood = 1
                if (numrec == 1) print ">"a[7];
                else print "\n>"a[7];
        } else {
            lastgood = 0;    
        }        
    } else {
        if (lastgood == 1) {
            gsub("U", "T", $0);
            printf $0;
        }
    }
}' SILVA_132_SSURef_tax_silva.fasta > SILVA_132_SSURef_tax_silva_species.fasta
```

This will give us a list of unique strains and species that have no strain designation. Extract only strains:

```sh
awk '$0 ~ "[^ ]+ [^ ]+ .*" { print $0 }' SILVA_132_SSURef_tax_silva_species.txt | sort | uniq > SILVA_132_SSURef_tax_silva_strains.txt
```

Also get a FASTA with strain names:

```sh
awk '{
    if ($0 ~ "^>") {
        n = split($0, a, " ");
        if (n > 2) {
            lastgood = 1;
            print $0;
        } else {
            lastgood = 0;
        }
    } else {
        if (lastgood == 1) print $0;
    }
}' SILVA_132_SSURef_tax_silva_species.fasta > SILVA_132_SSURef_tax_silva_strains.fasta
```

Now generate a BLAST database from this FASTA file. We will use this database to match strain sequences from HiMAP database to this and look for 100% hits to see which strains have an exact match:

```sh
makeblastdb -dbtype nucl -in SILVA_132_SSURef_tax_silva_strains.fasta -out SILVA_132_SSURef_tax_silva_strains
```

Count unique species assignments (replace multiple spaces from uniq with just one for easy manipulation):

```sh
cut -d" " -f1-2 SILVA_132_SSURef_tax_silva_species.txt | sed -E s/\'//g | grep -v "sp\.$" | sort | uniq -c > SILVA_132_SSURef_tax_silva_species_unique_counts.txt

sed -E "s/[ ]+/ /g" SILVA_132_SSURef_tax_silva_species_unique_counts.txt > SILVA_132_SSURef_tax_silva_species_unique_counts_table.txt
```

Count unique genera:

```sh
cut -d" " -f3 SILVA_132_SSURef_tax_silva_species_unique_counts_table.txt | sort | uniq -c | sed -E 's/[ ]+/ /g' > SILVA_132_SSURef_tax_silva_genus_unique_counts_table.txt
```

Count unique strains:

```sh
sort SILVA_132_SSURef_tax_silva_strains.txt | uniq -c | wc -l
```

Uniques: 3242 genera, 12734 species, 116543 strains.
Also generate a single-line FASTA file with the same species:

```sh
awk ' { 
    if ($0 ~ "^>") {
        split($0, a, ";"); 
        if (a[7] !~ "[Uu][nidentified|ncultured]" && a[7] !~ " [Bb]acterium" && 
            a[7] !~ "Bacterium" && a[7] !~ "^[Bb]acterium" && a[7] ~ " " && 
            a[7] !~ "^[a-z]") { 
                sub("Candidatus[ ]?", "", a[7]); 
                sub("\\[", "", a[7]); 
                sub("\\]", "", a[7]);
                numrec += 1
                lastgood = 1
                if (numrec == 1) print $0;
                else print "\n"$0;
        } else {
            lastgood = 0;    
        }        
    } else {
        if (lastgood == 1) {
            gsub("U", "T", $0);
            printf $0;
        }
    }
 }' SILVA_132_SSURef_tax_silva.fasta > SILVA_132_SSURef_tax_silva_species_wtax.fasta
```

What is the length distribution of these reference sequences? Let’s use this for determining BLAST word size in one of the next steps.

```sh
awk '$0 !~ "^>" { print length($0) }' SILVA_132_SSURef_tax_silva_species_wtax.fasta | sort -g | uniq -c | less
```

Most are between 1200 and 1500. Exact matching first between SILVA and GG was then used to find alignments between 100% hits. Generate a BLAST database:

```sh
makeblastdb -dbtype nucl -in SILVA_132_SSURef_tax_silva_species_wtax.fasta -out SILVA_132_SSURef_tax_silva_species_wtax
```

Then BLAST the SILVA fasta file against this database.
