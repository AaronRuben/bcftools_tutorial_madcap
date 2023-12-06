# A brief introduction to the Variant Call Format (VCF) and BCFtools

Genetic data is frequently distributed in the [variant call format (VCF)](https://samtools.github.io/hts-specs/VCFv4.4.pdf). VCF files follow standardized specifications that allow the easy analysis with common software tools, such as [BCFtools](https://samtools.github.io/bcftools/howtos/index.html). BCFtools is a powerful software suits providing tools for viewing, manipulating, and analyzing VCF files. In the following, I will first introduce the VCF file format, and then provide a brief introduction of BCFtools. This tutorial does not explore all the functionalities that come with BCFtools, but it is rather meant as an icebreaker. I recommend looking at the extensive [documentation](https://samtools.github.io/bcftools/bcftools.html) of BCFtools of learn about more powerful ways in which it can be used.

## Variant Call Format (VCF)
VCF is a tab-separated text file format that stores genetic variation data. It always starts with meta-information lines that are prefixed with "##" followed by a header line that is prefixed with "#" and finally data lines. The header line indicates the different fields in the data lines, which always contain information about the physical position (i.e., chromosome number and position) as well as genotype information for each sample at a given position. Missing data is represented with a dot (".") as zero length fields are not allowed. Let's take a look at an example header:
```
##fileformat=VCFv4.4
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS       ID	       REF    ALT      QUAL   FILTER   INFO                                FORMAT        NA00001	  NA00002	    NA00003
20      14370     rs6054257    G      A	       29     PASS     NS=3;DP=14;AF=0.5;DB;H2             GT:GQ:DP:HQ   0|0:48:1:51,51	  1|0:48:8:51,51    1/1:43:5:.,.
20	17330	  .	       T      A	       3      q10      NS=3;DP=11;AF=0.017                 GT:GQ:DP:HQ   0|0:49:3:58,50	  0|1:3:5:65,3	    0/0:41:3:.
20	1110696   rs6040355    A      G,T      67     PASS     NS=2;DP=10;AF=0.333,0.667;AA=T;DB   GT:GQ:DP:HQ   1|2:21:6:23,27	  2|1:2:0:18,2	    2/2:35:4:.
20	1230237	  .	       T      .	       47     PASS     NS=3;DP=13;AA=T	                   GT:GQ:DP:HQ   0|0:54:7:56,60	  0|0:48:4:51,51    0/0:61:2:.
20	1234567	  microsat1    GTC    G,GTCT   50     PASS     NS=3;DP=9;AA=G                      GT:GQ:DP      0/1:35:4	  0/2:17:2	    1/1:40:3
```
Here, we have SNPs and three samples. The data lines always follow the following structure: 
 1. CHROM - chromosome ID
 2. POS - position (1-based)
 3. ID - variant ID
 4. REF - reference allele
 5. ALT - alternative allele
 6. QUAL - PHRED-scaled quality score of the alternative allele
 7. FILTER - passed filters
 8. INFO - typically provides general annotation of the given variant. Descriptions of the annotations can be found in the corresponding meta-information lines.
If genotype data is present, there will also be
  9. FORMAT field - specifying the data and order of individual sample information. Again descriptions of the different data fields within the FORMAT field can be found in the meta-information lines.
  10. Individual sample genotype data.
  11. ...

Individual genotype data can either be unphased or phased. If the data is phased, we know which allele was inherited from the mother and father, respectively. Unphased genotypes are separated by **/**, while phased genotypes are separated by **|**. Furthermore, 0 always refers to the REF allele, 1 refers to the first ALT allele, and if multi-allelic, 2 refers to the second ALT allele, etc.

To reduce disk space, VCFs are often block compressed using bgzip (indicated by a **.gz** suffix). Compressed files are not human-readable per se and need to be uncompressed first or viewed with BCFtools (more on this later).

Furthermore, to enable efficient data lookups VCFs are often indexed with [tabix](http://www.htslib.org/doc/tabix.html), creating an additional file with the exact same name as the VCF augmented with **.tbi** suffix.

For a more detailed description of VCF, I refer you to the [manual](https://samtools.github.io/hts-specs/VCFv4.4.pdf).

## Introduction to BCFtools
As mentioned earlier BCFtools is a power software suite to view, manipulate, and analyze VCF files. In the following, I'll guide you through the installations (that's often the hardest part) and some simple examples. Note that I assume that you are using a Mac (or Linux distribution), although most of the code should also work on Windows computers. If you use Windows, I recommend installing [Ubuntu terminal](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview). If you run into problems, please, let me know.

### Installing BCFtools
To install BCFtools, open a terminal and run the following lines of code:
```
# Download htslib and BCFtools.
user@bash:~ $ git clone --recurse-submodules https://github.com/samtools/htslib.git
user@bash:~ $ git clone https://github.com/samtools/bcftools.git
# Do the actual installation
user@bash:~ $ cd bcftools
user@bash:bcftools $ make
```
Finally, to make the bcftools executable available from any directory (i.e., make it a global variable), we have to add the directory in which we installed BCFtools to our `PATH` variable:
```
user@bash:~ $ export PATH=$PATH:~/bcftools:~/bcftools/misc/
```
If you don't want to run the above command every time you start a new terminal, I recommend adding this line to the bottom of your `~/.bash_profile` file (just use your favorite editor).


To test if the installation worked, run:
```
user@bash:~ $ bcftools --help
```
Congratulations, you are all set, if you see something like this:
```
Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
Version: 1.18-25-g44deedcd (using htslib 1.18-52-g2140d03e)

Usage:   bcftools [--version|--version-only] [--help] <command> <argument>

Commands:

 -- Indexing
    index        index VCF/BCF files

 -- VCF/BCF manipulation
    annotate     annotate and edit VCF/BCF files
    concat       concatenate VCF/BCF files from the same set of samples
    convert      convert VCF/BCF files to different formats and back
    head         view VCF/BCF file headers
    isec         intersections of VCF/BCF files
    merge        merge VCF/BCF files files from non-overlapping sample sets
    norm         left-align and normalize indels
    plugin       user-defined plugins
    query        transform VCF/BCF into user-defined formats
    reheader     modify VCF/BCF header, change sample names
    sort         sort VCF/BCF file
    view         VCF/BCF conversion, view, subset and filter VCF/BCF files

 -- VCF/BCF analysis
    call         SNP/indel calling
    consensus    create consensus sequence by applying VCF variants
    cnv          HMM CNV calling
    csq          call variation consequences
    filter       filter VCF/BCF files using fixed thresholds
    gtcheck      check sample concordance, detect sample swaps and contamination
    mpileup      multi-way pileup producing genotype likelihoods
    roh          identify runs of autozygosity (HMM)
    stats        produce VCF/BCF stats

 -- Plugins (collection of programs for calling, file manipulation & analysis)
    0 plugins available, run "bcftools plugin -l" for help

 Most commands accept VCF, bgzipped VCF, and BCF with the file type detected
 automatically even when streaming from a pipe. Indexed VCF and BCF will work
 in all situations. Un-indexed VCF and BCF and streams will work in most but
 not all situations
```
This prints out the different programs implemented in BCFtools. To get more information on a specific program, e.g., `view`, run:
```
user@bash:~ $ bcftools view --help
```
### Download toy data set for this tutorial
To clone this repository and download the toy data set, run the following lines of code:
```
user@bash:~ $ git clone https://github.com/AaronRuben/bcftools_tutorial_madcap.git
user@bash:~ $ cd bcftools_tutorial_madcap
```

### View the data
In the `data` folder of this repository, you will find downsampled genotype data on chromosome 8 for 20 individuals from the 1000 genomes project (10 samples from YRI and 10 samples from GBR). Note, that all positions refer to hg19 coordinates.
Let's take a look at the data.
```
user@bash:bcftools_tutorial_madcap $ bcftools view data/example.vcf.gz
```
This is a **bad idea**, it will print a lot of lines (you can abort the command with CTRL+C).

<details>
  <summary>Try working with the documentation to print the first 10 data lines only! <b>Hint:</b> You want to use the "|" operator and "head" command.</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -H data/example.vcf.gz | head -n 10</code></pre>

  The `-H` flag suppresses bcftools from printing the meta-information lines. In bash, the `|` operator passes on (pipes) the out the following command. In this case, `head -n 10`, which only prints the first 10 lines.  
</details>

### Subsetting VCFs
In this section, we will subset a VCF in different ways.
#### Subset the VCF by sample
We can also subset a VCF by samples (e.g., case-vs-control or by population descriptor). For example, to extract the GBR samples only, we can run.
<details>
  
  <summary>Try working with the documentation and split the VCF into two VCFs by population. You can find YRI & GBR sample IDs in <code>data/yri_samples.txt</code> and <code>data/gbr_samples.txt</code>.</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -S data/yri_samples.txt -Oz -o data/yri_example.vcf.gz data/example.vcf.gz </code></pre>

This command will extract all samples contained in (`-S`) `data/yri_samples.txt` (one ID per line) from `data/example.vcf.gz` and write it to a new compressed VCF (`-Oz`) named (`-o`) `data/yri_example.vcf.gz`. 
</details>

<details>
  <summary>Now create a tabix index for the newly generated VCFs. This will allow efficient efficient look-ups in the future.</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools index -t data/gbr_example.vcf.gz</code></pre>

The `-t` flag tells bcftools to generate a tabix index rather than a CSI index, which is the default.
</details>

#### Subset by genomic region
<details>
  <summary>Assuming your favorite gene is <em>PCAT1</em>, try to extract variants that fall into the gene body. On hg19, the coordinates for <em>PCAT1</em> are: <b>chr8:128,025,399-128,033,259</b>. How many variants fall into this region? <b>Hint:</b> bcftools is sensitive to the chromosome notations (i.e., chr8 vs 8) and you can count lines with <code>wc -l</code></summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -H -r 8:128025399-128033259 data/example.vcf.gz | wc -l</code></pre>

The correct answer is **301**. The `-r` variable specifies a genomic region and only variants within this reason will be included.
</details>

#### Subset a set of variants by variant type, number of alleles present at a site, allele count, and allele frequency
<details>
  <summary>How many single-nucleotide polymorphisms (SNPs) fall into the <em>PCAT1</em> region?</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -H -r 8:128025399-128033259 -v snps data/example.vcf.gz | wc -l</code></pre>

The correct answer is **295**. The `-v` variable specifies which variant types to include. In this case, only SNPs are included.
</details>

<details>
  <summary>How many of those SNPs are biallelic?</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -H -r 8:128025399-128033259 -v snps --min-ac 2 --max-ac 2 data/example.vcf.gz | wc -l</code></pre>

The correct answer is **293**. The `--min-alleles` and `--max-alleles` variables specify how many alleles must be at least and at most present for a site to be included.
</details>

<details>
  <summary>How many of the biallelic SNPs have an allele count greater than or equal to 10?</summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools view -H -r 8:128025399-128033259 -v snps --min-ac 2 --max-ac 2 --min-ac 10 data/example.vcf.gz | wc -l</code></pre>

The correct answer is **8**. The `--min-ac` variables specify the minimum number of non-reference allele counts for a site to be included.
</details>

<details>
  <summary>How many of all variants in the <em>PCAT1</em> region are common in Africa (AF > 0.01) and rare in Europe (AF <=0.01)? <b>Hint:</b> you may want to use <code>bcftools query</code></summary>

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\t%GT]\n' -i 'AFR_AF > 0.01 && EUR_AF <= 0.01' -r 8:128025399-128033259 data/example.vcf.gz | wc -l</code></pre>

The correct answer is **26**. The `-f` variable specifies which fields to print (i.e., the output format), while the `-i` variable specifies conditions for a site to be included. Note that quotation marks are required!! *&&* indicates that both conditions must be met for a site to be included. If you use *||* instead, a site will be included if either condition is met. 
</details>
<details>
  <summary><h2>Bonus</h2></summary>
  
### Get summary statistics
Let's look at some of the statistics, describing the data:

<pre><code>user@bash:bcftools_tutorial_madcap $ bcftools stats data/example.vcf.gz > example_stats.txt</code></pre>

`bcftools stats` will compute some basic statistics, such as the allele frequency spectrum and read depth distribution. We will write this information to a file called `example_stats.txt`.

If you have `pdflatex` or `tectonic` installed, you can try to visualize the statistics, using:
<pre><code>plot-vcfstats -p plots example_stats.txt</code></pre>
This will create a series of figures in a new directory called `plots`.
</details>
