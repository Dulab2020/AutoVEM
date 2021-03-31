# Guides for AutoVEM

[toc]

## 0. Basic Messages

Program: AutoVEM (works on Linux machine)

Description: SARS-CoV-2 epidemic trends analysis tool

Version: V1.0

Authors: Xi Binbin

​				 Du Hongli

Contacts: 201766841276@mail.scut.edu.cn 

​					hldu@scut.edu.cn

Year: 2020

License: Released under GNU General Public License

## 1. Prerequisites

- Python (v3.8.6 or higher) with pandas, numpy and matplotlib packages (latest version).
- Java Runtime Environment (should meet the requirements of Haploview tool).

## 2. Dependencies

Bowtie2 (v2.4.2), SAMtools (v1.10), BCFtools (v1.10.2) and VCFtools (v0.1.16) are required. Please make sure you have installed these tools and add them to the environment variables globally before using AutoVEM. You can visit the following websites to install them:

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

SAMtools: http://samtools.sourceforge.net/

BCFtools: http://samtools.github.io/bcftools/bcftools.html

VCFtools: http://vcftools.sourceforge.net/

Haploview.jar (v4.2) has been provided. Please make sure the Java Runtime Environment of your computer meets the requirements of Haploview. For more details about Haploview, please visit: https://www.broadinstitute.org/haploview/tutorial

## 3. Installation

The tool does not need to be installed and can be used directly.

## 4. Format of Genome Sequences

Virus sequences should be whole genome sequence downloaded from GISAID. All genome sequences can be put in a fasta format file or several fasta format files. Put the file(s) into a directory. The directory will be the input of AutoVEM.

## 5. Usage

```py
cd path/AutoVEM
python main.py {auto|individual} --input [--sites] [--frequency] --output
```

You can use `auto` mode to run your work from A to Z. It will complete **SNVs calling**, **obtaining key mutations and finding haplotype of every sequence**, **visualizing the epidemic trends of these haplotypes**. The `--input` of this mode is a directory that stores the virus genome sequences.

You can also use `individual` mode to do more customizable analysis. This mode complete **obtaining key mutations and finding haplotype of every sequence, visualizing the epidemic trends of these haplotypes**. The `--input` of this mode is the `snp_merged.tsv` file produced by `auto` mode.

There are two optional arguments:  `--frequency` and `--sites`

`--frequency`: If provided, sites with mutation frequency less than the given frequency will be filtered out, default 0.05.

`--sites`: A list of space delimited integer numbers, such as `--sites 241 3037 14408 23403`. If provided, the sites will be as the specific sites which will used to find haplotypes.

**Note**: If both `--frequency` and `--sites`  arguments are not provided, will use `--frequency 0.05` by default. If both `--frequency` and `--sites` arguments are provided, intersection will be used. That is, if mutation frequency of a provided site less than the given frequency, it will be filtered out.

## 6.  Output Files

` snp_merged.tsv` 	Produced by ` auto` mode, ` \t` delimited. It stores information about all **SNV** mutations of all genome sequences. 

> **Id**: the identifier of the sequence 
>
> **Date**: collection date of the sequence 
>
> **Country**: country or region where the sequence was collected 
>
> **Position**: mutation position 
>
> **Ref**: reference base 
>
> **Alt**: mutation base

![image-20210314151850978](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314151850978.png)

**Note**: If  position is `0`, it means this sequence has no mutation. And the `Ref` and `Alt` will be `NA`.



`seq_info.tsv`	Produced by `auto` mode. It stores summary information about quality control.

![image-20210314151929583](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314151929583.png)



`snp_sites.tsv`	Produced by `auto` mode or `individual` mode, `\t` delimited. It stores information of the specific sites.

> Position: mutation position
>
> Ref: reference base
>
> Alt: mutation base
>
> Frequency: mutation frequency of the site

![image-20210314151955513](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314151955513.png)



`plot.LD.PNG`	Linkage analysis result of the specific sites. Produced by `auto` mode or `individual` mode.

![plot.LD]https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/plot.LD.PNG)

`haplotypes.tsv`	Produced by `auto` mode or `individual` mode, `\t` delimited. It stores the information of the haplotypes.

![image-20210314152039368](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314152039368.png)

`data_plot.tsv`	Produced by `auto` mode or `individual` mode, `\t` delimited. 

> Id: identifier of the sequence.
>
> Date: time when the sequence was collected.
>
> Country: country or region where the sequence was collected.
>
> Hap: haplotype sequence of the genome sequence.
>
> Name: haplotype name of the genome sequence.

![image-20210314152112776](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314152112776.png)



`hap_date_country.tsv`	Produced by `auto` mode or `individual` mode, `\t` delimited. It stores the number of different haplotypes of different countries or regions in different periods.

![image-20210314152434602](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/image-20210314152434602.png)



`hap_date.pdf`	Produced by `auto` mode or `individual` mode.  Visualization of the haplotypes epidemic trends.

![hap_date](https://github.com/Dulab2020/AutoVEM/blob/main/README.assets/hap_date.png)
