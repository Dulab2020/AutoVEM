PROGRAM: AutoVEM
DESCRIPTION: Virus epidemic trends analysis tool
AUTHOR: Xibinbin, Duhongli et al.
CONTACT: 201766841276@mail.scut.edu.cn, hldu@scut.edu.cn
YEAR: 2020
LICENSE: Released under GNU General Public License

1. Prerequisite
Python (v3.8.6 or higher) with pandas, numpy and matplotlib packages(latest version)
JAVA

2. Dependency
Bowtie2 (v2.4.2), samtools (v1.10), bcftools (v1.10.2) and vcftools (v0.1.16) are required. Please make sure you have installed these tools and add them to the environment path globally.
You can visit the following websites to install them:
bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
samtools: http://samtools.sourceforge.net/
bcftools: http://samtools.github.io/bcftools/bcftools.html
vcftools: http://vcftools.sourceforge.net/

Haploview.jar (v4.2) has been provided under the tools directory. And please make sure the Java Runtime Environment meets the requirements of Haploview. 
For more detail about Haploview, please visit:
https://www.broadinstitute.org/haploview/tutorial

3. INSTALLATION
The tool does not need to be installed and can be used directly.

4. USAGE
cd path/AutoVEM
python main.py auto|individual --input [--sites] [--frequency] --output