#!/env/bin/python

'''
author='xibinbin et al.'
'''

import argparse
import os
import sys
import re
import pandas as pd 
import matplotlib.pyplot as plt 
import datetime
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors

dirname, filename = os.path.split(os.path.abspath(sys.argv[0]))
indexDirectory = os.path.join(dirname, 'index')
refsequence = os.path.join(indexDirectory, 'NC_045512.2.fa')
toolsDirectory = os.path.join(dirname, 'tools')
haploview = os.path.join(toolsDirectory, 'Haploview.jar')

def arguments():
    dulab1_parser = argparse.ArgumentParser(prog='AutoVEM V1.1', 
                                            description='SARS-CoV-2 epidemic trends analysis tool using data from GISAID')

    dulab1_parser.add_argument('modes', choices=['auto', 'individual'], help='Analysis mode')
    dulab1_parser.add_argument('--input', type=str, required=True,
                                        help='In auto mode, this should be a directory contains fasta format genome sequences files.\nIn individual mode, this should be the snp_merged.tsv file produced by AutoVEM')
    dulab1_parser.add_argument('--sites', type=int, nargs='+', default=None,
                                          help = 'Mutation sites of interest, at least two sites. Space delimited list of integer numbers.')
    dulab1_parser.add_argument('--frequency', type=float, default=None,
                                              help = 'SNV sites with mutation frequency equal to or more than the given frenquency will be retained.\nIf the --sites parameter is also provided, this parameter will be ignored.')
    dulab1_parser.add_argument('--output', type=str, required=True, help='Output directory')
    
    args = dulab1_parser.parse_args()

    return args

def extract_sequence(directory, genomeDirectory):
    '''
    split fasta file(s)

    :param directory: input genome directory
    :param genomeDirectory: the absolute output genome files directory
    '''
    in_dir = os.path.abspath(directory)
    files = os.listdir(in_dir)
    files_list = []
    for file in files:
        file = os.path.join(in_dir, file)
        files_list.append(file)

    id_path = ''
    pointer = 1
    for file in files_list:
        with open(file, 'r') as fhand:
            for line in fhand.readlines():
                line = line.rstrip()
                if len(line)==0:
                    continue
                if line[0]=='>':
                    id = line.split('|')[1]
                    id = id.replace(" ","_")
                    id = id + '.fa'
                    id_path = os.path.join(genomeDirectory, id)
                    if os.path.exists(id_path):
                        pointer = 1
                        continue
                    else:
                        with open(id_path, 'a') as fhand_sequence:
                            pointer = 0
                            fhand_sequence.write(line+'\n')
                else:
                    if pointer == 1:
                        continue
                    else:
                        line = line.replace(" ", "")
                        line = line.replace("\t", "")
                        with open(id_path, 'a') as fhand_sequence:
                            fhand_sequence.write(line+'\n')

def get_file_list(directory):
    '''
    genome sequences

    :param directory: extracted genome files directory
    :returns filesPath(list, str):
    '''
    filesPath = []
    filesList = os.listdir(directory)
    for path in filesList:
        tmp_path = os.path.join(directory, path)
        filesPath.append(tmp_path)
    return filesPath

def genome_quality_control(file):
    '''
    quality control for the first time

    :param file: a fasta format file of a genome sequence
    :returns flag(int): (-1/0, not pass/pass)
    '''
    with open(file, 'r', encoding='utf-8') as fhand:
        sequence = ''
        number_N = 0
        number_db = 0
        pattern = {'A', 'T', 'G', 'C', 'N'}
        for line in fhand.readlines():
            if line[0] == ">":
                continue
            else:
                line = line.rstrip()
                if len(line) != 0:
                    sequence = sequence + line

        len_sequence = len(sequence)
        for item in sequence:
            if item == 'N':
                number_N = number_N + 1
                continue
            elif item not in pattern:
                number_db = number_db + 1
                continue
            else :
                pass

    if (len_sequence<29000 or number_N>15 or number_db>50):
        flag = -1
        return flag
    else:
        flag = 0
        return flag

def split_sequence(file, path):
    '''
    split genome sequence to reads(30-100nt)

    :param file: genome sequence
    :param path: output directory
    :returns bassicMessage(dict), splitFile(file), tempAnalysisDirectory(temp directory):
    '''
    with open(file,'r') as fhand:
        sequence = ""
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line) == 0:
                pass
            elif line[0] == ">":
                header = line
            else:
                line = line.upper()
                sequence = sequence + line 
    
    pattern_date = r'20\d{2}-\d{1,2}-\d{1,2}$'
    match_date = re.findall(pattern_date, header)
    len_match = len(match_date)
    if len_match == 0:
        return -1, -1, -1    
    else:
        date = match_date[0]

    pattern_id = r'EPI_ISL_\d*'
    match_id = re.findall(pattern_id, header)
    len_match = len(match_id)
    if len_match == 0:
        return -1, -1, -1    
    else:
        id = match_id[0]
    region = header.split("/")[1]

    basicMessage = {}
    basicMessage['Id'] = id
    basicMessage['Country'] = region
    basicMessage['Date'] = date

    splitFastaFileName = basicMessage['Id'] + '_' + 'split' +'.fa'
    tempAnalysisDirectory = os.path.join(path, id)
    os.mkdir(tempAnalysisDirectory)
    splitedFile = os.path.join(tempAnalysisDirectory, splitFastaFileName)

    length = len(sequence)
    subSequence = list()
    for i in list(range(0,length-30)):
        if (i+100) >= length:
            sub = sequence[i:]
        else:
            sub = sequence[i:(i+100)]
        subSequence.append(sub)

    with open(splitedFile, 'a') as fhand:
        for i,subsequence in enumerate(subSequence):
            readName = ">read" + str(i) + "." + header[1:] + "\n" + subsequence + "\n"
            fhand.write(readName)

    return basicMessage, splitedFile, tempAnalysisDirectory

def align(file, ref, directory):
    '''
    alignment

    :param file: XXX_split.fa
    :param ref: reference genome sequence
    :param directory: temp directory
    :returns samFile: temp.sam
    '''
    samFile = os.path.join(directory, 'temp.sam')
    os.system(f'bowtie2 -f -x {ref} -U {file} -S {samFile}')
    return samFile

def sort(file, directory):
    '''
    sort reads

    :param file: temp.sam
    :param directory: temp directory
    :returns bamFile: temp.bam
    '''
    bamFile = os.path.join(directory, 'temp.bam')
    os.system(f'samtools sort {file} > {bamFile}')
    os.system(f"samtools index {bamFile}")

    return bamFile

def mpileup(file, directory):
    '''
    add @RG to the bam file

    :param file: temp.bam
    :param directory: temp directory
    :returns vcfFile: temp.vcf
    '''
    addHeader = os.path.join(directory, "temp_addheader.bam")
    os.system(f"picard AddOrReplaceReadGroups -I {file} -O {addHeader} --RGID Sample --RGLB AMPLICON --RGPL ILLUMINA --RGPU unit1 --RGSM Sample")
    os.system(f"samtools index {addHeader}")

    return addHeader

def call(vcfFile, directory, refSeq):
    '''
    call SNPs

    :param vcfFile: temp_addheader.bam
    :param directory: temp directory
    :param refSeq: reference genome
    :returns flag, snpFile
    '''
    snp_indel_file_path = os.path.join(directory, 'snp_indel.vcf')
    snp_file_path = os.path.join(directory, 'snp.vcf')
    indel_file_path = os.path.join(directory, 'indel.vcf')

    # can add [--threads <int>] to use multithreading
    os.system(f'gatk HaplotypeCaller -I {vcfFile} -O {snp_indel_file_path} -R {refSeq} -ploidy 1')
    os.system('vcftools --vcf %s --recode --keep-only-indels --stdout > %s' % (snp_indel_file_path, indel_file_path))
    n_indels = 0
    with open(indel_file_path, 'r') as fhand:
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line)==0:
                continue
            if line[0]=='#':
                continue
            else:
                n_indels = n_indels + 1
    if n_indels > 3:
        return -1, -1
    else:
        os.system('vcftools --vcf %s --recode --remove-indels --stdout > %s' % (snp_indel_file_path, snp_file_path))
        return 0, snp_file_path

def snp_mutation_information(file):
    '''
    obtain SNPs information

    :param file: snp.vcf
    :returns mutationInformation(dict): key=['Position', 'Ref', 'Alt']
    '''

    mutation_lines = list()
    with open(file, 'r') as fhand:
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line)==0:
                pass
            elif line[0]=="#":
                pass
            else:
                mutation_lines.append(line)
    
    mutationInformation = list()
    if len(mutation_lines) == 0:
        record = {'Position':0, 'Ref':'NA', 'Alt':'NA'}
        mutationInformation.append(record)
    else:
        for item in mutation_lines:
            item = item.split()
            pos = int(item[1])
            ref = str(item[3])
            alt = str(item[4])
            snpMutationMessage = dict()
            snpMutationMessage['Position'] = pos
            snpMutationMessage['Ref'] = ref
            snpMutationMessage['Alt'] = alt
            mutationInformation.append(snpMutationMessage)

    return mutationInformation

def snp_filter(file, directory, sites=None, fre=None):
    '''
    filter SNP sites

    :param file: snp_merged.tsv
    :param directory: output directory
    :param sites: snp sites of interest
    :param fre: frequency of snp sites that more than fre will be retained
    :returns snp_pos(list), snp_ref_alt(dict)
    '''
    ## storage the retained SNP sites
    snp_sites = os.path.join(directory, 'snp_sites.tsv')
    os.system(f'touch {snp_sites}')
    ## read the SNV files
    df = pd.read_csv(file, sep='\t')
    ## count the number of sequences
    ids = df['Id'].unique().tolist()
    n_genome = len(ids)
    ## calculate the mutation frequency at each position
    counts = df['Position'].value_counts()
    frequency = counts/n_genome
    frequency = frequency.round(decimals=4)
    frequency = frequency.sort_index()
    frequency = frequency[frequency.index!=0]
    snp_dict = dict(zip(frequency.index.tolist(), frequency.values.tolist()))

    snp_pos = list()
    ## sites with mutation frequency bigger than the given frequency or default frequency (0.05) will be retained.
    if sites is None:
        if fre is None:
            snpSites = frequency[frequency.values>=0.05]
            snp_pos = snpSites.index.tolist()
        else:
            snpSites = frequency[frequency.values>=fre]
            snp_pos = snpSites.index.tolist()
    ## if provided sites, these sites will be retained
    else:
        snp_pos = sites
    ## whether the number of sites bigger than 1
    if len(snp_pos)<=1:
        print('There are no or too little sites that meet your requirements.')
        sys.exit()
    else:
        snp_pos.sort()
        ## print the site retained
        out_line = ""
        for site in snp_pos:
            out_line = out_line + " " + str(site)
        print(f"The following sites will be retained: {out_line}")

        ## get the mutation information of the retained sites
        snp_ref_alt = dict()
        for snp in snp_pos:
            df2 = df[df['Position']==snp]
            Ref = df2['Ref'].value_counts().index.tolist()[0]
            Alt = df2['Alt'].value_counts().index.tolist()[0]
            snp_ref_alt[snp] = (Ref, Alt, snp_dict[snp])
        ## write the information of retained sites to the record file
        header = 'Position\tRef\tAlt\tFrequency\n'
        with open(snp_sites, 'a') as fhand:
            fhand.write(header)
            for pos, (Ref, Alt, Fre) in snp_ref_alt.items():
                record = str(pos) + '\t' + Ref + '\t' + Alt + '\t' + str(Fre) + '\n'
                fhand.write(record)
            
        return snp_pos, snp_ref_alt

def ref_haplotype(positions=None):
    '''
    reference haplotype sequence

    :param position: snp position
    :returns referenceHaplotype: reference haplotype sequence
    '''
    with open(refsequence, 'r') as fhand:
        referenceGenomeSequence = ''
        for line in fhand.readlines():
            line = line.rstrip()
            if len(line) == 0:
                pass
            elif line[0] == ">":
                pass
            else :
                referenceGenomeSequence = referenceGenomeSequence + line
    referenceHaplotype = ''
    for pos in positions:
        referenceHaplotype = referenceHaplotype + referenceGenomeSequence[pos-1]

    return referenceHaplotype

def genome_haplotype(file, positions, referenceHaplotype, snp_alt_dict, directory):
    '''
    get haplotype sequence of genome

    :param file: snp_merged.tsv
    :param position
    :param referenceHaplotype: reference haplotype sequence
    :param snp_alt_dict: {pos:[ref,alt,fre]}
    :param directory: absolute path of output directory
    :returns filePath: data.tsv
    '''
    ## the file stores the haplotype sequence of each sequence
    filePath = os.path.join(directory, 'data.tsv')
    df = pd.read_csv(file, sep='\t')
    ids = df['Id'].unique().tolist()

    Date = []
    Country = []
    Case_id = []
    snp_positions = []

    for idx in ids:
        df1 = df[df['Id']==idx]
        date = df1['Date'].value_counts().index.tolist()[0]
        country = df1['Country'].value_counts().index.tolist()[0]
        Case_id.append(idx)
        Date.append(date)
        Country.append(country)
        mutation_positions = dict(zip(df1['Position'].tolist(), df1["Alt"].tolist()))
        snp_positions.append(mutation_positions)

    haplotypes = list()
    for record in snp_positions:
        hap_seq = ''
        for i, snp in enumerate(positions):
            if snp in record:
                hap_seq = hap_seq + snp_alt_dict[snp][1]
            else:
                hap_seq = hap_seq + referenceHaplotype[i]
        haplotypes.append(hap_seq)

    data_df = pd.DataFrame(data={'Id': Case_id,
                                   'Date': Date,
                                   'Country':Country,
                                   'Hap': haplotypes})
    data_df.to_csv(filePath, sep='\t', index=False)
    return filePath

def block_file(position, directory):
    '''
    block.txt

    :param position: snp_pos
    :param directory: output directory
    :returns blockFile: block.txt
    '''
    ## block file
    num = len(position)
    blockFile = os.path.join(directory, 'block.txt')
    with open(blockFile, 'a') as fhand:
        for i in list(range(num)):
            fhand.write(str(i+1))
            fhand.write('\t')
    
    return blockFile

def map_file(position, directory):
    '''
    get the snp.info file

    :param position: snp_ref_alt
    :param directory: the absolute path of the output directory
    :returns mapFile: snp.info
    '''
    ## info file
    mapFile = os.path.join(directory, 'snp.info')
    with open(mapFile, 'a') as fhand:
        for key, item in position.items():
            name = str(item[0])+str(key)+str(item[1])
            record = name + '\t' + str(key) + '\n'
            fhand.write(record)
    
    return mapFile

def ped_file(file, directory):
    '''
    生成ped格式文件

    :param file: data.tsv
    :param directory: output directory
    :returns pedFile: snp.ped
    '''
    ## the snp.ped file
    pedFile = os.path.join(directory, 'snp.ped')
    df = pd.read_table(file)
    Ids = df.Id.tolist()
    HaplotypeSequence = df.Hap.tolist()
    with open(pedFile, 'a') as f:
        record = ""
        constantString = '\t0\t0\t0\t0\t'
        for i,idx in enumerate(Ids):
            hap_seq = HaplotypeSequence[i]
            tmp = list(hap_seq)
            genotype = list()
            for base in tmp:
                genotype.append(base)
                genotype.append(base)
            hap = "\t".join(genotype)
            record = str(i) + "\t" + str(idx) + constantString + hap + "\n"
            f.write(record)
            record = ""

    return pedFile

def linkage_analysis(ped, mapf, block, directory):
    '''
    linkage analysis

    :param ped: snp.ped
    :param mapf: snp.info
    :param block: block.txt
    :param directory: output directory
    :returns haplotypesFile: plot.CUSTblocks
    '''

    temp = os.path.join(directory, 'plot')
    haplotypesFile = os.path.join(directory, 'plot.CUSTblocks')

    os.system(f'java -jar {haploview} -n -skipcheck -pedfile {ped} -info {mapf} -blocks {block} -png -out {temp}')
    if os.path.exists(ped):
        os.remove(ped)
    if os.path.exists(mapf):
        os.remove(mapf)
    if os.path.exists(block):
        os.remove(block)
        
    return haplotypesFile

def haplotyper(file, dataFile, directory):
    '''
    hyplotype

    :param file: plot.CUSTblocks
    :param dataFile: data.tsv
    :param directory: output directory
    :returns dataPlot, haplotypes: data_plot.tsv
    '''
    dataPlot = os.path.join(directory, 'data_plot.tsv')
    haplotypes = os.path.join(directory, 'haplotypes_temp.tsv')
    nt_dict = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
    ## obtain the haplotype sequence
    hap_dict = dict()
    with open(file, 'r') as fhand:
        for i, line in enumerate(fhand.readlines()):
            sequence = ''
            if i == 0:
                pass
            else:
                hap = 'H' + str(i)
                sequence_num = line.split()[0]
                proportion = float(line.split()[1][1:-1])
                ## if proportion less than 0.01, it means the variant population is too small
                ## so this variant population will be filtered out
                if proportion < 0.01:
                    break
                else:
                    ## get the haplotype sequence and name
                    for num in sequence_num:
                        sequence = sequence + nt_dict[num]
                    hap_dict[sequence] = hap

    df = pd.read_table(dataFile)
    Name = list()
    for hapSeq in df.Hap.tolist():
        if hapSeq in hap_dict:
            Name.append(hap_dict[hapSeq])
        else:
            Name.append("other")
    df["Name"] = Name
    df.to_csv(dataPlot, sep="\t", index=None)
    
    with open(haplotypes, 'a') as fhand:
        for key, value in hap_dict.items():
            tem = ''
            tem = value + '\t' + key + '\n'
            fhand.write(tem)
    if os.path.exists(dataFile):
        os.remove(dataFile)
    if os.path.exists(file):
        os.remove(file)

    return dataPlot, haplotypes

def plot(file, file2, directory):
    '''
    结果可视化

    :param file: data_plot.tsv
    :param file2: haplotypes_temp.tsv
    :param directory: plot directory
    '''
    pdfPath = os.path.join(directory, 'hap_date.pdf')
    hap_temp = os.path.join(directory, 'haplotypes.tsv')
    hap_date = os.path.join(directory, 'hap_date_country.tsv')

    df = pd.read_csv(file, sep='\t')
    df.index = pd.to_datetime(df["Date"])
    df = df.sort_index(ascending=True)
    
    n_genome = df.shape[0]
    HaplotypeCounts = df['Name'].value_counts()
    Fre = HaplotypeCounts/n_genome
    Fre = Fre.round(decimals=4)
    Fre_dict = dict(zip(Fre.index.tolist(), Fre.values.tolist()))

    header = 'Name\tSequence\tFrequency\n'
    with open(hap_temp, 'a') as f:
        with open(file2, 'r') as fhand:
            f.write(header)
            for line in fhand.readlines():
                line = line.rstrip()
                line = line.split('\t')
                name = line[0]
                Frequency = str(Fre_dict[name])
                record = ''
                record = name + '\t' + line[1] + '\t' + Frequency + '\n'
                f.write(record)
        record = 'other' + '\t' + 'NA' + '\t' + str(Fre_dict['other'])
        f.write(record)
    os.system('rm -rf %s' % file2)


    hap = df['Name'].unique()
    hap = hap.tolist()
    n_hap=len(hap)

    country = df["Country"].unique()
    country = country.tolist()
    n_country = len(country)
    
    n=7
    start = df.index[0]
    end = df.index[-1]
    detal = end - start
    m = math.ceil( detal.days/n )
    num = list(range(0,m))
    date_list=[]
    date_list.append(start.strftime('%Y/%m/%d'))
    for c in num:
        a = pd.to_datetime(start+datetime.timedelta(days = int(c*n)))
        b = pd.to_datetime(start+datetime.timedelta(days = int((c+1)*n)))
        date_list.append(b.strftime('%Y/%m/%d'))
    
    def name(n):
        haps = []
        for i in list(range(n-1)):
            s = 'H' + str(i+1)
            haps.append(s)
        haps.append('other')
        return haps

    def size(n):
        m=0
        if(n>=1000):
            m=1
        elif(n>=100):
            m=0.85
        elif(n>=10):
            m=0.7
        else:
            m=0.5
        return m

    dx = 0.01
    dy = 0.001
    a = x_origin = 0.18
    b = y_origin = 0.05
    c = x_len = 0.618
    d = y_len = 0.90
    x_size = 16
    lx = (x_len-dx)/(m + 1)
    y_size = x_size*(n_country+1)*0.608/((m+1)*0.899)
    ly = (y_len-dy)/(n_country+1)
    sub_lx = lx/x_len
    sub_ly = ly/y_len
    xticks = []
    yticks = []
    for x in list(range(m+1)):
        tick = dx/x_len + x*sub_lx
        xticks.append(tick)
    for y in list(range(n_country+1)):
        tick = dy/x_len + sub_ly/2 + y*sub_ly
        yticks.append(tick)

    xticklabels = date_list.copy()
    con = country[::-1]
    yticklabels = con.copy()
    yticklabels.append("")

    # 单倍型配色方案
    sort_hap = name(n_hap)
    colors=list(mcolors.CSS4_COLORS.keys())
    col=colors[10:(n_hap+10)]

    header = 'Date\tCountry'
    for item in sort_hap:
        header = header + '\t' + item
    header = header + '\n'
    with open(hap_date, 'a') as fhand:
        fhand.write(header)

    pdf = PdfPages(pdfPath)
    plt.figure(figsize=[x_size, y_size])
    ax1 = plt.axes([a, b, c, d])
    ax1.set_xticks(xticks)
    ax1.set_yticks(yticks)
    ax1.set_xticklabels(xticklabels)
    ax1.xaxis.set_tick_params(rotation=30, labelsize=12)
    ax1.set_yticklabels(yticklabels, size=14)

    
    tp = ly*(n_country-1) + y_origin + dy
    for c in num:
        a = pd.to_datetime(df.index[0]+datetime.timedelta(days = int(c*n)))
        b = pd.to_datetime(df.index[0]+datetime.timedelta(days = int((c+1)*n)))
        if b == end:
            df2 = df[(df.index >=a) & (df.index <= b)]
        else :
            df2 = df[(df.index >=a) & (df.index < b)]
        date_section = a.strftime('%Y/%m/%d') + '-' + b.strftime('%Y/%m/%d')
        for r, j in enumerate(country):
            record = ''
            df3 = df2[df2["Country"] == j]
            con = j

            haps = ''
            tmp=[]
            for k in sort_hap:
                df4 = df3[df3["Name"] == k]
                l = len(df4.index.tolist())
                haps = haps + '\t' + str(l)
                tmp.append(l)
            record = date_section + '\t' + con + haps + '\n'
            with open(hap_date, 'a') as fhand:
                fhand.write(record)
            
            n_sum = sum(tmp)
            if n_sum == 0:
                continue
            p1 = x_origin + dx + lx*c
            p2 = tp - ly*r
            p3 = lx
            p4 = ly
            ax = plt.axes([p1, p2, p3, p4])
            ax.pie(tmp, colors=col, radius=size(n_sum), normalize=True)
            
    for ci, cj in enumerate(sort_hap):
        st = 0.05 + n_hap*ly
        ax3=plt.axes([0.80, st-ci*ly, lx, ly])
        ax3.pie([1],labels=[""],colors=[col[ci]],rotatelabels=True, radius=0.5)
        ax3.annotate(cj,xy=(0.80, st-ci*ly), fontsize=12)

    k=ly
    p = 0.80
    en = 0.05 + (n_hap+4)*ly + 3*ly
    ax4=plt.axes([p, en, lx, ly])
    ax4.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=1.0)
    ax4.annotate("n>=1000",xy=(p, en), fontsize=12)
    ax4=plt.axes([p, en-ly, lx, ly])
    ax4.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.85)
    ax4.annotate("n>=100",xy=(p, en-ly), fontsize=12)
    ax4=plt.axes([p, en-2*ly, lx, ly])
    ax4.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.7)
    ax4.annotate("n>=10",xy=(p, en-2*ly), fontsize=12)
    ax4=plt.axes([p, en-3*ly, lx, ly])
    ax4.pie([1],labels=[""],colors=['white'],rotatelabels=True, wedgeprops={'linewidth': 1, 'edgecolor': "black"}, radius=0.5)
    ax4.annotate("n>=1",xy=(p, en-3*ly), fontsize=12)

    pdf.savefig()
    pdf.close()
    plt.close()


def module1(inputDirectory, outputDirectory):
    '''
    call SNVs

    :param inputDirectory: genome sequences directory
    :param outputDirectory
    :returns snpFile: snp_merged.tsv
    '''
    ## the absolute directory of the output filess
    path = os.path.abspath(outputDirectory)
    if not os.path.exists(path):
        os.mkdir(path)
    snpFile = os.path.join(path, 'snp_merged.tsv')
    seq_info = os.path.join(path, 'seq_info.tsv')
    ## genomeDirectory: the extracted genome sequences were stored in this directory temporarily
    genomeDirectory = os.path.join(path, 'genome')
    os.mkdir(genomeDirectory)
    ## extracted each genome sequence and stored as a single fasta file for a sequence
    extract_sequence(inputDirectory, genomeDirectory)
    ## get the undercalling variants files list
    files_path = get_file_list(genomeDirectory)
    ## get the length of the files list
    total = len(files_path)
    ## write the header line of the variant storage file, snp_merged.tsv
    header = '''Id\tDate\tCountry\tPosition\tRef\tAlt\n'''
    with open(snpFile, 'a') as fhand:
        fhand.write(header)
    ## pass_n: the number of sequence passing quality control
    pass_n = 0
    for file in files_path:
        flag = genome_quality_control(file)
        if flag == -1:
            continue
        else:
            basicmessage, splitfasta, tempdirectory = split_sequence(file, path)
            if basicmessage == -1:
                continue
            else:
                ## alig reads to the reference genome
                samFile = align(splitfasta, refsequence, tempdirectory)
                ## sort
                bamFile = sort(samFile, tempdirectory)
                ## add header
                vcfFile = mpileup(bamFile, tempdirectory)
                ## call variants
                filter, vcfSnpFile = call(vcfFile, tempdirectory, refsequence)
                if filter == -1:
                    os.system(f'rm -rf {tempdirectory}')
                    continue
                else:
                    pass_n = pass_n + 1
                    snpMutationInformation = snp_mutation_information(vcfSnpFile)
                    Id = basicmessage['Id']
                    Date = basicmessage['Date']
                    Country = basicmessage['Country']
                    with open(snpFile, 'a', encoding='utf-8') as fhand:
                        for snp in snpMutationInformation:
                            Position = snp['Position']
                            Ref = snp['Ref']
                            Alt = snp['Alt']
                            record = Id + "\t" + Date + "\t" + Country + "\t" + str(Position) + "\t" + Ref + "\t" + Alt +'\n'
                            fhand.write(record)
                    os.system('rm -rf %s' % tempdirectory)

    os.system(f'rm -rf {genomeDirectory}')

    not_pass_n = total - pass_n
    with open(seq_info, 'a') as fhand:
        line1 = f'Total genome sequences: {total}\n'
        line2 = f'Pass quality control: {pass_n}\n'
        line3 = f'Not pass quality control: {not_pass_n}\n'
        fhand.write(line1)
        fhand.write(line2)
        fhand.write(line3)
    print('Have successfully obtained SNV mutations information.')
    return snpFile

def module2(file, directory, sites=None, frequency=None):
    '''
    haplotype and plot

    :param file: snp_merged.tsv
    :param directory: output directory
    :param sites: sites of intrest
    :param frequency: sites with frequency less than the given frequency will be filtered out
    '''
    ## the absolute path of the output directory
    path = os.path.abspath(outputDirectory)
    if not os.path.exists(path):
        os.mkdir(path)
    directory = os.path.abspath(directory)
    ## filter snps
    print('Dealing with snp sites...')
    snp_position, snp_ref_alt = snp_filter(file, directory, sites=sites, fre=frequency)
    print('Done!')
    ## get the reference haplotype sequence and the haplotype sequence of each genome
    print('Obtaining haplotype sequences of each genome sequence...')
    ref_haplotype_sequence = ref_haplotype(positions=snp_position)
    dataFile = genome_haplotype(file, snp_position, ref_haplotype_sequence, snp_ref_alt, directory)
    print('Done!')
    ## get the block.txt file
    print('Obtaining block.txt file...')
    blockFile = block_file(snp_position, directory)
    print('Done!')
    ## get the snp.info file
    print('Obtaining map file...')
    mapFile = map_file(snp_ref_alt, directory)
    print('Done!')
    ## get the snp.ped file
    print('Obtaining snp.ped file...')
    pedFile = ped_file(dataFile, directory)
    print('Done!')
    ## linkage analysis
    print('Linkage analyzing...')
    haplotypesFile = linkage_analysis(pedFile, mapFile, blockFile, directory)
    os.system(f'rm -rf {mapFile} {pedFile} {blockFile}')
    print('Done!')
    ## haplotyping
    print('Haplotyping......')
    data_plot, haplotypes = haplotyper(haplotypesFile, dataFile, directory)
    print('Done!')
    ## visualization
    print('Plotting...')
    plot(data_plot, haplotypes, directory)
    print('Done!')

if __name__ == '__main__':
    ## parse the arguments
    args = arguments()
    ## func: analysis mode
    func = args.modes
    inp = os.path.abspath(args.input)
    outputDirectory = os.path.abspath(args.output)
    sites = args.sites
    fre = args.frequency
    ## check if the argument is suitable
    if not os.path.exists(inp):
        print(f"ERROR: {inp} doesn't exists.")
        sys.exit()
    if fre is None:
        pass
    else:
        if (fre<0 or fre>1):
            print('Frequency should be greater than 0 and less than or equal to 1.')
            sys.exit()

    if sites is None:
        pass
    else:
        if len(sites) <= 1:
            print('Too little sites.')
            sys.exit()

        if len(sites) > 1:
            sites.sort()
            under = sites[0]
            top = sites[-1]
            if (under<=0 or top>29903):
                print('Positions are out of range (should be bigger than 0 and no more than 29903).')
                sys.exit()

    if func == 'auto':
        snpFile = module1(inp, outputDirectory)
        module2(snpFile, outputDirectory, sites=sites, frequency=fre)
        print('Done!')
        sys.exit()
    elif func == 'individual':
        module2(inp, outputDirectory, sites=sites, frequency=fre)
        sys.exit()
    else:
        print("Invalid mode!")
        sys.exit()


