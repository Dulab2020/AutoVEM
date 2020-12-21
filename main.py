#!/env/bin/python

'''
author='xibinbin'
'''

import argparse
import os
import sys
import re
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import datetime
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors

dulab1Directory = os.path.split(os.path.abspath(sys.argv[0]))[0]
dataDirectory = os.path.join(dulab1Directory, 'data')
toolsDirectory = os.path.join(dulab1Directory, 'tools')

ploidy = os.path.join(dataDirectory, 'ploidy.txt')
refsequence = os.path.join(dataDirectory, 'reference.fa')
haploview = os.path.join(toolsDirectory, 'Haploview.jar')

# print(dulab1Directory)
# print(dataDirectory)
# print(toolsDirectory)
# print(ploidy)
# print(refsequence)
# print(haploview)

def arguments():
    dulab1_parser = argparse.ArgumentParser(prog='AutoVEM', 
                                            description='A tool for analyzing SARS-Cov-2 genome mutations and epidemic trends', 
                                            usage='python main.py modes --input [--sites] [--frequency] --output')

    dulab1_parser.add_argument('modes', choices=['auto', 'individual'], nargs=1, help='Analysis mode')

    dulab1_parser.add_argument('--input', type=str, nargs=1,
                                        help='In auto mode, this should be a directory contains genome sequences. In individual mode, this should be the snp_merged.tsv file produced by dulab1')

    dulab1_parser.add_argument('--sites', type=int, nargs='+',
                                          help = 'Sites you interest(1-29903), should more than one site')

    dulab1_parser.add_argument('--frequency', type=float, nargs=1,
                                              help = 'Analysis SNVs with mutation frequency equal to or more than given frenquency(0-1)')
    
    dulab1_parser.add_argument('--output', type=str, nargs=1,
                                         help='Output directory, empty or not exists')
    
    args = dulab1_parser.parse_args()


    return args

def extract_sequence(directory, genomeDirectory):
    '''
    split fasta file

    :param directory: input directory
    :param genomeDirectory: genome directory
    '''
    files = os.listdir(directory)
    files_list = []
    for file in files:
        file = os.path.join(directory, file)
        files_list.append(file)
    
    pattern_head = r'^>hCoV-19'
    pattern_id = r'EPI_ISL_\d*'
    id_path = ''
    pointer = 0
    for file in files_list:
        with open(file, 'r') as fhand:
            for line in fhand.readlines():
                line = line.rstrip()
                header = re.findall(pattern_head, line)
                len_header = len(header)
                if len_header == 1:
                    id = re.findall(pattern_id, line)[0]
                    id = id + '.fa'
                    id_path = os.path.join(genomeDirectory, id)
                    if os.path.exists(id_path):
                        # pointer：file already exists?
                        # yes/no：1/0
                        pointer = 1
                        continue
                    with open(id_path, 'a') as fhand_sequence:
                        pointer = 0
                        fhand_sequence.write(line+'\n')
                if len_header == 0:
                    if pointer == 1:
                        continue
                    with open(id_path, 'a') as fhand_sequence:
                        fhand_sequence.write(line+'\n')

def get_file_list(directory):
    '''
    genome sequences

    :param directory: genome directory
    :returns filesPath(list, str):
    '''
    # fileslist：待分析的genome文件名列表
    # filespath：待分析genome的绝对路径文件名列表
    filesPath = []
    filesList = os.listdir(directory)
    for path in filesList:
        path = os.path.join(directory, path)
        filesPath.append(path)
    return filesPath
 
def create_index(path):
    '''
    create index

    :param path: index directory
    '''

    os.system('bowtie2-build -f %s %s/index' % (refsequence, path))

def genome_quality_control(file):
    '''
    quality control for the first time

    :param file: a fasta format file of a genome sequence
    :returns flag(int): （-1/0，not pass/pass）
    '''
    with open(file, 'r', encoding='utf-8') as fhand:
        sequence = ''
        number_N = 0
        number_db = 0
        pattern = ['A', 'T', 'G', 'C', 'N']
        i = 0
        for line in fhand.readlines():
            i = i + 1
            if i == 1:
                continue
            line = line.rstrip()
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
                continue

    if len_sequence < 29000:
        flag = -1
        return flag
    if number_N > 15 or number_db > 50:
        flag = -1
        return flag
    else :
        flag = 0
        return flag

def split_sequence(file, path):
    '''
    split genome sequence to reads(30-100nt)

    :param file: genome sequence
    :param path: output directory
    :returns bassicMessage(dict), splitFile(file), tempAnalysisDirectory(temp directory):
    '''
    China = ['taiwan', 'guangzhou', 'fujian', 'sichuan','wuhan', 'hangzhou', 'jiangsu', 
             'shanghai','shandong', 'guangdong', 'foshan', 'nanchang', 'yingtan','hunan', 
             'shangrao', 'yichun', 'beijing','jingzhou','pingxiang', 'shenzhen', 'lishui', 
             'zhejiang', 'fuzhou', 'shaoxing', 'xinyu','yunnan', 'jiujiang', 'chongqing', 
             'henan', 'hefei', 'fuyang', 'changzhou', 'jiangxi', 'ganzhou', 'hong kong', 
             'macau', 'qingdao', 'liaoning', 'harbin', 'tianmen', 'jian']   
    delete = ['lion', 'cat', 'env', 'canine', 'tiger', 'mink']

    subSequence = []
    with open(file,'r') as fhand:
        for line in fhand.readlines():
            line = line.rstrip()
            subSequence.append(line) 
    
    if len(subSequence[-1]) < 30:
        subSequence[-2] = subSequence[-2] + subSequence[-1]
        subSequence.pop()

    # Id Region Date
    header = subSequence[0]
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
    pattern_region = r'(?<=>hCoV-19/)([a-zA-Z\s]*)/'
    match_region = re.findall(pattern_region, header)
    len_match = len(match_region)
    if len_match == 0:
        return -1, -1, -1 
    else :
        region = match_region[0]
    
    basicMessage = {}
    basicMessage['Id'] = id
    temp = region.casefold()
    if temp in China:
        region = 'China'
    if temp in delete:
        return -1, -1, -1
        
    basicMessage['Country'] = region
    basicMessage['Date'] = date
    splitFastaFileName = basicMessage['Id'] + '_' + 'split' +'.fa'
    tempAnalysisDirectory = os.path.join(path, id)
    os.mkdir(tempAnalysisDirectory)
    splitedFile = os.path.join(tempAnalysisDirectory, splitFastaFileName)
    with open(splitedFile, 'a') as fhand:
        for i, subsequence in enumerate(subSequence):
            if i == 0:
                pass
            else :
                readName = ">" + str(i) + "." + header[1:] + "\n" + subsequence
                fhand.write("\n" + readName)

    return basicMessage, splitedFile, tempAnalysisDirectory

def align(file, directory1, directory2):
    '''
    alignment

    :param file: XXX_split.fa
    :param directory1: index directory
    :param directory2: temp directory
    :returns samFile: temp.sam
    '''
    samFile = os.path.join(directory2, 'temp.sam')
    os.system('bowtie2 -f -x %s/index -U %s -S %s' % (directory1, file, samFile))
    return samFile

def sort(file, directory):
    '''
    sort reads

    :param file: temp.sam
    :param directory: temp directory
    :returns bamFile: temp.bam
    '''
    bamFile = os.path.join(directory, 'temp.bam')
    os.system('samtools sort %s > %s' % (file, bamFile))
    return bamFile

def mpileup(file, directory):
    '''
    bam format -> vcf format

    :param file: temp.sam
    :param directory: temp directory
    :returns vcfFile: temp.vcf
    '''
    vcfFile = os.path.join(directory, 'temp.vcf')
    os.system('bcftools mpileup %s --fasta-ref %s > %s' % (file, refsequence, vcfFile))
    return vcfFile

def call(vcfFile, directory):
    '''
    call SNPs

    :param id: genome sequence id
    :param vcfFile: temp.vcf
    :param directory: temp directory
    :returns flag, snpFile: -1|0 pass|not pass
    '''
    snp_indel_file_path = os.path.join(directory, 'snp_indel.vcf')
    snp_file_path = os.path.join(directory, 'snp.vcf')
    indel_file_path = os.path.join(directory, 'indel.vcf')

    os.system('bcftools call --ploidy-file %s -vm %s -o %s' % (ploidy, vcfFile, snp_indel_file_path))

    os.system('vcftools --vcf %s --recode --keep-only-indels --stdout > %s' % (snp_indel_file_path, indel_file_path))
    os.system('rm -rf out.log')
    with open(indel_file_path, 'r') as fhand:
        n_indels = 0
        for line in fhand.readlines():
            if line[0] == '#':
                continue
            else :
                n_indels = n_indels + 1
    if n_indels >= 3:
        return -1, -1
    else :
        os.system('vcftools --vcf %s --recode --remove-indels --stdout > %s' % (snp_indel_file_path, snp_file_path))
        os.system('rm -rf out.log')
        return 0, snp_file_path

def snp_mutation_information(file, directory):
    '''
    obtain SNPs information

    :param file: snp.vcf
    :param directory: temp directory
    :returns mutationInformation(dict): key=['Position', 'Ref', 'Alt']
    '''
    tempFile = os.path.join(directory, 'temp.txt')
    os.system('awk \'$1 ~ /^NC_045512.2/ {print $2, $4, $5}\' %s > %s' % (file, tempFile))
    size = os.path.getsize(tempFile)

    mutationInformation = []
    if size == 0:
        mutationInformation.append({'Position': 0, 'Ref':'NA', 'Alt':'NA'})
    else:    
        with open(tempFile, 'r') as fhand:
            for line in fhand.readlines():
                snpMutationMessage = dict()
                line = line.rstrip()
                line = line.split()
                pos = line[0]
                ref = line[1]
                alt = line[2]
                snpMutationMessage['Position'] = pos
                snpMutationMessage['Ref'] = ref
                snpMutationMessage['Alt'] = alt
                mutationInformation.append(snpMutationMessage)

    return mutationInformation

def snp_filter(file, directory, sites, fre):
    '''
    filter SNP sites

    :param file: snp_merged.tsv
    :param directory: output directory
    :param sites: snp sites thar we are interested
    :param fre: frequency of snp sites that more than fre will be obtained
    :returns snp_pos, snp_ref_alt(dict)
    '''
    print(sites, fre)
    
    snp_sites = os.path.join(directory, 'snp_sites.tsv')
    os.system('touch %s'%snp_sites)

    df = pd.read_csv(file, sep='\t')
    ids = df['Id'].unique().tolist()
    n_genome = len(ids)

    counts = df['Position'].value_counts()
    frequency = counts/n_genome
    frequency = frequency.round(decimals=4)
    snp_dict = dict(zip(frequency.index.tolist(), frequency.values.tolist()))
    snp_dict.pop(0)
    for item in list(range(1, 29904)):
        if item not in snp_dict:
            snp_dict[item] = 0

    snp_pos = []
    if len(sites)==0:
        if fre==0:
            for key,item in snp_dict.items():
                if item>=0.05:
                    snp_pos.append(key)
        else:
            for key,item in snp_dict.items():
                if item>=fre:
                    snp_pos.append(key)
    else:
        if fre==0:
            for item in sites:
                tmp = snp_dict[item]
                if tmp > 0:
                    snp_pos.append(item)
        else:
            for item in sites:
                tmp = snp_dict[item]
                if tmp>=fre:
                    snp_pos.append(item)
    if len(snp_pos)<=1:
        print('There are no or too little sites that meet your requirements.')
        sys.exit()

    snp_pos.sort()

    snp_ref_alt = dict()
    for item in snp_pos:
        print(item)
        df2 = df[df['Position']==item]
        Ref = df2['Ref'].value_counts().index.tolist()[0]
        Alt = df2['Alt'].value_counts().index.tolist()[0]
        snp_ref_alt[item] = (Ref, Alt, snp_dict[item])
    header = 'Position\tRef\tAlt\tFrequency\n'
    with open(snp_sites, 'a') as fhand:
        fhand.write(header)
        for key, item in snp_ref_alt.items():
            record = str(key)+'\t'+item[0]+'\t'+item[1]+'\t'+str(item[2])+'\n'
            fhand.write(record)
        
    return snp_pos, snp_ref_alt

def ref_haplotype(position):
    '''
    reference haplotype sequence

    :param position：snp position
    :returns referenceHaplotype：reference haplotype sequence
    '''
    with open(refsequence, 'r') as fhand:
        referenceGenomeSequence = ''
        i = 0
        for line in fhand.readlines():
            i = i + 1
            line = line.rstrip()
            if i == 1:
                continue
            else :
                referenceGenomeSequence = referenceGenomeSequence + line
    referenceHaplotype = ''
    for item in position:
        referenceHaplotype = referenceHaplotype + referenceGenomeSequence[item-1]
    return referenceHaplotype

def genome_haplotype(file, position, positionAlt, referenceHaplotype, directory):
    '''
    get haplotype sequence of genome

    :param file: snp_merged.tsv
    :param position
    :param positionAlt
    :param referenceHaplotype: reference haplotype sequence
    :param directory: output directory
    :returns filePath: data.tsv
    '''
    filePath = os.path.join(directory, 'data.tsv')

    df = pd.read_csv(file, sep='\t')
    id = df['Id'].unique().tolist()

    Date = []
    Country = []
    Case_id = []
    snp_positions = []

    for item in id:
        df1 = df[df['Id']==item]
        case = item
        date = df1['Date'].value_counts().index.tolist()[0]
        country = df1['Country'].value_counts().index.tolist()[0]
        Case_id.append(case)
        Date.append(date)
        Country.append(country)
        mutation_positions = df1['Position'].tolist()
        mutation_positions = [int(x) for x in mutation_positions]
        snp_positions.append(mutation_positions)

    haplotypes = []
    for item in snp_positions:
        sites = ''
        i = -1
        for site in position:
            i = i + 1
            site = int(site)
            if site in item:
                sites = sites + positionAlt[site][1]
            else :
                sites = sites + referenceHaplotype[i]
        haplotypes.append(sites)
    data = pd.DataFrame(data={'Id': Case_id,
                                   'Date': Date,
                                   'Country':Country,
                                   'Hap': haplotypes})
    data.to_csv(filePath, sep='\t', index=False)
    return filePath

def block_file(position, directory):
    '''
    block.txt

    :param position: snp_pos
    :param directory: output directory
    :returns blockFile：block.txt
    '''
    num = len(position)
    blockFile = os.path.join(directory, 'block.txt')
    with open(blockFile, 'a') as fhand:
        for i in list(range(num)):
            fhand.write(str(i+1))
            fhand.write('\t')
    
    return blockFile

def map_file(position, directory):
    '''
    snp.info

    :param position: snp_ref_alt
    :returns mapFile: snp.info
    '''
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
    pedFile = os.path.join(directory, 'snp.ped')
    with open(pedFile, 'a') as f:
        with open(file, 'r') as fhand:
            i = 0
            for line in fhand.readlines():
                i = i + 1
                if i == 1:
                    continue
                line = line.rstrip()
                line = line.split('\t')
                record = line[0]+'\t'+line[0]+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'
                for item in line[-1]:
                    record = record + '\t' + item + '\t' + item
                record = record + '\n'
                f.write(record)
    
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

    os.system('java -jar %s -n -skipcheck -pedfile %s -info %s -blocks %s -png -out %s' % (haploview, ped, mapf, block, temp))

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
    
    hap_dict = dict()
    with open(file, 'r') as fhand:
        i = -1
        for line in fhand.readlines():
            sequence = ''
            i = i + 1
            if i == 0:
                continue 
            hap = 'H' + str(i)
            sequence_list = line.split()
            sequence_num = sequence_list[0]
            for num in sequence_num:
                sequence = sequence + nt_dict[num]
            hap_dict[sequence] = hap

    header = 'Id\tDate\tCountry\tHap\tName\n'
    with open(dataPlot, 'a') as f:
        f.write(header)
        with open(dataFile, 'r') as fhand:
            i = 0
            for line in fhand.readlines():
                i = i + 1
                if i == 1:
                    continue
                line = line.rstrip()
                line = line.split()
                hap = line[-1]
                record = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\t'
                if hap in hap_dict:
                    name = hap_dict[hap]
                else:
                    name = 'other'
                record = record + name + '\n'
                f.write(record)
    
    with open(haplotypes, 'a') as fhand:
        for key, value in hap_dict.items():
            tem = ''
            tem = value + '\t' + key + '\n'
            fhand.write(tem)

    return dataPlot, haplotypes

def plot(file, file2, directory):
    '''
    结果可视化

    :param file: data_plot.tsv
    :param directory: plot directory
    :param file2: haplotypes_temp.tsv
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

    :param inputDirectory: genome sequences
    :param outputDirectory
    :returns snpFile: snp_merged.tsv
    '''
    path = os.path.abspath(outputDirectory)
    if not os.path.exists(path):
        os.mkdir(path)
    snpFile = os.path.join(path, 'snp_merged.tsv')
    seq_info = os.path.join(path, 'seq_info.tsv')
    indexDirectory = os.path.join(path, 'index')
    os.mkdir(indexDirectory)
    genomeDirectory = os.path.join(path, 'genome')
    os.mkdir(genomeDirectory)

    extract_sequence(inputDirectory, genomeDirectory)
    files_path = get_file_list(genomeDirectory)
    total = len(files_path)

    create_index(indexDirectory)

    header = '''Id\tDate\tCountry\tPosition\tRef\tAlt\n'''
    with open(snpFile, 'a') as fhand:
        fhand.write(header)

    pass_n = 0
    for file in files_path:
        flag = genome_quality_control(file)
        if flag == -1:
            continue
        else :
            basicmessage, splitfasta, tempdirectory = split_sequence(file, path)
            if basicmessage == -1:
                os.system('rm -rf %s' % tempdirectory)
                continue
            else :
                samFile = align(splitfasta, indexDirectory, tempdirectory)
                bamFile = sort(samFile, tempdirectory)
                vcfFile = mpileup(bamFile, tempdirectory)
                filter, vcfSnpFile = call(vcfFile, tempdirectory)
                if filter == -1:
                    os.system('rm -rf %s' % tempdirectory)
                    continue

                pass_n = pass_n + 1
                snpMutationInformation = snp_mutation_information(vcfSnpFile, tempdirectory)
                snp_mutation = {}
                snp_mutation['Id'] = basicmessage['Id']
                snp_mutation['Date'] = basicmessage['Date']
                snp_mutation['Country'] = basicmessage['Country']

                with open(snpFile, 'a', encoding='utf-8') as fhand:
                    for snp in snpMutationInformation:
                        snp_mutation['Position'] = snp['Position']
                        snp_mutation['Ref'] = snp['Ref']
                        snp_mutation['Alt'] = snp['Alt']
                        record = snp_mutation['Id'] + "\t" + snp_mutation['Date'] + "\t" + snp_mutation['Country'] + "\t" + str(snp_mutation['Position']) + "\t" + snp_mutation['Ref'] + "\t" + snp_mutation['Alt'] +'\n'
                        fhand.write(record)
        
        os.system('rm -rf %s' % tempdirectory)
    os.system('rm -rf %s %s'%(indexDirectory, genomeDirectory))

    not_pass_n = total - pass_n
    with open(seq_info, 'a') as fhand:
        line1 = 'Total genome sequences: %s' % total + '\n'
        line2 = 'Pass quality control: %s' % pass_n + '\n'
        line3 = 'Not pass quality control: %s' % not_pass_n + '\n'
        fhand.write(line1)
        fhand.write(line2)
        fhand.write(line3)

    print('Have successfully obtained SNV mutations information.')
    return snpFile

def module2(file, directory, sites, frequency):
    '''
    haplotype and plot

    :param file: snp_merged.tsv
    :param directory: output directory
    '''
    path = os.path.abspath(outputDirectory)
    if not os.path.exists(path):
        os.mkdir(path)
        
    print('Dealing with snp sites...')
    snp_position, snp_ref_alt = snp_filter(file, directory, sites=sites, fre=frequency)
    print('Done!')

    print('Obtaining data.tsv file...')
    ref_haplotype_sequence = ref_haplotype(snp_position)
    dataFile = genome_haplotype(file, snp_position, snp_ref_alt, ref_haplotype_sequence, directory)
    print('Done!')

    print('Obtaining block.txt file...')
    blockFile = block_file(snp_position, directory)
    print('Done!')

    print('Obtaining map file...')
    mapFile = map_file(snp_ref_alt, directory)
    print('Done!')

    print('Obtaining snp.ped file...')
    pedFile = ped_file(dataFile, directory)
    print('Done!')

    print('Linkage analyzing...')
    haplotypesFile = linkage_analysis(pedFile, mapFile, blockFile, directory)
    os.system('rm -rf %s %s %s'%(mapFile, pedFile, blockFile))
    print('Done!')
    
    print('Haplotyping......')
    data_plot, haplotypes = haplotyper(haplotypesFile, dataFile, directory)
    print('Done!')

    print('Plotting...')
    plot(data_plot, haplotypes, directory)
    print('Done!')


args = arguments()
func = args.modes[0]
inp = args.input[0]
inp = os.path.abspath(inp)
outputDirectory = args.output[0]
outputDirectory = os.path.abspath(outputDirectory)

if not args.sites:
    sites = []
else:
    sites = args.sites
if not args.frequency:
    frequency=0
else:
    frequency=args.frequency[0]

if (frequency<0) or (frequency>1):
    print('Frequency should be greater than 0 and less than or equal to 1.')
    sys.exit()

if len(sites)==1:
    print('Too little sites.')
    sys.exit()

if len(sites)>1:
    sites.sort()
    under = sites[0]
    top = sites[-1]
    if under<=0 or top>29903:
        print('Site(s) is(are) out of range.')
        sys.exit()

print(args)
print(func)
print(inp)
print(sites)
print(frequency)
print(outputDirectory)

if __name__ == '__main__':
    if func == 'auto':
        snpFile = module1(inp, outputDirectory)
        module2(snpFile, outputDirectory, sites, frequency)
        print('Done!')
    if func == 'individual':
        module2(inp, outputDirectory, sites, frequency)
        print('Done!')


