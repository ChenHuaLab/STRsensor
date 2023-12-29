#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

#*************************************************************************
#    > File Name: STRSimulatorV4.py
#    > Author: xlzh
#    > Mail: xiaolongzhang2015@163.com 
#    > Created Time: 2020年09月18日 星期五 14时40分20秒
#*************************************************************************

import sys
import os.path
import random
from random import randint
import codecs
import gzip

EPS = 1e-4
READLEN = 150

class STRLocus(object):
    def __init__(self, chrom, start, end, motif_len, exclude_base, is_haplotype, shift_base, fa_obj):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.motif_len = motif_len
        self.exclude_base = exclude_base
        self.is_haplotype = is_haplotype
        self.shift_base = shift_base
        self.ref_seq = fa_obj.get_seq(chrom, start, end)
        self.flank5 = fa_obj.get_seq(chrom, start-READLEN, start-1)
        self.flank3 = fa_obj.get_seq(chrom, end+1, end+READLEN)


class Fasta(object):
    def __init__(self, fa_file):
        self.fa_file = fa_file
        self.__fa_fp = self.__fa_init()
        self.__fai_dict = self.__fai_load()

    def __fa_init(self):
        # Func: open the fasta handle
        if not os.path.exists(self.fa_file):
            sys.stderr.write('Err: No such file: %s\n' %(self.fa_file))
            sys.exit(-1)

        return open(self.fa_file, 'r')

    def __fai_load(self):
        # Func: load the index of fasta file
        if not os.path.exists(self.fa_file + '.fai'):
            sys.stderr.write('Err: No such index file: %s\n' %(self.fa_file+'.fai'))
            sys.exit(-1)

        fai_dict = {} # {'chr1':(offset, b_len, l_len), 'chr2':(), ...}

        fai_fp = open(self.fa_file+'.fai', 'r')
        for line in fai_fp:
            l = line.split()
            fai_dict[l[0]] = (int(l[2]), int(l[3]), int(l[4]))

        return fai_dict

    def __read_base(self):
        # Func: read 1-base each time
        base = self.__fa_fp.read(1)
        if base == '': return 0

        return base

    def get_seq(self, chrom, start, end):
        # Func: get the fasta sequence from start to end
        seq = []

        if (chrom not in self.__fai_dict) or (end < start):
            sys.stderr.write('Warning: Invaild input: %s:%d-%d\n' %(chrom,start,end))
            return ''

        idx = self.__fai_dict[chrom]
        d_offset = idx[0] + (start-1) + (start-1)/idx[1]*(idx[2]-idx[1])

        self.__fa_fp.seek(d_offset) # put the file handle to the start position
        ch, seq_len = 1, end-start+1

        while ch and seq_len:
            ch = self.__read_base()
            if ch not in ['\r', '\n', 0]: seq.append(ch); seq_len -= 1

        return ''.join(seq)


class FqWrite(object):
    # mode = "gz": compressed fq.gz file
    def __init__(self, outBase, mode=''):
        self.out_base = outBase
        self.mode = mode

        if self.mode == "gz":
            self.Read1 = gzip.open(outBase + "_R1.fq.gz", "wb")
            self.Read1 = codecs.getwriter('utf-8')(self.Read1)
        else:
            self.Read1 = open(outBase + "_R1.fq", "w")

    def write(self, group):
        # group = [seqname1, seq1, '+', qual1]
        for i in range(4):
            self.Read1.write("%s\n" % group[i])

    def close(self):
        self.Read1.close()


def read_region(region_file, fa_file):
    # region_dict = {'DYS19':locus_obj, ...}
    region_dict = {}
    region_fp = open(region_file, 'r')

    fa_obj = Fasta(fa_file)
    line = region_fp.readline() # remove header

    for line in region_fp:
        l = line.rstrip().split()
        is_haplotype = True if l[6] == "Yes" else False
        region_dict[l[0]] = STRLocus(l[1], int(l[2]), int(l[3]), int(l[4]), int(l[5]), is_haplotype, int(l[7]), fa_obj)

    return region_dict


def read_frequency(freq_file, region_locus):
    ''' freq_dict = {
            'DYS19': {14:0.0022, 15:0.25, 16:0.28, 17:0.24, 17.3:0.2278}
            'DYS391':{11.3:0.003, 12:0.35, 12.3:0.0043}
        }
    '''
    freq_dict = {}

    for t in region_locus: # set locus in the target locus
        freq_dict[t] = {}

    if not freq_file: # frequency file is not given
        sys.stderr.write("[Error] frequency file of %s is not given!\n" % freq_file)
        sys.exit(-1)
    
    freq_fp = open(freq_file, 'r')

    for line in freq_fp:
        l = line.rstrip().split()
        if l[0] not in freq_dict: continue # only keep the freq occured in the region_locus
        for af in l[1:]: freq_dict[l[0]][float(af.split(':')[0])] = float(af.split(':')[1])

    return freq_dict


def read_stutter(stutter_file, region_locus):
    ''' stutter_dict = {
            'STRName': [normal, add_1, add_2, delete_2, delete_1]
            'DYS19': [0.920079, 0.005710, 0.000676, 0.004641, 0.068894]
        }
    '''
    stutter_dict = {}

    for t in region_locus: # set locus in the target locus
        stutter_dict[t] = []

    if not stutter_file: # stutter file is not given
        sys.stderr.write("[Error] stutter file of %s is not given!\n" % stutter_dict)
        sys.exit(-1)
    
    stutter_fp = open(stutter_file, 'r')

    for line in stutter_fp:
        l = line.rstrip().split()
        if l[0] not in stutter_dict: continue # only keep the stutter occured in the region_locus
        for lp in l[1:]: stutter_dict[l[0]].append(float(lp))

    return stutter_dict


def get_allele(locus_freq):
    ''' locus_freq = {14.0:0.0022, 15.0:0.25, 16.3:0.28, 17.0:0.24, 17.3:0.2278}
    '''
    freq_table = []
    freq_list = list(locus_freq.values())
    allele_list = list(locus_freq.keys())

    for i in range(1, len(freq_list)+1): 
        freq_table.append(sum(freq_list[:i]))

    freq_table[-1] = 1.0
    prob = random.random() # decide which allele block to choose

    for i in range(len(freq_list)):
        if freq_table[i] > prob: return allele_list[i]

    return allele_list[-1]


def true_allele(cur_allele, motif_len, exclude_base):
    ''' Func: get allele which has excluded the NOT COUNTED bases
        For allele 15.3 -> allele_pre: 15 & allele_rest: 3
    '''
    allele_pre = int(cur_allele)
    allele_rest = int(str(cur_allele).split('.')[1])

    cur_allele_base = allele_pre * motif_len + allele_rest
    true_allele_base = cur_allele_base - exclude_base

    true_pre = int(true_allele_base / motif_len)
    true_rest = true_allele_base % motif_len

    return true_pre + 0.1 * true_rest


def generate_allele(freq_dict, region_dict, person_id):
    ''' allele_dict = {'DYS391':15.0, 'DYS456':18.0, DS1356:(15.0, 17.3), ...}
    '''
    allele_dict = {}

    summary_fp = open("Summary.txt", "a")
    summary_fp.write("%s\t" % person_id) # write the person ID

    for str_name in region_dict:
        reg_obj = region_dict[str_name]
        m_len, e_base = reg_obj.motif_len, reg_obj.exclude_base

        if reg_obj.is_haplotype:
            allele = get_allele(freq_dict[str_name])
            allele_dict[str_name] = allele
            summary_fp.write("%s:%.1f\t" % (str_name, true_allele(allele, m_len, e_base)))
        else:
            allele1, allele2 = get_allele(freq_dict[str_name]), get_allele(freq_dict[str_name])
            true1, true2 = true_allele(allele1, m_len, e_base), true_allele(allele2, m_len, e_base)
            allele_dict[str_name] = (allele1, allele2)

            if true1 <= true2: summary_fp.write("%s:%.1f|%.1f\t" % (str_name, true1, true2))
            else: summary_fp.write("%s:%.1f|%.1f\t" % (str_name, true2, true1))

    summary_fp.write("\n"); summary_fp.close()

    return allele_dict


def gen_str_seq(str_obj, read_pre, read_rest):
    #                     start                           end
    #                      |         reference STR         |
    #     -----------------AGATAGAT AGATAGATAGATCCTGCCTGCCTG-----------------
    #              shift-bases(8bp) | <-- add/remove unit
    #
    whole_seq = ''
    m_len = str_obj.motif_len
    f5_seq, ref_seq, f3_seq = str_obj.flank5, str_obj.ref_seq, str_obj.flank3

    # get the reference and read allele
    ref_len = len(ref_seq)
    ref_pre, ref_rest = int(ref_len / m_len), ref_len % m_len
    ref_allele_count = ref_pre + 0.1 * ref_rest
    read_allele_count = read_pre + 0.1 * read_rest

    # generate read sequence by allele count and outframe shift base
    if read_allele_count == ref_allele_count:
        whole_seq = f5_seq + ref_seq + f3_seq

    elif read_allele_count < ref_allele_count: # remove unit from ref_seq
        remove_base = (ref_pre-read_pre) * m_len + (ref_rest-read_rest)
        core_seq = ref_seq[:str_obj.shift_base] + ref_seq[remove_base+str_obj.shift_base:]
        whole_seq = f5_seq + core_seq + f3_seq

    elif read_allele_count > ref_allele_count: # add unit from ref_seq
        add_base = (read_pre-ref_pre) * m_len + (read_rest-ref_rest)
        candidate_seq = 50 * ref_seq[str_obj.shift_base:][:m_len]
        core_seq = ref_seq[:str_obj.shift_base] + candidate_seq[:add_base] + ref_seq[str_obj.shift_base:]
        whole_seq = f5_seq + core_seq + f3_seq

    return whole_seq


def generate_read(str_obj, str_stutter, a_count):
    #                     start                       end
    #                      |          repeat           |
    #     -----------------AGATAGATAGATAGATAGATAGATAGAT-----------------
    #      |--- 150 bp ---|                           |--- 150 bp ---|
    #
    # str_stutter = [normal, add1, add2, remove2, remove1]

    slip_table = [0, 1, 2, -2, -1]
    prob_list = [sum(str_stutter[:i+1]) for i in range(5)] # normal -> add -> remove
    frag_start = randint(str_obj.start-READLEN, str_obj.start+str_obj.motif_len*(int(a_count-2)))

    a_pre, a_rest = int(a_count), int(str(a_count).split('.')[1])
    prob = random.random()
    slip_unit = 0

    for i in range(5):
        if prob_list[i] > prob: slip_unit = slip_table[i]; break

    whole_seq = gen_str_seq(str_obj, a_pre+slip_unit, a_rest)
    str_start = frag_start - (str_obj.start - READLEN)
    read_seq = whole_seq[str_start:str_start+READLEN]

    if len(read_seq) != READLEN:
        sys.stderr.write("[Error] the length of sequence and quality is different!\n")
        sys.exit(-1)

    # generate random base sequencing error with prob of 0.001 [p(read) = 150 * 0.001 = 0.15]
    prob = random.random()
    sub_base = ['A', 'T', 'G', 'C']

    if prob < 0.15:
        change_idx = randint(0, len(read_seq)-1)
        change_base = sub_base[randint(0, 3)]
        read_seq = read_seq[:change_idx] + change_base + read_seq[change_idx+1:]

    return frag_start, read_seq


def check_str_name(region_dict, freq_dict, stutter_dict):
    status = False

    for k in region_dict:
        if (not freq_dict[k]) or (not stutter_dict[k]): 
            return (status, k)

    return (True, None)


def main():
    ''' Coverage Calculation: C = L * N / G
        C: stands for coverage
        L: stands for read length
        N: number of reads
        G: genome length for the reads covered
    '''
    args = sys.argv
    if len(args) != 7:
        sys.stderr.write("Usage: python STRSimulator.py <str.locus> <freq.txt> <stutter.txt> <genome.fa> <out_prefix> <depth>\n")
        sys.exit(-1)

    region_dict = read_region(args[1], args[4])
    freq_dict = read_frequency(args[2], region_dict)
    stutter_dict = read_stutter(args[3], region_dict)

    status = check_str_name(region_dict, freq_dict, stutter_dict)
    if status[0] is False:
        sys.stderr.write("[Error]: Locus of %s is not occured in frequnecy or stutter file!\n" % status[1])
        sys.exit(-1)

    # generate allele for each locus by their frequency
    fq_obj = FqWrite(args[5], "gz")
    allele_dict = generate_allele(freq_dict, region_dict, os.path.basename(args[5]))
    sys.stderr.write("[*] start to generate simulate STR reads for %s\n" % os.path.basename(args[5]))

    for str_name in allele_dict:
        coverage = int(int(args[6]) * 1.5) # N = 30 * (50 + 150) / 150 
        for i in range(coverage): 
            if region_dict[str_name].is_haplotype: # X or Y
                allele = allele_dict[str_name]
            else:
                allels = allele_dict[str_name]
                allele = allels[0] if random.random() > 0.5 else allels[1]

            f_start, f_seq = generate_read(region_dict[str_name], stutter_dict[str_name], allele)
            fq_obj.write(['@%s:%.1f_hg19:%d_rand:%d' % (str_name, allele, f_start, randint(1,10000)), f_seq, '+', 'I'*150])
            
    fq_obj.close()


if __name__ == '__main__':
    main()
