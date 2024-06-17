#!/usr/bin/env python

'''
find_genes_nearby_loc.py

Purpose:
    Find nearby gene names from a gene table with location information ([chr,] [strand,] start, end) by its location ([chr,] [strand,] start [, end])

License: GNU General Public License v2.0 (http://www.gnu.org/licenses/gpl-2.0.html)
Author: Dr. Xiao-Qin Xia
Email: xqxia70@gmail.com

Usage:
    find_genes_nearby_loc.py  [options]  parameters

Options:
    -h, --help:     Display this message, then exit
    --version:      Print version information
    -i, --input=:   Input tabular file (ITF). If not provided, pop the first value in
                    {parameters}, or STDIN
    -o, --output=:  Output tabular file. If not provided, pop the first value in
                    {parameters}, or STDOUT
    --gene-file=:   Mandatory. Gene tabular file (GTF)
    --gene-name=:   Mandatory. The column for gene name in the GTF
    --gene-chr=:    Optional. The column for chromosome name in the GTF
    --gene-strand=: Optional. The column for strand in the GTF
    --gene-start=:  Mandatory. The column for start location in the GTF
    --gene-end=:    Mandatory. The column for end location in the GTF
    --loc-chr=:     Optional. The column for chromosome name in the ITF
    --loc-strand=:  Optional. The column for strand in the ITF
    --loc-start=:   Mandatory. The column for start location in the ITF
    --loc-end=:     Optional. The column for end location in the GTF
    --adjacent=:    Optional. The default is 100. If other values provided as a
                    number XXX, in case that on gene found at the location, the
                    program will search around in adjacent area within XXX
                    bases. Search for nearest downstream and upstream genes.

parameters:
    Input file, Output file in order, in case that these augments are not provided as options.

'''

import os, sys, re
from xqlibs.base.basetools import parseOpts, exitMsg, NoOptionErr, NoColumnErr, TableFile, OutputFile, TableCols

import banyan

__version__ = '0.0.1'

def read_gene_loc_with_strand(gene_file, cols, i_chr=None, i_chain=None, byupdator=banyan.RankUpdator):
    '''
    cols:    [chr,] [strand,] start, [end,] name
    gene_file: A tabular file that should contain columns in cols.
    i_chr, i_chain: The index of chr and strand in cols. 
    return a dict: {'+': [a dict of] banyan.SortedDict((start, end): name)}
    '''
    def get_val_with_strand(line, i_chr, i_chain, key_idx=[], sort_loc_by_strand=True):
        '''
        line:    [chr,] [chain,] start, [end,], name
        return chr, strand, loc, gene
        '''
        chr = i_chr is not None and line[i_chr] or None
        strand = i_chain is not None and line[i_chain] or None 
        for i in key_idx:  # remove chr, strand from line
            line.pop(i)
        i = int(line[0])  # start
        if len(line) > 2:  # start, end, name
            j = int(line[1])
            if strand:
                if sort_loc_by_strand:
                    i = strand=='+' and (i < j and (i, j) or (j, i)) or (i < j and (j, i) or (i, j))
                else:  # make sure start < end
                    i = i > j and (j, i) or (i, j) 
            else:
                strand = i < j and '+' or '-'
                i = sort_loc_by_strand and (i, j) or i > j and (j, i) or (i, j)  
        return(chr, strand or '+', i, line[-1])
    key_idx = [i for i in [i_chr, i_chain] if i is not None]
    key_idx.sort(reverse=True)
    if i_chr is not None: 
        gene_dic = chr_dic = {'+': {}, '-': {}}
        for line in TableCols(gene_file, cols=cols):
            chr, strand, loc, gene = get_val_with_strand(line, i_chr, i_chain, key_idx)
            loc_dic = chr in chr_dic[strand] and chr_dic[strand][chr] or chr_dic[strand].setdefault(chr, banyan.SortedDict(updator=byupdator))
            loc_dic[loc] = gene
    else:
        gene_dic = loc_dic = {'+': banyan.SortedDict(updator=byupdator), '-': banyan.SortedDict(updator=byupdator)}
        for line in TableCols(gene_file, cols=cols):
            chr, strand, loc, gene = get_val_with_strand(line, i_chr, i_chain, key_idx)
            loc_dic[strand][loc] = gene
    return(gene_dic)

def tail_ahead(gene_dic, updator):
    if type(gene_dic) is dict:
        new_dic = {}
        for k, v in gene_dic.items():
            new_dic[k] = tail_ahead(v, updator)
    else:
        new_dic = banyan.SortedDict(updator=updator)
        for k, v in gene_dic.items():
            new_dic[(k[1], k[0])] = v
    return(new_dic)

def Main(fsrc, fobj, gene_file, col_gene_name=None, col_gene_chr=None, col_gene_strand=None, col_gene_start=None, col_gene_end=None, col_loc_chr=None, col_loc_strand=None, col_loc_start=None, col_loc_end=None, adjacent=100, sep='\t', joint_str=' / ', linesep=os.linesep):
    # determine what columns to read
    if not col_gene_name or not col_gene_start or not col_gene_end or not col_loc_start:
        raise NoOptionErr('Not all mandatory augments have been provided: gene_name, gene_start, gene_end, loc_start')
    strand_specified = bool(col_loc_strand) # find gene only on the same strand
    need_strand = bool(adjacent)  # need to know the real Start site, which is bigger then End on the negative strand.
    # if col_gene_strand is not specified, it will be determined by the compare the value of Start and End
    cols_gene, cols_input = [], []
    i_chr = i_chain_gene = i_chain_loc = None
    i = 0
    has_chr = col_gene_chr and col_loc_chr
    if has_chr:
        cols_gene.append(col_gene_chr)
        cols_input.append(col_loc_chr)
        i_chr = i
        i += 1
    if col_gene_strand:
        cols_gene.append(col_gene_strand)
        i_chain_gene = i
    if col_loc_strand:
        cols_input.append(col_loc_strand)
        i_chain_loc = i
    cols_gene.extend([col_gene_start, col_gene_end, col_gene_name])  # chr, [strand,] start, end, gene
    cols_input.append(col_loc_start)
    byupdator = banyan.RankUpdator
    by_empty = banyan.SortedDict(updator=byupdator)
    if col_loc_end:  
        cols_input.append(col_loc_end)
    # read gene info
    gene_dic_5 = read_gene_loc_with_strand(gene_file, cols_gene, i_chr, i_chain_gene, byupdator)
    gene_dic_3 = tail_ahead(gene_dic_5, byupdator)
    # convert colname to idx
    fsrc = TableFile(fsrc)
    cols_src = fsrc.next()
    key_idx = [i for i in [i_chr, i_chain_loc] if i is not None]
    key_idx.sort(reverse=True)
    try:
        cols_idx = map(cols_src.index, cols_input)
        if key_idx:  # convert to index in input file 
            # key_idx = map(cols_idx.__getitem__, key_idx)
            if i_chr is not None:
                i_chr = cols_idx[i_chr]
            if i_chain_loc is not None:
                i_chain_loc = cols_idx[i_chain_loc]
            for i in key_idx:
                cols_idx.pop(i)
    except:
        raise NoColumnErr('No necessary columns in input ' + repr(cols_input))
    # prepare output file
    cols_src.extend(["5'_"+col_gene_name, "5'_Shift", "3'_"+col_gene_name, "3'_Shift"])
    fobj = OutputFile(fobj)
    fobj.write(sep.join(cols_src)+linesep)
    gene_dic_5_pos, gene_dic_3_pos = gene_dic_5.get('+', by_empty), gene_dic_3.get('+', by_empty)
    gene_dic_5_neg, gene_dic_3_neg = gene_dic_5.get('-', by_empty), gene_dic_3.get('-', by_empty)
    if i_chr is None:
        loc_dic_5_pos, loc_dic_3_pos = gene_dic_5_pos, gene_dic_3_pos
        loc_dic_5_neg, loc_dic_3_neg = gene_dic_5_neg, gene_dic_3_neg

    # deal with input
    for line in fsrc:
        strand = i_chain_loc is not None and line[i_chain_loc] or '+'
        if i_chr is not None:
            chr = line[i_chr] #tuple(map(line.__getitem__, key_idx))
            loc_dic_5_pos, loc_dic_3_pos = gene_dic_5_pos.get(chr, by_empty), gene_dic_3_pos.get(chr, by_empty)
            loc_dic_5_neg, loc_dic_3_neg = gene_dic_5_neg.get(chr, by_empty), gene_dic_3_neg.get(chr, by_empty)
        loc = map(lambda i:int(line[i]), cols_idx)
        loc.sort()
        loc = tuple(loc)
        genes_5, shifts_5, genes_3, shifts_3 = [], [], [], []
        if not i_chain_loc or strand == '+':
            # for 5' regulation (+)
            loc_point = col_loc_end is None and loc or (loc[1],)
            k_loc = loc_dic_5_pos.order(loc_point)
            if k_loc < len(loc_dic_5_pos):
                key = loc_dic_5_pos.kth(k_loc)
                dif = key[0] - loc_point[0]
                if dif <= adjacent:  # satisfied
                    genes_5.append(loc_dic_5_pos[key])
                    shifts_5.append('-'+str(dif))
            # for 3' regulation (+)
            loc_point = col_loc_end is None and loc or (loc[0],)
            k_loc = loc_dic_3_pos.order(loc_point)
            if k_loc > 0:
                key = loc_dic_3_pos.kth(k_loc-1)
                dif = loc_point[0] - key[0]
                if dif <= adjacent:  # satisfied
                    genes_3.append(loc_dic_3_pos[key])
                    shifts_3.append(str(dif))
        if not i_chain_loc or strand == '-':
            # for 5' regulation (-)
            loc_point = col_loc_end is None and loc or (loc[0],)
            k_loc = loc_dic_5_neg.order(loc_point)
            if k_loc > 0:
                key = loc_dic_5_neg.kth(k_loc-1)
                dif = loc_point[0] - key[0]
                if dif <= adjacent:  # satisfied
                    genes_5.append(loc_dic_5_neg[key])
                    shifts_5.append(str(dif))
            # for 3' regulation (-)
            loc_point = col_loc_end is None and loc or (loc[1],)
            k_loc = loc_dic_3_neg.order(loc_point)
            if k_loc < len(loc_dic_3_neg):
                key = loc_dic_3_neg.kth(k_loc)
                dif = key[0] - loc_point[0]
                if dif <= adjacent:  # satisfied
                    genes_3.append(loc_dic_3_neg[key])
                    shifts_3.append('-'+str(dif))

        line.extend(map(joint_str.join, [genes_5, shifts_5, genes_3, shifts_3]))
        fobj.write(sep.join(line)+linesep)

if __name__ == '__main__':
    from getopt import getopt
    optlst, args = getopt(sys.argv[1:], 'hi:o:', ['help', 'version', 'input=', 'output=', 
        'gene-file=', 'gene-name=', 'gene-chr=', 'gene-strand=', 'gene-start=', 'gene-end=', 
        'loc-chr=', 'loc-strand=', 'loc-start=', 'loc-end=', 'adjacent='])

    # print help information
    optdic = dict(optlst)
    if '-h' in optdic or '--help' in optdic:
        doc = re.sub(r'^\n\S+|(?<=\nUsage:\n\t)\S+\b', os.path.basename(sys.argv[0]), __doc__)
        exitMsg(doc)
    if '--version' in optdic:
        exitMsg(__file__ + ': ' + __version__)

    err = []
    # options that can have multiple values, e.g., '--rep-opt'
    repopts = [] 
    # mandatory options that musted be supplied, e.g., ['--mand-opt', ('-c', '--col'), ...]
    mandopts = ['--gene-file', '--gene-name', '--gene-start', '--gene-end', '--loc-start'] 
    # options that will be automatically dealt with. e.g. [(int, ('-n', '--num'), 'num', 3), ...] -- (type, (option names), target variable names, default values). type can be int, float, str, bool, or any self-defined functions, e.g. lambda s:set(re.split(r'[,;|]', s)).
    autopts = [(str, '--gene-file', 'gene_file', None), (str, '--gene-name', 'col_gene_name', None), 
            (str, '--gene-chr', 'col_gene_chr', None), (str, '--gene-strand', 'col_gene_strand', None), 
            (str, '--gene-start', 'col_gene_start', None), (str, '--gene-end', 'col_gene_end', None),
            (str, '--loc-chr', 'col_loc_chr', None), (str, '--loc-strand', 'col_loc_strand', None), 
            (str, '--loc-start', 'col_loc_start', None), (str, '--loc-end', 'col_loc_end', None), 
            (int, '--adjacent', 'adjacent', 100)] 
    # tuples for options (in repopts) that should have the same length, e.g., ('gene_chr', 'gene_strand', 'gene_file')
    matchvars = [] #('col_gene_chr', 'col_loc_chr'), ('col_gene_strand', 'col_loc_strand'), ('col_gene_start', 'col_loc_start')] 
    # variables with fixed location, will be passed to Main using *fixopts, the rest will be passed to Main using **dictopts
    fixvars = [] # positional variables

    optlist, optdict = parseOpts(optlst, repopts=repopts, mandopts=mandopts, autopts=autopts, matchopts=matchvars, fixopts=fixvars, err=err)

    # specific options can be added over here
    try:
        fsrc = optdic.get('-i', optdic.get('--input', None)) or args.pop(0)
    except: fsrc = sys.stdin
    try:
        fobj = optdic.get('-o', optdic.get('--output', None)) or args.pop(0)
    except: fobj = sys.stdout

    if args:
        err.append('Error - unrecognized parameters:\n\t%s' % ', '.join(args))
    # quit on error
    if err:
        err.append('\n\nPlease type "%s -h " for help.\n' % sys.argv[0])
        exitMsg(err, out=sys.stderr)

    # start job
    Main(fsrc, fobj, *optlist, **optdict)
    
