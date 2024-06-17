#!/usr/bin/env python

'''
find_genes_by_loc.py

Purpose:
    Find gene names from a gene table with location information ([chr,] [strand,] start, end) by its location ([chr,] [strand,] start [, end]). Multiple matches will be output in multiple rows. To merge these rows, refer to command "mergerows".

License: GNU General Public License v2.0 (http://www.gnu.org/licenses/gpl-2.0.html)
Author: Dr. Xiao-Qin Xia
Email: xqxia70@gmail.com

Usage:
    find_genes_by_loc.py  [options]  parameters

Options:
    -h, --help:     Display this message, then exit
    --version:      Print version information
    -i, --input=:   Input tabular file (ITF). If not provided, pop the first value in
                    {parameters}, or STDIN
    -o, --output=:  Output tabular file. If not provided, pop the first value in
                    {parameters}, or STDOUT
    --gene-file=:   Mandatory. Gene tabular file (GTF)
    --gene-name=:   Optional. The column for gene name in the GTF. If provided, this
                    column will be added to output, otherwise all columns in
                    GTF will be outputed.
    --gene-chr=:    Optional. The column for chromosome name in the GTF
    --gene-strand=: Optional. The column for strand in the GTF
    --gene-start=:  Mandatory. The column for start location in the GTF
    --gene-end=:    Mandatory. The column for end location in the GTF
    --loc-chr=:     Optional. The column for chromosome name in the ITF
    --loc-strand=:  Optional. The column for strand in the ITF
    --loc-start=:   Mandatory. The column for start location in the ITF
    --loc-end=:     Optional. The column for end location in the GTF

    --no-title-GTF: Optional. No title row in the GTF. No title in output too
    --no-title-ITF: Optional. No title row in the ITF. No title in output too
    --by-index-GTF: Optional. Columns are specified by their indice in GTF
    --by-index-ITF: Optional. Columns are specified by their indice in ITF

parameters:
    Input file, Output file in order, in case that these augments are not provided as options

'''

import os, sys, re
from xqlibs.base.basetools import parseOpts, exitMsg, NoOptionErr, NoColumnErr, TableFile, OutputFile, TableCols, setdiff

import banyan

__version__ = '0.0.1'

def read_gene_loc_old(gene_file, cols, key_idx=None):
    '''
    cols:    [chr,] [strand,] start, end, name
    gene_file: A tabular file that should contain columns in cols.
    key_idx: The index of chr (and/or strand) in cols. If key_idx is provided,
            this function will return a dict with banyan.SortedDict as values
            and with (chr, strand) as keys, otherwise just a banyan.SortedDict
            will be returned.
    '''
    if key_idx: 
        gene_dic = chr_dic = {}
        for line in TableCols(gene_file, cols=cols):
            key = tuple(map(line.__getitem__, key_idx))
            loc_dic = key in chr_dic and chr_dic[key] or chr_dic.setdefault(key, banyan.SortedDict(updator=banyan.OverlappingIntervalsUpdator))
            for v in key_idx: # remove chr, strand from line
                line.pop(0)
            loc_dic[(int(line[0]), int(line[1]))] = line[2]
    else:
        gene_dic = loc_dic = banyan.SortedDict(updator=banyan.OverlappingIntervalsUpdator)
        for line in TableCols(gene_file, cols=cols):
            loc_dic[(int(line[0]), int(line[1]))] = line[2]
    return(gene_dic)

def read_gene_loc(gene_file, cols, cols_info, key_idx=None, no_title=False, by_index=False, sep='\t'):
    '''
    cols:    [chr,] [strand,] start [, end]
    cols_info: name or other columns (a list) will be used by the caller function.
    gene_file: A tabular file that should contain columns in cols.
    key_idx: The index of chr (and/or strand) in cols. If key_idx is provided,
            this function will return a dict with banyan.SortedDict as values
            and with (chr, strand) as keys, otherwise just a banyan.SortedDict
            will be returned.
    '''
    def set_val(loc_dic, line, strand=None, has_end=True, sort_loc=True, sep=sep):
        '''
        loc_dic: a dict or a banyan.SortedDcit - { (start[, end]) : name, ... }
        line:    start, [end,], name
        '''
        try:
            i = int(line[0])
        except:
            print(line)
            print(loc_dic)
        i = int(line[0])
        if has_end:  # start, end
            j = int(line[1])
            i = sort_loc and i > j and (j, i) or (i, j)  # make sure start < end
        loc_dic[i] = sep.join(line[1+(has_end and 1 or 0):])
    has_end = len(cols) - (key_idx and len(key_idx) or 0) > 1
    if cols_info:
        fsrc = TableCols(gene_file, cols=cols+cols_info, no_title=no_title, by_index=by_index)
    else:  # read all columns
        fsrc = TableCols(gene_file)
        if no_title: # must by_index, cols_info will not be filled since no_title
            i_cols = map(lambda a:int(a)-1, cols)
            def next_line(fsrc, idx=i_cols):
                i_ord = idx[:]
                i_ord.sort()
                for line in fsrc:
                    i_row = list(range(len(line)))
                    map(i_row.remove, idx)
                    yield map(line.__getitem__, idx+i_row)
        else:
            head = fsrc.next()
            i_locs = by_index and map(lambda a:int(a)-1, cols) or map(head.index, cols)
            idx = list(range(len(head)))
            # move i_locs to the front in idx
            map(idx.remove, i_locs)
            idx = i_locs + idx
            # find the title names for output
            idx_title = idx[1+(has_end and 1 or 0)+(key_idx and len(key_idx) or 0):]
            cols_info.extend(map(head.__getitem__, idx_title))
            def next_line(fsrc, idx=idx):
                for line in fsrc:
                    yield map(line.__getitem__, idx)
        fsrc = next_line(fsrc)
    if key_idx: 
        gene_dic = chr_dic = {}
        for line in fsrc:  # TableCols(gene_file, cols=cols):
            key = tuple(map(line.__getitem__, key_idx))
            loc_dic = key in chr_dic and chr_dic[key] or chr_dic.setdefault(key, banyan.SortedDict(updator=banyan.OverlappingIntervalsUpdator))
            for v in key_idx:  # remove chr, strand from line
                line.pop(0)
            set_val(loc_dic, line, has_end=has_end)
    else:
        gene_dic = loc_dic = banyan.SortedDict(updator=banyan.OverlappingIntervalsUpdator)
        for line in fsrc:  # TableCols(gene_file, cols=cols):
            set_val(loc_dic, line, has_end=has_end)
    return(gene_dic)

def Main(fsrc, fobj, gene_file, col_gene_name=None, col_gene_chr=None, col_gene_strand=None, col_gene_start=None, col_gene_end=None, col_loc_chr=None, col_loc_strand=None, col_loc_start=None, col_loc_end=None, no_title_gtf=False, no_title_itf=False, by_index_gtf=False, by_index_itf=False, sep='\t', joint_str=' / ', linesep=os.linesep):
    # determine what columns to read
    if not col_gene_start or not col_gene_end or not col_loc_start:
        raise NoOptionErr('Not all mandatory augments have been provided: gene_name, gene_start, gene_end, loc_start')
    cols_gene, cols_input, key_idx = [], [], []
    i = 0
    has_chr = col_gene_chr and col_loc_chr
    if has_chr:
        cols_gene.append(col_gene_chr)
        cols_input.append(col_loc_chr)
        key_idx.append(i)
        i += 1
    has_strand = col_gene_strand and col_loc_strand
    if has_strand:
        cols_gene.append(col_gene_strand)
        cols_input.append(col_loc_strand)
        key_idx.append(i)
        i += 1
    #cols_gene.extend([col_gene_start, col_gene_end, col_gene_name]) # chr, strand, start, end, gene
    cols_gene.extend([col_gene_start, col_gene_end]) # chr, strand, start, end
    cols_input.append(col_loc_start)
    if col_loc_end: 
        cols_input.append(col_loc_end)
    cols_info = []
    if col_gene_name:
        cols_info.append(col_gene_name)
    # read gene info, and cols_info!
    gene_dic = read_gene_loc(gene_file, cols_gene, cols_info, key_idx, no_title=no_title_gtf, by_index=by_index_gtf, sep=sep)
    # convert colname to idx
    fsrc = TableFile(fsrc)
    cols_src = not no_title_itf and fsrc.next() or []
    try:
        cols_idx = (no_title_itf or by_index_itf) and map(lambda a:int(a)-1, cols_input) or map(cols_src.index, cols_input)
    except:
        dif = setdiff(cols_input, cols_src)
        raise NoColumnErr('No necessary column%s "%s" in input.' % (len(dif)-1 and 's' or '', '", "'.join(cols_input)))
    if key_idx: # convert to index in input file 
        key_idx = map(cols_idx.__getitem__, key_idx)
        for v in key_idx:
            cols_idx.pop(0)

    # prepare output file
    fobj = OutputFile(fobj)
    if not no_title_gtf and not no_title_itf:
        #cols_src.append(col_gene_name)
        cols_src.extend(cols_info)
        fobj.write(sep.join(cols_src)+linesep)
    # deal with input
    if not key_idx:
        loc_dic = gene_dic
    for line in fsrc:
        if key_idx:
            key = tuple(map(line.__getitem__, key_idx))
            loc_dic = key in gene_dic and gene_dic[key] or banyan.SortedDict(updator=banyan.OverlappingIntervalsUpdator)
        loc = tuple(map(lambda i:int(line[i]), cols_idx))
        if col_loc_end:
            match_locs = loc_dic.overlap(loc)
        else:
            match_locs = loc_dic.overlap_point(loc[0])
        # genes = joint_str.join(map(loc_dic.__getitem__, match_locs))
        # line.append(genes)
        # fobj.write(sep.join(line)+linesep)
        line = sep.join(line)
        genes = map(loc_dic.__getitem__, match_locs)
        if not genes:
            fobj.write(line+linesep)
        else:
            for gene in genes:  # output in multiple rows.
                fobj.write(line+sep+gene+linesep)


if __name__ == '__main__':
    from getopt import getopt
    optlst, args = getopt(sys.argv[1:], 'hi:o:', ['help', 'version', 'input=', 'output=', 
        'gene-file=', 'gene-name=', 'gene-chr=', 'gene-strand=', 'gene-start=', 'gene-end=', 
        'loc-chr=', 'loc-strand=', 'loc-start=', 'loc-end=', 'no-title-GTF', 'no-title-ITF',
        'by-index-GTF', 'by-index-ITF'])

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
    mandopts = ['--gene-file', '--gene-start', '--gene-end', '--loc-start'] 
    # options that will be automatically dealt with. e.g. [(int, ('-n', '--num'), 'num', 3), ...] -- (type, (option names), target variable names, default values). type can be int, float, str, bool, or any self-defined functions, e.g. lambda s:set(re.split(r'[,;|]', s)).
    autopts = [(str, '--gene-file', 'gene_file', None), (str, '--gene-name', 'col_gene_name', None), 
            (str, '--gene-chr', 'col_gene_chr', None), (str, '--gene-strand', 'col_gene_strand', None), 
            (str, '--gene-start', 'col_gene_start', None), (str, '--gene-end', 'col_gene_end', None),
            (str, '--loc-chr', 'col_loc_chr', None), (str, '--loc-strand', 'col_loc_strand', None), 
            (str, '--loc-start', 'col_loc_start', None), (str, '--loc-end', 'col_loc_end', None),
            (bool, '--no-title-GTF', 'no_title_gtf', False), (bool, '--no-title-ITF', 'no_title_itf', False),
            (bool, '--by-index-GTF', 'by_index_gtf', False), (bool, '--by-index-GTF', 'by_index_itf', False)] 
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
    
