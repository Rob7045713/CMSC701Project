# match.py
# Copyright (c) 2013 Rob Argue (rargue@cs.umd.edu) and Tommy Pensyl (tpensyl@cs.umd.edu)

__author__ = 'Rob Argue'
__author__ = 'Tommy Pensyl'

import time
import numpy
import cProfile
import pickle
import sys

# alphabet for DNA
dna_alphabet = {'A':0,'C':1,'G':2,'T':3}

def parseFASTA(fasta):
    """ Parse a FASTA file into a list.

    Returns:
        List of tuples of the form (pattern_name, pattern)

    Arguments:
        fasta - Filename as a string
    """

    f = open(fasta, 'r')

    strList = []
    name = ''
    str = ''

    for line in f:
        if line.startswith('>'):
            if name != '':
                strList.append((name, str))

            name = line.split()[0].strip('>')
            str = ''

        else:
            str += line.strip()

    if name != '':
        strList.append((name, str))
        
    return strList



def get_freqs(string, alphabet, max_length = None, freqs = None):
    """ Calculate the letter frequencies of a string

    Returns:
        List of frequencies of each letter, based on the alphabet

    Arguments:
        string      - String to get the frequencies for

        alphabet    - Dictionary mapping letters in the alphabet to indicies

        max_length  - Maximum length of the string to consider
                      Default = None (uses full string)

        freqs       - Array to use instead of creating a new one. Optional
                      parameter used for optimization.
                      Default = None (creates new array)
    """

    if max_length == None:
        max_length = len(string)

    if freqs == None:
        freqs = numpy.array([0] * len(alphabet))
    else:
        for i in range(len(alphabet)):
            freqs[i] = 0

    for i in range(max_length):
        freqs[alphabet[string[i]]] += 1

    return freqs


    
def get_data_freqs(data, alphabet, max_length = sys.maxint):
    """ Caclulate the letter frequencies of all strings in the data

    Returns:
        List of arrays of letter frequencies, based on the alphabet

    Arguments:
        data        - Data set in the form of a list of (string_name, string)

        alphabet    - Dictionary mapping letters in the alphabet to indicies
        
        max_length  - Maximum length of the string to consider
                      Default = sys.maxint
    """

    freqs = []
    for tup in data:
        string = tup[1]
        freq = get_freqs(string, alphabet, min(max_length, string))
        freqs.append(freq)
    return freqs



def cull(pat_freq, str_freq, max_dist):
    """ Determine whether to cull based on letter freqencies

    Returns:
        Boolean which is true if the pair of strings can be culled

    Arguments:
        pat_freq    - Array of letter frequencies in the pattern

        str_freq    - Array of letter frequencies in the string

        cutoff      - Maximum edit distance allowed between the strings
    """

    sum = 0

    for let in range(len(pat_freq)):
        pat_let = pat_freq[let]
        str_let = str_freq[let]
        
        if pat_let > str_let:
            sum += pat_let - str_let
        else:
            sum += str_let - pat_let

    return sum > 2 * max_dist



def edit_distance(pattern, string, max_dist, top_row = None, bot_row = None):
    """ Computes the edit distance of a pair of strings (ignoring insertions
        at the end of the pattern) within a specified maximum distance. Uses
        dynamic programming optimized by not considering anything outside of
        the maximum edit distance. Only stores two rows of the dynamic
        programming table at any given time.

    Returns:
        Minimum of the edit distance between the strings and the maximum 
        distance plus 1

    Arguments:
        pattern     - Pattern string

        string      - String being matched to

        max_dist    - Maximum edit distance to be considered

        top_row     - First array to use for storage. Optional parameter used
                      for optimization.
                      Default = None (creates a new array)

        bot_row     - Second array to use for storage. Optional parameter used
                      for optimization.
                      Default = None (creates a new array)
    """

    row_len = len(string) + 1
    
    if (top_row == None):
        top_row = numpy.array([0] * row_len)
    if (bot_row == None):
        bot_row = numpy.array([0] * row_len)
    
    min_in_row = max_dist + 1
    min_idx = 1
    max_idx = max_dist + 1
    
    for col in range(max_idx + 1):
        top_row[col] = col  
    
    for row in range(1, len(pattern) + 1):
        min_in_row = max_dist + 1
        bot_row[0] = top_row[0] + 1
                
        for col in range(min_idx, max_idx + 1):
            
            # compute 3 possible vals
            a = top_row[col - 1]
            if pattern[row - 1] != string[col - 1]:
                a += 1
            b = top_row[col] + 1
            c = bot_row[col - 1] + 1
            
            # find min of 3
            if a < b and a < c:
                val = a
            elif b < c:
                val = b
            else:
                val = c
            
            bot_row[col] = val
            
            # min / max update
            if  val > max_dist:
                min_idx = min_idx + 1
                
            if col == max_idx and val <= max_dist:
                max_idx = max_idx + 1
                if max_idx > row_len - 1:
                    max_idx = row_len - 1
                else:
                    bot_row[max_idx] = max_dist + 1
            
            # update min in row
            if val < min_in_row:
                min_in_row = val
                
        # early return
        if min_in_row > max_dist:
            return max_dist + 1
                    
        # swap rows
        top_row, bot_row = bot_row, top_row

    return min_in_row



def run(reference, query, max_dist = 5, freq_len = 100, out_file = None, verbose = False, profile = False):
    
    if out_file == None:
        out_file = time.strftime('%y_%m_%d_%H_%M_%S') 
        out_file += '_max_dist_' + str(max_dist)
        out_file += '_freq_len_' + str(freq_len)

    out = open(out_file + '.out', 'w')

    if profile:
        prof = cProfile.Profile()
        prof.enable()

    if verbose:         
        print 'Run info:'
        print '  Time : ' + time.asctime()
        if isinstance(reference, str):
            print '  Reference file : ' + reference
        if isinstance(query, str):
            print '  Query file : ' + query
        print '  Output file : ' + out_file
        print '  Max edit distance : ' + str(max_dist)
        print '  Frequency length cutoff : ' + str(freq_len)
        print ''


    time_start = time.time()



    ### \/ Loading data \/ ###

    time_loading_start = time.time()

    if verbose:
        print 'Loading data'

    if isinstance(reference, str):
        if verbose:
            print '  Loading ' + reference

        ref = parseFASTA(reference)
    
    else:
        ref = reference
    
    ################
    #ref = ref[0:5]#                                     # Hack here
    ################
    
    if isinstance(query, str):
        if verbose:
            print '  Loading ' + query
    
        qry = parseFASTA(query)

    else:
        qry = query

    #####################
    #qry = qry[0:100000]#                                # and here
    #####################

    if verbose:
        print '  Load time : %.3f s' % (time.time() - time_loading_start)
        print ''

    ### /\ Loading data /\ ###



    ### \/ Query preprocessing \/ ###

    time_prep_start = time.time()
    if verbose:
        print 'Preprocessing queries'
        print '  - Num queries'

    num_queries = len(qry)

    if verbose:
        print '  - Max length'

    max_len = max([len(q[1]) for q in qry])
    #freq_len = min([len(x[1]) for x in qry])       # TODO : this

    if verbose:
        print '  - Query frequencies'

    q_freq = get_data_freqs(qry, dna_alphabet, freq_len)

    if verbose:
        print '  Time : %.3f s' % (time.time() - time_prep_start)
        print ''

    ### /\ Query preprocessing /\ ###



    ### \/ Matching \/ ###

    time_match_start = time.time()

    if verbose:
        print 'Matching'
        print ''

    count = 0
    arr1 = numpy.array([0] * (max_len + 1))
    arr2 = numpy.array([0] * (max_len + 1))
    pat_freq = numpy.array([0] * len(dna_alphabet))


    for pat in ref:
    
        time_single_match_start = time.time()

        line = pat[0]

        pat_freq = get_freqs(pat[1], dna_alphabet, freq_len, pat_freq)

        num_culled = 0    
        num_matched = 0
        for i in range(len(qry)):
            q = qry[i]
            if not cull(pat_freq, q_freq[i], max_dist):
                edit_dist = edit_distance(q[1], pat[1], max_dist, arr1, arr2)
                if edit_dist <= max_dist:
                    num_matched += 1
                    line = line + ' ' + q[0]
            else:
                num_culled+=1
            
        line = line + '\n'
        out.write(line)

        if verbose:
            count += 1
            percent_culled = 100.0 * float(num_culled) / float(num_queries)

            print 'Reference ' + str(count) + ' (' + pat[0] + '):' 
            print '  Time : %.3f s' % (time.time() - time_single_match_start)
            print '  Queries culled : ' + str(num_culled) + ' / ' + str(num_queries) + ' (%.2f%%)' % percent_culled
            print '  Queries matched : ' + str(num_matched)
    
    if verbose:
        print ''
        print 'Matching time : %.3f s' % (time.time() - time_match_start)
        print ''

    ### /\ Matching /\ ###



    out.close()

    if verbose:
        print 'Total run time : %.3f s' % (time.time() - time_start)
    
    if profile:
        prof.create_stats()
        prof.print_stats()

        ## apparently you can't do this?
        #out = open(out_file + '.prof', 'w')
        #pickle.dump(prof, out)
        #out.close()
