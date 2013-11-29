__author__ = 'Rob Argue'
__author__ = 'Tommy Pensyl'

import profile
import time
import numpy

alphabet = {'A':0,'C':1,'G':2,'T':3}

def parseFASTA(fasta):
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
    
def get_q_freqs(qry, search_range):
    q_freq = []
    for q in qry:
        string = q[1]
        freq_arr = [0]*4
        for i in range(min(search_range,len(string))):
            freq_arr[alphabet[string[i]]] += 1
        q_freq.append(freq_arr)
    return q_freq

def dynamic_opt(pat, str, cutoff, top_row = None, bot_row = None):
    row_len = len(str) + 1
    
    if (top_row == None):
        top_row = numpy.array([0] * row_len)
    if (bot_row == None):
        bot_row = numpy.array([0] * row_len)
    
    min_in_row = cutoff + 1
    min_idx = 1
    max_idx = cutoff + 1
    
    for col in range(max_idx + 1):
        top_row[col] = col  
    
    for row in range(1, len(pat) + 1):
        min_in_row = cutoff + 1
        bot_row[0] = top_row[0] + 1
                
        for col in range(min_idx, max_idx + 1):
            
            # compute 3 possible vals
            a = top_row[col - 1]
            if pat[row - 1] != str[col - 1]:
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
            if  val > cutoff:
                min_idx = min_idx + 1
                
            if col == max_idx and val <= cutoff:
                max_idx = max_idx + 1
                if max_idx > row_len - 1:
                    max_idx = row_len - 1
                else:
                    bot_row[max_idx] = cutoff + 1
            
            # update min in row
            if val < min_in_row:
                min_in_row = val

        ###
        #print bot_row
        #debug_marker = numpy.array([0] * (len(str) + 1))
        #debug_marker[min_idx] = 1
        #debug_marker[max_idx] = 2
        #print debug_marker
        #print max_idx
        #print row_len
        ###
                
        # early return
        if min_in_row > cutoff:
            return min_in_row
                    
        # swap rows
        top_row, bot_row = bot_row, top_row

    return min_in_row


# gets letter counts for a string
def get_freqs(string, length = None, freqs = None):
    if length == None:
        length = len(string)

    if freqs == None:
        freqs = numpy.array([0] * len(alphabet))
    else:
        for i in range(len(alphabet)):
            freqs[i] = 0

    for i in range(length):
        freqs[alphabet[string[i]]] += 1

    return freqs

# optimimzed culling function
def cull(pat_freq, str_freq, cutoff):
    sum = 0

    for let in range(len(alphabet)):
        pat_let = pat_freq[let]
        str_let = str_freq[let]
        
        if pat_let > str_let:
            sum += pat_let - str_let
        else:
            sum += str_let - pat_let

    return sum > 2 * cutoff

def run(reference, query, cutoff = 5, search_range = 100, out_file = None, verbose = False):
    
    if out_file == None:
        out_file = time.strftime('%y_%m_%d_%H_%M_%S') + '_cutoff_' + str(cutoff) + '_range_' + str(search_range) + '.txt'

    out = open(out_file, 'w')

    if verbose:         
        print 'Run info:'
        print '  Time : ' + time.asctime()
        if isinstance(reference, str):
            print '  Reference file : ' + reference
        if isinstance(query, str):
            print '  Query file : ' + query
        print '  Output file : ' + out_file
        print '  Cutoff : ' + str(cutoff)
        print '  Search range : ' + str(search_range)
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
    ref = ref[0:5]#                                     # Hack here
    ################
    
    if isinstance(query, str):
        if verbose:
            print '  Loading ' + query
    
        qry = parseFASTA(query)

    else:
        qry = query

    #####################
    qry = qry[0:100000]#                                # and here
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
    #search_range = min([len(x[1]) for x in qry])       # TODO : this

    if verbose:
        print '  - Query frequencies'

    q_freq = get_q_freqs(qry, search_range)

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
    pat_freq = numpy.array([0] * len(alphabet))


    for pat in ref:
    
        time_single_match_start = time.time()

        line = pat[0]

        pat_freq = get_freqs(pat[1], search_range, pat_freq)

        num_culled = 0    
        num_matched = 0
        for i in range(len(qry)):
            q = qry[i]
            if not cull(pat_freq, q_freq[i], cutoff):
                edit_dist = dynamic_opt(pat[1], q[1], cutoff, arr1, arr2)
                if edit_dist <= cutoff:
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
    
