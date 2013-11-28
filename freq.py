__author__ = 'Rob Argue'
__author__ = 'Tommy Pensyl'

from project import *
import profile
import time
import numpy

#global? I just kinda guessed
map1 = {'A':0,'C':1,'G':2,'T':3}
    
def get_q_freqs(qry, search_range):
    q_freq = []
    for q in qry:
        string = q[1]
        freq_arr = [0]*4
        for i in range(min(search_range,len(string))):
            freq_arr[map1[string[i]]] += 1
        q_freq.append(freq_arr)
    return q_freq

# gets letter count for all prefixes
def get_cumulative_freqs(string):
    freqs = numpy.matrix([[0]*4]*len(string))
    freq_arr = [0]*4
    for i in range(0,len(string)):
        freq_arr[map1[string[i]]] += 1
        freqs[i] = freq_arr
    return freqs

def run(reference, query, cutoff, out_file, search_range = 100, verbose = False):
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
    
    ref = ref[0:5]                                      # Hack here

    
    if isinstance(query, str):
        if verbose:
            print '  Loading ' + query
    
        qry = parseFASTA(query)

    else:
        qry = query

    qry = qry[0:100000]                                 # and here
   
    if verbose:
        print '  Load time : %.3f s' % (time.time() - time_loading_start)
        print ''

    ### /\ Loading data /\ ###



    num_queries = len(qry)

    out = open(out_file, 'w')

    time_start = time.time()

    max_len = max([len(q[1]) for q in qry])
    arr1 = numpy.array([0] * (max_len + 1))
    arr2 = numpy.array([0] * (max_len + 1))

    #search_range = min([len(x[1]) for x in qry])       # TODO : this
    


   ### \/ Query frequencies \/ ###

    time_freq_start = time.time()

    if verbose:
        print 'Calculating query letter frequencies'

    q_freq = get_q_freqs(qry, search_range)

    if verbose:
        print '  Calc time : %.3f s' % (time.time() - time_freq_start)
        print ''

    ### /\ Query frequencies /\ ###



    ### \/ Matching \/ ###

    time_match_start = time.time()

    if verbose:
        print 'Matching'
        print ''

    count = 0

    for pat in ref:
    
        time_single_match_start = time.time()

        line = pat[0]

        pat_freq = get_cumulative_freqs(pat[1])

        num_culled = 0    
        matched = 0
        for i in range(len(qry)):
            q = qry[i]
            if sum([abs(pat_freq[search_range-1,j]-q_freq[i][j]) for j in [0,1,2,3]]) <= 2*cutoff:
                edit_dist = dynamic_opt(pat[1], q[1], cutoff, arr1, arr2)
                if edit_dist <= cutoff:
                    matched += 1
                    line = line + ' ' + q[0]
                    line = line + '\n'
            else:
                num_culled+=1
        out.write(line)

        if verbose:
            count += 1
            percent_culled = 100.0 * float(num_culled) / float(num_queries)

            print 'Reference ' + str(count) + ' (' + pat[0] + '):' 
            print '  Time : %.3f s' % (time.time() - time_single_match_start)
            print '  Queries culled : ' + str(num_culled) + ' / ' + str(num_queries) + ' (%.2f%%)' % percent_culled
            print '  Queries matched : ' + str(matched)
    
    if verbose:
        print ''
        print 'Matching time : %.3f s' % (time.time() - time_match_start)
        print ''

    ### /\ Matching /\ ###

    out.close

    if verbose:
        print 'Total run time : %.3f s' % (time.time() - time_start)
    

# make sure to manually edit file so that it ends on a proper line.
def get_small_query_file(file, n):
    in_file = open('query.fna','r')
    out = open('query2.fna','w')
    i = 0
    for line in in_file:
        i+=1
        if (i >= n):
            break
        out.write(line)
    in_file.close()
    out.close()

