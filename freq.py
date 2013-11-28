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

def run(reference, query, cutoff, out_file):
    ref = parseFASTA(reference)
    ref = ref[0:5]                                      # Hack here
    qry = parseFASTA('query2.fna')
    qry = qry[0:20000]                                  # and here
    out = open(out_file, 'w')

    time0 = time.time()
	
    max_len = max([len(q[1]) for q in qry])
    arr1 = numpy.array([0] * (max_len + 1))
    arr2 = numpy.array([0] * (max_len + 1))

    #search_range = min([len(x[1]) for x in qry])
    search_range = 100                                  # and here too
    
    q_freq = get_q_freqs(qry, search_range)

    for pat in ref:
        line = pat[0]

        pat_freq = get_cumulative_freqs(pat[1])

        num_culled = 0    
        for i in range(len(qry)):
            q = qry[i]
            if sum([abs(pat_freq[search_range-1,j]-q_freq[i][j]) for j in [0,1,2,3]]) <= 2*cutoff:
                edit_dist = dynamic_opt(pat[1], q[1], cutoff, arr1, arr2)
                if edit_dist <= cutoff:
                    line = line + ' ' + q[0]
                    line = line + '\n'
            else:
                num_culled+=1
        out.write(line)
        print 'culled '+str(num_culled)
		
    out.close

    print "Total Time: " + str(time.time()-time0)

profile.run("run('reference.fna','query2.fna',5,'out1.txt')")



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

