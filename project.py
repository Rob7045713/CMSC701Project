__author__ = 'Rob Argue'
__author__ = 'Tommy Pensyl'

import numpy

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
	
def dynamic_naive(pat, str):
	table = numpy.matrix([[0] * (len(pat) + 1)] * (len(str) + 1))
	
	for row in range(len(pat) + 1):
		table[0, row] = 0
		
	for col in range(len(str) + 1):
		table[col, 0] = col
		
	for col in range(1, len(str) + 1):
		for row in range(1, len(pat) + 1):
			a = table[col - 1, row - 1]
			if pat[row - 1] != str[col - 1]:
				a += 1
			b = table[col, row - 1] + 1
			c = table[col - 1, row] + 1
			table[col, row] = min(a, min(b,c))
				
	return table
	
def dynamic_better(pat, str, cutoff):
	top_row = numpy.array([0] * (len(str) + 1))
	bot_row = numpy.array([0] * (len(str) + 1))
	min_in_row = cutoff + 1
	min_idx = 1
	max_idx = cutoff + 1
	
	for col in range(max_idx + 1):
		top_row[col] = col	
	
	###	
	#print top_row	
	#debug_marker = numpy.array([0] * (len(str) + 1))
	#debug_marker[min_idx] = 1
	#debug_marker[max_idx] = 2
	#print debug_marker
	###
	
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
			val = min(a, min(b,c))
			bot_row[col] = val
			
			# min / max update
			if  val > cutoff:
				min_idx = min_idx + 1
				
			if col == max_idx and val <= cutoff:
				max_idx = max_idx + 1
				bot_row[max_idx] = cutoff + 1
			
			# update min in row
			min_in_row = min(min_in_row, val)
			
		# early return
		if min_in_row > cutoff:
			return min_in_row
								
		###
		#print bot_row
		#debug_marker = numpy.array([0] * (len(str) + 1))
		#debug_marker[min_idx] = 1
		#debug_marker[max_idx] = 2
		#print debug_marker
		###

		# swap rows
		temp = top_row
		top_row = bot_row
		bot_row = temp

	return min_in_row
	
def dynamic_opt(pat, str, cutoff):
	top_row = numpy.array([0] * (len(str) + 1))
	bot_row = numpy.array([0] * (len(str) + 1))
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
				bot_row[max_idx] = cutoff + 1
			
			# update min in row
			if val < min_in_row:
				min_in_row = val

		# early return
		if min_in_row > cutoff:
			return min_in_row
					
		# swap rows
		temp = top_row
		top_row = bot_row
		bot_row = temp

	return min_in_row

def test(ref, qry, cutoff, its):
	for i in range(its):
		print dynamic_opt(qry[i][1], ref[0][1], cutoff)
	
def run(reference, query, cutoff, out_file):
	ref = parseFasta(reference)
	qry = parseFasta(query)
	out = open(out_file, 'w')
	
	for pat in ref:
		out.write(pat[0])
		for str in qry:
			table = dynamic_naive(pat[1], str[1])
			if table[-1,-1] <= cutoff:
				out.write(' ' + str[0])
			
			