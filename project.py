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

def test(ref, qry, cutoff, start, end):
	max_len = max([len(q[1]) for q in qry])
	arr1 = numpy.array([0] * (max_len + 1))
	arr2 = numpy.array([0] * (max_len + 1))
	ret = 0
	for i in range(start, end):
		ret = dynamic_opt(qry[i][1], ref[0][1], cutoff, arr1, arr2)
		if (ret <= cutoff):
			print qry[i][0]
	
def run(reference, query, cutoff, out_file):
	ref = parseFASTA(reference)
	qry = parseFASTA(query)
	out = open(out_file, 'w')
	
	max_len = max([len(q[1]) for q in qry])
	arr1 = numpy.array([0] * (max_len + 1))
	arr2 = numpy.array([0] * (max_len + 1))
	
	for pat in ref:
		line = pat[0]
		
		for str in qry:
			edit_dist = dynamic_opt(pat[1], str[1], cutoff, arr1, arr2)
			if edit_dist <= cutoff:
				line = line + ' ' + str[0]
		
		line = line + '\n'
		out.write(line)	
		
	out.close
			