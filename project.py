__author__ = 'Rob'

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