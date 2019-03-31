def window(seq, n=2):
	""" Returns a sliding window (of width n) over data from the iterable
	    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... """
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield ''.join(result)
	for elem in it:
		result = result[1:] + (elem,)
		yield ''.join(result)

#eg.
# seq = 'TGAGGGATCACTCCTGTTTGTGAGGCACAGAGAAGGCACCCAGACCTGTACTGGCCA'
# kmers = window(seq)
# kmers:
#	TGA, GAG, AGG, .....
