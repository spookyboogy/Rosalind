#!/usr/bin/env python3
## http://rosalind.info/problems/list-view/

import os
from math import factorial as fact
from urllib.request import urlopen


def nucleotide_count(string):

	"""
	string -> sequence of DNA (GACT) or RNA (GACU) nucleotides

	Returns the count of the number of appearances of each nucleotide.
	"""
	
	s = string = string.upper()

	dna = ['A', 'C', 'G', 'T']; rna = ['A', 'C', 'G', 'U']

	if list(set(string)) == list(set(dna)):
		a, c = s.count("A"), s.count("C")
		g, t = s.count("G"), s.count("T")
		return a, c, g, t
	elif list(set(string)) == list(set(rna)):
		a, c = s.count("A"), s.count("C") 
		g, u = s.count("G"), s.count("U")
		return a, c, g, u
	else:
		raise ValueError("You must enter a DNA or RNA sequence.")


def dna_to_rna(rna_string):

	"""
	Takes a sequence of RNA nucleotides and returns a transcribed
	sequence of DNA nucleotides.
	"""

	r = rna_string = rna_string.upper()
	d = ''

	for nt in r:
		if nt == "T":
			d += "U"
		else:
			d += nt
	return d


def reverse_compliment(dna_string):

	"""
	Returns the reverse complement of a string of DNA nucleotides.
	"""

	s = dna_string = dna_string.upper()[::-1]

	compliments = {"A":"T", "C":"G", "G":"C", "T":"A"}

	c = ''
	for i in s:
		c += compliments[i]
	return c


def fatsa_read(fatsa_file, names = True):

	"""
	Takes a file in FATSA format and returns list of
	(ID, genetic_string) tuples.
	"""
	
	f = fatsa_file
	if not os.path.isfile(f):
		raise ValueError("Invalid file")
		return
	else:
		fopen = open(f, 'r')
		l = fopen.read().split('>')
		l = [i for i in l if i != '']

	for i in range(len(l)):
		entry = l[i]
		name = entry[:entry.index('\n')]
		dna = entry[entry.index('\n') + 1:]
		l[i] = tuple([name, dna.replace('\n', '')])
	
	return l	
	

def gc_content(fatsa_file):

	"""
	Takes a fatsa-formatted file containing 1 or more DNA strings
	and returns the ID and GC-content of the string with the
	greatest GC-content, ie the percentage of a given DNA string
	consiting of either G or C.
	"""
	
	l = fatsa_read(fatsa_file)
	
	max_gc = 0
	max_entry = None

	for tup in l:
			
		s = tup[1]
		gc = 100 * ((s.count("C") + s.count("G")) / len(s))
		if gc > max_gc:
			max_gc = gc
			max_entry = tup[0] 
	
	print("\n{}\n{}\n".format(max_entry, max_gc))


def hamm(s, t):

	"""
	Returns the hamming distance between two strings s and t;
	ie, the number of corresponding symbols that differ
	between s and t.
	"""

	d = 0

	if len(s) != len(t):
		raise ValueError("Strings should be the same length")
	else:
		for i in range(len(s)):
			if s[i] != t[i]:
				d += 1
	return d			
	

def rabbits(n, m = 0, k = 1):

	"""
	n -> number of months for breeding
	k -> number of rabbit-pair (M|F) offspring per generation

	Returns the number of rabbit pairs at the end of n months given
	that we start with 1 pair, k rabbit pairs are produced
	each generation and each rabbit pair may reproduce after 1 month
	of being born.
	"""	

	def real_rabbits(n, k = 1):

		"""
		Could be extended to account for death or to 
		simulate population growth when rabbits don't stick in
		pairs or whatever. Too memory intesive for large populations.
		"""

		class rabbit_pair:
			
			def __init__(self):
				
				self.age = 0
				self.mature = 0
				#self.gender = gender //for single-rabbit class
			
			def envejecer(self):
				
				self.age += 1
				if self.age > 0:
					self.mature = 1
			
			def reproduce(self, pairs = 3):

				if self.mature:
					return [rabbit_pair() for i in range(pairs)]
				else:
					return [None]


		population = [rabbit_pair()]
		
		for i in range(1, n+1):
			
			print("\nMonth = {}, pop = {}".format(i, len(population)))
			offspring = []
			
			for rabbit in population:
					
				litter = rabbit.reproduce(pairs = k)
				offspring += [i for i in litter if i != None] 
				rabbit.envejecer()	
				if rabbit.age == m:
					population.remove(rabbit)

			population += offspring
			
		return len(population)


	def rabbits_easy_mode(n, k):

		"""
		Treating the rabbits less like rabbits.
		"""

		population = ['r']

		for i in range(1, n):
			print("i = {}, pop = {}".format(i, len(population)))	
			offspring = []
		
			for i in range(len(population)):
				if population[i] == 'r':
					population[i] = 'R'
				else:
					offspring += k * ['r']
			
			population += offspring
		
		return len(population)


	def math_rabbits(n, m, k = 1):

		if not m:
			s = [0, 1, 1, 1 + k]
			i = 3
			while i < n:
				s += [s[i] + (k*s[i-1])]
				i += 1
			return s[-1]
		else:
			print("m = {}".format(m))
			s = [0, 1, 1, 2] #2, 3, 4] # m =3
			i = 3
			while i < n:
				if len(s) < m + 1:
					s += [s[i] + s[i-1]]
				elif len(s) == m + 1:
					s += [s[i] + s[i-1] - 1]
				else:
					s += [s[i] + s[i-1] - s[i-m]]
				i += 1
			return s[-1]
	
	return math_rabbits(n, m, k)


def dominance_prb(k, m, n):

	"""
	k -> Number of homozygous dominant individuals 
	m -> Number heterozygous individuals
	n -> Number homozygous recessive individuals

	Returns the probability that two randomly selected mates will 
	produce an individual possessing a dominant allele (and thus 
	displaying the dominant phenotype). Assume that any two organisms
	can mate.
	"""
	
	def dom_prob(g1, g2):

		if g1 in k or g2 in k:
			return 1
		elif g1 in m:
			if g2 in m:
				return .75
			else: 
				return .5
		else: 
			if g2 in m:
				return .5
			else:
				return 0

	k = ["AA" for i in range(k)]
	m = ["Aa" for i in range(m)]
	n = ["aa" for i in range(n)]	
	
	tot = len(k) + len(m) + len(n)
	
	prbmm = (len(m)/tot)*((len(m)-1)/(tot-1))
	prbmn = (len(m)/tot)*(len(n)/(tot-1))
	prbnn = (len(n)/tot)*((len(n)-1)/(tot-1))
	
	prb_nodom_mm = prbmm * (1-(dom_prob(m[0], m[0])))
	prb_nodom_mn = prbmn * (1-(dom_prob(m[0], n[0])))
	prb_nodom_nn = prbnn * (1-(dom_prob(n[0], n[0])))
	
	return 1 - (prb_nodom_mm + 2*(prb_nodom_mn) + prb_nodom_nn)


rna_codons = {

	'UUU':'F'   , 'CUU':'L', 'AUU':'I', 'GUU':'V',
	'UUC':'F'   , 'CUC':'L', 'AUC':'I', 'GUC':'V',
	'UUA':'L'   , 'CUA':'L', 'AUA':'I', 'GUA':'V',
	'UUG':'L'   , 'CUG':'L', 'AUG':'M', 'GUG':'V',
	'UCU':'S'   , 'CCU':'P', 'ACU':'T', 'GCU':'A',
	'UCC':'S'   , 'CCC':'P', 'ACC':'T', 'GCC':'A',
	'UCA':'S'   , 'CCA':'P', 'ACA':'T', 'GCA':'A',
	'UCG':'S'   , 'CCG':'P', 'ACG':'T', 'GCG':'A',
	'UAU':'Y'   , 'CAU':'H', 'AAU':'N', 'GAU':'D',
	'UAC':'Y'   , 'CAC':'H', 'AAC':'N', 'GAC':'D',
	'UAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E',
	'UAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E',
	'UGU':'C'   , 'CGU':'R', 'AGU':'S', 'GGU':'G',
	'UGC':'C'   , 'CGC':'R', 'AGC':'S', 'GGC':'G',
	'UGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G',
	'UGG':'W'   , 'CGG':'R', 'AGG':'R', 'GGG':'G'}


def rna_to_prot(rna_string):
	
	"""
	Takes a string of RNA nt's and returns the 
	protein for which it encodes.
	"""
	
	r = rna_string = rna_string.upper()
	p = ''

	for i in range(0, len(r), 3):
		s = r[i:i+3]
		if rna_codons[s] == 'Stop':
			p += ' '
			continue
		p += rna_codons[s]
	
	print(len(r))
	return p


def subs(string, substring):

	"""
	Returns the starting positions of substring in string.
	"""

	s = string = string.upper()
	ss = substring = substring.upper()

	indexes = []

	for i in range(len(s) - len(ss) + 1):
		if s[i : i + len(ss)] == ss:
			indexes += [i+1]

	for i in indexes: 
		print(i, end = ' ')
	print()
	return indexes 

def cons(fatsa_file):

	"""
	Essential Goal: Average a set of genetic strings
	Returns: A Consensus string and Profile Matrix P for a fatsa_file
		     containing n genetic strings of length m.

	Profile Matrix P:
		An IxJ matrix where:
		  I = # of distinct [Nucleotide] symbols used within
			  the set of given genetic strings (eg 4 for ACGT)
		  J = m, the length of the given strings
		  P_ij = the total count of occurences of the i_th distinct
		  	     symbol (eg A) amongst the j_th index of the given
		  	     strings. 
	Consensus String:
		A genetic string g where g[i] is the nucleotide which appeared
		most frequently at the i_th indexes of the given strings.
	"""

	string_entries = fatsa_read(fatsa_file)
	data_matrix = []
	for i in range(len(string_entries)):
		data_matrix += [string_entries[i][1]]	
	
	Sym_set = []
	for i in range(len(data_matrix)):
		Sym_set += list(set(data_matrix[i]))
	Sym_set = set(Sym_set)

	p = {i : [] for i in Sym_set}
	consensus = ''
	
	length = min((len(i) for i in data_matrix))
	if length != max((len(i) for i in data_matrix)):
		print("Strings aren't of same length as they should be.")

	for j in range(length):	
		col = [i[j] for i in data_matrix]
		for key in p.keys():
			p[key] += [col.count(key)]
	
	for j in range(length):
		cons_sym = ''
		max_freq = 0
		for i in p.items():
			if i[1][j] > max_freq:
				max_freq = i[1][j]
				cons_sym = i[0]
		consensus += cons_sym

	p_syms_ordered = sorted(p, key = lambda entry: entry[0])	
	
	with open("output_{}".format(fatsa_file), 'w') as fout: 
		print(consensus)
		fout.write("{}\n".format(consensus))
		for sym in p_syms_ordered:
			print("{}: ".format(sym), end = ' ')
			fout.write("{}: ".format(sym))
			for i in p[sym]:
				print(i, end = ' ')
				fout.write("{} ".format(i))
			print()
			fout.write("\n")
		

def overlap_graph(fatsa_file):

	"""
	Takes a fatsa_file containing DNA strings and prints the
	adjacency list corresponding to O_3, ie the overlap graph
	of the strings where overlaps are at least 3 nucleotides long. 
	(Also writes an output file for easy keeping)
	"""

	data = fatsa_read(fatsa_file)
	k = 3
	directions = {}

	for tail in data: 
		suffix = tail[1][-k:] 
		overlaps = []
		for head in data:		
			prefix = head[1][:k]
			if suffix == prefix: 
				if head[0] != tail[0]:
					overlaps += [head[0]]
		if len(overlaps) != 0:
			directions[tail[0]] = overlaps

	with open("output_" + fatsa_file, 'w') as fout:
		for tail in directions.keys():
			for head in directions[tail]:
				print("{} {}".format(tail, head))
				fout.write("{} {}\n".format(tail, head)) #WOWWWEWAHHHH! :D


def expected_offspring(AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa):

	"""
	Takes 6 positive integers representing the number of couples
	in a population possessing each genotype pairing for a given factor
	Returns the expected number of offspring displaying the dominant
	phenotype in the next generation, under the assumption that every
	couple has exactly two offspring.
	"""
	
	pairings = [AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa]
	dom_prob = [1, 1, 1, .75, .5, 0]
	EX = 0

	for i in range(len(pairings)):
		EX += 2*(pairings[i]*dom_prob[i])

	return EX

def motif(fatsa_file):

	"""
	Returns the longest common substring of the genetic strings
	contained in fatsa_file.
	"""

	data = fatsa_read(fatsa_file)
	data = [i[1] for i in data]
	longest = ''

	for i in range(len(data)):
		a = data[i]
		lngst = ''
		for c in range(len(a)-1):
			testlen = max(2, len(longest))
			while testlen <= len(a) - c:
				teststr = a[c:c+testlen]
				common = True
				for j in data:
					if teststr in j:
						continue
					else:
						common = False
				if common:
					if len(teststr) > len(lngst):	
						lngst = teststr
					testlen += 1
				else:
					break
		if len(lngst) > len(longest):
			longest = lngst
	with open("output_" + fatsa_file, 'w') as fout:
		fout.write(longest)
	return longest	


def choose(n, k):

	return fact(n)/(fact(k)*fact(n-k))

	
def ind_alleles(k, n):

	"""
	Given: A population starts at generation 0 with 
		   1 organism of genotype AaBb. This organism and 
		   all offspring mate with type AaBb, and have two 
		   offspring.
	Returns: The probability that n organisms of genotype
			 AaBb will belong to the k_th generation. 
	"""

	tot = 2**k
	prb = 0
	for i in range(n, tot+1):
		prb += ((1/4)**i)*((3/4)**(tot-i))*choose(tot,i)
	return prb


def motif_location(prot, motif = "N{P}[ST]{P}"):

	"""
	Takes a protein sequence and a protein sequence and
	returns a list of the motif's locations within the 
	protein, if any.
	
	Note: A motif (at least with respect to this program)
		  is a definite-length protein sequence. Protein
		  units enclosed in a '[]' indicate that anyone one of 
		  the enclosed units can appear at that point. Units
		  enclosed in a '{}' indicate that any protein unit
		  NOT enclosed can appear at that point.
	"""	

	def is_motif(string, motif = motif, get_motif_len = False):
		
		"""
		Returns True if string accords with motif,
		False otherwise. Currently does not check 
		input-validity.
		"""

		motif = motif.upper(); 
		m = {}; i = 0; char_count = 0
		while i < len(motif):
			if motif[i] == '[':
				end = i + motif[i:].index(']')
				m[char_count] = [c for c in motif[i+1:end]]
				i = end + 1
			elif motif[i] == '{':
				end = i + motif[i:].index('}')
				m[char_count] = []
				for c in set(rna_codons.values()):
					if c not in motif[i+1:end]:
						m[char_count] += [c]
				i = end + 1
			else:
				m[char_count] = [motif[i]]
				i += 1
			char_count += 1

		if get_motif_len:
			return len(m.keys())
		else:
			s = string = string.upper()
			if len(s) != len(m.keys()):
				print("Motif and string not of equal length.") 
				return			
			for i in range(len(s)):
				if s[i] not in m[i]:
					return False
			return True
	
	locations = []	
	m_len = is_motif(None, motif, get_motif_len = True)
	for i in range(len(prot) - m_len + 1):
		if is_motif(prot[i:i+m_len], motif):
			locations += [i+1]
	return locations		


def uniprot_get(uniprot_ID):

	"""
	uniprot_ID -> (string) Protein ID on uniprot.org

	Fetches the amino acid sequence of the given protein from the
	uniprot database. Returns a string. 
	"""

	url = "http://www.uniprot.org/uniprot/{}.fasta".format(uniprot_ID)
	data = urlopen(url).read().decode()

	string = data[data.index('\n')+1:].replace('\n', '')

	return (uniprot_ID, string)


def motif_in_proteins(protein_IDs, motif = "N{P}[ST]{P}", fout = True):

	"""
	Proteins -> Takes 3 types of values:
					- Single protein ID
					- List of protein IDs
					- File containing newline-separated IDs
	Motif -> A valid protein motif string

	Output:
		- Prints the locations at which the given motif appears in
		  each given protein.
		- Returns a list of tuples corresponding to the printed output.
		- If fout is true, output is written to an output file.
		  If a file '<fname>.txt'' is provided, output is written to
		  'output_<fname>.txt'; otherwise, writes to 'output_mprt.txt'.
	Requires access to uniprot.org.
	"""

	results = []
	fout_name = 'output_mprt.txt'

	if os.path.isfile(protein_IDs):
		prot_IDs = open(protein_IDs, 'r').read().split('\n')
		prot_IDs = [prot for prot in prot_IDs if prot not in ['', ' ']]
		prots = [uniprot_get(prot) for prot in prot_IDs]
		fout_name = "output_{}".format(protein_IDs)
	elif type(protein_IDs) == list:
		prots = [uniprot_get(prot) for prot in protein_IDs]
	else:
		try:
			prots = [uniprot_get(protein_IDs)]
		except Exception as ex:
			print(ex)
			raise ValueError('Invalid input type or nonexistent protein ID.')

	locations = [motif_location(prot[1], motif = motif) for prot in prots]
	if len(locations) != len(prots):
		raise Exception("Something went wrong... Exiting...")
	for i in range(len(prots)):
		results += [(prots[i][0], locations[i])] 
	results = [i for i in results if i[1] not in [[], '', None]]
	
	with open(fout_name, 'w') as fout:
		for tup in results:
			print(tup[0]); fout.write("{}\n".format(tup[0]))
			for i in tup[1]:
				print(i, end = ' '); fout.write("{} ".format(i)) 
			print(); fout.write('\n')

	return results















