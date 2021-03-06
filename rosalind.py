#!/usr/bin/env python3
## http://rosalind.info/problems/list-view/

"""
Solutions to Rosalind bioinformatics problems.
"""

import os
import operator
from itertools import product
from urllib.request import urlopen
from math import log10, floor, ceil, factorial as fact
from collections import OrderedDict as Dict #dict is non-ordered
from scipy.misc import comb

rna_codons = {
    #'AUG' : 'Start',
    'UUU' : 'F',    'CUU' : 'L', 'AUU' : 'I', 'GUU' : 'V',
    'UUC' : 'F',    'CUC' : 'L', 'AUC' : 'I', 'GUC' : 'V',
    'UUA' : 'L',    'CUA' : 'L', 'AUA' : 'I', 'GUA' : 'V',
    'UUG' : 'L',    'CUG' : 'L', 'AUG' : 'M', 'GUG' : 'V',
    'UCU' : 'S',    'CCU' : 'P', 'ACU' : 'T', 'GCU' : 'A',
    'UCC' : 'S',    'CCC' : 'P', 'ACC' : 'T', 'GCC' : 'A',
    'UCA' : 'S',    'CCA' : 'P', 'ACA' : 'T', 'GCA' : 'A',
    'UCG' : 'S',    'CCG' : 'P', 'ACG' : 'T', 'GCG' : 'A',
    'UAU' : 'Y',    'CAU' : 'H', 'AAU' : 'N', 'GAU' : 'D',
    'UAC' : 'Y',    'CAC' : 'H', 'AAC' : 'N', 'GAC' : 'D',
    'UAA' : 'Stop', 'CAA' : 'Q', 'AAA' : 'K', 'GAA' : 'E',
    'UAG' : 'Stop', 'CAG' : 'Q', 'AAG' : 'K', 'GAG' : 'E',
    'UGU' : 'C'   , 'CGU' : 'R', 'AGU' : 'S', 'GGU' : 'G',
    'UGC' : 'C'   , 'CGC' : 'R', 'AGC' : 'S', 'GGC' : 'G',
    'UGA' : 'Stop', 'CGA' : 'R', 'AGA' : 'R', 'GGA' : 'G',
    'UGG' : 'W'   , 'CGG' : 'R', 'AGG' : 'R', 'GGG' : 'G'}

dna_codons = {
    #'ATG' : 'Start',
    'TTT' : 'F',    'CTT' : 'L', 'ATT' : 'I', 'GTT' : 'V',
    'TTC' : 'F',    'CTC' : 'L', 'ATC' : 'I', 'GTC' : 'V',
    'TTA' : 'L',    'CTA' : 'L', 'ATA' : 'I', 'GTA' : 'V',
    'TTG' : 'L',    'CTG' : 'L', 'ATG' : 'M', 'GTG' : 'V',
    'TCT' : 'S',    'CCT' : 'P', 'ACT' : 'T', 'GCT' : 'A',
    'TCC' : 'S',    'CCC' : 'P', 'ACC' : 'T', 'GCC' : 'A',
    'TCA' : 'S',    'CCA' : 'P', 'ACA' : 'T', 'GCA' : 'A',
    'TCG' : 'S',    'CCG' : 'P', 'ACG' : 'T', 'GCG' : 'A',
    'TAT' : 'Y',    'CAT' : 'H', 'AAT' : 'N', 'GAT' : 'D',
    'TAC' : 'Y',    'CAC' : 'H', 'AAC' : 'N', 'GAC' : 'D',
    'TAA' : 'Stop', 'CAA' : 'Q', 'AAA' : 'K', 'GAA' : 'E',
    'TAG' : 'Stop', 'CAG' : 'Q', 'AAG' : 'K', 'GAG' : 'E',
    'TGT' : 'C'   , 'CGT' : 'R', 'AGT' : 'S', 'GGT' : 'G',
    'TGC' : 'C'   , 'CGC' : 'R', 'AGC' : 'S', 'GGC' : 'G',
    'TGA' : 'Stop', 'CGA' : 'R', 'AGA' : 'R', 'GGA' : 'G',
    'TGG' : 'W'   , 'CGG' : 'R', 'AGG' : 'R', 'GGG' : 'G'}

prot_mass_table = {
    'water' : 18.01056,
    'A' : 71.03711 , 'C' : 103.00919, 'D' : 115.02694,
    'E' : 129.04259, 'F' : 147.06841, 'G' : 57.02146 ,
    'H' : 137.05891, 'I' : 113.08406, 'K' : 128.09496,
    'L' : 113.08406, 'M' : 131.04049, 'N' : 114.04293,
    'P' : 97.05276 , 'Q' : 128.05858, 'R' : 156.10111,
    'S' : 87.03203 , 'T' : 101.04768, 'V' : 99.06841 ,
    'W' : 186.07931, 'Y' : 163.06333 }


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


def dna_to_rna(dna_string):

    """
    dna_string -> A DNA string or a file containing one.

    Returns the RNA string transcribed from dna_string.
    If a file is provided, output is written to 'output_<fname>'.
    """

    if os.path.isfile(dna_string):
        f_out = True
        d = open(dna_string, 'r').read().replace('\n', '').upper()
    else:
        f_out = False
        d = dna_string = dna_string.upper()
    r = ''

    for nt in d:
        if nt == "T":
            r += "U"
        else:
            r += nt

    if f_out:
        with open("output_{}".format(dna_string), 'w') as fout:
            fout.write(r)
    return r


def reverse_complement(dna_string):

    """
    Returns the reverse complement of a string of DNA nucleotides.
    """

    s = dna_string = dna_string.upper()[::-1]

    complements = {"A":"T", "C":"G", "G":"C", "T":"A"}
    c = ''.join(complements[i] for i in s)
    return c


def fasta_read(fasta_file, names=True):

    """
    Takes a file in fasta format and returns list of
    (ID, genetic_string) tuples.
    """

    #names option not implemented yet

    f = fasta_file
    if not os.path.isfile(f):
        raise ValueError("Invalid file")
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


def gc_content(fasta_file):

    """
    Takes a fasta-formatted file containing 1 or more DNA strings
    and returns the ID and GC-content of the string with the
    greatest GC-content, ie the percentage of a given DNA string
    consiting of either G or C.
    """

    l = fasta_read(fasta_file)

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

    if len(s) != len(t):
        raise ValueError("Strings should be the same length")

    return sum((s[i] != t[i]) for i in range(len(s)))


def rabbits(n, m=0, k=1):

    """
    n -> number of months for breeding
    k -> number of rabbit-pair (M|F) offspring per generation

    Returns the number of rabbit pairs at the end of n months given
    that we start with 1 pair, k rabbit pairs are produced
    each generation and each rabbit pair may reproduce after 1 month
    of being born.
    """

    def real_rabbits(n, k=1):

        """
        Could be extended to account for death or to
        simulate population growth when rabbits don't stick in
        pairs or whatever. Too memory intesive for large populations.
        """

        class rabbit_pair:
            'A real rabbit.'

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

                litter = rabbit.reproduce(pairs=k)
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


    def math_rabbits(n, m, k=1):
        'Rabbits are numbers.'

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
        """
        g1, g2 -> genotypes belonging to k, m, or n

        Returns probability of dominant genotype in punnet square
        g1 x g2.
        """

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


def subs(string, substring, zero_based=True):

    """
    Returns the starting positions of substring in string.
    """

    s = string = string.upper()
    ss = substring = substring.upper()

    indexes = []

    for i in range(len(s) - len(ss) + 1):
        if s[i : i + len(ss)] == ss:
            if zero_based:
                indexes += [i]
            else:
                indexes += [i+1]

    return indexes


def cons(fasta_file):

    """
    Essential Goal: Average a set of genetic strings
    Returns: A Consensus string and Profile Matrix P for a fasta_file
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

    string_entries = fasta_read(fasta_file)
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

    p_syms_ordered = sorted(p, key=lambda entry: entry[0])

    with open("output_{}".format(fasta_file), 'w') as fout:
        print(consensus)
        fout.write("{}\n".format(consensus))
        for sym in p_syms_ordered:
            print("{}: ".format(sym), end=' ')
            fout.write("{}: ".format(sym))
            for i in p[sym]:
                print(i, end=' ')
                fout.write("{} ".format(i))
            print()
            fout.write("\n")


def overlap_graph(fasta_file):

    """
    Takes a fasta_file containing DNA strings and prints the
    adjacency list corresponding to O_3, ie the overlap graph
    of the strings where overlaps are at least 3 nucleotides long.
    (Also writes an output file for easy keeping)
    """

    data = fasta_read(fasta_file)
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

    with open("output_" + fasta_file, 'w') as fout:
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


def motif(fasta_file):

    """
    Returns the longest common substring of the genetic strings
    contained in fasta_file.
    """

    data = fasta_read(fasta_file)
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
    with open("output_" + fasta_file, 'w') as fout:
        fout.write(longest)
    return longest


def choose(n, k):

    return int(fact(n)/(fact(k)*fact(n-k)))

def permute(n, k):

    return int(fact(n)/fact(n-k))


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
          If a file '<fname>.txt' is provided, output is written to
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
            raise ValueError('Invalid input type or nonexistent',
                             'protein ID.')

    locations = [motif_location(prot[1], motif=motif) for prot in prots]
    if len(locations) != len(prots):
        raise Exception("Something went wrong... Exiting...")

    for i in range(len(prots)):
        results += [(prots[i][0], locations[i])]
    #results += [(prots[i][0], locations[i]) for i in range(len(prots))]
    results = [i for i in results if i[1] not in [[], '', None]]

    with open(fout_name, 'w') as fout:
        for tup in results:
            print(tup[0]); fout.write("{}\n".format(tup[0]))
            for i in tup[1]:
                print(i, end = ' '); fout.write("{} ".format(i))
            print(); fout.write('\n')

    return results


def num_mrna_from_prot(prot):

    """
    prot -> Protein string or file containing protein

    Returns n modulo 1,000,000, where n is the number of
    possible mrna strings which could encode the given protein.
    """

    if os.path.isfile(prot):
        f = open(prot, 'r')
        prot = f.read().replace('\n', '').upper()
    else:
        prot = prot.upper()

    count = 1
    inferences = [0 for i in range(len(prot))]

    for i in range(len(prot)):
        for key in rna_codons.keys():
            if rna_codons[key] == prot[i]:
                inferences[i] += 1
        count *= inferences[i]

    count *= 3 # Don't forget about the 3 possible stop codons
    return count % int(1E6)


def reading_frames(dna_string, simple = True):

    """
    dna_string -> A DNA string or a file containing one.

    Returns: If simple = True:
                A list containing the reading 6 reading frames
                pertaining to dna_string.
             If simple = False:
                A list of two 3-tuples; the first containing the
                5'-to-3' reading frames of dna_string, and the second
                containing the 3'-to-5' reading frames of dna_string
                (or vice-versa, depending on dna_string directionality.
    """

    if os.path.isfile(dna_string):
        f = open(dna_string, 'r').read()
        d = f.replace('\n', '')
    else:
        d = dna_string.upper()

    reading_frames = []
    d_c = reverse_complement(d)

    for s in [d, d_c]:
        frames = []
        for shift in range(3):
            frame = ''
            for i in range(shift, len(s), 3):
                codon = s[i:i+3]
                if len(codon) == 3:
                    frame += codon
            frames += [frame]
        if simple:
            reading_frames += frames
        else:
            reading_frames += tuple(frame for frame in frames)

    return reading_frames


def dna_to_proteins(dna_string):

    """
    dna_string -> A DNA string or a fasta file containing a
                  DNA string.

    Returns every distinct candidate protein that can be translated
    from the Open Reading Frames of dna_string.
    If a fasta file is provided, output is written to 'output_<filename>'.
    """

    if os.path.isfile(dna_string):
        d = fasta_read(dna_string)[0][1]
        f_out = True
    else:
        d = dna_string = dna_string.upper()
        f_out = False

    d_frames = reading_frames(d)
    proteins = []

    for frame in d_frames:

        starts = [i for i in subs(frame, "ATG") if (i+3)%3 == 0]
        if len(starts) == 0:
            continue

        for start in starts:
            s = frame[start:]
            prot = ''
            stop_found = False
            for i in range(0, len(s), 3):
                if dna_codons[s[i:i+3]] == 'Stop':
                    stop_found = True
                    break
                prot += dna_codons[s[i:i+3]]
            if prot and stop_found:
                proteins += [prot]

    proteins = list(set(proteins))

    if f_out:
        with open('output_{}'.format(dna_string), 'w') as fout:
            for prot in proteins:
                fout.write("{}\n".format(prot))

    return proteins


def permute(l):

    """
    l -> list of set of (preferably) ordered numbers
         (or objects, I guess)

    Returns the list of all permutations of l.
    """

    perms = []

    if len(l) == 0:
        return []
    elif len(l) == 1:
        return [l]
    else:
        for i in l:
            i_perms = permute([j for j in l if j != i])
            for perm in i_perms:
                perms += [[i] + perm]
        return perms


def permutations(n, f_out=True, give_total=False):

    """
    n -> Positive integer

    Returns the total number of permuations of length n, ie, the
    permutations of the first n positive integers, followed by a list
    of all such permuations.

    If f_out is true, output is written to 'output_permuations.txt'.
    If give_total is true, a tuple is returned, with the first element
    being the total number of permuations and the second element being
    the list of permutations.
    """

    total = fact(n)
    perms = permute(list(range(1, n+1)))

    if f_out:
        with open('output_permuations.txt', 'w') as fout:
            fout.write("{}\n".format(total))
            for perm in perms:
                for i in perm:
                    fout.write('{} '.format(i))
                fout.write('\n')

    if give_total:
        return total, perms
    else:
        return perms


def protein_mass(protein_string):

    """
    protein_string -> A protein sequence or a file containing one

    Returns the monoisotopic mass of the protein, assuming that all
    amino acids are residues.
    """

    if os.path.isfile(protein_string):
        prot = open(protein_string, 'r').read().replace('\n', '').upper()
    else:
        prot = protein_string.upper()

    mass = 0
    for i in prot:
        mass += prot_mass_table[i]
    return mass


def reverse_palindromes(dna_string, zero_based=True):

    """
    dna_string -> A DNA string or a fasta file containing one.

    Returns the position and length of every reverse palindrome in the
    string having length between 4 and 12. (Related to restriction
    enzymes)

    If a fasta file is provided, output is written to 'output_<fname>'.
    """

    if os.path.isfile(dna_string):
        d = fasta_read(dna_string)[0][1]
        f_out = True
    else:
        d = dna_string.upper()
        f_out = False

    pals = []

    for i in range(len(d) - 3):
        for testlen in range(4, 13, 2):
            a = d[i : i + testlen]
            if len(a) >= testlen:
                b = reverse_complement(a)
                if a == b:
                    if zero_based:
                        pals += [(i, testlen)]
                    else:
                        pals += [(i+1, testlen)]

    if f_out:
        with open('output_ {}'.format(dna_string), 'w') as fout:
            for i in pals:
                fout.write('{} {}\n'.format(i[0], i[1]))
    return pals


def dna_and_introns_to_protein(fasta_file):

    """
    fasta_file -> A fasta file with the first entry being a DNA string d
                  and the remaining entries being substrings of d
                  acting as introns of d (which implies that the
                  substrings should not overlap).

    Returns the protein transcribed by the mRNA translated from the exons
    of the given DNA string. Output is written to 'output_<fname>'.
    """
    ## To Do:
    #### Write this in a more cohesive and less retarded way.
    #### Check that the intron positions and lengths don't overlap and
    #### therefore produce valid exons.
    #### Check that the exon string is evenly divisble by 3 and therefore
    #### nicely transcribable.
    #### Somehow allow for multiple dna strings to be processed.
    #### Somehow account for UTR's.

    dna_and_introns = fasta_read(fasta_file)

    dna = dna_and_introns[0][1]
    introns = [i[1] for i in dna_and_introns[1:]]

    coding_region, protein = str(), str()
    intron_intervals = []

    for i in introns:
        positions = subs(dna, i)
        for pos in positions:
            intron_intervals += [range(pos, pos + len(i))]

    for i in range(len(dna)):
        is_exon_part = True
        for rng in intron_intervals:
            if i in rng:
                is_exon_part = False
        if is_exon_part:
            coding_region += dna[i]

    for i in range(0, len(coding_region) - 2, 3):
        s = coding_region[i:i+3]
        try:
            if dna_codons[s] == 'Stop':
                protein += ' '
            else:
                protein += dna_codons[s]
        except Exception as ex:
            print(ex)

    with open('output_{}'.format(fasta_file), 'w') as fout:
        fout.write(protein)
    return protein


def combinations(items, k, rep=False):

    """
    Returns a list of all combinations of length n of items.
    If rep = True, reptititions of items are allowed.
    Returns list in lexicographic order.
    """

    #if type(items) != list:
    #   raise ValueError("Items must be a list.")
    if type(k) != int:
        raise ValueError("n must be an integer.")

    if k == 0:
        return []
    elif k == 1:
        return [[i] for i in items]
    elif k > len(list(items)) and not rep:
        return []

    combos = []
    if not rep:
        for i in range(len(items) - k + 1):
            head = [items[i]]
            for tail in combinations(items[i+1:], k - 1):
                combos += [head + tail]
        return combos
    else:
        for i in range(len(items)):
            head = [items[i]]
            for tail in combinations(items[i:], k-1, rep=1):
                combos += [head + tail]
        return combos


## Should be extended to allow for repetition.
def signed_combinations(n, k):

    """
    n -> Positive integer
    k -> Positive integer => n

    Returns a list of all k-combinations of the first n positive integers
    including each integer negated. Does not include combinations like
    -1, 1 or 2, -2 nor repetitive combinations such as -1, -1 or 2, 2.
    """

    def rec_combine(l, k):

        """
        l -> List of range of integers

        Recursively generates the signed k-combinations of l.
        """

        combinations = []
        if len(l) == 0:
            return []
        elif k == 0:
            return []
        elif k == 1:
            return [[i] for i in l] + [[-i] for i in l]
        else:
            for i in range(len(l)):
                head = l[i]
                tails = rec_combine(l[i+1:], k-1)
                for tail in tails:
                    combinations += [[head] + tail] + [[-head] + tail]
        return combinations

    if type(n) != int or type(k) != int:
        raise ValueError("n and k must be integers.")
    elif n < k:
        return []
    else:
        combinations = rec_combine(range(1, n+1), k)
        return combinations


def signed_permutations(n):

    """
    n -> Positive integer

    Returns the total number of signed permutations of length
    n, followed by a list of all such permutations.

    Writes output to 'output_<n>_sign_perms'.
    """

    total = 0
    signed_perms = []

    for combination in signed_combinations(n, n):
        for permutation in permute(combination):
            total += 1
            signed_perms += [permutation]

    with open('output_{}_sign_perms'.format(n), 'w') as fout:
        fout.write(str(total) + '\n')
        for i in signed_perms:
            for j in i:
                fout.write('{} '.format(j))
            fout.write('\n')
    return total, signed_perms


def lex_perms(ordered_alphabet, n):

    """
    ordered_alphabet -> A collection of symbols defining and ordered
                        alphabet.
    n -> If n >= len(alphabet), all permutations
         of ordered_alphabet are given in lexicographic order, as
         defined by ordered_alphabet.
         If n < len(alphabet), all permutations of length n are given
         in lexicographic order.

    Note: Includes permutations with repetitions. I have yet to (but
          should) implement non-repetition-including lex_perms.

    Output is written to 'output_lex_perms.txt'.
    """

    perms = []
    alph = ordered_alphabet

    if n == 0:
        return []
    elif n == 1:
        return [[i] for i in alph]
    else:
        for i in range(len(alph)):
            head = [alph[i]]
            for tail in lex_perms(alph, n-1):
                perms += [head + tail]

    #with open('output_lex_perms.txt', 'w') as fout:
    #   for i in perms:
    #       for j in i:
    #           fout.write('{}'.format(j))
    #       fout.write('\n')
    return perms


def has_perfect_matching(rna_string):

    """
    rna_string -> An rna string or a fasta file containing one.

    Returns true if there exists a perfect match in rna_string,
    else returns False.
    """

    if os.path.isfile(rna_string):
        r = fasta_read(rna_string)[0][1]
    else:
        r = rna_string.upper()

    if len(r)%2 != 0:
        return False
    else:
        return r.count('A') == r.count('U') and r.count('G') == r.count('C')


def perfect_match_count(rna_string):

    """
    rna_string -> An rna string containing as many occurences of
                  'A' as 'U' and as many occurences of 'G' as 'C'
                  or a fasta file containing such a string.

    Returns the number of possible perfect matchings in the bonding
    graph of rna_string.
    """

    if os.path.isfile(rna_string):
        rna = fasta_read(rna_string)[0][1].upper()
    else:
        rna = rna_string.upper()

    if not has_perfect_matching(rna_string):
        raise Exception("Unequal number of occurences of A-U or G-C.")
    else:
        a_count, g_count = rna.count('A'), rna.count('G')
        return fact(a_count)*fact(g_count)


def partial_perms_count(n, k):

    """
    n -> Positive integer
    k -> Positive integer <= n

    Returns the total number of partial permutations of length k
    of the first n positive integers modulo 1,000,000.
    """

    a = choose(n, k)
    b = fact(k)
    return (a*b)%int(1E6)


def recursive_LIS(sequence):

    """
    Takes an integer sequence or a file containing one and returns
    the longest increasing subsequence and the longest decreasing
    subsequence. If a file is given, it should be formatted as follows:

    <Length of sequence>
    <S e q u e n c e>

    Output is written to 'output_<fname>'
    """

    if type(sequence) == str:
        if os.path.isfile(sequence):
            f_out = True
            with open(sequence, 'r') as f:
                s = f.readlines()[1].replace(' ', '')
        else: f_out = False
    else:
        f_out = False
        try:
            s = [int(i) for i in sequence]
        except:
            raise ValueError("First arg must be an integer sequence.")

    def rec_subseqs(seq):

        """
        Takes a list of non-repeated integers and returns all
        increasing subsequences in the list if inc = True, all
        decreasing subsequences in the list if inc = False.
        """
        ## Way too slow to work on large sequences u_u
        subseqs = []
        s_len = len(seq)
        if s_len == 0:
            return []
        elif s_len == 1:
            return [seq]
        elif s_len == 2:
            if seq[0] > seq[1]:
                return [[seq[0]]]
            else:
                return [seq]
        else:
            for i in range(s_len):
                head = seq[i]
                tails = list((n for n in seq[i+1:] if n > head))
                #tails = list(filter(lambda x: x > h, seq[i+1:]))
                #tails = [j for j in seq[i+1:] if j > head]
                for j in range(len(tails)):
                    for subseq in rec_subseqs(tails[j:]):
                        subseqs += [[head] + subseq]
        return subseqs

    max_inc, max_dec = [], []
    for seq in rec_subseqs(s):
        if len(seq) > len(max_inc):
            max_inc = seq
    for seq in rec_subseqs([-i for i in s]):
        if len(seq) > len(max_dec):
            max_dec = seq
    max_dec = [-i for i in max_dec]

    if f_out:
        with open('output_{}'.format(sequence), 'w') as fout:
            for seq in [max_inc, max_dec]:
                for i in seq:
                    fout.write('{} '.format(i))
                fout.write('\n')
    return max_inc, max_dec


def LIS(sequence):

    """
    Takes an integer sequence or a file containing one and returns
    the longest increasing subsequence and the longest decreasing
    subsequence. If a file is given, it should be formatted as follows:

    <Length of sequence>
    <S e q u e n c e>

    Output is written to 'output_<fname>'
    """

    if type(sequence) == str:
        if os.path.isfile(sequence):
            f_out = True
            with open(sequence, 'r') as f:
                s = [int(i) for i in f.readlines()[1].split()]
        else:
            s = [int(i) for i in sequence]
            f_out = False
    else:
        f_out = False
        try:
            s = [int(i) for i in sequence]
        except:
            raise ValueError("First arg must be an integer sequence.")

    def get_lis(s):
        'Returns increasing subsequences of sequence s.'

        lis = [[s[0]]]
        for i in range(1, len(s)):
            l_i = []
            max_l_j = 0
            for j in range(0, i):
                if s[i] > lis[j][-1]:
                    if len(lis[j]) > max_l_j:
                        max_l_j = len(lis[j])
                        l_i = lis[j]
            lis += [l_i + [s[i]]]
        return lis

    LIS_ = []
    for seq in get_lis(s):
        if len(seq) > len(LIS_):
            LIS_ = seq
    LDS = []
    for seq in get_lis([-i for i in s]):
        if len(seq) > len(LDS):
            LDS = seq
    LDS = [-i for i in LDS]

    if f_out:
        with open('output_{}'.format(sequence), 'w') as fout:
            for seq in [LIS_, LDS]:
                for i in seq:
                    fout.write('{} '.format(i))
                fout.write('\n')
    return LIS_, LDS


def edges_to_form_tree(graph_file):

    """
    graph_file -> A file formatted as follows:
                      n
                      <adjacency list>
                  Where n is the number of nodes in a graph containing
                  no cycles. The adjacencies should be listed as follows:
                      1 2
                      3 4
                      ...
    Returns the minimum number of edges required to form a tree from
    the given graph.
    """

    if os.path.isfile(graph_file):
        with open(graph_file, 'r') as f:
            f = f.readlines()
            n = int(f[0])
            adj_list = [[int(j) for j in i.split()] for i in f[1:]]
    else:
        raise ValueError('Input must be a file. See docstring.')

    return n - len(adj_list) - 1


def strprob(input_file):

    """
    inpWut_file -> A file formatted as follows:
                      <DNA-string>
                      <Arbitrary, Space-separated GC-contents>
    Returns an array of the same length as the GC-contents list representing
    the probability that a random string constructed with each GC-content will
    match DNA-string exactly.
    Output is written to 'output_<input_file>'.
    """

    if os.path.isfile(input_file):
        with open(input_file, 'r') as f:
            f = f.readlines()
            d = f[0].strip('\n').upper()
            gc_contents = [float(i) for i in f[1].split()]
    else:
        raise ValueError('Input must be a file. \n{}'.format(strprob.__doc__))

    def get_prob(gc_content):

        prob = 1
        C = G = gc_content/2
        A = T = (1 - gc_content)/2
        for i in d:
            prob *= eval(i)
        return prob

    probs = [log10(get_prob(gc)) for gc in gc_contents]
    with open('output_{}'.format(input_file), 'w') as fout:
        for i in probs:
            fout.write('{} '.format(i))
    return probs


def shortest_superstring(fasta_file):

    """
    Takes a fasta_file containing DNA reads and returns the shortest
    superstring containing all the given strings. For practical purposes,
    this function is written on the assumption that there exists a unique
    way to reconstruct the entire superstring from the reads by "gluing
    together" pairs of reads that overlap by more than half their length.
    """

    if os.path.isfile(fasta_file):
        reads = [i[1] for i in fasta_read(fasta_file)]
    else:
        raise ValueError("Input must be a fasta file.")


    def merge_longest_overlap(reads):

        """
        Takes a list of reads and returns the same list with the two
        longest-overlapping strings having been joined.
        """

        longest = str()
        for i in range(len(reads)):
            for j in range(len(reads)):
                if i == j:
                    continue

                if len(reads[i]) <= len(reads[j]):
                    s1, s2 = reads[i], reads[j]
                else:
                    s1, s2 = reads[j], reads[i]

                overlap = str()
                testlen = len(s1)
                while testlen > floor(len(s1)/2):
                    if s1[:testlen] == s2[len(s2)-testlen:]:
                        overlap = s2[:len(s2)-testlen] + s1
                        break
                    elif s2[:testlen] == s1[len(s1)-testlen:]:
                        overlap = s1[:len(s1)-testlen] + s2
                        break
                    testlen -= 1

                if len(overlap) > len(longest):
                    longest = overlap
                    indices = (i, j)
        if longest:
            reads = [reads[i] for i in range(len(reads)) if i not in indices]
            reads += [longest]
        return reads


    # This would never stop executing if there were a non-overlapping string.
    # If this module were robust, this would be taken into consideration.
    while len(reads) > 1:
        reads = merge_longest_overlap(reads)

    with open('output_{}'.format(fasta_file), 'w') as fout:
        fout.write(reads[0])
    return reads[0]



## Implement the input_type parameter to skip type checking to speed up
## repetitive calls to the function. All the type checking and subsequent
## control flows make this more laggy
## than it ought to be.

def spliced_motif(fasta_file, zero_based=True):

    """
    fasta_file -> A fasta-formatted file containing two strings; the first
                  being a DNA string dna_string and the second being
                  a motif string. Also accepts an indexable container
                  of two such strings.

    Returns the indices of dna_string where motif occurs as a
    subsequence of dna_string. Does not return more than one occurrence.

    If zero_based is True, indices returned are 0-based, else 1-based.
    """

    if type(fasta_file) == str:
        f_out = True
        if os.path.isfile(fasta_file):
            with open(fasta_file, 'r') as f:
                data = [t[1] for t in fasta_read(fasta_file)]
                d = data[0]
                motif = data[1]
        else:
            raise ValueError("Invalid input. See .__doc__")
    else:
        f_out = False
        try:
            d = fasta_file[0]
            motif = fasta_file[1]
        except:
            raise ValueError("Invalid input parameter. See docstring.")

    locations = []
    dna_index, motif_index = 0, 0
    while dna_index < len(d) and motif_index < len(motif):
        if d[dna_index] == motif[motif_index]:
            locations += [dna_index]
            motif_index += 1
        dna_index += 1

    if not zero_based:
        locations = [i+1 for i in locations]

    if f_out:
        with open('output_{}'.format(fasta_file), 'w') as fout:
            for i in locations:
                fout.write('{} '.format(i))
    return locations


def is_subsequence(seq, subseq):

    """
    seq -> A DNA string
    subseq -> A DNA string of length less then seq

    Returns True is subseq is a subsequence of seq, else False.
    """

    ss_len = len(subseq)
    s_index, ss_index = 0, 0
    while s_index < len(seq) and ss_index < ss_len:
        if seq[s_index] == subseq[ss_index]:
            ss_index += 1
        s_index += 1

    return ss_index == ss_len


def is_purine(nucleobase):
    "Returns True if nucleobase is a purine, else False."

    return True if nucleobase.upper() in ['A', 'G'] else False

def is_pyrimadine(nucleobase):
    "Returns True if nucleobase is a pyrimadine, else False."

    return True if nucleobase.upper() in ['C', 'T', 'U'] else False


## If this were a robust function it would give accurate ratios
## for both DNA and RNA.

def trans_tranv_ratio(fasta_file):

    """
    fasta_file -> A fasta-formatted file containing two DNA strings
                  of equal length.

    Returns the translation to transversion ratio between the
    two strings.
    """

    if os.path.isfile(fasta_file):
        with open(fasta_file, 'r') as f:
            data = [i[1] for i in fasta_read(fasta_file)]
            s1, s2 = data[0], data[1]
    else:
        raise ValueError("Input must be a fasta file. See docstring.")
    if len(s1) != len(s2):
        raise ValueError("Input strings must be of equal length.")

    translations, transversions = int(), int()
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            if sum([is_purine(s1[i]), is_purine(s2[i])]) in [0, 2]:
                translations += 1
            else:
                transversions += 1
    R = translations/transversions

    with open('output_{}'.format(fasta_file), 'w') as fout:
        fout.write(str(R))
    return R


def max_matchings_count(rna_string):

    """
    rna_stirng -> An RNA string or a fasta file containing one.

    Returns the total possible number of maximum matchings of
    basepair edges in the bonding graph of rna_string.
    """

    if os.path.isfile(rna_string):
        r = fasta_read(rna_string)[0][1].upper()
    else:
        r = rna_string.upper()

    a_count, u_count = r.count('A'), r.count('U')
    g_count, c_count = r.count('G'), r.count('C')

    pur_min, pur_max = min(a_count, u_count), max(a_count, u_count)
    pyr_min, pyr_max = min(g_count, c_count), max(g_count, c_count)

    matching_count = fact(pur_min) * choose(pur_max, pur_min)
    matching_count *= fact(pyr_min) * choose(pyr_max, pyr_min)

    return matching_count


def catalan_number(n):

    "Returns the nth Catalan number."

    catalans = [1, 1]
    for i in range(2, n+1):
        c = int()
        for k in range(1, i+1):
            c += catalans[k-1] * catalans[i-k]
        catalans += [c]
    return catalans[-1]


def pncm_count(rna_string):

    """
    rna_string -> An RNA string or a fasta file containing one,
                  where the number of occurences of A equals that
                  of U and the number of occurences of G equals that
                  of C.

    Returns the number of perfect non-crossing matchings in the
    bonding graph of rna_string modulo 1,000,000.
    """

    if os.path.isfile(rna_string):
        rna = fasta_read(rna_string)[0][1]
    else:
        rna = rna_string.upper()

    if not has_perfect_matching(rna_string):
        return 0

    bondings = {'A':'U', 'U':'A', 'G':'C', 'C':'G'}
    encountered = {}

    def rec_pncm_count(rna):

        "Counts the number of perfect noncrossing matchings on rna."

        if len(rna) <= 2: #Function is only passed perfect matchings
            return 1
        else:
            if rna in encountered:
                return encountered[rna]

            pncm = 0
            for i in range(1, len(rna), 2):
                if bondings[rna[0]] == rna[i]:
                    l, r = rna[1:i], rna[i+1:]
                    if has_perfect_matching(l): #then so does r
                        l_count = rec_pncm_count(l)
                        r_count = rec_pncm_count(r)
                        pncm += l_count * r_count

            encountered[rna] = pncm_count
            return pncm

    return rec_pncm_count(rna) % int(1E6)


def inodes(n):

    """
    Returns the number of internal nodes in unrooted binary
    tree with n leaves.
    """

    return n - 2


def corr(fasta_file):

    """
    fasta_file -> A fasta-formatted file containing DNA reads of
                  equal length where some contain single-nucleotide
                  errors. For each read, one of the following applies:

                     - The read is correct and appears at least twice
                       in the set (possibly as its reverse complement).
                     - The read is incorrect and appears only once,
                       and has a hamming distance of 1 from exactly
                       one correct read in the set (or its reverse
                       complement).

    Returns a list of all read corrections (with 1 symbol substitution)
    in the form: [old read]->[new read]. Output is written to
    'output_<fasta_file>'.
    """

    if os.path.isfile(fasta_file):
        reads = [i[1] for i in fasta_read(fasta_file)]
    else:
        raise ValueError("Input must be a proper fasta file. See .__doc__")

    corrections, correct_reads, incorrect_reads = [], [], []

    for i in reads:
        if reads.count(i) + reads.count(reverse_complement(i)) > 1:
            if i not in correct_reads:
                correct_reads += [i]
        else:
            incorrect_reads += [i]

    for i in incorrect_reads:
        for j in correct_reads:
            if hamm(i,j) == 1:
                corrections += [(i, j)]
            elif hamm(i, reverse_complement(j)) == 1:
                corrections += [(i, reverse_complement(j))]

    with open('output_{}'.format(fasta_file), 'w') as fout:
        for i in set(corrections):
            fout.write('{}->{}\n'.format(i[0], i[1]))
    return list(set(corrections))


## Needs a small refactoring to account for values of k other than 4.
## Written in this dumb way using ordered dict because of failure of
## 'AAAA'.count('AA'), for example. Should be 3 but gives 2.

def k_mer_comp(dna_string):

    """
    dna_string -> A DNA string or a fasta file containing one.

    Returns the 4-mer compositoin of dna_string.
    """

    if os.path.isfile(dna_string):
        dna = fasta_read(dna_string)[0][1]
        f_out = True
    else:
        dna = dna_string.upper()
        f_out = False

    k_mers = [''.join(i for i in j) for j in lex_perms(['ACGT'], 4)]
    composition = Dict([[i, 0] for i in k_mers])

    for i in range(0, len(dna)-3):
        composition[dna[i:i+4]] += 1

    if f_out:
        with open('output_{}'.format(dna_string), 'w') as fout:
            for i in composition.keys():
                fout.write('{} '.format(composition[i]))

    return [composition[key] for key in composition]


def failure_array(string):

    """
    string -> A [genetic] string of a fasta_file containing one.

    Returns the failure array of string.
    Output is written to 'f_array.txt' or 'output_<string>' if string
    is a file.
    """

    if os.path.isfile(string):
        s = fasta_read(string)[0][1]
        f_out = True
    else:
        s = string.upper()
        f_out = False

    f_arr = [0]*len(s)
    length = 0
    for i in range(1, len(s)):
        while length > 0 and s[length] != s[i]:
            length = f_arr[length-1]
        if s[length] == s[i]:
            length += 1
        f_arr[i] = length

    if f_out:
        with open('output_{}'.format(string), 'w') as fout:
            for i in f_arr:
                fout.write('{} '.format(i))
    return f_arr


## The following ranking algorithm can solve the worst case input in
## 5 minutes or so, but there is a yet-quicker way about this abusing the
## nature of itertools.product. It's worth learning how to construct from
## scratch the fastest algorithm at solving this.
def lex_sort_strings(alphabet, n):

    """
    alphabet -> A list or string containing/defining an ordered alphabet.
    n -> A non-negative integer less than or equal to the length
         of alphabet.

    Returns all strings of length at most n formed from alphabet, ordered
    lexicographically.
    Output is written to 'ouput_<n>_lex_sort.txt'.
    """

    if not type(alphabet) in [str, list]:
        raise ValueError("Alphabet must be of type string or list.")
    elif not (type(n) == int and n > 0  and n <= len(alphabet)):
        raise ValueError(print(lex_sort_strings.__doc__))

    alph = [' '] + [i for i in alphabet]
    pperms = [list(product(alphabet, repeat=i) for i in range(1, n+1))]
    perms = []
    for i in pperms:
        for j in i:
            for k in j:
                s = ''.join(c for c in k)
                while len(s) < n:
                    s += ' '
                perms += [s]

    ranks = [0]*len(perms)
    ordered = [0]*len(ranks)
    for i in range(len(perms)):
        rank = 0
        for j in range(len(perms)):
            for c in range(n):
                if perms[i][c] != perms[j][c]:
                    if alph.index(perms[i][c]) > alph.index(perms[j][c]):
                        rank += 1
                    break
        ranks[i] = rank

    for i in range(len(ranks)):
        place = ranks[i]
        ordered[place] = perms[i].strip(' ')

    with open('output_{}_lex_perms'.format(n), 'w') as fout:
        for i in ordered:
            fout.write('{}\n'.format(i))
    return ordered


def distance_matrix(fasta_file):

    """
    fasta_file -> A fasta-formatted file containing dna strings
                  of equal length.

    Returns a matrix P where Pij is the p distance between the
    i_th and j_th given string.
    Output is written to 'output_<fasta_file>'.
    """

    if os.path.isfile(fasta_file):
        data = [i[1] for i in fasta_read(fasta_file)]
    else:
        raise ValueError(print(distance_matrix.__doc__))

    P = [[0 for j in range(len(data))] for i in range(len(data))]

    for i in range(len(data)):
        s1 = data[i]
        for j in range(i):
            s2 = data[j]
            p_dist = sum((s1[c] != s2[c]) for c in range(len(s1)))
            p_dist = p_dist / len(s1)
            P[i][j], P[j][i] = p_dist, p_dist

    with open('output_{}'.format(fasta_file), 'w') as fout:
        for i in P:
            for j in i:
                fout.write('{} '.format(j))
            fout.write('\n')
    return P


def sizeof_powerset(n):

    """
    Returns the size of the powerset of the first n positive integers
    modulo 1,000,000.
    """

    return (2 ** n) % int(1E6)


def reversal_distance(seq_file, sortfirst=False):

    """
    seq_file -> A file containing eqaul-length sequence pairs formatted
                as follows:
                    s e q 1
                    s e q 2

                    s e q 3
                    s e q 4
                    ...

    Returns the reversal distances of each seqeuence pair.
    Output is written to 'output_<seq_file>'.

    If sortfirst = True, only the first pair is considered and the output is
    returned and written to the output file as follows:

        n
        i0 i1
        i2 i3
        ...

    where n is the reversal distance and the i-pairs are reversal indices
    which result in the reversal sort.
    """

    if os.path.isfile(seq_file):
        f = open(seq_file).readlines()
        f = [i.strip('\n').split() for i in f if i not in ['', '\n']]
        f = [[int(i) for i in j] for j in f]
        if not sortfirst:
            pairs = [(f[i], f[i+1]) for i in range(0, len(f), 2)]
        else:
            pairs = [(f[0], f[1])]
    else:
        raise ValueError('Argument must be a valid file.\n{}'
                         .format(reversal_distance.__doc__))


    def breakpoints(a):
        'Returns a list of breakpoint indices of a.'

        b = []
        a = [0] + a + [len(a) + 1]
        for i in range(len(a)-1):
            if abs(a[i] - a[i+1]) != 1:
                b += [i]
        return b


    def inverse_perm(a):
        'Returns the inverse permutation of a.'

        perm = [0]*len(a)
        for i in range(len(a)):
            perm[a[a[i] - 1] - 1] = a[i]
        return perm


    def apply_permutation(a, b):
        'Applies permutation b to a.'

        return [a[i-1] for i in b]


    def rev_dist(a, b, sortfirst=sortfirst):
        """
        Greedy breakpoint sort algorithm for finding reversal distance.
        Finding the reversal distance between a and b is the same as finding
        the distance between b^-1 * a and the identity permutation.
        The identity permutation has 0 break points, so the goal of the
        algorithm is to maximally reduce the number of breakpoints at each
        branching point.
        """

        def breaks_to_strips(breaks):
            'Converts breakpoints to a set of candidate strips.'

            strips = [i for i in product(breaks, repeat=2)]
            strips = [i for i in strips if i[1] - i[0] > 1]
            return strips

        rev = lambda p, i, j: p[:i] + p[i:j+1][::-1] + p[j+1:]

        inv_b = inverse_perm(b)
        p = apply_permutation(inv_b, a)

        breaks = breakpoints(p)
        break_count = len(breaks)
        if break_count == 0:
            return 0

        strips = breaks_to_strips(breaks)
        reversals = [rev(p, i[0], i[1]) for i in strips]

        d = 0
        while True:
            d += 1

            # Future candidate permutations will have highest break_delta
            candidates = {delta : [] for delta in [0, 1, 2]}
            break_delta = 0
            for branch in reversals:
                temp_breaks = breakpoints(branch)
                temp_break_len = len(temp_breaks)
                # This test or the other break_count test might be redundant.
                if temp_break_len == 0:
                    return d
                if break_count - temp_break_len > break_delta:
                    break_delta = break_count - temp_break_len
                    candidates[break_delta] += [branch]
                elif break_count - temp_break_len == break_delta:
                    if branch not in candidates[break_delta]:
                        candidates[break_delta] += [branch]

            break_count -= break_delta
            if break_count == 0:
                return d

            new_revs = []
            for branch in candidates[break_delta]:
                temp_breaks = breakpoints(branch)
                temp_strips = breaks_to_strips(temp_breaks)
                new_revs += [rev(branch, i[0], i[1] - 1) for i in temp_strips]
            reversals = new_revs

    distances = []
    for i in pairs:
        # For some reason, one order may not be optimal.
        distances += [min(rev_dist(i[0], i[1]), rev_dist(i[1], i[0]))]

    with open('output_{}'.format(seq_file), 'w') as fout:
        for distance in distances:
            fout.write('{} '.format(distance))
    return distances


def reversal_sort(input_file):

    """
    pass
    """

    return reversal_distance(input_file, sortfirst=True)


def rstr(input_file):

    """
    input_file -> A list containing (in order) a sample size N,
                  GC-content, and DNA string or such a file, formatted as:

                      N GC-Content
                      DNA string

    Returns the probability that if N random strings (of equal length
    to DNA string) are generated with GC-Content, at least 1 of the strings
    will equal DNA string.
    """

    try:
        try:
            with open(input_file, 'r') as f:
                input_file = f.read().split()
        except:
            pass
        n = int(input_file[0])
        gc = float(input_file[1])
        s = str(input_file[2])
    except:
        raise ValueError('Invalid input.\n{}'.format(rstr.__doc__))

    p_cg = gc/2
    p_at = (1-gc)/2
    s_prob = p_cg ** len([1 for i in s if i in ['G', 'C']])
    s_prob *= p_at ** len([1 for i in s if i in ['A', 'T']])

    return 1 - (1- s_prob) ** n


def nCk(n, k):

    'n >= k >= 0'

    return (fact(n) / (fact(n-k)*fact(k))) % int(1E6)


def binomial_sum(n, m):

    'Returns the sum of nCk for m <= k <= n modulo 1,000,000.'

    s = sum([comb(n, k, exact=True) % int(1E6) for k in range(m, n+1)])
    return s % int(1E6)


def newick_distance(input_file, quiet=True):

    """
    input_file -> A file formatted as follows:

                      Newick-formatted tree 1;
                      x1 y1

                      Newick-formatted tree 2;
                      x2 y2
                      ...

    input_file: A collection of n trees in Newick format, with each tree
        followed by the labels of two nodes in the tree.

    Returns: A collection of n positive integers, where n1 represents
        the distance between x1 and y1 in newick_tree 1.

    Output is written to 'output_<input_file>'.
    """

    if os.path.isfile(input_file):
        with open(input_file, 'r') as f:
            f = f.readlines()
            f = [i.strip(';\n') for i in f if i != '\n']
            D = [[f[i], f[i+1]] for i in range(0, len(f), 2)]
    else:
        raise ValueError("Invalid argument.\n{}"
                         .format(newick_distance.__doc__))


    class Node(object):

        """
        Tree node.
        """

        def __init__(self, number, parent, name=None):

            self.number = number
            self.name = [name, 'Node_{}'.format(number)][name is None]
            self.parent = parent
            self.children = []

        def report(self, indent=1):
            'Utilize me!'

            print('{}Node : {}'.format(indent*' ', self.name))
            if self.children != []:
                print('{}Children = ['.format((indent+2)*' '))
                for child in self.children:
                    child.report(indent=indent+4)
                print('{}]'.format((indent+2)*' '))
            else:
                print('{}[is leaf]'.format((indent + 2)*' '))

        def add_child(self, child):
            'Add a child node.'
            self.children += [child]


    class Newick(object):

        """
        Newick tree.
        """

        def __init__(self, tree):

            self.rawtree = tree
            self.tree = []
            self.nodes = []
            self.parse(tree)
            self.node_refs = {node.name: node for node in self.nodes}


        def __repr__(self):

            return str(self.tree)


        def report(self, indent=2):
            'Prints information about the tree\'s nodes'

            print('\n\n\nRaw Tree = {}\nsplit tree = {}\ntree report:\n'
                  .format(self.rawtree, self.tree))
            roots = []
            for n in self.nodes:
                if n.parent.number == -1:
                    roots += [n]
            print('Base Nodes = ')
            for n in roots:
                if n.children == []:
                    print('{}{} [is leaf]'.format(indent*' ', n.name))
                else:
                    print('{}{}'.format(indent*' ', n.name))
            print()
            for node in roots:
                if node.children != []:
                    node.report(indent=indent)
            print()


        def parse(self, tree):

            """
            Parses the tree, fills self.nodes and creates relations
            among nodes.
            """

            tree = tree.replace('(', '( ').replace(')', ' )')
            tree = tree.replace(',', " ").split()
            self.tree = tree

            cur_parent = Node(-1, None)
            for i in range(len(tree)):

                if tree[i] == '(':
                    cur_parent = Node(len(self.nodes), cur_parent)
                    self.nodes += [cur_parent]
                    if len(self.nodes) > 1:
                        cur_parent.parent.add_child(cur_parent)

                elif tree[i][0] == ')':
                    if len(tree[i]) > 1:
                        cur_parent.name = tree[i][1:]
                    cur_parent = cur_parent.parent

                elif tree[i] == '':
                    n = Node(len(self.nodes), cur_parent)
                    cur_parent.add_child(n)
                    self.nodes += [n]

                else:
                    n = Node(len(self.nodes), cur_parent, name=tree[i])
                    cur_parent.add_child(n)
                    self.nodes += [n]


        def distance(self, node1, node2):
            """
            node1, node2 -> Names of nodes in tree.
            Returns the distance between node1 and node2 in tree.
            """

            if node1 not in self.node_refs or node2 not in self.node_refs:
                return None
            else:
                n1, n2 = self.node_refs[node1], self.node_refs[node2]

            n1_path, n2_path = [], []

            cur_node = n1
            while cur_node.number != -1:
                n1_path += [cur_node.name]
                cur_node = cur_node.parent
            cur_node = n2
            while cur_node.number != -1:
                n2_path += [cur_node.name]
                cur_node = cur_node.parent
            n1_path, n2_path = n1_path[::-1], n2_path[::-1]

            overlap = []
            for i in range(len(min(n1_path, n2_path, key=len))):
                if n1_path[i] == n2_path[i]:
                    overlap += [n1_path[i]]

            d = len(n1_path) + len(n2_path) - 2*len(overlap)
            if not quiet:
                print('search nodes = {}'.format([node1, node2]))
                print('\nn1_path = {}\nn2_path = {}\n'
                      .format(n1_path, n2_path))
                print('overlap = {}'.format(overlap))
                print('distance = {}'.format(d))
            return d


    distances = []
    for i in D:
        tree = i[0]
        nodes = i[1].split()
        t = Newick(tree)
        if not quiet:
            t.report()
        d = t.distance(nodes[0], nodes[1])
        distances += [d]

    with open('output_{}'.format(input_file), 'w') as fout:
        fout.write(' '.join(str(d) for d in distances))
    return distances


def LCS(fasta_file, indexed=False):

    """
    fasta_file -> A fasta-formatted file containing two DNA strings.
                  or an indexable container of two such strings.

    Returns the longest common subsequence of the two strings.
    Output is written to 'output_<fasta_file>'.
    """

    if os.path.isfile(str(fasta_file)):
        f_out = True
        strings = [i[1] for i in fasta_read(fasta_file)[:2]]
        s1, s2 = strings[0], strings[1]
    else:
        f_out = False
        try:
            s1, s2 = fasta_file[0].upper(), fasta_file[1].upper()
        except:
            raise ValueError("Invalid input.\n{}".format(LCS.__doc__))

    if indexed:
        index = []

    if len(s1) > len(s2):
        s, t = s1, s2
    else:
        s, t = s2, s1

    Arr = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    for i in range(len(s)):
        for j in range(len(t)):
            if s[i] == t[j]:
                Arr[i+1][j+1] = Arr[i][j] + 1
            else:
                Arr[i+1][j+1] = max(Arr[i+1][j], Arr[i][j+1])


    seq = ''
    i, j = len(s), len(t)
    while i*j != 0:
        if Arr[i][j] == Arr[i-1][j]:
            i -= 1
        elif Arr[i][j] == Arr[i][j-1]:
            j -= 1
        else:
            seq += s[i-1]
            if indexed:
                index += [[i-1, j-1]]
            i -= 1
            j -= 1
    seq = seq[::-1]

    if f_out:
        with open('ouput_{}'.format(fasta_file), 'w') as fout:
            fout.write(seq)
    if indexed:
        index = index[::-1]
        return [[seq[i], index[i]] for i in range(len(seq))]
    else:
        return seq


def edit_distance(s, t, quiet=True):

    'Returns the edit distance between strings s and t.'

    s, t = str(s).upper(), str(t).upper()

    if len(s) > len(t):
        s, t = t, s
    lens, lent = len(s), len(t)

    if not quiet:
        print('\ns = {}\nt = {}\n'.format(s, t))

    M = [[i for i in range(lens + 1)]]
    M += [[j + 1] + [0 for i in range(lens)] for j in range(lent)]

    for i in range(1, lent + 1):
        for j in range(1, lens + 1):
            if s[j-1] == t[i-1]:
                d = min(M[i-1][j-1], M[i-1][j] + 1, M[i][j-1] + 1)
            else:
                d = min(M[i-1][j-1], M[i-1][j], M[i][j-1]) + 1
            M[i][j] = d

    if not quiet:
        print(' ' * 3 + "\' \' " + '  '.join(i for i in s))
        temp_t = ' ' + t
        for i in range(len(M)):
            if i == 0:
                print("\' \'", end='')
            else:
                print(' ' + temp_t[i], end=' ')
            print(M[i])
        print('\n>> Edit distance = {}\n'.format(M[-1][-1]))

    return M[-1][-1]


def edit(fasta_file, quiet=True):

    """
    fasta_file -> A fasta-formatted file containing two strings.

    Returns the edit distance between the strings in fasta_file.
    """

    if os.path.isfile(fasta_file):
        s, t = [i[1] for i in fasta_read(fasta_file)][:2]
        return edit_distance(s, t, quiet=quiet)
    else:
        raise ValueError("Invalid input.\n{}".format(edit.__doc__))


def string_prob(string, GC_content):

    """
    Returns the probability of string occuring when constructed with
    GC-content GC_content.
    """

    s = string.upper()

    p_gc = GC_content / 2
    p_at = (1 - GC_content) / 2

    p = p_at ** (s.count('A') + s.count('T'))
    p *= p_gc ** (s.count('G') + s.count('C'))

    return p


def EX_occurences(input_file, n=int, s=str, gc=float):

    """
    input_file -> None if keywords are to be used, else a file formatted as
                  follows:

                      n   # Length of a DNA string
                      s   # A potential DNA substring of length <= n
                      gc1 gc2 ... gcM   # M GC-content values

    Returns a list of length M where M[i] is the expected number of
    occurences of s in a DNA string of length n constructed with given
    GC-content gc[i].
    If you only want the result of once such case, pass input_file as None
    and use the keywords.
    """

    if not input_file:
        return (n - (len(s) - 1)) * string_prob(s, gc)
    else:
        with open(input_file, 'r') as f:
            f = f.readlines()[:3]
            n = int(f[0].strip())
            s = f[1].strip().upper()
            M = [float(i) for i in f[2].split()]

    fits = n - (len(s) - 1)
    EX = lambda gc: fits * string_prob(s, gc)
    P = [EX(gc) for gc in M]

    with open('output_{}'.format(input_file), 'w') as fout:
        fout.write(' '.join(str(i) for i in P))
    return P


def motzkin_matchings(fasta_file):

    """
    fasta_file -> A fasta-formatted file containing an RNA string or, simply,
                  an RNA string.

    Returns the total number of noncrossing matchings of basepair edges
    in the bonding graph of the given RNA string, modulo 1,000,000.
    """

    try:
        if os.path.isfile(fasta_file):
            rna = fasta_read(fasta_file)[0][1]
        else:
            rna = fasta_file.upper()
    except:
        raise ValueError("Invalid input.\n{}"
                         .format(motzkin_matchings.__doc__))


    bondings = {'A':'U', 'U':'A', 'G':'C', 'C':'G'}
    memo = {}

    def rec_motz(rna):
        'Recursively counts all noncrossing matchings of rna.'

        if len(rna) <= 1:
            return 1
        elif rna in memo:
            return memo[rna]

        parts = []
        for i in range(1, len(rna)):
            if bondings[rna[0]] == rna[i]:
                parts += [(rna[1:i], rna[i+1:])]

        partsum = sum(rec_motz(part[0]) * rec_motz(part[1]) for part in parts)
        memo[rna] = partsum + rec_motz(rna[1:])
        return memo[rna]

    return rec_motz(rna) % int(1E6)



def scsp(input_file):

    """
    input file -> A file containing to strings formatted as follows:

                      string_s
                      string_t

    Returns the shortest common supersequence of s and t.
    If an input file is given, output is written to 'output_<input_file>'.
    """

    try:
        if os.path.isfile(input_file):
            f_out = True
            with open(input_file, 'r') as f:
                f = [i.strip() for i in f.readlines() if i not in ['', ' ']]
                s, t = f[0], f[1]
        else:
            f_out = False
            s, t = input_file[0], input_file[1]
    except:
        raise ValueError("Invalid input.\n{}".format(scsp.__doc__))

    lcs = LCS([s, t], indexed=True)

    spsq = ''
    pos0 = (0, 0)
    for i in range(len(lcs)):

        pos1 = (lcs[i][1][0], lcs[i][1][1])
        s_block = s[pos0[0]:pos1[0]]
        t_block = t[pos0[1]:pos1[1]]

        spsq += s_block + t_block
        spsq += lcs[i][0]

        pos0 = (pos1[0] + 1, pos1[1] + 1)

    spsq += s[lcs[-1][1][0] + 1:] + t[lcs[-1][1][1] + 1:]

    if f_out:
        with open('output_{}'.format(input_file), 'w') as fout:
            fout.write(spsq)
    return spsq


def set_operations(input_file, quiet=True):

    """
    input_file -> A file containing a positive integer n and two subsets
                  A and B of {1, 2, ... n}, formatted as follows:

                      n
                      {A}
                      {B}

    Returns A ∪ B, A ∩ B, A - B, B - A, A¬ and B¬ where the complements are
    based off of {1, 2, ... n}.
    Output is written to 'output_<input_file>' and printed if quiet = False.
    """

    if type(input_file) == str:
        if os.path.isfile(input_file):
            f = open(input_file, 'r').readlines()
            f = [i.strip('\n') for i in f]
            n = int(f[0])
            A = [i for i in sorted(eval(f[1]))]
            B = [i for i in sorted(eval(f[2]))]
    else:
        raise ValueError("Invalid input...\n{}".format(set_operations.__doc__))

    union = sorted(list(set(A + B)))
    intersection = [i for i in A if i in B]
    difA = [i for i in A if i not in B]
    difB = [i for i in B if i not in A]
    Ac = [i for i in range(1, n + 1) if i not in A]
    Bc = [i for i in range(1, n + 1) if i not in B]

    with open('output_{}'.format(input_file), 'w') as fout:
        for i in [union, intersection, difA, difB, Ac, Bc]:
            fout.write('{' + ', '.join(str(i[j]) for j in range(len(i) - 1)))
            fout.write(', ' + str(i[-1]) + '}\n')
            if not quiet:
                print('{', end='')
                print(', '.join(str(i[j]) for j in range(len(i)-1)), end=', ')
                print(str(i[-1]) + '}')


def de_bruijn(input_file):

    """
    input_file -> A file containing DNA strings of equal length (not exceeding
                  50 bp) separated by newlines corresponding to a set
                  S of (k+1)-mers of some genetic string.

    Returns the adjacency list as a set of tuples corrseponding to the de
    Bruijn graph of S ∪ Src where Src is the set of the reverse complements
    of S.

    De Bruijn graph Bk of order k, a digraph on S ∪ Src is defined as follows:

        - The nodes of Bk correspond to all k-mers that are present as a
          substring of a (k+1)-mer from S ∪ Src.
        - The edges of Bk are enconded by the (k+1)-mers of S ∪ Src such that
          for each (k+1)-mer r in S ∪ Src, the directed edge
          (r[0:k-1], r[1:k]) is formed.

    Output is written to 'output_<input_file>' as the directed edges (tuples)
    of the  De Bruijn graph, separated by newlines.
    """

    try:
        with open(input_file, 'r') as f:
            data = [i.strip() for i in f.readlines() if i not in ['\n', '']]
            data = list(set(data))
    except:
        raise ValueError("Invalid input.\n{}".format(de_bruijn.__doc__))

    data += list(set([reverse_complement(i) for i in data]))
    B = [(r[:-1], r[1:]) for r in data]

    with open('output_{}'.format(input_file), 'w') as fout:
        fout.write('\n'.join('({}, {})'.format(i[0], i[1]) for i in B))
    return B






###################################
############ To Do ################
###################################
##
##  Learn how itertools.product is
##+ so fast.
##
##  Do the same with scipy.comb.
##
##  Optimize, refactor, etc.
##
###################################
