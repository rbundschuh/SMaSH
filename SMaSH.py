#
# COPYRIGHT (c)2016 THE OHIO STATE UNIVERSITY
# ALL RIGHTS RESERVED
#
# PERMISSION IS GRANTED TO USE, COPY, CREATE DERIVATIVE WORKS AND
# REDISTRIBUTE THIS SOFTWARE AND SUCH DERIVATIVE WORKS FOR NONCOMMERCIAL
# EDUCATION AND RESEARCH PURPOSES, SO LONG AS NO FEE IS CHARGED, AND SO
# LONG AS THE COPYRIGHT NOTICE ABOVE, THIS GRANT OF PERMISSION, AND THE
# DISCLAIMER BELOW APPEAR IN ALL COPIES MADE; AND SO LONG AS THE NAME OF
# THE OHIO STATE UNIVERSITY IS NOT USED IN ANY ADVERTISING OR PUBLICITY
# PERTAINING TO THE USE OR DISTRIBUTION OF THIS SOFTWARE WITHOUT
# SPECIFIC, WRITTEN PRIOR AUTHORIZATION.
#
# THIS SOFTWARE IS PROVIDED AS IS, WITHOUT REPRESENTATION FROM THE OHIO
# STATE UNIVERSITY AS TO ITS FITNESS FOR ANY PURPOSE, AND WITHOUT
# WARRANTY BY THE OHIO STATE UNIVERSITY OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE OHIO STATE
# UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES, INCLUDING SPECIAL,
# INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, WITH RESPECT TO ANY
# CLAIM ARISING OUT OF OR IN CONNECTION WITH THE USE OF THE SOFTWARE,
# EVEN IF IT HAS BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.
#
# SMaSH - Sample Matching using SNPs in Humans
#
# Maximilian Westpha, David Frankhouser, Carmine Sonzone, Peter G. Shields,
# Pearlly Yan, and Ralf Bundschuh
#
from collections import defaultdict
import math
import numpy as np
from scipy.misc import comb
from scipy import stats
import argparse
from time import gmtime, strftime
import pysam
from decimal import Decimal
import pickle
import os
import glob

#calculate beta function
def B(alpha, beta_val):
	return Decimal(math.factorial(alpha - 1)) * Decimal(math.factorial(beta_val - 1)) / Decimal(math.factorial(alpha + beta_val - 1))

#Alpha/Beta for Beta functions
homozygous = 30
heterozygous = 2

#priors for being different samples or the same sample
p_s_0 = 0.99
p_s_1 = 0.01

#n_tot = total reads, n_A = alternate reads
def Q(n_tot, n_A, alpha, beta_val):
	return comb(n_tot, n_A, exact=True) * B(alpha + n_A, beta_val + n_tot - n_A) / B(alpha, beta_val)

#q = allele frequency
def p_s(q, h, s):
	if s == 1:
		if h == 0: return (1-q)**2 #WW-WW
		elif h == 1: return (2*q)*(1-q) #WA-WA
		elif h == 2: return q**2 #AA-AA
		elif h >= 3 and h <= 8: return 0. #WW-WA, WA-WW, WW-AA, AA-WW, AA-WA, WA-AA
	elif s == 0:
		if h == 0: return (1-q)**4 #WW-WW
		elif h == 1: return 4 * q**2 * (1-q)**2 #WA-WA
		elif h == 2: return q ** 4 #AA-AA
		elif h == 3 or h == 4: return 2*q*(1-q) ** 3 #WW-WA, WA-WW
		elif h == 5 or h == 6: return q**2 * (1-q) ** 2 #WW-AA, AA-WW
		elif h == 7 or h == 8: return 2 * q**3 * (1-q) #AA-WA, WA-AA
	else: print "h or s out of range"


#set s1 = WW, s2 = AA, s3 = WA
def calc_matrix(n_tot, n_A):
	Q_WW = Q(n_tot, n_A, 1, homozygous)# * Decimal(p_s(p, "WW"))
	Q_AA = Q(n_tot, n_A, homozygous, 1)# * Decimal(p_s(p, "AA"))
	Q_WA = Q(n_tot, n_A, 2, 2)# * Decimal(p_s(p, "WA"))
	matrix = [Q_WW, Q_AA, Q_WA]
	return matrix

def comb_matrix(mA, mB, q, s):
	h_sum = 0.
	if s == 0: maximum = 9
	elif s == 1: maximum = 3
	for h in range(0, maximum):
		p = p_s(q,h,s)
		if h == 0:
			product = mA[0] * mB[0]
		elif h == 1:
			product = mA[2] * mB[2]
		elif h == 2:
			product = mA[1] * mB[1]
		elif h == 3:
			product = mA[0] * mB[2]
		elif h == 4:
			product = mA[2] * mB[0]
		elif h == 5:
			product = mA[0] * mB[1]
		elif h == 6:
			product = mA[1] * mB[0]
		elif h == 7:
			product = mA[1] * mB[2]
		elif h == 8:
			product = mA[2] * mB[1]
		tmp_sum = float(product) * p
		h_sum += tmp_sum
	return h_sum


def write_list(inList, outFile): outFile.write("\t".join( map( str, inList ) ) + "\n" )

parser = argparse.ArgumentParser()
parser.add_argument('-bam', action='store', dest='bam', required=False, default='ALL', 
	help='input bam/sam files to be tested (comma separated list) or ALL to use all bam in current dir. Note bam files must end in .bam and be indexed')
parser.add_argument('-i', action='store', dest='infile', required=False,
	default='/fs/lustre/osu7905/references/snps/snps_snp_id.vcf', help='input locations file')
parser.add_argument('-o', action='store', dest='outname', required=False, default='pval_out.txt',
	help='output file name [pval_out.txt]')
parser.add_argument('-chr_index', action='store', dest='chr_index', required=False, default=0, type=int,
	help='index of chromosome column in locations file [0]')
parser.add_argument('-pos_index', action='store', dest='pos_index', required=False, default=1, type=int,
	help='index of position column in locations file [1]')
parser.add_argument('-ref_index', action='store', dest='ref_index', type=int, required=False, default=3,
	help='reference allele index')
parser.add_argument('-alt_index', action='store', dest='alt_index', type=int, required=False, default=4,
	help='alternate allele index')
parser.add_argument('-regenerate', action='store', dest='regenerate', type=str, required=False, choices=['True', 'False'],
	default='False', help='regenerate SNP read counts and all calculations? [False]')

args = parser.parse_args()

bams = args.bam
infile = args.infile
outname = args.outname
chr_index = args.chr_index
pos_index = args.pos_index
ref_index = args.ref_index
alt_index = args.alt_index
regenerate = args.regenerate

if bams == "ALL": bams = glob.glob("*.bam")
else: bams = bams.split(',')

def dd():
	return defaultdict(list)

locs_name = infile.split('/')[-1]
if regenerate == "True":
	data = defaultdict(dd)
	matrix_dic = defaultdict(dd)
else:
	if os.path.isfile(locs_name + ".mat.p") == True:
		matrix_dic = pickle.load( open(locs_name + ".mat.p", "rb"))
	else: 
		matrix_dic = defaultdict(dd)

	if os.path.isfile(locs_name + ".loc.p") == True:
		data = pickle.load(open(locs_name + ".loc.p", "rb"))
	else:
		data = defaultdict(dd)

#data = defaultdict(lambda: defaultdict(list))
#read bam/sam file type and open bam file
new_entry = False

for bam in bams:
	if bam.split('.')[-1] == 'bam': samfile = pysam.Samfile(bam, 'rb')
	elif bam.split('.')[-1] == 'sam': samfile = pysam.Samfile(bam, 'r')
	else: sys.exit('Cannot tell file type of input bam/sam file. Make sure file extension is bam or sam.')
	print 'Reading sample', bam, 'for variant read counts.', strftime("%Y-%m-%d %H:%M:%S")
	file_name = '.'.join(bam.split('.')[0:-1])
	try: #checks the bam file for the chromosome annotation type and selects chr1 or 1 for example
		for alignedread in samfile.fetch("chr1", 1, 1):
			chr_in_annot = True
		chrom_refname = "chr"
	except(ValueError):
		chrom_refname = ""
	with open(infile, 'r') as readfile:

		for lines in readfile:
			if lines[0] == "#": continue
			cols = lines.strip().split('\t')
			chrom = cols[chr_index]
			#checks if the bam file's chromosome annotation matches the vcf file annotation
			#and edits the chromosome accordingly to be able to properly retrieve reads from the bam file
			if "chr" not in chrom and chrom_refname == "chr": chrom = "chr" + chrom
			elif "chr" in chrom and chrom_refname == "":
				if chrom != "chrX" and chrom != "chrY" and chrom != "chrM":
					tmp_chrom = []
					for i in list(chrom):
						if i.isdigit() == True: tmp_chrom.append(i)
					chrom = ''.join(tmp_chrom)
				else: chrom = chrom[-1]
			loc = cols[chr_index] + ":" + cols[pos_index]
			pos = int(cols[pos_index])
			ref = cols[ref_index]
			alt = cols[alt_index]
			pysam_pos = pos - 1
			reads = []
			if data[loc][bam] == []:
				new_entry = True
				for alignedread in samfile.fetch(chrom, pos -1, pos): #pysam fetches with standard coordinates
					try:
						index = alignedread.positions.index(pos - 1) #pysam lists with python 0-based coordinates
						reads.append(alignedread.query[index])
					except(ValueError):
						continue #ValueError occurs when the read covers the requested base's position via splicing
				nts = []
				nts.append(reads.count(ref))
				non_ref = len(reads) - reads.count(ref)
				nts.append(non_ref)
				data[loc][bam] = nts
				INFO = cols[7].split(';')
				#looks for AF in any field of the info column
				#for field in INFO:
				#	if field.split('=')[0] == 'AF': 
				#		AF=float(field.split('=')[1])
						#break
				AF = float(INFO[1].split('=')[1])
				data[loc]['AF'] = AF
			if matrix_dic[loc][bam] == []:
				new_entry = True
				n_A = data[loc][bam][1]
				n_tot = sum(data[loc][bam])
				matrix = calc_matrix(n_tot, n_A)
				matrix_dic[loc][bam] = matrix
	readfile.close()

if new_entry == True:
	pickle.dump(matrix_dic, open(locs_name + ".mat.p", "wb"))
	pickle.dump(data, open(locs_name + ".loc.p", "wb"))


header = ['.'] + bams
outfile = open(outname, 'w')
write_list(header, outfile)
comp_data = defaultdict(lambda: defaultdict(float))
comps = set([])
p_out = open(outname + '.pvals.test.txt', 'w')
p_header = ['#sampA', 'sampB', 'loc', 'sampA [ref, alt]', 'sampB [ref, alt]', 'pval s=0', 'pval s=1']
write_list(p_header, p_out)
#matrix_dic = defaultdict(lambda: defaultdict(list))
print 'Beginning p-value calculations.', strftime("%Y-%m-%d %H:%M:%S")
for file_row in bams:
	#outline = [file_row]
	for file_col in bams:
		compA = file_row + '-vs-' + file_col
		compB = file_col + '-vs-' + file_row
		if compA in comps or compB in comps: continue
		print 'Calculating p-value for', file_row, 'vs', file_col + ".", strftime("%Y-%m-%d %H:%M:%S")
		comps.add(compA)
		comps.add(compB)
		s_0 = Decimal(1.)
		s_1 = Decimal(1.)
		for loc in sorted(data.keys()):
			#added to remove instances of no reads
			q_snp = data[loc]['AF']
			matrixA = matrix_dic[loc][file_row]
			matrixB = matrix_dic[loc][file_col]
			prob_0 = comb_matrix(matrixA, matrixB, q_snp, 0)
			prob_1 = comb_matrix(matrixA, matrixB, q_snp, 1)
			s_0 *= Decimal(prob_0)
			s_1 *= Decimal(prob_1)
			#print loc, matrix_dic[loc][file_row], matrix_dic[loc][file_col], prob_0, prob_1, s_0, s_1
			outline = [file_row, file_col, loc, data[loc][file_row], data[loc][file_col], prob_0, prob_1, s_0, s_1]
			write_list(outline, p_out)
		prob = s_0 * Decimal(p_s_0) / (s_0 * Decimal(p_s_0) + s_1 * Decimal(p_s_1))
		comp_data[file_row][file_col] = prob
		comp_data[file_col][file_row] = prob

reformat_output = open('best_guesses.' + outname, 'w')
for file_row in bams:
	outline = [file_row]
	output = [file_row]
	for file_col in bams:
		outline.append(comp_data[file_row][file_col])
		if file_col == file_row: continue
		elif comp_data[file_row][file_col] == 'NOTEST': continue
		elif comp_data[file_row][file_col] <= 0.05: output.append([file_col, comp_data[file_row][file_col]])
	write_list(outline, outfile)
	write_list(output, reformat_output)

print 'finished', strftime("%Y-%m-%d %H:%M:%S")
