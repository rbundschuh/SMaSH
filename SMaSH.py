#! /usr/bin/env python

from __future__ import print_function

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
# Maximilian Westphal, David Frankhouser, Carmine Sonzone, Peter G. Shields,
# Pearlly Yan, and Ralf Bundschuh
#
# Usability and Python3 update by Alan Hoyle, UNC Lineberger Bioinformatics

from collections import defaultdict
import math
import numpy
from scipy.special import comb
from scipy import stats
import argparse
from time import gmtime, strftime
import pysam
from decimal import Decimal
import pickle
import os
import sys
import glob

def eprint (*args, **kwargs):
    # function that easily prints to stderr
    print (*args, file=sys.stderr,**kwargs)

# calculate beta function
def B(alpha, beta_val):
	return Decimal(math.factorial(alpha - 1)) * Decimal(math.factorial(beta_val - 1)) / Decimal(math.factorial(alpha + beta_val - 1))

# Alpha/Beta for Beta functions
homozygous = 30
heterozygous = 2

# priors for being different samples or the same sample
p_s_0 = 0.99
p_s_1 = 0.01

# n_tot = total reads, n_A = alternate reads
def Q(n_tot, n_A, alpha, beta_val):
	return comb(n_tot, n_A, exact=True) * B(alpha + n_A, beta_val + n_tot - n_A) / B(alpha, beta_val)

# q = allele frequency
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
	print ("h or s out of range")


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

parser.add_argument('-bam', '--bam',action='store', dest='bam_comma_list', required=False, default='', 
	help="[deprecated] input BAM/SAM files to be tested (comma separated list) or 'ALL' to use all BAMs in current dir. ")
parser.add_argument('-i', '--positions', action='store', dest='infile', required=False,
	default= os.path.join(os.path.dirname(os.path.realpath(__file__)),'snps_hg19.vcf'), help='input locations file')
parser.add_argument('-o', '--output_file', action='store', dest='outname', required=False, default='pval_out.txt',
	help='output file name [pval_out.txt]')
parser.add_argument('-chr_index', '--chr_index', action='store', dest='chr_index', required=False, default=0, type=int,
	help='index of chromosome column in locations file [0]')
parser.add_argument('-pos_index', '--pos_index', action='store', dest='pos_index', required=False, default=1, type=int,
	help='index of position column in locations file [1]')
parser.add_argument('-ref_index', '--ref_index', action='store', dest='ref_index', type=int, required=False, default=3,
	help='reference allele index')
parser.add_argument('-alt_index', '--alt_index', action='store', dest='alt_index', type=int, required=False, default=4,
	help='alternate allele index')
parser.add_argument('-regenerate', '--regenerate', action='store_true', dest='regenerate', required=False, default='False', 
	help='regenerate SNP read counts and all calculations')
parser.add_argument('-output_dir', '--output_dir', action='store', required=False, default='.', 
	help='The directory to save output files.  [default: ./]')
parser.add_argument('-include_rgid', '--include_rgid', action='store_true',dest='include_rgid', required=False,
	default='False', 
	help="include BAM's Read Group ID value in output and only print the BAM's basename, not full path" )
parser.add_argument('bam',nargs='*', help = 'BAM/SAM files to check.  Note BAMs must end in .bam and be indexed')

args = parser.parse_args()

bams = args.bam
bam_comma_list = args.bam_comma_list
infile = args.infile
outname = args.outname
chr_index = args.chr_index
pos_index = args.pos_index
ref_index = args.ref_index
alt_index = args.alt_index
regenerate = args.regenerate
output_dir = args.output_dir
include_rgid = args.include_rgid


if  bams == ['ALL'] or bams == ['*']: 
	print (strftime("[%Y-%m-%d %H:%M:%S]"), "Selecting *.bam from current directory")
	bams = glob.glob("*.bam")

elif bam_comma_list == "ALL" or bam_comma_list == "*" or bam_comma_list == "*.bam" :
	print (strftime("[%Y-%m-%d %H:%M:%S]"), "Adding *.bam from current directory...")
	bams += glob.glob("*.bam")
elif len(bam_comma_list)>4:  # has to be greater than four for the file to be x.bam
	bams += bam_comma_list.split(',')

elif len(bam_comma_list)>0: 
	eprint("-bam[",bam_comma_list, "]is not valid.  use 'ALL' or a comma-separated list of BAM/SAM files.")
	sys.exit (1)

if len(bams) < 2:  
	eprint ('Not enough bams specified. ')
	sys.exit (1)


def dd():
	return defaultdict(list)

print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Reading variants in this VCF:', infile)

f = open (infile, 'r')
variants = f.read().splitlines()
f.close()

print (strftime("[%Y-%m-%d %H:%M:%S]"), 'num variants:',len(variants))

locs_name = infile.split('/')[-1]
if regenerate == True:
	data = defaultdict(dd)
	matrix_dic = defaultdict(dd)
else:
	if os.path.isfile(output_dir + os.path.sep + locs_name + ".mat.p") == True:
		print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Reading cached values from', output_dir + os.path.sep + locs_name + ".mat.p")

		matrix_dic = pickle.load( open(output_dir + os.path.sep + locs_name + ".mat.p", "rb"))
	else: 
		matrix_dic = defaultdict(dd)

	if os.path.isfile(output_dir + os.path.sep + locs_name + ".loc.p") == True:
		print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Reading cached values from', output_dir + os.path.sep + locs_name + ".loc.p")
		data = pickle.load(open(output_dir + os.path.sep + locs_name + ".loc.p", "rb"))
	else:
		data = defaultdict(dd)

# data = defaultdict(lambda: defaultdict(list))
# read bam/sam file type and open bam file
new_entry = False

rgids = {}
printable_header = ['.']

if bams ==['']:
	eprint ('ERROR No bams specified')
	parser.print_help()
	sys.exit()

for bam in bams:
	if bam.split('.')[-1] == 'bam': samfile = pysam.Samfile(bam, 'rb')
	elif bam.split('.')[-1] == 'sam': samfile = pysam.Samfile(bam, 'r')
	else: sys.exit('Cannot tell file type of '+bam+'. Make sure file extension is bam or sam.')

	if include_rgid == True:
		rgids[bam] = samfile.header['RG'][0]['ID']

	printable_bam = bam
	if include_rgid == True:
		printable_bam = os.path.basename(bam) + ' ('+rgids[bam]+')'
	printable_header = printable_header + [printable_bam]


	print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Reading sample variant read counts from', printable_bam)
	file_name = '.'.join(bam.split('.')[0:-1])
	chr_in_annot = False
	try: # checks the bam file for the chromosome annotation type and selects chr1 or 1 for example
		samfile.fetch("chr1", 1, 1) # just calling it to see if it generates an error
		chr_in_annot = True
		chrom_refname = "chr"
	except(ValueError):
		chrom_refname = ""
	if not chr_in_annot:
		try: # checks reversely if bam file contains 1 as chromosome; if not something is wrong with the bam or bai files
			samfile.fetch("1", 1, 1) # just calling it to see if it generates an error
		except(ValueError):
			eprint ('ERROR bam files neither contain "chr1" nor "1" as chromosomes; you may have forgotten to provide bam index files')
			sys.exit()
		

	for lines in variants:
		if lines[0] == "#": continue
		cols = lines.strip().split('\t')
		chrom = cols[chr_index]
		# checks if the bam file's chromosome annotation matches the vcf file annotation
		# and edits the chromosome accordingly to be able to properly retrieve reads from the bam file
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

if new_entry == True:
	pickle.dump(matrix_dic, open(output_dir + os.path.sep + locs_name + ".mat.p", "wb"))
	pickle.dump(data, open(output_dir + os.path.sep + locs_name + ".loc.p", "wb"))


header = ['.'] + bams
outfile = open(output_dir + os.path.sep + outname, 'w')
write_list(printable_header, outfile)
revoutfile = open(output_dir + os.path.sep + 'issameindividual_' + outname, 'w')
revoutfile.write('# probabilities to be from the same individual' + "\n");
write_list(printable_header, revoutfile)
comp_data = defaultdict(lambda: defaultdict(float))
rev_comp_data = defaultdict(lambda: defaultdict(float))
comps = set([])
p_out = open(output_dir + os.path.sep + outname + '.pvals.test.txt', 'w')
p_header = ['#sampA', 'sampB', 'loc', 'sampA [ref, alt]', 'sampB [ref, alt]', 'pval s=0', 'pval s=1']
write_list(p_header, p_out)

print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Beginning p-value calculations...')
for file_row in bams:
	printable_row = file_row
	if include_rgid == True:
		printable_row = os.path.basename(file_row)+'_('+rgids[file_row]+')'
	#outline = [file_row]
	for file_col in bams:
		printable_col = file_col
		if include_rgid == True:
			printable_col = os.path.basename(file_col)+ '_('+rgids[file_col]+')'

		compA = file_row + '-vs-' + file_col
		compB = file_col + '-vs-' + file_row
		if compA in comps or compB in comps: continue
		print (strftime("[%Y-%m-%d %H:%M:%S]"), 'Calculating p-value for', printable_row, 'vs', printable_col + ".")
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
			#print (loc, matrix_dic[loc][file_row], matrix_dic[loc][file_col], prob_0, prob_1, s_0, s_1)
			outline = [printable_row, printable_col, loc, data[loc][file_row], data[loc][file_col], prob_0, prob_1, s_0, s_1]
			write_list(outline, p_out)
		prob = s_0 * Decimal(p_s_0) / (s_0 * Decimal(p_s_0) + s_1 * Decimal(p_s_1))
		comp_data[file_row][file_col] = prob
		comp_data[file_col][file_row] = prob
		revprob = s_1 * Decimal(p_s_1) / (s_0 * Decimal(p_s_0) + s_1 * Decimal(p_s_1))
		rev_comp_data[file_row][file_col] = revprob
		rev_comp_data[file_col][file_row] = revprob

reformat_output = open(output_dir + os.path.sep + 'best_guesses.' + outname, 'w')
for file_row in bams:
	printable_row = file_row
	if include_rgid == True:
		printable_row = os.path.basename(file_row) + '_('+rgids[file_row]+')'
	outline = [printable_row]
	revoutline = [printable_row]
	output = [printable_row]
	for file_col in bams:
		printable_col = file_col
		if include_rgid == True:
			printable_col = os.path.basename(file_col) + '_('+rgids[file_col]+')'

		outline.append(comp_data[file_row][file_col])
		revoutline.append(rev_comp_data[file_row][file_col])
		if file_col == file_row: continue
		elif comp_data[file_row][file_col] == 'NOTEST': continue
		elif comp_data[file_row][file_col] <= 0.05: output.append([printable_col, comp_data[file_row][file_col]])
	write_list(outline, outfile)
	write_list(revoutline, revoutfile)
	write_list(output, reformat_output)

print (strftime("[%Y-%m-%d %H:%M:%S]"), 'SMaSH finished. Look for results in directory',output_dir + os.path.sep)
