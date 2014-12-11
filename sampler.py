## about the data:
## 1. safe mode: fix length (4 or 5), and fix coverage (1 or 2), and fix coverage
##	because the true case may be complicated, and we can use this to guarantee a converged result?
## 2. ...

## two things to demonstrate in the report:
## 1. trend of likelihood changes
## 2. trend of error rate changes


## questions:
## 1. read index start from 1 or 0?


import os
import sys
import numpy as np
import math


## make every parameter and drawer global
epsilon = 0.01
omega = 0.002
N_site = 0
N_ref = 0
N_reads = 0
position = []
TRUE = []
H = []
reference = []
S = []
reads = []
R = []
reads_reshaped = []
N_round = 10000



def like():
	global omega
	global epsilon
	global reads
	global R
	global H
	global reference
	global position

	sum = 0

	## sum for R-H related potentials, and S-H related potentials
	for i in range(N_reads):
		read = reads[i]
		for site in read:  # (index, allele)
			index = site[0]
			allele = site[1]
			if allele == 1:
				if R[i] == 1:  # favor above
					if H[0][index] == 1:
						sum += math.log(1 - epsilon)
					else:
						sum += math.log(epsilon)
				if R[i] == -1:  # favor below
					if H[1][index] == 1:
						sum += math.log(1 - epsilon)
					else:
						sum += math.log(epsilon)
			else:
				if R[i] == 1:  # favor above
					if H[0][index] == 1:
						sum += math.log(epsilon)
					else:
						sum += math.log(1 - epsilon)
				if R[i] == -1:  # favor below
					if H[1][index] == 1:
						sum += math.log(epsilon)
					else:
						sum += math.log(1 - epsilon)

	for i in range(N_site - 1):
		# zeta
		if H[0][i] == reference[S[0][i]][i]:
			sum += math.log(1 - omega)
		else:
			sum += math.log(omega)

		if H[1][i] == reference[S[1][i]][i]:
			sum += math.log(1 -omega)
		else:
			sum += math.log(omega)

		# tau
		tune1 = 1 # the ratio
		tune2 = 1000 # N in the literature
		r = (position[i+1] - position[i]) * 0.00000001
		if S[0][i] == S[0][i + 1]:
			sum += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
		else:
			sum += (1 - math.exp( - tune1 * r )) / tune2

		if S[1][i] == S[1][i + 1]:
			sum += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
		else:
			sum += (1 - math.exp( - tune1 * r )) / tune2

	# last zeta
	if H[0][N_site - 1] == reference[S[0][N_site - 1]][N_site - 1]:
		sum += math.log(1 - omega)
	else:
		sum += math.log(omega)

	if H[1][N_site - 1] == reference[S[1][N_site - 1]][N_site - 1]:
		sum += math.log(1 -omega)
	else:
		sum += math.log(omega)

	likelihood = math.exp(sum)
	return likelihood



def error():  # error rate: errors / (2 * haplotype length)
	global TRUE
	global H
	global N_site

	error = 0
	for i in range(N_site):
		if TRUE[0][i] != H[0][i]:
			error += 1
	for i in range(N_site):
		if TRUE[1][i] != H[1][i]:
			error += 1
	return error * 1.0 / (2 * N_site)



if __name__ == '__main__':


	print "This is the haplotype phasing program..."

	###================================ read the file to generate the following: initialization ===============================
	###========================================================================================================================
	## physical positions
	position = []  # the physical positions for the inferred haplotypes (only heterozygous sites in)
	file = open("xxx.legend", 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		pos = int((line.strip()).split(" ")[1])
		position.append(pos)
	N_site = len(position)


	## haplotype
	# real one
	file = open("xxx.phase", 'r')
	lines = file.readlines()
	file.close()
	h1 = map(lambda x: int(x), (lines[0].strip()).split(" ") )
	h2 = map(lambda x: int(x), (lines[1].strip()).split(" ") )
	TRUE = [h1, h2]  # two haplotype (only heterozygous sites in)
	# random generated one
	h1 = np.random.binomial(1, 0.5, N_site) * 2 - 1
	h2 = - h1
	H = [h1, h2]


	## reference; removed all homozygous sites according to the above testing sample already
	lines = lines[2:]
	reference = []
	for line in lines:
		ref = map(lambda x: int(x), (line.strip()).split(" ") )
		reference.append(ref)
	N_ref = len(reference)
	s1 = [0] * N_site  # without loss of generality, we can choose 0 as the start reference panel
	s2 = [0] * N_site
	S = [s1, s2]



	## reads, and reshaped reads
	file = open("xxx.reads", 'r')
	lines = file.readlines()
	file.close()
	N_reads = len(lines)
	reads = [] # each list is one short read; each tuple is a specific site, index and allele value
	for line in lines:
		line = map(lambda x: int(x), (line.strip()).split(" ")[1:] )
		read = []
		i = 0
		while i < len(line):
			index = line[i] - 1  # assume the index start from 0
			allele = line[i + 1]
			i += 2
			read.append((index, allele))
		reads.append(read)
	R = np.random.binomial(1, 0.5, N_reads) * 2 - 1  ## NOTE: if r = 1, favor the first haplotype; otherwise the second

	reads_reshaped = [[ ]] * N_site  # each list is one site; each tuple has index of read and corresponding allele value
	for i in range(len(reads)):
		read = reads[i]
		for site in read:  # site: (index, allele)
			index = site[0]
			allele = site[1]
			reads_reshaped[index].append((i, allele))
	###=========================================================================================================================
	###=========================================================================================================================



	round = 0
	for round in range(N_round):  ## or until equilibrium
		print "sampling round#",
		print round,
		print ":"



		###======================== sample the R with fixed H (S now is conditionally independent)
		for i in range(N_reads):
			read = reads[i]

			sum1 = 0  # favoring score for haplotype#1
			sum2 = 0  # favoring score for haplotype#2
			for site in read:
				index = site[0]
				allele = site[1]

				if allele * H[0][index] == 1:
					theta1 = math.log(1 - epsilon)
					theta2 = math.log(epsilon)
				else:
					theta1 = math.log(epsilon)
					theta2 = math.log(1 - epsilon)

				sum1 += theta1
				sum2 += theta2

			p1 = math.exp(sum1) / ( math.exp(sum1) + math.exp(sum2) )

			## sample r of this read
			r = np.random.binomial(1, p1) * 2 - 1
			R[i] = r

		###======================== sample the S with fixed H (R now is conditionally independent)
		### pending












		###======================== sample the H with fixed R and S (NOTE: to be checked)
		## sample h1 -- H[0]
		for i in range(N_site):
			sum1 = 0 # favor h1 = 1
			sum2 = 0 # favor h1 = -1
			sum3 = 0 # favor h2 = 1
			sum4 = 0 # favor h2 = -1

			for read in reads_reshaped[i]:  # all the reads under a specific site
				allele = read[1]
				if allele = 1:
					if R[i] == 1:  # for h1
						sum1 += math.log(1 - epsilon)
						sum2 += math.log(epsilon)
					else:  # for h2
						sum3 += math.log(1 - epsilon)
						sum4 += math.log(epsilon)
				else:
					if R[i] == 1:
						sum1 += math.log(epsilon)
						sum2 += math.log(1 - epsilon)
					else:
						sum3 += math.log(epsilon)
						sum4 += math.log(1 - epsilon)

			if reference[S[0][i]][i] == 1: # favor h1 = 1
				sum1 += math.log(1 - omega)
				sum2 += math.log(omega)
			else:  # favor h1 = -1
				sum1 += math.log(omega)
				sum2 += math.log(1 - omega)

			if reference[S[1][i]][i] == 1: # favor h2 = 1
				sum3 += math.log(1 - omega)
				sum4 += math.log(omega)
			else:  # favor h2 = -1
				sum3 += math.log(omega)
				sum4 += math.log(1 - omega)

			p1 = math.exp(sum1) / ( math.exp(sum1) + math.exp(sum2) )
			p2 = math.exp(sum3) / ( math.exp(sum3) + math.exp(sum4) )

			## sample h of this site in both haplotypes
			h1 = np.random.binomial(1, p1) * 2 - 1
			H[0][i] = h1
			h2 = np.random.binomial(1, p2) * 2 - 1
			#H[1][i] = h2
			H[1][i] = -h1  # they should be composite


		###======================== calculate (or summarize) the present likelihood value; and phasing errors
		## likelihood; NOTE: maybe every 10 rounds?
		if round % 10 = 0:
			print "present likelihood value is",
			print like()

		## error rate
		print "present error rate is:",
		print error()



	print "sampling done..."
	print "the final error rate is:",
	print error()