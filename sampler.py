## about the data:
## 1. safe mode: fix length (4 or 5), and fix coverage (1 or 2), and fix coverage
##	because the true case may be complicated, and we can use this to guarantee a converged result?
## 2. ...

## two things to demonstrate in the report:
## 1. trend of likelihood changes
## 2. trend of error rate changes


################################################################################################
import os
import sys
import numpy as np
import math


## make every parameter and drawer global
##=======================================
epsilon = 0.01
omega = 0.002
N_site = 0
N_ref = 0
N_reads = 0
position = []
TRUE = []
H = []
reference = []
emission = []
S = []
reads = []
R = []
reads_reshaped = []
N_round = 10000
N_test = 100  # the number of tested allele sites
N_test2 = 50  # the number of reference panels used here

# the following two are used in \tau, they should be tuned
tune1 = 0.0001 # the ratio
tune2 = 1000 # N in the literature
##=======================================



def like():  # global likelihood under the present parameter setting
	global omega
	global epsilon
	global reads
	global R
	global H
	global reference
	global position
	global tune1
	global tune2

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
			sum += math.log(1 - omega)
		else:
			sum += math.log(omega)

		# tau
		r = (position[i + 1] - position[i])
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
		sum += math.log(1 - omega)
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
		if TRUE[1][i] != H[1][i]:
			error += 1
	return error * 1.0 / (2 * N_site)



def sampler_R():  # sampling R with fixed H (S is conditionally independent)
	global N_reads
	global reads
	global H
	global R
	global epsilon

	for i in range(N_reads):
		read = reads[i]

		sum1 = 0  # favoring score for haplotype#1
		sum2 = 0  # favoring score for haplotype#2
		for site in read:
			index = site[0]
			allele = site[1]
			theta1 = 0
			theta2 = 0

			## check all sites, and allocate score based on the evidence from allele and h
			if allele * H[0][index] == 1:
				theta1 = math.log(1 - epsilon)
			else:
				theta1 = math.log(epsilon)
			## because H[0][i] may be different from H[1][i]
			if allele * H[1][index] == 1:
				theta2 = math.log(1 - epsilon)
			else:
				theta2 = math.log(epsilon)


			sum1 += theta1
			sum2 += theta2

		p1 = math.exp(sum1) / ( math.exp(sum1) + math.exp(sum2) )

		## sample r of this read
		r = np.random.binomial(1, p1) * 2 - 1
		R[i] = r

	return



def sampler_S():  # we use the sampling procedure described in the paper
	global emission
	global N_ref
	global N_site
	global S
	global H
	global position
	global tune1
	global tune2

	## emission matrix: [  [[1-omega, omega, 1-omega, ...],  [1-omega, omega, 1-omega, ...]],  []  ]
	## transition matrix: called when necessary

	F1 = {}
	F2 = {}
	for i in range(N_ref):
		F1[i] = [0] * N_site
		F2[i] = [0] * N_site


	#####============================== sampling S[0] =============================
	#####======================== forward marginalization =========================
	for i in range(N_site):
		if i == 0:  # no content
			continue

		if i == 1:  # V1(s2)
			for j in range(N_ref):
				## fill F1[j][i] for each j
				sum = 0
				for k in range(N_ref):
					## calculate exp( \xi() + \tau() + \xi() )
					temp = 0

					## \xi
					if H[0][i - 1] == 1:  # favor S = 1, use emission[0]
						temp += emission[0][k][i - 1]
					else:
						temp += emission[1][k][i - 1]

					## \tau
					r = (position[i] - position[i - 1])
					if k == j:
						temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
					else:
						temp += (1 - math.exp( - tune1 * r )) / tune2

					## \xi
					if H[0][i] == 1:  # favor S = 1, use emission[0]
						temp += emission[0][j][i]
					else:
						temp += emission[1][j][i]

					sum += math.exp(temp)
				F1[j][i] = sum
			continue

		### otherwise i >= 2
		for j in range(N_ref):
			## fill F1[j][i] for each j
			sum = 0
			for k in range(N_ref):
				## calculate V * exp( \tau() + \xi() )
				temp = 0

				## \tau
				r = (position[i] - position[i - 1])
				if k == j:
					temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
				else:
					temp += (1 - math.exp( - tune1 * r )) / tune2

				## \xi
				if H[0][i] == 1:  # favor S = 1, use emission[0]
					temp += emission[0][j][i]
				else:
					temp += emission[1][j][i]

				sum += F1[k][i - 1] * math.exp(temp)
			F1[j][i] = sum

	## now we have F1 ready (there are no contents in the first column)

	#####======================== backward sampling =========================
	sum = 0
	for i in range(N_ref):
		sum += F1[i][N_site - 1]
	for i in range(N_ref):
		F2[i][N_site - 1] = F1[i][N_site - 1] / sum

	## sample P(s_L)
	p = []
	for i in range(N_ref):
		p.append(F2[i][N_site - 1])
	l = np.random.multinomial(1, p, size=1)[0]
	for i in range(N_ref):
		if l[i] == 1:
			S[0][N_site - 1] = i
			break

	for i in reversed(range(N_site)):
		if i == (N_site - 1) or i == 0:
			continue
		else:
			for j in range(N_ref):
				tau = 0
				## \tau
				r = (position[i + 1] - position[i])
				if j == S[0][i + 1]:
					tau = (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
				else:
					tau = (1 - math.exp( - tune1 * r )) / tune2

				F2[j][i] = F1[j][i] * math.exp(  tau  )
			sum = 0
			for j in range(N_ref):
				sum += F2[j][i]
			for j in range(N_ref):
				F2[j][i] = F2[j][i] / sum

			## sample P(s_i | s_{i+1}), i = 2:(L-1)
			p = []
			for j in range(N_ref):
				p.append(F2[j][i])
			l = np.random.multinomial(1, p, size=1)[0]
			for j in range(N_ref):
				if l[j] == 1:
					S[0][i] = j
					break

	for j in range(N_ref):
		temp = 0

		## \xi
		if H[0][0] == 1:  # favor S = 1, use emission[0]
			temp += emission[0][j][0]
		else:
			temp += emission[1][j][0]

		## \tau
		r = (position[1] - position[0])
		if j == S[0][1]:
			temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
		else:
			temp += (1 - math.exp( - tune1 * r )) / tune2

		F2[j][0] = math.exp(temp)

	sum = 0
	for j in range(N_ref):
		sum += F2[j][0]
	for j in range(N_ref):
		F2[j][0] = F2[j][0] / sum

	## sample P(s_1 | s_{2})
	p = []
	for j in range(N_ref):
		p.append(F2[j][0])
	l = np.random.multinomial(1, p, size=1)[0]
	for j in range(N_ref):
		if l[j] == 1:
			S[0][0] = j
			break
	#####============================== sampling S[0] done =============================





	#####============================== sampling S[1] =================================
	#####======================== forward marginalization =============================
	for i in range(N_site):
		if i == 0:  # no content
			continue

		if i == 1:  # V1(s2)
			for j in range(N_ref):
				## fill F1[j][i] for each j
				sum = 0
				for k in range(N_ref):
					## calculate exp( \xi() + \tau() + \xi() )
					temp = 0

					## \xi
					if H[1][i - 1] == 1:  # favor S = 1, use emission[0]
						temp += emission[0][k][i - 1]
					else:
						temp += emission[1][k][i - 1]

					## \tau
					r = (position[i] - position[i - 1])
					if k == j:
						temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
					else:
						temp += (1 - math.exp( - tune1 * r )) / tune2

					## \xi
					if H[1][i] == 1:  # favor S = 1, use emission[0]
						temp += emission[0][j][i]
					else:
						temp += emission[1][j][i]

					sum += math.exp(temp)
				F1[j][i] = sum
			continue

		### otherwise i >= 2
		for j in range(N_ref):
			## fill F1[j][i] for each j
			sum = 0
			for k in range(N_ref):
				## calculate V * exp( \tau() + \xi() )
				temp = 0

				## \tau
				r = (position[i] - position[i - 1])
				if k == j:
					temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
				else:
					temp += (1 - math.exp( - tune1 * r )) / tune2

				## \xi
				if H[1][i] == 1:  # favor S = 1, use emission[0]
					temp += emission[0][j][i]
				else:
					temp += emission[1][j][i]

				sum += F1[k][i - 1] * math.exp(temp)
			F1[j][i] = sum

	## now we have F1 ready (there are no contents in the first column)

	#####======================== backward sampling =========================
	sum = 0
	for i in range(N_ref):
		sum += F1[i][N_site - 1]
	for i in range(N_ref):
		F2[i][N_site - 1] = F1[i][N_site - 1] / sum

	## sample P(s_L)
	p = []
	for i in range(N_ref):
		p.append(F2[i][N_site - 1])
	l = np.random.multinomial(1, p, size=1)[0]
	for i in range(N_ref):
		if l[i] == 1:
			S[1][N_site - 1] = i
			break

	for i in reversed(range(N_site)):
		if i == (N_site - 1) or i == 0:
			continue
		else:
			for j in range(N_ref):
				tau = 0
				## \tau
				r = (position[i + 1] - position[i])
				if j == S[1][i + 1]:
					tau = (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
				else:
					tau = (1 - math.exp( - tune1 * r )) / tune2

				F2[j][i] = F1[j][i] * math.exp(  tau  )
			sum = 0
			for j in range(N_ref):
				sum += F2[j][i]
			for j in range(N_ref):
				F2[j][i] = F2[j][i] / sum

			## sample P(s_i | s_{i+1}), i = 2:(L-1)
			p = []
			for j in range(N_ref):
				p.append(F2[j][i])
			l = np.random.multinomial(1, p, size=1)[0]
			for j in range(N_ref):
				if l[j] == 1:
					S[1][i] = j
					break

	for j in range(N_ref):
		temp = 0

		## \xi
		if H[1][0] == 1:  # favor S = 1, use emission[0]
			temp += emission[0][j][0]
		else:
			temp += emission[1][j][0]

		## \tau
		r = (position[1] - position[0])
		if j == S[1][1]:
			temp += (math.exp( - tune1 * r ) + (1 - math.exp( - tune1 * r ))) / tune2
		else:
			temp += (1 - math.exp( - tune1 * r )) / tune2

		F2[j][0] = math.exp(temp)

	sum = 0
	for j in range(N_ref):
		sum += F2[j][0]
	for j in range(N_ref):
		F2[j][0] = F2[j][0] / sum

	## sample P(s_1 | s_{2})
	p = []
	for j in range(N_ref):
		p.append(F2[j][0])
	l = np.random.multinomial(1, p, size=1)[0]
	for i in range(N_ref):
		if l[i] == 1:
			S[1][0] = i
			break
	#####============================== sampling S[1] done ============================

	return



def sampler_H():  # sampling H conditional on R and S
	global N_site
	global reads_reshaped
	global R
	global reference
	global S
	global H
	global epsilon
	global omega

	for i in range(N_site):
		sum1 = 0 # favor h1 = 1
		sum2 = 0 # favor h1 = -1
		sum3 = 0 # favor h2 = 1
		sum4 = 0 # favor h2 = -1

		for read in reads_reshaped[i]:  # all the reads under a specific site
			allele = read[1]
			if allele == 1:
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
		H[1][i] = h2
		#H[1][i] = -h1  # they should be composite; NOTE: should they? --> I don't think so. This will break the whole inference

	return





if __name__ == '__main__':


	## manually set the index first
	index = 0


	###================================ read the file to generate the following: initialization ===============================
	###========================================================================================================================
	## physical positions
	position = []  # the physical positions for the inferred haplotypes (only heterozygous sites in)
	file = open("legend" + str(index), 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		pos = int((line.strip()).split("\t")[1])
		position.append(pos)
	position = position[0: N_test]
	N_site = len(position)


	## haplotype
	# real one
	file = open("phase" + str(index), 'r')
	lines = file.readlines()
	file.close()
	h1 = (map(lambda x: int(x) * 2 - 1, (lines[2 * index].strip()).split(" ") ))[0: N_test]
	h2 = (map(lambda x: int(x) * 2 - 1, (lines[2 * index + 1].strip()).split(" ") ))[0: N_test]
	TRUE = [h1, h2]  # two haplotype (only heterozygous sites in)
	del lines[2 * index : 2 * index + 2]  # the left are the reference
	# random generated one
	h1 = np.random.binomial(1, 0.5, N_site) * 2 - 1
	h2 = - h1
	H = [h1, h2]


	## reference; removed all homozygous sites according to the above testing sample already
	## should be lines
	reference = []
	for line in lines:
		ref = (map(lambda x: int(x) * 2 - 1, (line.strip()).split(" ") ))[0: N_test]
		reference.append(ref)
	reference = reference[0: N_test2]
	N_ref = len(reference)
	s1 = [0] * N_site  # without loss of generality, we can choose 0 as the start reference panel
	s2 = [0] * N_site
	S = [s1, s2]


	# also, once we have the reference read into memory, we can set the emission matrix
	temp1 = []
	temp2 = []
	for i in range(N_ref):
		ref1 = reference[i][:]
		ref2 = reference[i][:]
		for j in range(N_site):
			ref2[j] = -ref2[j]
		temp1.append(ref1)
		temp2.append(ref2)
	emission = []
	emission.append(temp1)
	emission.append(temp2)
	for i in range(2):
		for j in range(N_ref):
			for k in range(N_site):
				if emission[i][j][k] == 1:
					emission[i][j][k] = 1 - omega
				else:
					emission[i][j][k] = omega


	## reads, and reshaped reads
	file = open("read" + str(index), 'r')
	lines = file.readlines()
	file.close()
	N_reads = len(lines)
	reads = [] # each list is one short read; each tuple is a specific site, containing index and allele value
	for line in lines:
		line = map(lambda x: int(x), (line.strip()).split(" ") )
		read = []
		i = 0
		while i < len(line):
			index = line[i]
			allele = line[i + 1] * 2 - 1
			read.append((index, allele))
			i += 2
		reads.append(read)
	R = np.random.binomial(1, 0.5, N_reads) * 2 - 1  ## NOTE: if r = 1, favor the first haplotype; otherwise the second

	reads_reshaped = {}  # each hashing value is one site; each tuple has index of read and corresponding allele value
	for i in range(N_site):
		reads_reshaped[i] = []
	for i in range(N_reads):
		read = reads[i]
		for site in read:  # site: (index, allele)
			index = site[0]
			allele = site[1]
			reads_reshaped[index].append((i, allele))
	###=========================================================================================================================
	###=========================================================================================================================




	for round in range(N_round):  ## or until equilibrium
		print "sampling round#",
		print round,
		print ":"


		###======================== sample the R with fixed H (S now is conditionally independent)
		print "sampling R..."
		sampler_R()


		###======================== sample the S with fixed H (R now is conditionally independent)
		print "sampling S..."
		sampler_S()


		###======================== sample the H with fixed R and S (NOTE: to be checked)
		print "sampling H..."
		sampler_H()


		###======================== calculate (or summarize) the present likelihood value; and phasing errors
		## likelihood; NOTE: maybe every 10 rounds?
		if round % 10 == 0:
			print "present likelihood value is",
			print like()

		## error rate
		print "present error rate is:",
		print error()


	print "sampling done..."
	print "the final error rate is:",
	print error()
