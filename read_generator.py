import os
import sys



if __name__ == '__main__':

	for index in range(60):



		print "now generating read for index number",
		print index


		file = open("phase" + str(index), 'r')
		lines = file.readlines()
		file.close()
		for i in range(len(lines)):
			lines[i] = map(lambda x: int(x), (lines[i].strip()).split(" "))

		test1 = lines[2 * index]
		test2 = lines[2 * index + 1]
		del lines[2 * index : 2 * index + 2]
		reference = lines

		print len(reference)
		

		###===================== generating the reads =======================
		print "length of tested haplotype1:",
		print len(test1)
		print "length of tested haplotype2:",
		print len(test2)

		## 5/4 units length reads every 3 position for both h1 and h2
		file = open("read" + str(index), 'w')

		## h1: 5
		i = 0
		while i < len(test1):
			## a new read
			j = i
			while j < len(test1) and (j - i) <= 4:
				file.write(str(j) + " " + str(test1[j]) + " ")
				j += 1
			file.write("\n")
			i += 3
		## h2: 5
		i = 0
		while i < len(test2):
			## a new read
			j = i
			while j < len(test2) and (j - i) <= 4:
				file.write(str(j) + " " + str(test2[j]) + " ")
				j += 1
			file.write("\n")
			i += 3
		## h1: 4
		i = 0
		while i < len(test1):
			## a new read
			j = i
			while j < len(test1) and (j - i) <= 3:
				file.write(str(j) + " " + str(test1[j]) + " ")
				j += 1
			file.write("\n")
			i += 3
		## h2: 4
		i = 0
		while i < len(test2):
			## a new read
			j = i
			while j < len(test2) and (j - i) <= 3:
				file.write(str(j) + " " + str(test2[j]) + " ")
				j += 1
			file.write("\n")
			i += 3

		print "now read generation finished for read index number",
		print index


		file.close()
