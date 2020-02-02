######################
#
# Knotty 
#
######################
import argparse
import numpy as np
import re

class operations():

	# Intersection Function based on the Non-Culling Moller Trumbore Algorithm
	def intersection(A, B, C, O, E):

		E1 = B - A # edge 1
		E2 = C - A # edge 2

		D = (E - O) - A

		P = np.cross(D, E2) # P vector for determinant 
		det = np.dot(P, E1) # determinant

		if det < 0.000001 and det > -0.000001: # determines if system has solution
			return 0 

		inv_det = 1 / det

		# vector representing distance from s2 to p1
		T = O - A

		###
		# In order for point to be inside triangle, u and v must
		# be positive and their sum must be less than one.
		### 

		# solving system for u parameter
		u = np.dot(T, P) * inv_det

		# check if u is positive
		if u < 0 or u > 1:
		 	return 0 

		# solving system for v parameter
		Q = np.cross(T, E1)
		v = np.dot(D, Q) * inv_det

		# check if v is positive that sum of u and v is no greater than 1
		if v < 0 or (u + v) > 1:
			return 0 

		# solving for parameter t
		# Segment is repsented by S = O + tD
		# t is the scalar that multiples direction vector so that it intersects with triangle 
		t = np.dot(E2, Q) * inv_det # distance from origin to point A

		# scalar must be no bigger than one, cannot extend passed segment
		if t < 0 or t > 1:
			return 0 
		else:
			return 1


	def colinear(points):

		for x in range(2, len(points)):
			mat = np.array([points[x-2], points[x-1], points[x]])
			det = np.linalg.det(mat)

			if det * 0.5 != 0:
				return 0

		return 1 


class knotFinder():

	def __init__(self):
		pass

	def trefoil(p1, p2, p3, s1, s2):

		return operations.intersection(p1, p2, p3, s1, s2)


	def scan(sequence):

		iterations = 1

		while iterations <= 500:

			if (iterations % 10)  == 0:
				print("*** Iterations: ", (iterations)," ***", sep="")

			coline_check = operations.colinear(sequence)
			if coline_check == 1:
				print(coline_check)
				print("No Knots")
				return "No Knots"

			for i in range(1, len(sequence) - 1):

				seq_list = sequence[i-1:i+2]
				i_prime = np.array(list(map(lambda x : x/3, np.sum(seq_list, axis=0))))

				# terminal i is essentiall an i prime
				prev_i = sequence[i-1]		

				prev_tri = [prev_i, sequence[i], i_prime]
				next_tri = [sequence[i], i_prime, sequence[i+1]]

				knot_presence = False

				if i < len(sequence) - 3:
					for j in range(i + 2, len(sequence) - 1):

						line_seg = sequence[j:j+2]

						prev_list = tuple(prev_tri + line_seg)
						next_list = tuple(next_tri + line_seg)
						prev_knot = knotFinder.trefoil(*prev_list)
						next_knot = knotFinder.trefoil(*next_list)

						# potential point for not knots
						if prev_knot == 1 or next_knot == 1:
							print("Knot")
							knot_presence = True
							break
							#print("Knots")
							#return "Knots"


				if i > 2 and knot_presence != True: # if preceding segments are available
					for k in range(0, i-2): # iterates over previous line segments

						line_seg = sequence[k:k+2]

						prev_list = tuple(prev_tri + line_seg)
						next_list = tuple(next_tri + line_seg)
						prev_knot = knotFinder.trefoil(*prev_list)
						next_knot = knotFinder.trefoil(*next_list)

						if prev_knot == 1 or next_knot == 1:
							print("Knot")
							knot_presence = True
							break
				
				if knot_presence == False:
					sequence[i] = i_prime
				else:
					knot_presence = False

			iterations += 1


		if operations.colinear(sequence) == 1:
			print("No Knots")
			return "No Knots"
		else:
			print("Knots")
			return "Knots"



class protein():

	def __init__(self, InputFile):   

		self.backbone = []

		with open(InputFile, "r") as fi:
			lines = fi.readlines()

		self.backbone = np.array([list(map(float, line.split())) for line in lines])

		
############################################
# Argument Parser 


def main():
	parser = argparse.ArgumentParser(description='Protein Knot Detection')
	parser.add_argument('-i','--input', help='file of amino acid sequence coordinates', required=True)
	parser.add_argument('-o','--output', help='name of output file. If not specified, default is given with overwrite capabilities', required=True)
	args = vars(parser.parse_args())

	InputFile = args["input"]
	OutputFile = args["output"]

	pro_sequence = protein(InputFile)

	aa_chain = pro_sequence.backbone

	result = knotFinder.scan(sequence = aa_chain)

	with open(OutputFile, "w") as fo:
		fo.write(result)


if __name__ == '__main__':
	main()
