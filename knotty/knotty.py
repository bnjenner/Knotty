######################
#
# Knotty 
#
######################
import argparse
import numpy as np
import re

class operations():

	def triNormal(p1, p2, p3):

		return np.cross((p2 - p1), (p3 - p1))


	def triSideNormals(p1, p2, p3, N):

		N_12 = np.cross((p2 - p1), N) 
		N_23 = np.cross((p3 - p2), N)
		N_31 = np.cross((p1 - p3), N)

		return [N_12, N_23, N_31]


	def plane_dist(p1, N):

		return ( -1 * np.dot(N, p1) )


	def tri_distance(xp, tp, N_side):

		return ( np.dot((xp - tp), N_side) / np.linalg.norm(N_side) )

	def tri_intersection(ex_point, points, side_normals):

		R = ex_point

		for i in range(len(side_normals)):
			p = points[i]
			n = side_normals[i]
			dist = operations.tri_distance(R, p, n)

			if dist > 0:
				return 0
			elif dist < 0 and i == 2:
				return 1

	def linear_distance(p1, p2):

		return np.sqrt(np.dot(p1, p2))

	def linear(points, length):

		total_length = 0 

		for x in range(len(points)-1):
			inner_list = tuple(points[x:x+2])
			inner_dist = operations.linear_distance(*inner_list)
			total_length += inner_dist

			if total_length > length:
				return 0 
			else:
				pass

		if length == total_length:
			return 1 


class knotFinder():

	def __init__(self):
		pass

	def trefoil(p1, p2, p3, s1, s2):


		N = operations.triNormal(p1, p2, p3)

		points = [p1, p2, p3]

		D = operations.plane_dist(p1, N)

		Unit_N = N / np.linalg.norm(N)

		t = ( np.dot((D + Unit_N), s1) ) / ( np.dot(Unit_N, (s2 - s1)) )

		if t < 0 or t > 1:
			return 0 

		R = s1 + t * (s2 - s1)

		Tri_Side_Normals = operations.triSideNormals(p1, p2, p3, N)

		intersect = operations.tri_intersection(R, points, Tri_Side_Normals)

		return intersect



	def figeight(self, protein):

		pass

	def quinqefoil(self, protein):
		
		pass 

	def pentaknot(self, protein):
		
		pass 

	def scan(sequence):

		iterations = 1

		straight_length = operations.linear_distance(sequence[1], sequence[-1])

		while iterations < 501:

			if ( ( (iterations / 500) * 100) % 10)  == 0:
				print("*** ", (iterations / 500 * 100),"% Complete ***", sep="")

			if operations.linear(sequence, straight_length) == 1:
				return "No Knots"

			for i in range(1, len(sequence) - 3):
				seq_list = sequence[i-1:i+2]
				aa_prime = np.array(list(map(lambda x : x/3, np.sum(seq_list, axis=0))))

				prev_list = [sequence[i-1], sequence[i], aa_prime]
				next_list = [sequence[i], aa_prime, sequence[i+1]]

				knot_presence = False

				for j in range(i + 2, len(sequence) - 1):
					prev_tri = tuple(prev_list + sequence[j-1:j+1])
					next_tri = tuple(next_list + sequence[j:j+2])
					prev_knot = knotFinder.trefoil(*prev_tri)
					next_knot = knotFinder.trefoil(*next_tri)

					if prev_knot == 1 or next_knot == 1:
						knot_presence = True
						break

				if knot_presence == False:
					sequence[i] = aa_prime
				else:
					return "Knots"

			iterations += 1

		if operations.linear(sequence, straight_length) == 1:
				return "No Knots"
		else:
			return "Knots"



class protein():

	def __init__(self, InputFile):

		self.backbone = []

		with open(InputFile, "r") as fi:
			lines = fi.readlines()

		temp = list(map(lambda x : re.sub("\s+", ",", x.strip()), lines))
		temp_lst = list(map(lambda x : x.split(','), temp))

		for point in temp_lst:
			temp = np.array([float(x) for x in point])
			self.backbone.append(temp)

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
