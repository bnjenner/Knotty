######################
#
# Knotty 
#
######################
import argparse
import numpy as np
import re
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import time

class dll_object():

	def __init__(self, data):
		self.coord = data
		self.before = None
		self.after = None

	def get_data(self):
		return self.coord

	def set_next(self, point):
		self.after = point

	def get_next(self):
		return self.after

	def set_before(self, point):
		self.before = point

	def get_before(self):
		return self.before


class operations():

	def visualize(head):

		X = []
		Y = []
		Z = [[], []]

		curr_node = head

		while curr_node != None:
			point = curr_node.get_data()
			X.append(point[0])
			Y.append(point[1])
			Z[0].append(point[2])
			Z[1].append(point[2])

			curr_node = curr_node.get_next()

		X = np.array(X)
		Y = np.array(Y)
		Z = np.array([np.array(Z[0]), np.array(Z[1])])

		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		ax.plot_wireframe(X, Y, Z)

		plt.show() 

	# Intersection Function based on the Non-Culling Moller Trumbore Algorithm
	def intersection(A, B, C, O, E):


		E1 = B - A # edge 1
		E2 = C - A # edge 2

		D = (E - O)

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

		# scalar must be no bigger than one, cannot extend past segment
		if t < 0 or t > 1:
			return 0 
		else:
			return 1


	def tube_eval(A, B, C):

		E1 = B - C
		E2 = A - C


		proj = ( np.dot(E1, E2) / ( np.linalg.norm(E2) * np.linalg.norm(E2)) ) * E2

		d = np.linalg.norm(( E1 - proj )) 

		if d <= 0.5:
			return 1
		else:
		 	return 0


	def forward_search(llist):

		pass

	def backward_search(llist):

		pass 


class knotFinder():

	def __init__(self):
		pass

	def trefoil(p1, p2, p3, s1, s2):

		return operations.intersection(p1, p2, p3, s1, s2)


	def scan(head):

		prev_node = head
		curr_node = head.get_next()
		next_node = curr_node.get_next()

		smoothed = True

		while True:

			if next_node == None :
				break

			colinear_check = operations.tube_eval(prev_node.get_data(),
												  curr_node.get_data(), 
												  next_node.get_data())
			
			if colinear_check == 1:
				prev_node.set_next(next_node) # node deletion
				next_node.set_before(prev_node)
				prev_node = next_node
				curr_node = next_node.get_next()

				if curr_node != None:
					next_node = curr_node.get_next()
				else:
					next_node = None

				continue

			smoothed = False

			try:
				seq_list = [prev_node.get_data(), curr_node.get_data(), next_node.get_data()] # list of i-1, i, and i+1
			except:
				print(colinear_check)
				print(prev_node, curr_node, next_node)
				exit()

			i_prime = np.array(list(map(lambda x : x/3, np.sum(seq_list, axis=0)))) # i' calculation 
			prev_i = prev_node.get_data() # previous i 

			prev_tri = [prev_i, seq_list[1], i_prime] # first triangle
			next_tri = [seq_list[1], i_prime, seq_list[2]] # second triangle

			knot_presence = False

			line_start = next_node.get_next() 
			if line_start != None:
				line_end = line_start.get_next()

				while line_end != None:

					line_seg = [line_start.get_data(), line_end.get_data()]

					prev_list = tuple(prev_tri + line_seg)
					next_list = tuple(next_tri + line_seg)
					prev_knot = knotFinder.trefoil(*prev_list)
					next_knot = knotFinder.trefoil(*next_list)

					# potential point for not knots
					if prev_knot == 1: 
						knot_presence = True
						break


					elif next_knot == 1:
						knot_presence = True
						break

					line_start = line_end
					line_end = line_end.get_next()


			if prev_node != head:
				line_start = prev_node.get_before()
				line_end = line_start.get_before()

				while line_end != None and knot_presence == False:

					line_seg = [line_start.get_data(), line_end.get_data()]

					prev_list = tuple(prev_tri + line_seg)
					next_list = tuple(next_tri + line_seg)
					prev_knot = knotFinder.trefoil(*prev_list)
					next_knot = knotFinder.trefoil(*next_list)

					# potential point for not knots
					if prev_knot == 1: 
						knot_presence = True
						break


					elif next_knot == 1:
						knot_presence = True
						break

					if line_end.get_before() != None:
						line_end = line_end.get_before()
						line_start = line_end.get_next()
					else:
						break		


			if knot_presence == True:	
				prev_node = curr_node
				curr_node = prev_node.get_next()
				next_node = curr_node.get_next()
			else:
				new_node = dll_object(i_prime)
				prev_node.set_next(new_node) 
				new_node.set_next(next_node) 
				new_node.set_before(prev_node)
				next_node.set_before(new_node)
				prev_node = new_node
				curr_node = new_node.get_next()
				next_node = curr_node.get_next()


		if smoothed == True:
			return 0
		else:
			return 1


	def scan_aux(head):

		iterations = 1

		while True:

			if (iterations % 10)  == 0:
				print("*** Iterations: ", (iterations)," ***", sep="")

			knots = knotFinder.scan(head)

			iterations += 1

			if knots == 0:
				print("--- %s seconds ---" % (time.time() - start_time))
				operations.visualize(head)
				return "No Knots"

			if iterations > 500:
				print("--- %s seconds ---" % (time.time() - start_time))
				operations.visualize(head)
				return "Knots"
			

class protein():

	def __init__(self, InputFile):

    # creates a doubly linked list with each node being of Class "dll_object"
    # returns the head of the list aka one of the termini of the protein

		with open(InputFile, "r") as fi:
				lines = fi.readlines()

		node, node_before, head = None, None, None

		for line in lines:
			data = np.array(list(map(float, line.split())))

			if head is None:
				head = dll_object(data)
				node = head

			else:
				node_before, node = node, dll_object(data)
				node.set_before(node_before)
				node_before.set_next(node)

		self.backbone = head





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

	result = knotFinder.scan_aux(aa_chain)

	# with open(OutputFile, "w") as fo:
	# 	fo.write(result)


if __name__ == '__main__':
	start_time = time.time()
	main()
