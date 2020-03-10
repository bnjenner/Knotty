#!/usr/bin/env python3

######################
#
# Knotty 
#
######################

import argparse
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import time

class dll_object():

	def __init__(self, data):
		self.coord = data
		self.before = None
		self.after = None
		self.knot = False
		self.pos = None

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

	def set_status(self, status):
		self.knot = status

	def get_status(self):
		return self.knot

	def set_pos(self, pos):
		self.pos = pos

	def get_pos(self):
		return self.pos


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

		proj = ( np.dot(E1, E2) / ( np.linalg.norm(E2) * np.linalg.norm(E2)) ) * E1

		d = np.linalg.norm(( E1 - proj )) 

		if d < Epsilon:
			return 1
		else:
		 	return 0

	def forward_search(prev_tri, next_tri, next_node):
        
		line_start = next_node.get_next() 
		if line_start != None:
			line_end = line_start.get_next()

			while line_end != None:

				line_seg = [line_start.get_data(), line_end.get_data()]

				next_list = tuple(next_tri + line_seg)
				next_knot = knotFinder.trefoil(*next_list)

				if next_knot == 1:
					return 1

				prev_list = tuple(prev_tri + line_seg)
				prev_knot = knotFinder.trefoil(*prev_list)

				if prev_knot == 1:
					return 1

				line_start = line_end
				line_end = line_end.get_next()

		return 0

	def backward_search(prev_tri, next_tri, prev_node):

		line_start = prev_node.get_before()
		line_end = line_start.get_before()
		
		while line_end != None:

			line_seg = [line_start.get_data(), line_end.get_data()]


			next_list = tuple(next_tri + line_seg)
			next_knot = knotFinder.trefoil(*next_list)

			if next_knot == 1:
				return 1

			prev_list = tuple(prev_tri + line_seg)
			prev_knot = knotFinder.trefoil(*prev_list)

			if prev_knot == 1:
				return 1

			line_start = line_end
			line_end = line_start.get_before()	

		return 0


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

		while next_node != None:

			if curr_node.get_status() != True:
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

			seq_list = [prev_node.get_data(), curr_node.get_data(), next_node.get_data()] # list of i-1, i, and i+1

			i_prime = np.array(list(map(lambda x : x/3, np.sum(seq_list, axis=0)))) # i' calculation 

			prev_tri = [seq_list[0], seq_list[1], i_prime] # first triangle
			next_tri = [seq_list[1], i_prime, seq_list[2]] # second triangle

			knot_presence = False
		
			if prev_node != head:
				backward_result = operations.backward_search(prev_tri, next_tri, prev_node)
			else:
				backward_result = 0

			if backward_result == 0:
				forward_result = operations.forward_search(prev_tri, next_tri, next_node)
			else:
				forward_result = 1

			if forward_result == 1 or backward_result == 1:
				if curr_node.get_status() != True:
					prev_node.set_status(True)
					curr_node.set_status(True)
					next_node.set_status(True)

				knot_presence = True
						
			if knot_presence == True:
				prev_node = curr_node
				curr_node = prev_node.get_next()
				next_node = curr_node.get_next()
			else:
				new_node = dll_object(i_prime)
				new_node.set_pos(curr_node.get_pos())
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

			knots = knotFinder.scan(head)

			iterations += 1

			if (iterations % 50) == 0:
				print("*** Iterations: ", (iterations)," ***", sep="")


			if iterations >= MaxIterations:

				knot_pos = []
				current = head

				while current != None:

					if current.get_status() == True:
						knot_pos.append(str(current.get_pos()))
					
					current = current.get_next()


				if len(knot_pos) == 0:

					print("*** Max iterations reached, but no knot could be detected. ***")
					print("*** If complete smoothing is required, an increase in Epsilon or ***")
					print("*** 	the number of iterations is recommended. ***")

				else:
					knot_list = ", ".join(knot_pos)
					print("*** Max iterations reached and a knot is likely. ***")
					print("*** Knot possible at amino acids", knot_list, "***")

				print("*** Time Elapsed: %s seconds ***" % round((time.time() - start_time), 2))
				return

			if knots == 0:
				print("*** Smoothing was able to be completed. No knot is likely present ***")
				print("*** Time Elapsed: %s seconds ***" % round((time.time() - start_time), 2))
				return

			

class protein():

	def __init__(self, InputFile):

    # creates a doubly linked list with each node being of Class "dll_object"
    # returns the head of the list aka one of the termini of the protein

		with open(InputFile, "r") as fi:
				lines = fi.readlines()

		node, node_before, head = None, None, None

		pos = 0

		for line in lines:
			data = np.array(list(map(float, line.split())))

			if head is None:
				head = dll_object(data)
				head.set_pos(pos)
				node = head
				pos += 1
				

			else:
				node_before, node = node, dll_object(data)
				node.set_before(node_before)
				node.set_pos(pos)
				node_before.set_next(node)
				pos += 1

		self.backbone = head


	def crd_write(self, OutputFile):

		with open(OutputFile, "w") as fo:

			current = self.backbone

			while current != None:
				temp = list(map(str, list(current.get_data())))	
				fo.write("\t".join(temp) + "\n")
				current = current.get_next()

############################################

class findApp():

	def __init__(self):
		self.verbose = False

	def start(self, InputFile, OutputFile, max_iterations, epsilon):

		global MaxIterations
		global Epsilon 

		MaxIterations = max_iterations
		Epsilon = epsilon


		pro_sequence = protein(InputFile)
		aa_chain = pro_sequence.backbone

		print("*** Beginning Knot Detection... ***")
		print("*** Max Iterations:", MaxIterations, "Epsilon:", Epsilon, "***")

		knotFinder.scan_aux(aa_chain)

		pro_sequence.crd_write(OutputFile)


class visualizeApp():

	def __init__(self):
		self.verbose = False

	def start(self, InputFile):

		pro_sequence = protein(InputFile)
		aa_chain = pro_sequence.backbone

		operations.visualize(aa_chain)


############################################
# find 

def findParser(subparsers):
	find_parser = subparsers.add_parser('find', help='Protein knot detection tool')
	find_parser.add_argument('-i','--input', help='file of amino acid sequence coordinates', dest="InputFile", required=True)
	find_parser.add_argument('-o','--output', help='name of output file. If not specified, default is given with overwrite capabilities', dest="OutputFile", required=True)
	find_parser.add_argument('-m','--max-iterations', help='maximum number of iterations for smoothing algorithm (default: 250).', dest="max_iterations", type=int, default=250)
	find_parser.add_argument('-e','--epsilon', help='threshold (in Angstroms) for approximating colinearity of alpha carbons (default: 0.25).', dest="epsilon", type=float, default=0.25)
	
	return find_parser

class findCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = findApp()
   		return app.start(args.InputFile, args.OutputFile, args.max_iterations, args.epsilon) 

############################################
# visualize

def visualizeParser(subparsers):
	visualize_parser = subparsers.add_parser('visualize', help='Interactive 3D protein coordinate visualization tool')
	visualize_parser.add_argument('-i','--input', help='file of amino acid sequence coordinates', dest="InputFile", required=True)
	
	return visualize_parser

class visualizeCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = visualizeApp()
   		return app.start(args.InputFile) 


#####################################################################################

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Knotty is a program for detecting and visualizing knots in protein structures.', add_help=True,
        epilog="For questions or comments, contact somebody else.")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    findParser(subparsers)
    visualizeParser(subparsers)
    args = parser.parse_args()
    return args


def main():

	find = findCMD()
	visualize = visualizeCMD()

	global MaxIterations 
	global Epsilon 

	commands = {'find':find, 'visualize':visualize}
	args = parseArgs()
	commands[args.command].execute(args)


if __name__ == '__main__':
	start_time = time.time()
	main()
