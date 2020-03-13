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

#####################################################################################
# data structure classes

############################################
# doubly-linked list

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


############################################
# protein class

class protein():

	def __init__(self, InputFile):

    # creates a doubly linked list with each node being of Class "dll_object"
    # returns the head of the list aka one of the termini of the protein

		with open(InputFile, "r") as fi:
				lines = fi.readlines()

		node, node_before, head = None, None, None

		pos = 0 # index for position in polypeptide chain.

		for line in lines:
			data = np.array(list(map(float, line.split())))
			
			if len(data) != 0: #remove empty list entry from extra newline characters

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


#####################################################################################
# operations class 

class operations():


	############################################
	# user interface 

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

		X = np.array(X) # x coord for each node 
		Y = np.array(Y) # y coord for each node 
		Z = np.array([np.array(Z[0]), np.array(Z[1])]) # z coord for wire frame 
		Z_2 = np.array(np.array(Z[0])) # z coord for each node 
		fig = plt.figure()
		
		ax = fig.add_subplot(111, projection='3d')
		
		ax.scatter(X, Y, Z_2, c="red")
		ax.plot_wireframe(X, Y, Z)

		plt.show()


	def crd_write(head, OutputFile):

		# creates a crd file from protein / linked list object 

		with open(OutputFile, "w") as fo:

			current = head

			while current != None:
				temp = list(map(str, list(current.get_data())))	
				fo.write("\t".join(temp) + "\n")
				current = current.get_next()



	############################################
	# mathematical operations 

	def intersection(A, B, C, O, E): # Intersection Function based on the Non-Culling Moller Trumbore Algorithm


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

	def tube_eval(A, B, C): # Projection based collinearity approximation tool

		# transforma triangle so edges stem from origin 
		E1 = B - C 
		E2 = A - C

		# project vector E1 onto E2
		proj = ( np.dot(E1, E2) / ( np.linalg.norm(E2) * np.linalg.norm(E2)) ) * E1

		# height of vector perpendicular to projection
		d = np.linalg.norm(( E1 - proj )) 

		if d < Epsilon:
			return 1
		else:
		 	return 0


	############################################
	# knot finding traversals

	def forward_search(prev_tri, next_tri, next_node):
        

		line_start = next_node.get_next() 
		if line_start != None:
			line_end = line_start.get_next()

			while line_end != None:

				line_seg = [line_start.get_data(), line_end.get_data()]

				next_list = tuple(next_tri + line_seg)
				next_knot = knotFinder.trefoil(*next_list) # calls trefoil function which calls intersection function 

				if next_knot == 1:
					return 1

				prev_list = tuple(prev_tri + line_seg)
				prev_knot = knotFinder.trefoil(*prev_list) # calls trefoil function which calls intersection function 

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
			next_knot = knotFinder.trefoil(*next_list) # calls trefoil function which calls intersection function 

			if next_knot == 1:
				return 1

			prev_list = tuple(prev_tri + line_seg)
			prev_knot = knotFinder.trefoil(*prev_list) # calls trefoil function which calls intersection function 

			if prev_knot == 1:
				return 1

			line_start = line_end
			line_end = line_start.get_before()	

		return 0


#####################################################################################
# knot finder 

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

		# iterates over every node in linked list (minus termini)
		while next_node != None:

			# restricts application of colinear approximatory to
			# if current node is not involved in a currently knotted section 
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

					continue # restarts loop

			# section of code is only run if program fails colinearity check

			smoothed = False

			seq_list = [prev_node.get_data(), curr_node.get_data(), next_node.get_data()] # list of i-1, i, and i+1

			i_prime = np.array(list(map(lambda x : x/3, np.sum(seq_list, axis=0)))) # i' calculation 

			prev_tri = [seq_list[0], seq_list[1], i_prime] # first triangle
			next_tri = [seq_list[1], i_prime, seq_list[2]] # second triangle

			knot_presence = False
		
			# if backwards traversal is possible (n => 4)
			# checks backwards first cause traversal is usually shorter
			if prev_node != head:
				backward_result = operations.backward_search(prev_tri, next_tri, prev_node)
			else:
				backward_result = 0

			# if backwards result yielded no knot, check downstream of triangle
			if backward_result == 0:
				forward_result = operations.forward_search(prev_tri, next_tri, next_node)
			else:
				forward_result = 1

			# if knot is detected, parameters in nodes is marked as possibly being 
			# involved in the detection of a knot. 
			if forward_result == 1 or backward_result == 1:
				if curr_node.get_status() != True:
					prev_node.set_status(True)
					curr_node.set_status(True)
					next_node.set_status(True)

				knot_presence = True

			else:
				if curr_node.get_status() != True:
					prev_node.set_status(False)
					curr_node.set_status(False)
					next_node.set_status(False)
		

			# section for traversing the linked list. If knot was detected, next node is set to current and so on
			# if a knot is not detected, current node is replaced with i prime			
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

		# iteration accumulator variable
		iterations = 1

		while True:

			# call knot finder scanning function 
			knots = knotFinder.scan(head)

			iterations += 1

			if (iterations % 50) == 0:
				print("*** Iterations: ", (iterations)," ***", sep="")


			if iterations >= MaxIterations:

				knot_pos = []
				current = head

				# checks current linked list for any proteins potentially involved in knots.
				# returns their positions if true 
				while current != None:

					if current.get_status() == True:
						knot_pos.append(str(current.get_pos()))
					
					current = current.get_next()


				if len(knot_pos) == 0: # if no knots detected

					print("*** Max iterations reached and no knot detected. Visualize to confirm. ***")
					print("*** If complete smoothing is required, an increase in Epsilon or ***")
					print("*** 	the number of iterations is recommended. ***")

				else:
					knot_list = ", ".join(knot_pos)
					print("*** Max iterations reached and a knot is likely. ***")
					print("*** Knot possible at amino acids", knot_list, "***")

				print("*** Time Elapsed: %s seconds ***" % round((time.time() - start_time), 2))
				return

			 # if true collinearity was achieved.
			if knots == 0:
				print("*** Smoothing was able to be completed. No knot is likely present ***")
				print("*** Time Elapsed: %s seconds ***" % round((time.time() - start_time), 2))
				return


#####################################################################################
# app mains

############################################
# find 

class findApp():

	def __init__(self):
		self.verbose = False

	def start(self, InputFile, OutputFile, max_iterations, epsilon):

		global MaxIterations
		global Epsilon 

		# set parameters for find function
		MaxIterations = max_iterations
		Epsilon = epsilon

		# instantiate protein class and set pointer to AA chain (linked list).
		pro_sequence = protein(InputFile)
		aa_chain = pro_sequence.backbone

		print("*** Beginning Knot Detection... ***")
		print("*** Max Iterations:", MaxIterations, "Epsilon:", Epsilon, "***")

		# begin finding knots
		knotFinder.scan_aux(aa_chain)

		# write output of smoothed file
		operations.crd_write(aa_chain, OutputFile)


############################################
# visualize

class visualizeApp():

	def __init__(self):
		self.verbose = False

	def start(self, InputFile):

		# instantiate protein class and set pointer to AA chain (linked list).
		pro_sequence = protein(InputFile)
		aa_chain = pro_sequence.backbone

		# visualize chain
		operations.visualize(aa_chain)


#####################################################################################
# argument parser

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


############################################
# functions

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Knotty is a program for detecting and visualizing knots in protein structures.', add_help=True,
        epilog="Last Updated: March 13th, 2020.")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    findParser(subparsers)
    visualizeParser(subparsers)
    args = parser.parse_args()
    return args


#####################################################################################
# main

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
