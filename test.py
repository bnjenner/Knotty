# import numpy as np

# def triNormal(p1, p2, p3):

# 	return np.cross((p2 - p1), (p3 - p1))


# def triSideNormals(p1, p2, p3, N):

# 	N_12 = np.cross((p2 - p1), N) 
# 	N_23 = np.cross((p3 - p2), N)
# 	N_31 = np.cross((p1 - p3), N)

# 	return [N_12, N_23, N_31]


# def plane_dist(p1, N):

# 	return ( -1 * np.dot(N, p1) )


# def distance(xp, tp, N_side):

# 	return ( np.dot((xp - tp), N_side) / np.linalg.norm(N_side) )




# p1 = [1, 0, 0]
# p2 = [0, 1, 0]
# p3 = [0, 0, 1]

# s1 = [2, 2, 2]
# s2 = [1, 1, 1]

# p1 = np.array(p1)
# p2 = np.array(p2)
# p3 = np.array(p3)

# points = [p1, p2, p3]

# print(list(map(lambda x : x/3, np.sum(points, axis=0))))

# s1 = np.array(s1)
# s2 = np.array(s2)

# points = [p1, p2, p3]

# N = triNormal(p1, p2, p3)

# D = plane_dist(p1, N)

# Unit_N = N / np.linalg.norm(N)

# t = ( np.dot((D + Unit_N), s1) ) / ( np.dot(Unit_N, (s2 - s1)) )

# if t < 0 or t > 1:
# 	print("No Intersection")
# 	exit()

# R = s1 + t * (s2 - s1)

# Tri_Side_Normals = triSideNormals(p1, p2, p3, N)

# for i in range(len(Tri_Side_Normals)):
# 	p = points[i]
# 	n = Tri_Side_Normals[i]
# 	dist = distance(R, p, n)

# 	if dist > 0:
# 		print("No Intersection")
# 	elif dist < 0 and i == 2:
# 		print("Intersection")

# def func(a, b, c, d, e):
# 	print(a, b, c, d, e)

# sequence = "ABCDEFGHIJK"
# for i in range(1, len(sequence) - 3):
# 	seq_list = sequence[i-1:i+2]
# 	for j in range(i + 2, len(sequence) - 1):
# 		seq = seq_list + sequence[j:j+2]
# 		seq = tuple(seq_list + sequence[j:j+2])
# 		func(*seq)

def linear(points):

	for x in range(len(points)-1):
		print(points[x:x+2])


linear("ABCDEFGHIJKLMN")



