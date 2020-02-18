import numpy as np
import matplotlib.pyplot as mp
import statistics

# np.vstack((arr_a, arr_b)) --> stacks two arrays vertically
# np.hstack((arr_a, arr_b)) --> stacks two array inputs horizontally
# np.delete(arr, row_index, axis=0) --> deletes an entire row in array

class dll_object():
    def __init__(self, data):
        self.x = data[0]
        self.y = data[1]
        self.z = data[2]
        self.before = None
        self.after = None

    def get_data(self):
        return (self.x, self.y, self.z)

    def set_next(self, point):
        self.after = point

    def get_next(self):
        return self.after

    def set_before(self, point):
        self.before = point

    def get_before(self):
        return self.before


def clean_array(in_file):
    # map(function, iterable input) produces a list of float values for each line in *.crd
    # returns an outer list (for the entire file) containing inner lists (for each line) made of float values
    # np.array() turns an input list of lists into a 1 column array

    f = open(in_file, "r+")      # open *.crd file as a file_object
    n_dim = np.array([list(map(float, line.split())) for line in f.readlines()])
    f.close()
    return n_dim


def clean_list(in_file):
    # creates a doubly linked list with each node being of Class "dll_object"
    # returns the head of the list aka one of the termini of the protein

    f = open(in_file, "r+")
    node, node_before, list_head = None, None, None
    for line in f.readlines():
        data = list(map(float, line.split()))
        if list_head is None:
            list_head = dll_object(data)
            node = list_head
        else:
            node_before, node = node, dll_object(data)
            node.set_before(node_before)
            node_before.set_next(node)
    return list_head


def point_dist(pointA, pointB):
    # takes in two points and calculates the 3 dimensional distance between them

    a = np.array(pointA)
    b = np.array(pointB)

    dist = np.linalg.norm(b - a)
    return dist


def tube_eval(pointA, pointC, pointB):
    x = pointA - pointC
    t = (np.dot(pointB-pointC, x)/np.dot(x, x))
    distance = np.linalg.norm(t*(pointA-pointC)+(pointC-pointB))

    if distance <= 0.5:
        return True
    else:
        return distance


def visualize():
    pass
#     zdata = arr[:, [2]].flatten()
#     xdata = arr[:, [0]].flatten()
#     ydata = arr[:, [1]].flatten()
#     ax.scatter3D(xdata, ydata, zdata, cmap='Greens')
#     ax.plot(xdata, ydata, zdata, color='r')


def max_dist(array):
    list = []
    for i in range(2,len(array)):
        list.append(point_dist(array[i], array[i-1]))
    return (max(list), statistics.mean(list), min(list))


def main():

    file = "/Users/jormungandr/Downloads/knot1.crd"
    pos_arr = clean_array(file)   # produces an array of position values based on output of clean_array
    pos_list = None


if __name__ == "__main__":
    main()
