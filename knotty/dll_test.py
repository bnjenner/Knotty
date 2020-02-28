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


head = dll_object(0)
node1 = dll_object(1)
node2 = dll_object(2)
node3 = dll_object(3)
node4 = dll_object(4)
node5 = dll_object(5)
node6 = dll_object(6)
node7 = dll_object(7)

head.set_next(node1)
node1.set_before(head)
node1.set_next(node2)
node2.set_before(node1)
node2.set_next(node3)
node3.set_before(node2)
node3.set_next(node4)
node4.set_before(node3)
node4.set_next(node5)
node5.set_before(node4)
node5.set_next(node6)
node6.set_before(node5)
node6.set_next(node7)
node7.set_before(node6)

prev_node = head
curr_node = head.get_next()
next_node = curr_node.get_next()

while next_node != None:

	line_start = next_node.get_next() 

	if line_start != None:
		line_end = line_start.get_next()

		while line_end != None:

			print(prev_node.get_data(), curr_node.get_data(), next_node.get_data(),
				  line_start.get_data(), line_end.get_data())

			line_start = line_end
			line_end = line_end.get_next()

	if prev_node != head:
		line_start = prev_node.get_before()
		line_end = line_start.get_before()

		while line_end != None:

			print(prev_node.get_data(), curr_node.get_data(), next_node.get_data(),
				  line_start.get_data(), line_end.get_data())

			if line_end.get_before() != None:
				line_end = line_end.get_before()
				line_start = line_end.get_next()
			else:
				break

	prev_node = curr_node
	curr_node = curr_node.get_next()
	next_node = curr_node.get_next()


