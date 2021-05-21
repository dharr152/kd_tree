# Update this code as needed
import math

class kDTreeNode(object):
    """
    Class object that stores properties that allow for a functional kDTree to be implemented on point data.
    """

    def __init__(self, data, axis, left=None, right=None, depth=None):
        """initializes kDTreeNode object with data,axis,left,right, and depth properties"""
        self.data = data
        self.left = left
        self.right = right
        self.axis = axis
        self.depth = depth

    def __repr__(self):
        """returns a string indicating the point data and the object type"""
        return f'kDTreeNode{self.data}'

    def __getitem__(self, i):
        """returns coordinate at the i-th index"""
        return self.data[i]

    def __le__(self, other):
        """returns True if a node is less than or equal to its comparator along the axis of the node and False if
        otherwise"""
        if isinstance(other, kDTreeNode):
            if self.data[self.axis] <= other.data[self.axis]:
                return True
        else:
            if self.data[self.axis] <= other[self.axis]:
                return True
        return False

    def __ge__(self, other):
        """returns True if a node is greater than or equal to its comparator along the axis of the node and False if
        otherwise"""
        if isinstance(other, kDTreeNode):
            if self.data[self.axis] >= other.data[self.axis]:
                return True
        else:
            if self.data[self.axis] >= other[self.axis]:
                return True
        return False

    def __lt__(self, other):
        """returns True if a node is less than its comparator along the axis of the node and False if otherwise"""
        if isinstance(other, kDTreeNode):
            if other.data[self.axis] > self.data[self.axis]:
                return True
        else:
            if other[self.axis] > self.data[self.axis]:
                return True
        return False

    def __gt__(self, other):
        """returns True if a node is greater than its comparator along the axis of the node and False if otherwise"""
        if isinstance(other, kDTreeNode):
            if other.data[self.axis] < self.data[self.axis]:
                return True
        else:
            if other[self.axis] < self.data[self.axis]:
                return True
        return False

    def __eq__(self, other):
        """returns True if a node and its comparator are equal along the axis of the node and False if otherwise"""
        if isinstance(other, kDTreeNode):
            if other.data[self.axis] == self.data[self.axis]:
                return True
        else:
            if other[self.axis] == self.data[self.axis]:
                return True
        return False

    def point_distance(self, other):
        """Calculates linear distance between two points. Supports kDTreeNode objects and data with the same number of
        indices as point data on the node assuming indices for 0 and 1 are x and y coordinates"""
        if isinstance(other, kDTreeNode):
            other = (other.data[0], other.data[1])
        return math.sqrt((self.data[0] - other[0]) ** 2 + (self.data[1] - other[1]) ** 2)


class kDTree(object):
    """
    class object storing a collection of spatially indexed point or point like objects.  Supports O(log n) search,
    nearest neighbor, and circular range queries.
    Initialization time complexity is O(n^2 log n) but guarentees a balanced tree once the tree is constructed.
    Insertions are O(log n) once when the tree is first constructed, however, after many insertions, the
    balanced property of the tree is not guarenteed.  A rebalance method is offered, but this would run in the same
    O(n^2 log n) of intitialization, and essentially just rebuild the tree. Deletions are not supported.

    """

    def __init__(self, points=None, k=2):
        """
        points: an array of point or point like objects that would be initialized into a kDTree
        k: the dimensions of each individual point. Dimensions must match dimensions of individual point objects
        """
        self.k = k
        self.root = None
        if points is not None:
            self.initialize(points)

    def __len__(self):
        """
        returns the number of nodes in the tree
        """
        return len(self.nodes())

    def initialize(self, points):
        """
        public method to initialize a kDTree with an array like input of point or point like objects
        """
        sorted_data = sorted(points, key=lambda x: x[0])
        self.root = kDTreeNode(sorted_data[len(points) // 2], 0, depth=0)
        if len(sorted_data) > 2:
            self._add_subtree(self.root, sorted_data[:len(points) // 2], (self.root.axis) % self.k, left=True)
            self._add_subtree(self.root, sorted_data[len(points) // 2 + 1:], (self.root.axis) % self.k, left=False)
        elif len(sorted_data) == 2:  # if data is only two items, root is larger item , and smaller item is left child
            self.root.left = kDTreeNode(sorted_data[0], (self.root.axis) % self.k, depth=self.root.depth + 1)

    def _add_subtree(self, node, data, axis, left=True):
        """
        private recursive method to inititalize a balanced subtree starting at the input node. Converts the input
        data into kdTreeNodes.
        """
        if len(data[
                   0]) != self.k:  # if a point's dimensions are greater or less than the tree's dimentions, raise error
            raise TypeError('The dimensions of point object(s) do not match specified kDTree dimensions')
        data = sorted(data, key=lambda x: x[(axis + 1) % self.k])  # sort the data along the axis of node's children
        if len(data) > 1:
            new_node = kDTreeNode(data[len(data) // 2], (axis + 1) % self.k, depth=node.depth + 1)  # create child node
            if left:  # if the boolean parameter indicates creating a left child, store new node as left child
                node.left = new_node
            else:  # otherwise make the new node the right child
                node.right = new_node
            left_data = data[:len(data) // 2]  # copy smaller half of data
            right_data = data[len(data) // 2 + 1:]  # copy larger half of data
            if left_data:  # if there is any data with values less than new neode, create a left subtree
                self._add_subtree(new_node, left_data, new_node.axis, left=True)
            if right_data:  # if there is any data with values more than new neode, create a right subtree
                self._add_subtree(new_node, right_data, new_node.axis, left=False)
        else:  # if length of data is only one item, determine if item is less than or greater than node, store accordingly
            if data[0][axis] < node.data[axis]:
                node.left = kDTreeNode(data[0], (axis + 1) % self.k, depth=node.depth + 1)
            else:
                node.right = kDTreeNode(data[0], (axis + 1) % self.k, depth=node.depth + 1)

    def rebalance(self):
        """
        public method to rebalance a tree.  Lazy algorithm to re-initialize a balanced tree.
        """
        self.initialize(self.items())
        return

    def _subtree_preorder(self, node):
        """
        private method to return an a preorder generator of kDTreeNodes starting at the input node.
        """
        if node is not None:
            yield node
            for child in [node.left, node.right]:
                for node in self._subtree_preorder(child):
                    yield node

    def items(self):
        """
        public method to return a preorder of the data items stored at the kDTreeNodes in the tree as a list.
        """
        return [i.data for i in self._subtree_preorder(self.root)]

    def nodes(self):
        """
        public method to return a preorder of the kDTreeNodes in the tree as a list.
        """
        return [i for i in self._subtree_preorder(self.root)]

    def height(self, node=None):
        """
        public method to return the depth or height of a the kDTree subtree rooted at any node.  If no node is specified,
        the entire tree height is returned.
        """
        if node is None:
            return max([node.depth for node in self._subtree_preorder(self.root)])
        else:
            return max([node.depth for node in self._subtree_preorder(node)]) - node.depth

    def _subtree_search(self, node, data):
        """
        private search utility to return the node storing the input data within a subtree. If the data is not in the
        tree, the node holding the nearest value is returned
        """
        if isinstance(data, kDTreeNode):  # simplify data of a kDTreeNode if offered as a parameter
            data = (data.data[i] for i in range(len(data.data)))
        if node.data == data:
            return node  # if the item is found, return the node storing the item
        elif node.data[node.axis] >= data[node.axis]:
            if node.left is not None:  # if item is smaller than node and node has a left child, recur
                return self._subtree_search(node.left, data)
        elif node.data[node.axis] < data[node.axis]:
            if node.right is not None:  # if item is greater than node and node has a right child, recur
                return self._subtree_search(node.right, data)
        return node

    def search(self, data):
        """
        private search utility to return the node storing the input data within a kDTree. If the data is not in the
        tree, the node holding the nearest value is returned
        """
        return self._subtree_search(self.root, data)

    def insert(self, data):
        """
        public method to insert a point or point like data into the kDTree.  If the point data is already in the tree,
        no data is replaced and the method returns False. If data is inserted it returns True
        """
        if isinstance(data, kDTreeNode):  # simplify data if kDTreeNode is offered as a parameter
            data = (data.data[i] for i in range(len(data.data)))
        parent = self.search(data)  # find closest node
        if parent.data == data:  # if closest node has the same data as the data to be inserted, return
            return False
        if data[parent.axis] >= parent.data[
            parent.axis]:  # otherwise, insert data into right or left child appropriately
            parent.right = kDTreeNode(data, (parent.axis + 1) % self.k, depth=parent.depth + 1)
        if data[parent.axis] < parent.data[parent.axis]:
            parent.left = kDTreeNode(data, (parent.axis + 1) % self.k, depth=parent.depth + 1)
        return True

    def _visualize_from_node(self, node, margin, depth):
        """
        private method to visualize a subtree of a kDTree
        """
        space = margin * depth
        print(f'{space}data: {node.data}\n{space}axis:{node.axis}\n{space}left: {node.left}\n{space}right:{node.right}\
        \n{space}depth: {node.depth}')
        print('')
        if node.left is not None:
            self._visualize_from_node(node.left, margin, depth + 1)
        if node.right is not None:
            self._visualize_from_node(node.right, margin, depth + 1)

    def visualize(self):
        """
        public method to view the nodes of kDTree
        """
        return self._visualize_from_node(self.root, '     ', 0)

    def _update_neighbors_and_max_distance(self, node, point, neighbors, n):
        """
        private function that assess how and if a node is a neareast neighbor to an input point
        """
        distance = node.point_distance(point)
        for idx, data in enumerate(neighbors):
            if distance < data[1]:  # if distance is less than existing neighbor, add node before that neighbor
                neighbors.insert(idx, (node, distance))
                if len(neighbors) < n:  # if neighbors hold fewer then desired output, return no maximum distance
                    return float('inf')
                return neighbors[n - 1][1]  # otherwise return distance of nth neighbor
        neighbors.append((node, distance))  # add evaluated neighbor if there are no existing neighbors
        return float('inf')  # return no maximum distance f there is only one neighbor

    def _subtree_nearest_neighbor(self, node, point, neighbors, n, maxdist):
        """
        private recursive method to evaluate a subtree in the kDTree for nearest neighbors.  Assesses which subtree's
        must be visited based on the difference between the value of the nodes axis and the point axis
        """
        if node is None:
            return
        if node.left is None and node.right is None:  # base case-- closest neighbor is reached
            self._update_neighbors_and_max_distance(node, point, neighbors, n)
            return
        if node < point:  # define closer and farther subtrees based on comparison of node and pt
            closer, farther = node.right, node.left
        else:
            closer, farther = node.left, node.right
        self._subtree_nearest_neighbor(closer, point, neighbors, n, maxdist)  # recur on closest subtree
        maxdist = self._update_neighbors_and_max_distance(node, point, neighbors,
                                                          n)  # add neighbor to neighbors if needed
        if abs(node.data[node.axis] - point[node.axis]) < maxdist:  # if farther neighbor is still close enough, recur
            self._subtree_nearest_neighbor(farther, point, neighbors, n, maxdist)

    def nearest_neighbor(self, pt, n=1):
        """
        Finds the n number of closest neighbors to a given point
        """
        neighbors = []
        self._subtree_nearest_neighbor(self.root, pt, neighbors=neighbors, n=n, maxdist=float('inf'))
        return neighbors[:n]

    def _subtree_range_query(self, node, pt, radius, neighbors):
        """
        private method to return all kDTreeNodes within the radius of the input point in a subtree starting at input node
        """
        if node is None:
            return
        if isinstance(pt, kDTreeNode):
            pt = (pt.data[i] for i in range(len(pt.data)))
        elif len(pt) > self.k:
            raise ValueError('Input point dimensions must match kDTree dimensions.')
        if node >= [pt[i] + radius for i in
                    range(len(pt))]:  # if the node is greater than point, recur only on left child
            return self._subtree_range_query(node.left, pt, radius, neighbors)
        if node < [pt[i] - radius for i in
                   range(len(pt))]:  # if the node is less than point, recur only on right child
            return self._subtree_range_query(node.right, pt, radius, neighbors)
        dist = node.point_distance(pt)  # if node is possibly in the range, calculate distance
        if dist <= radius:
            neighbors.append((node, dist))  # recur on both children
        self._subtree_range_query(node.left, pt, radius, neighbors)
        self._subtree_range_query(node.right, pt, radius, neighbors)

    def range_query(self, pt, radius=5):
        """
        public method to return all kDTreeNodes within the radius of the input point. returns
        """
        neighbors = []
        self._subtree_range_query(self.root, pt, radius, neighbors)
        neighbors.sort(key=lambda x: x[1])
        return neighbors
