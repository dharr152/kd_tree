import arcpy
import math


class ArcArcKDTreeNode(object):
    """
    Class object that stores properties that allow for a functional ArcKDTree to be implemented on any acpy shape based
    on centroid.
    """

    def __init__(self, axis, left=None, right=None, depth=None, shape=None, x=None, y=None, data=None):
        """initializes ArcArcKDTreeNode object with data,axis,left,right, and depth properties"""
        self.left = left
        self.right = right
        self.axis = axis
        self.depth = depth
        self.shape = shape
        self.X = x
        self.Y = y
        self.data = data

    def point_distance(self, other):
        """Calculates linear distance between two points. Supports ArcArcKDTreeNode objects and ArcGIS Points
         as point data"""

        return math.sqrt((self.X - other.X) ** 2 + (self.Y - other.Y) ** 2)

    def __ge__(self, other):
        """returns True if a node is greater than or equal to its comparator along the axis of the node and False if
        otherwise"""
        if self.axis == 0:
            if self.X >= other.X:
                return True
            else:
                return False
        if self.axis == 1:
            if self.Y >= other.Y:
                return True
            else:
                return False

    def __lt__(self, other):
        """returns True if a node is less than its comparator along the axis of the node and False if otherwise"""
        if self.axis == 0:
            if self.X < other.X:
                return True
            else:
                return False
        if self.axis == 1:
            if self.Y < other.Y:
                return True
            else:
                return False


class ArcKDTree(object):
    """
    class object storing a collection of spatially indexed point or point like objects.  Supports O(log n) circular
    range queries. Initialization time complexity is O(n^2 log n) but guarentees a balanced tree once
    the tree is constructed.
    """

    def __init__(self, points=None):
        """
        points: an array of point or point like objects that would be initialized into a ArcKDTree
        k: the dimensions of each individual point. Dimensions must match dimensions of individual point objects
        """
        self.k = 2
        self.root = None
        if points is not None:
            self.initialize(points)

    def initialize(self, shapes):
        """
        public method to initialize a ArcKDTree with an array like input of point or point like objects
        """
        sorted_data = sorted(shapes, key=lambda x: x[0].centroid.X)
        shape = 0
        fields = 1
        self.root = ArcArcKDTreeNode(axis=0,
                               depth=0,
                               shape=sorted_data[len(shapes) // 2][shape],
                               x=sorted_data[len(shapes) // 2][shape].centroid.X,
                               y=sorted_data[len(shapes) // 2][shape].centroid.Y,
                               data=sorted_data[len(shapes) // 2][fields])
        if len(sorted_data) > 2:
            self._add_subtree(self.root, sorted_data[:len(shapes) // 2], (self.root.axis) % self.k, left=True)
            self._add_subtree(self.root, sorted_data[len(shapes) // 2 + 1:], (self.root.axis) % self.k, left=False)
        elif len(sorted_data) == 2:  # if data is only two items, root is larger item , and smaller item is left child
            self.root.left = ArcArcKDTreeNode(axis=(self.root.axis) % self.k, depth=self.root.depth + 1,
                                        shape=sorted_data[0][0], x=sorted_data[0][0].centroid.X,
                                        y=sorted_data[0][0].centroid.Y, data=sorted_data[1])

    def _add_subtree(self, node, data, axis, left=True):
        """
        private recursive method to inititalize a balanced subtree starting at the input node. Converts the input
        data into ArcArcKDTreeNodes.
        """
        shape = 0
        fields = 1
        data = sorted(data, key=lambda x: x[shape].centroid.Y if axis is 0 else x[
            shape].centroid.X)  # sort the data along the axis of node's children

        while True:  # loop ensures that only records with a shape and fields are loaded into the kd
            if len(data) > 0:
                if len([i for i in data[len(data) // 2][fields] if i is not None]) < 1:
                    data.pop(len(data) // 2)
                else:
                    break
            else:
                break

        if len(data) > 1:
            new_node = ArcArcKDTreeNode(axis=1 if axis is 0 else 0,
                                  shape=data[len(data) // 2][shape],
                                  depth=node.depth + 1,
                                  x=data[len(data) // 2][shape].centroid.X,
                                  y=data[len(data) // 2][shape].centroid.Y,
                                  data=data[len(data) // 2][fields])  # create child node
            if left:  # if the boolean parameter indicates creating a left child, store new node as left child
                node.left = new_node
            else:  # otherwise make the new node the right child
                node.right = new_node
            left_data = data[:len(data) // 2]  # copy smaller half of data
            right_data = data[len(data) // 2 + 1:]  # copy larger half of data
            if left_data:  # if there is any data with values less than new node, create a left subtree
                self._add_subtree(new_node, left_data, new_node.axis, left=True)
            if right_data:  # if there is any data with values more than new node, create a right subtree
                self._add_subtree(new_node, right_data, new_node.axis, left=False)
        elif len(
                data) > 0:  # if length of data is only one item, determine if item is less than or greater than node,
                            # store accordingly
            if node >= data[0][shape].centroid:
                node.left = ArcArcKDTreeNode(axis=1 if axis is 0 else 0,
                                       depth=node.depth + 1,
                                       shape=data[0][shape],
                                       x=data[0][shape].centroid.X,
                                       y=data[0][shape].centroid.Y,
                                       data=data[0][fields])
            else:
                node.right = ArcArcKDTreeNode(axis=1 if axis is 0 else 0,
                                        depth=node.depth + 1,
                                        shape=data[0][shape],
                                        x=data[0][shape].centroid.X,
                                        y=data[0][shape].centroid.Y,
                                        data=data[0][fields])

    def _subtree_range_query(self, node, pt, radius, neighbors):
        """
        private method to return all ArcArcKDTreeNodes within the radius of the input point in a subtree starting at
        input node
        """
        if node is None:
            return
        if node >= arcpy.Point(pt.X + radius,
                               pt.Y + radius):  # if the node is greater than point, recur only on left child
            return self._subtree_range_query(node.left, pt, radius, neighbors)
        if node < arcpy.Point(pt.X - radius,
                              pt.Y - radius):  # if the node is less than point, recur only on right child
            return self._subtree_range_query(node.right, pt, radius, neighbors)
        dist = node.point_distance(pt)  # if node is possibly in the range, calculate distance
        if dist <= radius:
            neighbors.append((node, dist))  # recur on both children
        self._subtree_range_query(node.left, pt, radius, neighbors)
        self._subtree_range_query(node.right, pt, radius, neighbors)

    def range_query(self, pt, radius_param=5):
        """
        public method to return all ArcArcKDTreeNodes within the radius of the input point. returns
        """
        neighbors = []
        self._subtree_range_query(self.root, pt, radius_param, neighbors)
        neighbors.sort(key=lambda x: x[1])
        return neighbors

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
        private recursive method to evaluate a subtree in the ArcKDTree for nearest neighbors.  Assesses which subtree's
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
        node_coord = node.X if node.axis == 0 else node.Y
        point_coord = point.X if node.axis == 0 else node.Y
        if abs(node_coord - point_coord) < maxdist:  # if farther neighbor is still close enough, recur
            self._subtree_nearest_neighbor(farther, point, neighbors, n, maxdist)

    def nearest_neighbor(self, pt, n=1):
        """
        Finds the n number of closest neighbors to a given point
        """
        neighbors = []
        self._subtree_nearest_neighbor(self.root, pt, neighbors=neighbors, n=n, maxdist=float('inf'))
        return neighbors[:n]


def load_kd(fc, fields):
    """
    Returns an ArcKDTree object with the fields parameter loaded as data for each feature in a feature class of
    shapefile.
    param fc: a feature class to be indexed based in its centroid
    param fields: fields form the feature class that will be stored as the 'data' property at each node in the
        kd tree. 
    """
    try:
        features = []
        with arcpy.da.SearchCursor(fc, ['SHAPE@'] + fields) as scur:
            for row in scur:
                num = len(fields)
                features.append([row[0], [row[i + 1] for i in range(num)]])
    except Exception as E:
        print(str(E))
        arcpy.AddError(str(E))

    return ArcKDTree(features)

