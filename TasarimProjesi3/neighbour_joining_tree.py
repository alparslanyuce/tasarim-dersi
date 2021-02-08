
from random import randint


class NeighbourJoiningTree:
    def __init__(self):
        self.max = 0
        self.max_matrix = []
        self.distances = []

    def calculate_max(self, matrix):
        max_value = 0
        for i in range(1, len(matrix)):
            max_value = max(max_value, len(matrix[i]))

        self.max = max_value

    def normalize(self, matrix):
        self.calculate_max(matrix)
        for row in range(0, len(matrix)):
            while len(matrix[row]) < self.max:
                random_int = randint(0, len(matrix[row]) - 1)
                matrix[row] = matrix[row][0:random_int] + "*" + matrix[row][random_int::]
        self.max_matrix = matrix

    def calculate_matrix_distances(self):
        self.distances = [[0] * len(self.max_matrix) for i in range(len(self.max_matrix))]

        length = len(self.max_matrix[0])
        for i in range(len(self.max_matrix) - 1):
            for j in range(i + 1, len(self.max_matrix)):
                distance = 0
                for k in range(length):
                    if self.max_matrix[i][k] != self.max_matrix[j][k]:
                        distance += 1

                self.distances[j][i] = distance

