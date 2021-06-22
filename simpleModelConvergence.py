import random
import numpy as np
from time import sleep
import expectedMatrix as em



class Model():

    def __init__(self, N, P):
        self.N = N
        self.P = P
        self.terminated = np.full(self.N, False)
        self.expected_a, indices = em.expectedMatrix(self.N, P)
        self.expected_a = 1/(P+1)*self.expected_a

    def initialOpinions(self):
        opinions = np.full(self.N, -1)
        count = 0
        for i in range(self.N):
            opinions[i] = random.randint(0, 1)
            count += opinions[i]
        print('There are {} many nodes starting with opinion 1'.format(count))
        print(np.mean(opinions))
        print('\n{}\n'.format(opinions))
        return opinions

    def randomRowStochasticMatrix(self, terminated):
        a = np.zeros((self.N,self.N))
        np.fill_diagonal(a, 1/(self.P+1))
        indices_in_rows = []
        for i in range(self.N):
            if terminated[i] == False:
                # print('Computing raw: {}'.format(i))
                j = 0
                indices = [i]
                while j < self.P:
                    index = random.randint(0, self.N-1)
                    if not index in indices:
                        j += 1
                        indices.append(index)
                for j in indices:
                    a[i][j] = 1/(self.P+1)
                # print('Raw {}: {}'.format(i, a[i]))
                indices_in_rows.append(indices)
            else:
                a[i] = np.zeros(self.N)
                a[i][i] = 1
                indices_in_rows.append([])
        return a, indices_in_rows

    def nextIteration(self, opinions, terminated):
        a, indices = self.randomRowStochasticMatrix(terminated)
        next_opinions = np.matmul(a, opinions)
        return next_opinions, indices

    def nextExpectedIteration(self, opinions):
        next_opinions = np.matmul(self.expected_a, opinions)
        return next_opinions

    def simpleModel(self):
        opinions = self.initialOpinions()
        next_opinions, indices = self.nextIteration(opinions, self.terminated)
        print('\n{}'.format(next_opinions))
        # next_opinions = np.around(next_opinions)
        # print('{}\n'.format(next_opinions))
        sleep(1)
        while True:
            next_opinions, indices = self.nextIteration(next_opinions, self.terminated)
            print('\n{}'.format(next_opinions))
            # next_opinions = np.around(next_opinions)
            # print('{}\n'.format(next_opinions))
            sleep(1)

    # simpleModel()