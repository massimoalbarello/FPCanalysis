import random
import numpy as np
from time import sleep

N = 10
P = 2
terminated = np.full(N, False)

def initialOpinions():
    opinions = np.full(N, -1)
    count = 0
    for i in range(N):
        opinions[i] = random.randint(0, 1)
        count += opinions[i]
    print('There are {} many nodes starting with opinion 1'.format(count))
    print(np.mean(opinions))
    print('\n{}\n'.format(opinions))
    return opinions

def randomRowStochasticMatrix(terminated):
    a = np.zeros((N,N))
    np.fill_diagonal(a, 1/(P+1))
    indices_in_rows = []
    for i in range(N):
        if terminated[i] == False:
            # print('Computing raw: {}'.format(i))
            j = 0
            indices = [i]
            while j < P:
                index = random.randint(0, N-1)
                if not index in indices:
                    j += 1
                    indices.append(index)
            for j in indices:
                a[i][j] = 1/(P+1)
            # print('Raw {}: {}'.format(i, a[i]))
            indices_in_rows.append(indices)
        else:
            a[i] = np.zeros(N)
            a[i][i] = 1
            indices_in_rows.append([])
    return a, indices_in_rows

def nextIteration(opinions, terminated):
    a, indices = randomRowStochasticMatrix(terminated)
    next_opinions = np.matmul(a, opinions)
    return next_opinions, indices


def simpleModel():
    opinions = initialOpinions()
    next_opinions, indices = nextIteration(opinions, terminated)
    print('\n{}'.format(next_opinions))
    # next_opinions = np.around(next_opinions)
    # print('{}\n'.format(next_opinions))
    sleep(1)
    while True:
        next_opinions, indices = nextIteration(next_opinions, terminated)
        print('\n{}'.format(next_opinions))
        # next_opinions = np.around(next_opinions)
        # print('{}\n'.format(next_opinions))
        sleep(1)

# simpleModel()