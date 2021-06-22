import numpy as np
import random

iterations = 10000

def randomMatrix(N, P):
    a = np.zeros((N,N))
    np.fill_diagonal(a, 1)
    for i in range(N):
        # print('Computing raw: {}'.format(i))
        j = 0
        indices = [i]
        while j < P:
            index = random.randint(0, N-1)
            if not index in indices:
                j += 1
                indices.append(index)
        for j in indices:
            a[i][j] = 1
        # print('Raw {}: {}'.format(i, a[i]))
    return a

def expectedMatrix(N, P):
    expected = np.zeros((N,N))
    for i in range(iterations):
        expected += 1/iterations * randomMatrix(N, P)
    # print(expected)
    return expected, []

# expectedMatrix(10, 3)