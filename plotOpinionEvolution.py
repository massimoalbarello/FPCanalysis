import simpleModelConvergence as model
import expectedMatrix as em
import numpy as np
from time import sleep
import matplotlib.pyplot as plt
from scipy.linalg import expm

def computeLaplacian(a):
    dOut = np.zeros(a.shape)
    laplacian = np.zeros(a.shape)
    for i in range(a.shape[0]):
        dOut[i][i] = np.sum(a[i])
    # print(dOut)
    laplacian = dOut - a
    print(laplacian)
    return laplacian


N = 10
P = 2
iterations = 30
terminated = np.full(N, False)
opinions_evolution = []
node_opinion_ev = []
expected_node_opinion_ev = []
time = np.linspace(0, iterations+2, iterations*10)


opinions = model.initialOpinions()
opinions_evolution.append(opinions)
next_opinions, indices = model.nextIteration(opinions, terminated)
opinions_evolution.append(next_opinions)
# print('\n{}'.format(next_opinions))
# next_opinions = np.around(next_opinions)
# print('{}\n'.format(next_opinions))
# sleep(1)
for i in range(iterations):
    next_opinions, indices = model.nextIteration(next_opinions, terminated)
    opinions_evolution.append(next_opinions)
    # print('\n{}'.format(next_opinions))
    # next_opinions = np.around(next_opinions)
    # print('{}\n'.format(next_opinions))
    # sleep(1)

for ops in opinions_evolution:
    node_opinion_ev.append(ops[0])

plt.plot(node_opinion_ev)
# plt.show()

expected = em.expectedMatrix(N, P)

exp_laplacian = computeLaplacian(expected)

for t in time:
    expected_node_opinion_ev.append(np.matmul(expm(-exp_laplacian*t), opinions)[0])

plt.plot(time, expected_node_opinion_ev)

plt.show()
