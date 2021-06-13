import numpy as np
import simpleModelConvergence as model
from time import sleep

def checkTerminatedNodes():
    for i in range(len(indices_in_rows)):
        # print(indices[i])
        queryer_opinion = next_opinions[i]
        if np.all([queryer_opinion == next_opinions[index] for index in indices_in_rows[i][1:]]):
            # node i has the same opinion of all the nodes it queried and so it can terminate
            print('Node {} terminated with opinion {}'.format(i, next_opinions[i]))
            terminated[i] = True

N = 100
terminated = np.full(N, False)

opinions = model.initialOpinions()
next_opinions, indices_in_rows = model.nextIteration(opinions, terminated)
# print('\n{}'.format(next_opinions))
next_opinions = [1 if val>=0.5 else 0 for val in np.around(next_opinions)]
print('\n{}'.format(next_opinions))
checkTerminatedNodes()
sleep(1)
while not np.all(terminated):
    next_opinions, indices_in_rows = model.nextIteration(next_opinions, terminated)
    # print('\n{}'.format(next_opinions))
    next_opinions = [1 if val>=0.5 else 0 for val in np.around(next_opinions)]
    print('\n{}'.format(next_opinions))
    checkTerminatedNodes()
    sleep(1)
