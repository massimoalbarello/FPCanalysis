import simpleModelConvergence
import numpy as np
from time import sleep
import matplotlib.pyplot as plt

def initialOpinions():
    opinions0 = np.full(int(N*1/4), 0)
    opinions1 = np.full(int(N*3/4), 1)
    # print(opinions0)
    # print(opinions1)
    opinions = np.concatenate((opinions0, opinions1))
    count = 0
    for i in range(N):
        # opinions[i] = random.randint(0, 1)
        count += opinions[i]
    print('There are {} many nodes starting with opinion 1'.format(count))
    print(np.mean(opinions))
    print('\n{}\n'.format(opinions))
    return opinions

def maxErrorInIteration(iteration):
    max = 0
    for realization in realizations_node_opinion_evolution:
        if abs(realization[iteration] - np.array(expected_opinions_evolution)[iteration,0]) > max:
            max = abs(realization[iteration] - np.array(expected_opinions_evolution)[iteration,0])
    return max

def setPlot(row, col, title, xlabel, ylabel):
    axes[row, col].set_title(title)
    for label in (axes[row, col].get_xticklabels() + axes[row, col].get_yticklabels()):
	    label.set_fontsize(15)
    axes[row, col].set_xlabel(xlabel, fontsize=15)
    axes[row, col].set_ylabel(ylabel, fontsize=15)
    if not(row == 1 and col == 2):
        axes[row, col].label_outer()

N = 100
P_VAL = [1, 5, 10, 20, 50]
opinions = initialOpinions()
figure, axes = plt.subplots(2, int((len(P_VAL)+1)/2), sharex=True, sharey=True)
err_labels = ["P = 1", "P = 5", "P = 10", "P = 20", "P = 50"]
plt.rcParams['font.size'] = '15'

row = 0
col = 0
for iter, P in enumerate(P_VAL):
    model = simpleModelConvergence.Model(N, P)
    iterations = 10
    realizations = 5
    terminated = np.full(N, False)
    realizations_node_opinion_evolution = []
    expected_opinions_evolution = []
    errors = []


    for i in range(realizations):

        opinions_evolution = []
        opinions_evolution.append(opinions)
        next_opinions, indices = model.nextIteration(opinions, terminated)
        opinions_evolution.append(next_opinions)

        for j in range(iterations):
            next_opinions, indices = model.nextIteration(next_opinions, terminated)
            opinions_evolution.append(next_opinions)

        axes[row, col].plot(np.array(opinions_evolution)[:,0])

        realizations_node_opinion_evolution.append(np.array(opinions_evolution)[:,0])



    expected_opinions_evolution.append(opinions)
    expected_next_opinions = model.nextExpectedIteration(opinions)
    expected_opinions_evolution.append(expected_next_opinions)

    for i in range(iterations):
        expected_next_opinions = model.nextExpectedIteration(expected_next_opinions)
        expected_opinions_evolution.append(expected_next_opinions)

    exp_plt, = axes[row, col].plot(np.array(expected_opinions_evolution)[:,0], label="Expected evolution")
    legend = axes[row, col].legend(handles=[exp_plt,])



    for iteration in range(iterations+2):
        errors.append(maxErrorInIteration(iteration))

    err_plt, = axes[1, 2].plot(errors, label=err_labels[iter])


    col += 1
    if col == int((len(P_VAL)+1)/2):
        col = 0
        row += 1

setPlot(0, 0, "Node opinion for P = 1", "iterations", "opinion")
setPlot(0, 1, "Node opinion for P = 5", "iterations", "opinion")
setPlot(0, 2, "Node opinion for P = 10", "iterations", "opinion")
setPlot(1, 0, "Node opinion for P = 20", "iterations", "opinion")
setPlot(1, 1, "Node opinion for P = 50", "iterations", "opinion")
setPlot(1, 2, "Max error for each iteration", "iterations", "error")

axes[1,2].legend()
plt.show()
