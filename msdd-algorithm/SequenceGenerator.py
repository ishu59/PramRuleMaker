from math import gcd
from itertools import combinations
from functools import reduce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint
class SequenceGenerator(object):
    def __init__(self, transition_prob):
        self.transition_prob = transition_prob
        self.states = list(transition_prob.keys())

    def next_state(self, current_state):
        return np.random.choice(
            self.states, p=[self.transition_prob[current_state][next_state]
                            for next_state in self.states])
    def generate_states(self, current_state, nt=10):
        future_states = []
        future_states.append(current_state)
        for i in range(nt):
            next_state = self.next_state(current_state)
            future_states.append(next_state)
            current_state = next_state
        return future_states

def test_seq_generator():
    transition_prob = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.8, 'R': 0.2}, 'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    # [[0.9,0.1,0.0],[0,0.8,0.2],[0.1,0.0,0.9]]
    SIR = SequenceGenerator(transition_prob=transition_prob)
    SIR.next_state(current_state='S')
    SIR.next_state(current_state='R')
    print(SIR.generate_states(current_state='S', nt=30))

def test_seq_generator_2d():
    transition_prob = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.8, 'R': 0.2}, 'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    transition_prob = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                  'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    # [[0.9,0.1,0.0],[0,0.8,0.2],[0.1,0.0,0.9]]
    SIR = SequenceGenerator(transition_prob=transition_prob)
    iterations = 30
    SIR.generate_states(current_state='S', nt=iterations)
    data_list = []
    for i in range(1000):
        data_list.append(SIR.generate_states(current_state='S', nt=iterations))
    dl = np.array(data_list)
    print(dl[0])
    freq_count_at_t = {'S':[],'I':[],'R':[]}
    for x in range(iterations):
        uniqueValues, occurCount  =  np.unique(dl[:,x],return_counts = True)
        for i in range(len(uniqueValues)):
            # print(b)
            if len(freq_count_at_t[uniqueValues[i]]) > 25:
                continue
            freq_count_at_t[uniqueValues[i]].append(occurCount[i])
    print()
    df = pd.DataFrame(freq_count_at_t)
    print(df)
    df.plot()
    plt.show()

if __name__ == '__main__':
    test_seq_generator()
    # test_seq_generator_2d()

 # 1  flu: (0.39 0.22 0.40 )   (  386.0   216.5   397.5 )   [1000.0]