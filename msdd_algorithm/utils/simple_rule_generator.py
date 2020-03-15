from collections import OrderedDict
import numpy as np
import pandas as pd
from SequenceGenerator import SequenceGenerator
# Generate transition Matrix
def simple_rule_generator(state_sequence, states=None):
    if states is None:
        states = list(np.unique(state_sequence))
    # states = ['S', 'I', 'R']
    num_notation = zip(states, range(1, len(states)+1))
    numeric_transition = OrderedDict(num_notation)
    transitions = []
    for i, item in enumerate(state_sequence):
        transitions.append(numeric_transition[item])
    numeric_state_translation = {val:key for key, val in numeric_transition.items()}
    n = 1+ max(transitions)
    M = [[0]*n for _ in range(n)]
    for (i,j) in zip(transitions,transitions[1:]):
        M[i][j] += 1
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    predicted_rule = {}
    for i in range(1, len(M)):
        child_dict = {}
        for j in range(1, len(M[i])):
            child_dict[numeric_state_translation[j]] = M[i][j]
        predicted_rule[numeric_state_translation[i]] = child_dict
    return predicted_rule


def test_simple_rule():
    transition_prob = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.8, 'R': 0.2},
                       'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    transition_prob = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.7, 'R': 0.3},
                       'R': {'S': 0.1, 'I': 0.2, 'R': 0.7}}
    SIR = SequenceGenerator(transition_prob=transition_prob)
    seq = SIR.generate_states(current_state='S', nt=2000)
    print(seq)
    m = simple_rule_generator(seq)
    print(m)

if __name__ == '__main__':
    # test
    test_simple_rule()