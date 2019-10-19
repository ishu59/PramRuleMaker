import numpy as np
from SequenceGenerator import SequenceGenerator

class MarkovChain:
    def __init__(self, transition_matrix, states):
        np.random.seed(0)
        self.transition_matrix = np.atleast_2d(transition_matrix)
        self.states = states
        self.index_dict = {self.states[index]: index for index in
                           range(len(self.states))}
        self.state_dict = {index: self.states[index] for index in
                           range(len(self.states))}

    def next_state(self, current_state):
        return np.random.choice(
            self.states,
            p=self.transition_matrix[self.index_dict[current_state], :])

    def generate_states(self, current_state, no=10):
        future_states = []
        for i in range(no):
            next_state = self.next_state(current_state)
            future_states.append(next_state)
            current_state = next_state
        return future_states

class HiddenModel:
    def __init__(self, hidden_states, observation_states, prior_probabilities,
                 transition_matrix, emission_probabilities):
        self.hidden_markov_chain = MarkovChain(
            transition_matrix=transition_matrix,
            states=hidden_states,
        )
        self.observation_states = observation_states
        self.prior_probabilities = np.atleast_1d(prior_probabilities)
        self.transition_matrix = np.atleast_2d(transition_matrix)
        self.emission_probabilities = np.atleast_2d(emission_probabilities)


    def observation_from_state(self, state):
        state_index = self.hidden_markov_chain.index_dict[state]
        return np.random.choice(self.observation_states,
                                p=self.emission_probabilities[state_index, :])


    def generate_samples(self, no=10):
        observations = []
        state_sequence = []
        initial_state = np.random.choice(
            self.hidden_markov_chain.states,
            p=self.prior_probabilities)
        state_sequence.append(initial_state)
        observations.append(self.observation_from_state(initial_state))
        current_state = initial_state
        for i in range(2, no):
            next_state = self.hidden_markov_chain.next_state(current_state)
            state_sequence.append(next_state)
            observations.append(self.observation_from_state(next_state))
            current_state = next_state
        return observations, state_sequence

def test_seq_generator():
    transition_prob = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.8, 'R': 0.2}, 'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    # [[0.9,0.1,0.0],[0,0.8,0.2],[0.1,0.0,0.9]]
    SIR = SequenceGenerator(transition_prob=transition_prob)
    SIR.next_state(current_state='S')
    SIR.next_state(current_state='R')
    print(SIR.generate_states(current_state='S', nt=30))

def test_mul_seq_generator():
    test_hmm = HiddenModel(hidden_states=['S','I','R'], observation_states=['school', 'home'], prior_probabilities=[1, 0, 0],
                           transition_matrix=[[0.9, 0.1, 0],[0, 0.8, 0.2], [0.1, 0.0, 0.9]],
                           emission_probabilities=[[0.8, 0.2], [0.5, 0.5],[0.1, 0.9 ]])
    data = test_hmm.generate_samples(no = 20)
    print(data)

if __name__ == '__main__':
    # test_seq_generator()
    test_mul_seq_generator()