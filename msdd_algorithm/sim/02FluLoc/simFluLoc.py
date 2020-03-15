import os,sys
from algorithm import generate_class_from_rule_dict_large, msdd_algorithm_simple, clean_data
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from collections import defaultdict
from pprint import pprint
from utils.DependentSequenceGenerator import HiddenModel

def algorithm_multi_sequence_data(sequence_data=None,observed_data=None,  precursor_width=None, sucessor_width=None, lagtime_between=None, dependency_evaluator_fn = None, num_rules = None):
    rule_dict = msdd_algorithm_simple(sequence_data, precursor_width= 1, sucessor_width= 1, lagtime_between=1, dependency_evaluator_fn = None, num_rules = None)
    emission_count_dict = {}
    emission_rule_dict = {}
    for hidden in rule_dict.keys():
        emission_count_dict[hidden] = defaultdict(float)
        emission_rule_dict[hidden] = defaultdict(float)
    for index, hidden in enumerate(sequence_data):
        emission_count_dict[(hidden,)][observed_data[index]] += 1
    # pprint(emission_count_dict)
    for hidden_key, hidden_val in emission_count_dict.items():
        total = sum(hidden_val.values())
        # print(total)
        for observed_key, observed_val in hidden_val.items():
            emission_rule_dict[hidden_key][observed_key] = observed_val/total
    # pprint(emission_rule_dict)
    return [rule_dict, emission_rule_dict]

def generate_multiseq_data():
    hidden_states = ['S', 'I', 'R']
    observation_states = ['school', 'home']
    emission_probabilities = [[0.8, 0.2], [0.5, 0.5], [0.1, 0.9]]
    transition_matrix = [[0.9, 0.1, 0], [0, 0.8, 0.2], [0.1, 0.0, 0.9]]
    initial_condition = [1, 0, 0]
    test_hmm = HiddenModel(hidden_states=hidden_states, observation_states=observation_states,
                           prior_probabilities=initial_condition,
                           transition_matrix=transition_matrix,
                           emission_probabilities=emission_probabilities)
    observed_data, sequence_data = test_hmm.generate_samples(no=10000)
    with open('data/multiseq.txt','w') as file:
        for item in sequence_data:
            file.writelines('{},'.format(item))
        file.writelines('\n')
        for item in observed_data:
            file.writelines('{},'.format(item))

def read_multiseq_data():
    with open('data/multiseq.txt', 'r') as file:
        l = file.readlines()
    # print(l[0].split(','))
    hid = l[0].split(',')
    obs = l[1].split(',')
    hid = clean_data(hid)
    obs = clean_data(obs)
    return  obs, hid

def test_algorithm_multi_seq_data():
    observed_data, sequence_data = read_multiseq_data()
    output = algorithm_multi_sequence_data(sequence_data, observed_data)
    pprint(output)
    has_attr_list = ['flu', 'flu']
    set_attr_list = ['flu', 'loc']
    for index,rule in enumerate(output):
        generate_class_from_rule_dict_large(rule,index, has_attr_list[index], set_attr_list[index])

if __name__ == '__main__':
    test_algorithm_multi_seq_data()
