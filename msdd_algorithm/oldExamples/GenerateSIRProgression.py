import sys, os
from algorithm import msdd_algorithm_simple, generate_class_from_rule_dict_large, clean_data

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from SequenceGenerator import SequenceGenerator
from pprint import pprint


def create_sir_data(transition=None, sample_size = 5000):
    transition = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                  'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    SIR = SequenceGenerator(transition_prob=transition)
    sequence_data = SIR.generate_states(current_state='S', nt=2000)
    with open('data/sir_data.txt','w') as file:
        for item in sequence_data:
            file.writelines('{},'.format(item))

def read_data_lower_middle_SIR():
    with open('data/sir_data.txt', 'r') as file:
        l = file.readlines()
    # print(l[0].split(','))
    l = l[0].split(',')
    l = clean_data(l)
    return l


def test_msdd_simple2():
    sequence_data = read_data_lower_middle_SIR()
    rule_dict = msdd_algorithm_simple(sequence_data, precursor_width= 1, sucessor_width= 1, lagtime_between=1, dependency_evaluator_fn = None, num_rules = None)
    pprint(rule_dict)
    generate_class_from_rule_dict_large(rule_dict)

if __name__ == '__main__':
    # create_sir_data()
    # print(read_data_lower_middle_SIR())
    # print(len(read_data_lower_middle_SIR()))
    test_msdd_simple2()