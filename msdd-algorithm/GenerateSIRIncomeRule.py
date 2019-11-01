import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from typing import List, Dict
from collections import defaultdict, OrderedDict
import numpy as np
import scipy as sp
from SequenceGenerator import SequenceGenerator
from pprint import pprint
from DependentSequenceGenerator import HiddenModel
from algorithm import msdd_algorithm_simple, clean_data
from PramRulesArributes import PramRulesArributes

def create_income_sir_data(transition_lower=None, transition_upper=None, sample_size = 5000):
    # A Dynamic Bayesian net.
    #
    # Initial state distribution:
    #     flu: (1 0 0)
    #     income: (0.5 0.5)
    #
    # Transition model:
    #                 L                    M
    #           S     I     R        S     I     R
    #     S  0.90  0.00  0.20     0.95  0.00  0.10
    #     I  0.10  0.75  0        0.05  0.50  0
    #     R  0     0.25  0.80     0     0.50  0.90

    # transition_lower = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
    #               'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    transition_lower = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.75, 'R': 0.25},
                        'R': {'S': 0.2, 'I': 0, 'R': 0.8}}
    SIR = SequenceGenerator(transition_prob=transition_lower)
    sequence_data = SIR.generate_states(current_state='S', nt=sample_size)
    with open('data/sir_lower_income.txt','w') as file:
        for item in sequence_data:
            file.writelines('{},'.format(item))

    transition_upper = {'S': {'S': 0.95, 'I': 0.05, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                        'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    SIR = SequenceGenerator(transition_prob=transition_upper)
    sequence_data = SIR.generate_states(current_state='S', nt=sample_size)
    with open('data/sir_middle_income.txt', 'w') as file:
        for item in sequence_data:
            file.writelines('{},'.format(item))

def read_data_lower_middle_SIR():
    with open('sim/01FluProgIncome/data/sir_lower_income.txt', 'r') as file:
        l = file.readlines()
    # print(l[0].split(','))
    l = l[0].split(',')
    l = clean_data(l)
    with open('sim/01FluProgIncome/data/sir_middle_income.txt', 'r') as file:
        m = file.readlines()
    # print(m[0].split(','))
    m = m[0].split(',')
    m = clean_data(m)
    return l,m

def generate_flu_income_rule_from_data():
    lower_seq_data, middle_seq_data = read_data_lower_middle_SIR()
    lower_rule_dict = msdd_algorithm_simple(lower_seq_data, precursor_width= 1, sucessor_width= 1, lagtime_between=1,
                                            dependency_evaluator_fn = None, num_rules = None)
    middle_rule_dict = msdd_algorithm_simple(middle_seq_data, precursor_width=1, sucessor_width=1, lagtime_between=1,
                                            dependency_evaluator_fn=None, num_rules=None)
    pram_rules = []
    prd = PramRulesArributes()
    inter_pram_rules = PramRulesArributes()
    for key, val in lower_rule_dict.items():
        prd = PramRulesArributes()
        prd.attribute = ['flu']
        prd.precursor = list(*key)
        for child_key, child_val in val.items():
            prd.add_new_successor(value=list(*child_key),
                                  successorAttribute=['flu'], probability=child_val)
        pram_rules.append(prd)
    pprint(pram_rules)
    rule_dict = dict(l = lower_rule_dict, m = middle_rule_dict)

    # pprint(rule_dict)
    # for rule_key, rule_val in rule_dict.items():
    #     generate_class_from_rule_dict_large(rule_val,rule_key)
    # generate_class_from_rule_dict_large(rule_dict)

if __name__ == '__main__':
    # create_income_sir_data()
    # lower, middle = read_data_lower_middle_SIR()
    # print(lower)
    # print(middle)
    # pass
    generate_flu_income_rule_from_data()