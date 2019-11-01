import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from typing import List, Dict, Union, Optional
from collections import defaultdict, OrderedDict
import numpy as np
import scipy as sp
from SequenceGenerator import SequenceGenerator
from pprint import pprint
from DependentSequenceGenerator import HiddenModel
from PramRulesArributes import PramRulesArributes

def get_sample_seq():
    transition = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                       'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    SIR = SequenceGenerator(transition_prob=transition)
    sequence_data = SIR.generate_states(current_state='S', nt=2000)
    return sequence_data

def compute_g_score(n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y):
    total = n1_x_y + n2_x_not_y + n3_not_x_y + n4_not_x_not_y
    if total == 0:
        return 0
    n1_hat = (n1_x_y + n2_x_not_y)*(n1_x_y+n3_not_x_y) / total
    n2_hat = (n1_x_y + n2_x_not_y)*(n2_x_not_y + n4_not_x_not_y) / total
    n3_hat = (n3_not_x_y + n4_not_x_not_y) * (n1_x_y + n3_not_x_y) / total
    n4_hat = (n3_not_x_y + n4_not_x_not_y) * (n2_x_not_y + n4_not_x_not_y) / total
    g_score = 2 * (n1_x_y * np.log(n1_x_y / n1_hat) +
                   (n2_x_not_y * np.log(n2_x_not_y / n2_hat)) +
                   (n3_not_x_y * np.log(n3_not_x_y / n3_hat))
                   + (n4_not_x_not_y * np.log(n4_not_x_not_y / n4_hat)))
    return g_score

def generate_n_tokens(seq:List, width):
    n_token = list(zip(*[seq[i:] for i in range(width)]))
    return list((set(n_token)))

def msdd_algorithm_simple(sequence_data=None, precursor_width=None, sucessor_width=None, lagtime_between=None, dependency_evaluator_fn = None, num_rules = None):
    '''
    :param sequence_data: Stream of data (sequence of categorical features)
    :param precursor_width: width of the Rule head
    :param sucessor_width: width
    :param lagtime_between: time between rule and effect
    :param dependency_evaluator_fn: fundtion to evaluate best rule
    :param num_rules: total Number of expected rules
    :return: rules list as dictionary {}
    '''
    if sequence_data is not None:
        seq = sequence_data
    else:
        # sucessor_width = 2
        # precursor_width = 3
        # time_lag_between = 6
        sucessor_width = 1
        precursor_width = 1
        lagtime_between = 1
        seq = ['S', 'S', 'R', 'I', 'I', 'I', 'R', 'S', 'S', 'S', 'S', 'I', 'R', 'S', 'I', 'I', 'I', 'I', 'R', 'I', 'I', 'I',
               'I', 'R', 'R', 'R', 'R', 'R', 'S', 'I', 'R', 'R', 'R', 'I', 'I', 'I', 'I', 'I', 'R', 'R', 'R', 'R', 'R', 'R',
               'R', 'R', 'R', 'I', ]
        seq = get_sample_seq()
        num_rules = None

    # Rule generation Code Begins
    token_width = lagtime_between if lagtime_between > (
            sucessor_width + precursor_width) else sucessor_width + precursor_width
    p_s_token_dict = defaultdict(int)
    for i in range(len(seq)):
        if i + token_width >= len(seq):
            break
        p_tok = tuple(seq[i:(i + precursor_width)])
        s_tok = tuple(seq[i + token_width - sucessor_width: i + token_width])
        p_s_token_dict[(p_tok, s_tok, token_width)] += 1
    s_token_dict = defaultdict(int)
    for i in range(len(seq)):
        if i + sucessor_width >= len(seq):
            break
        s_tok = tuple(seq[i:i + sucessor_width])
        s_token_dict[s_tok] += 1
    p_token_dict = defaultdict(int)
    for i in range(len(seq)):
        if i + precursor_width >= len(seq):
            break
        p_tok = tuple(seq[i:i + precursor_width])
        p_token_dict[p_tok] += 1
    p_s_token_prob_dict = defaultdict(int)
    for p_s_key in p_s_token_dict.keys():
        p, s, width = p_s_key
        p_s_token_prob_dict[p_s_key] = p_s_token_dict[p_s_key] / p_token_dict[p]
    p_s_not_occur_dict = defaultdict(int)
    for tok_key in p_s_token_dict.keys():
        p_token_key, s_token_key, _ = tok_key
        for i in range(len(seq)):
            p_tok_seq = tuple(seq[i:(i + precursor_width)])
            s_tok_seq = tuple(seq[i + token_width - sucessor_width: i + token_width])
            if p_token_key != p_tok_seq and s_token_key != s_tok_seq:
                p_s_not_occur_dict[tok_key] += 1
    g_score_p_s_dict = {}
    for token_key, token_value in p_s_token_dict.items():
        p_token, s_token, _ = token_key
        n1_x_y = token_value
        n2_x_not_y = p_token_dict[p_token] - token_value
        n3_not_x_y = s_token_dict[s_token] - token_value
        n4_not_x_not_y = p_s_not_occur_dict[token_key]
        g_score_p_s_dict[token_key] = compute_g_score(n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y)
    g_score_p_s_dict_ordered = OrderedDict(sorted(g_score_p_s_dict.items(), key=lambda x: x[1], reverse=True))
    # pprint(g_score_p_s_dict)
    # print('GScores')
    # pprint(g_score_p_s_dict_ordered)
    if num_rules is None:
        num_rules = len(p_s_token_prob_dict)
    unstructured_rule_dict = {}
    for i, token_key in enumerate(g_score_p_s_dict_ordered.keys()):
        if i > num_rules:
            break
        else:
            rule_left, rule_right, _ = token_key
            rule_prob = p_s_token_prob_dict[token_key]
            unstructured_rule_dict[(rule_left, rule_right)] = rule_prob
    rule_dict = {}
    for token_key, prob in unstructured_rule_dict.items():
        rule_left, rule_right = token_key
        if rule_dict.get(rule_left) is None:
            rule_dict[rule_left] = {}
        rule_dict[rule_left][rule_right] = prob
    # print(rule_dict)
    return rule_dict

def get_sample_seq2():
    transition = {'S': {'S': 0.9, 'I': 0.1, 'R': 0}, 'I': {'S': 0, 'I': 0.8, 'R': 0.2},
                       'R': {'S': 0.1, 'I': 0.0, 'R': 0.9}}
    SIR = SequenceGenerator(transition_prob=transition)
    sequence_data = SIR.generate_states(current_state='S', nt=2000)
    return sequence_data

def generate_class_from_rule_dict(rule_dict):
    import_str = "from pram.data   import GroupSizeProbe, ProbeMsgMode\nfrom pram.entity import Group, GroupQry, GroupSplitSpec, Site\nfrom pram.rule   import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule\nfrom pram.sim    import Simulation\n\n\n"
    rule_string = "class Autogenerated(Rule):\n\tdef apply(self, pop, group, iter, t):\n\t\t"
    rule_string = import_str + rule_string
    condition_string = "\t\tif group.has_attr({'flu': '"
    condition_string_2 = "'}): return ["
    g_spec_str_1 = "GroupSplitSpec(p= "
    g_spec_str_2 = ", attr_set={ 'flu': '"
    g_spec_str_3 = "' }),"
    condition_string_3 = "]\n"

    for upper_index, (key, val) in enumerate(rule_dict.items()):
        rs = condition_string + key[0] + condition_string_2
        for i, (child_key, child_val) in enumerate(val.items()):
            p_minus = 0
            if i > len(val) - 1:
                p = p_minus
            else:
                p = child_val
                p_minus += p
            rss = g_spec_str_1 + str(p) + g_spec_str_2 + child_key[0] + g_spec_str_3
            rs += rss
        rs += condition_string_3
        rule_string = rule_string + '\n' + rs

    print(rule_string)
    with open('Autogenerated.py', 'w') as file:
        file.writelines(rule_string)

def generate_class_from_rule_dict_large(rule_dict,class_index = 1, has_attr = 'flu', set_attr = 'flu'):
    import_str = "from pram.data   import GroupSizeProbe, ProbeMsgMode\nfrom pram.entity import Group, GroupQry, GroupSplitSpec, Site\nfrom pram.rule   import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule\nfrom pram.sim    import Simulation\n\n\n"
    rule_string = "class Autogenerated{}(Rule):\n\tdef apply(self, pop, group, iter, t):\n\t\t".format(class_index)
    rule_string = import_str + rule_string
    condition_string = "\t\tif group.has_attr({'"+has_attr +"': "
    condition_string_2 = "}): return ["
    g_spec_str_1 = "GroupSplitSpec(p= "
    g_spec_str_2 = ", attr_set={ '"+set_attr+"': "
    g_spec_str_3 = " }),"
    condition_string_3 = "]\n"

    for upper_index, (key, val) in enumerate(rule_dict.items()):
        if len(key) < 2:
            rs = condition_string +"'" + key[0] + "'" + condition_string_2
        else:
            rs = condition_string + str(key) + condition_string_2

        for i, (child_key, child_val) in enumerate(val.items()):
            p_minus = 0
            if i > len(val) - 1:
                p = p_minus
            else:
                p = child_val
                p_minus += p
            print(type(child_key))
            if len(child_key) < 2:
                rss = g_spec_str_1 + str(p) + g_spec_str_2 +"'"+ child_key[0] +"'"+ g_spec_str_3
            elif type(child_key) == 'str' or type(child_key == 'numpy.str_'):
                rss = g_spec_str_1 + str(p) + g_spec_str_2 + "'" + child_key + "'" + g_spec_str_3
            else:
                rss = g_spec_str_1 + str(p) + g_spec_str_2 + str(child_key) + g_spec_str_3
            rs += rss
        rs += condition_string_3
        rule_string = rule_string + '\n' + rs

    print(rule_string)
    new_file_name = 'Autogenerated{}.py'.format(class_index)
    with open(new_file_name, 'w') as file:
        file.writelines(rule_string)

def test_msdd_simple():
    rule_dict = msdd_algorithm_simple()
    generate_class_from_rule_dict(rule_dict)

def get_pram_attribute_from_rule(rules_dict:Dict):
    pram_rule_attr = PramRulesArributes()


def test_msdd_simple2():
    transition = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                  'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    SIR = SequenceGenerator(transition_prob=transition)
    sequence_data = SIR.generate_states(current_state='S', nt=2000)
    rule_dict = msdd_algorithm_simple(sequence_data, precursor_width= 1, sucessor_width= 1, lagtime_between=1, dependency_evaluator_fn = None, num_rules = None)
    pprint(rule_dict)
    # generate_class_from_rule_dict_large(rule_dict)

def clean_data(l:List):
    return [x for x in l if x not in ['',' ',',','\n']]

if __name__ == '__main__':
    # test_msdd_simple()
    test_msdd_simple2()
    # algorithm_multi_sequence_data()
    # test_algorithm_multi_seq_data()
    # create_income_sir_data()
    # read_data_lower_middle_SIR()