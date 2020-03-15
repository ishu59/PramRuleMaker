import os, sys
from itertools import combinations
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from typing import List, Dict, Union, Tuple
from collections import defaultdict, OrderedDict, Set
import numpy as np
from SequenceGenerator import SequenceGenerator
from pprint import pprint
from PramRulesArributes import PramRulesArributes


class RuleGenBase:

    def __init__(self, precursor:Union[List, Tuple, str], successor:Union[List, Tuple, str],
                 only_precursor_count: int, both_precursor_successor_count: int, only_successor_count: int,
                 total_data_count: int, attribute = None):
        self.attribute = attribute
        if self.attribute is None:
            self.attribute = 'default_attribute'
        self.total_data_count = total_data_count
        self.only_successor_count = only_successor_count
        self.both_precursor_successor_count = both_precursor_successor_count
        self.only_precursor_count = only_precursor_count
        self.precursor = precursor
        self.successor = successor
        self.precursor_str = ' '.join(precursor).strip()
        self.successor_str = ' '.join(successor).strip()
        self.probability = self.both_precursor_successor_count/ self.only_precursor_count
        self.g_score = self.get_g_score()

    def get_g_score(self):
        n1_x_y = self.both_precursor_successor_count
        n2_x_not_y = self.only_precursor_count - self.both_precursor_successor_count if self.only_precursor_count - self.both_precursor_successor_count > 0 else 0
        n3_not_x_y = self.only_successor_count - self.both_precursor_successor_count if self.only_successor_count - self.both_precursor_successor_count > 0 else 0
        n4_not_x_not_y = self.total_data_count + self.both_precursor_successor_count - (
                    self.only_successor_count + self.only_precursor_count)
        g_score = self._compute_g_score(n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y)
        return g_score

    def _compute_g_score(self, n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y):
        total = n1_x_y + n2_x_not_y + n3_not_x_y + n4_not_x_not_y
        if total == 0:
            return 0
        n1_hat = (n1_x_y + n2_x_not_y) * (n1_x_y + n3_not_x_y) / total
        n2_hat = (n1_x_y + n2_x_not_y) * (n2_x_not_y + n4_not_x_not_y) / total
        n3_hat = (n3_not_x_y + n4_not_x_not_y) * (n1_x_y + n3_not_x_y) / total
        n4_hat = (n3_not_x_y + n4_not_x_not_y) * (n2_x_not_y + n4_not_x_not_y) / total
        g_score = 2 * (n1_x_y * np.log(n1_x_y / n1_hat) +
                       (n2_x_not_y * np.log(n2_x_not_y / n2_hat)) +
                       (n3_not_x_y * np.log(n3_not_x_y / n3_hat))
                       + (n4_not_x_not_y * np.log(n4_not_x_not_y / n4_hat)))
        return g_score

    def __repr__(self):
        all_val = str(self.__dict__)
        return all_val


def test_RuleGenBase():
    r = RuleGenBase(successor='s',precursor='p',only_precursor_count=10,only_successor_count=12, both_precursor_successor_count=6,total_data_count=100)
    print(r)

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
    n1_hat = (n1_x_y + n2_x_not_y) * (n1_x_y + n3_not_x_y) / total
    n2_hat = (n1_x_y + n2_x_not_y) * (n2_x_not_y + n4_not_x_not_y) / total
    n3_hat = (n3_not_x_y + n4_not_x_not_y) * (n1_x_y + n3_not_x_y) / total
    n4_hat = (n3_not_x_y + n4_not_x_not_y) * (n2_x_not_y + n4_not_x_not_y) / total
    g_score = 2 * (n1_x_y * np.log(n1_x_y / n1_hat) +
                   (n2_x_not_y * np.log(n2_x_not_y / n2_hat)) +
                   (n3_not_x_y * np.log(n3_not_x_y / n3_hat))
                   + (n4_not_x_not_y * np.log(n4_not_x_not_y / n4_hat)))
    return g_score


def generate_n_tokens(seq: List, width):
    n_token = list(zip(*[seq[i:] for i in range(width)]))
    return list((set(n_token)))


def msdd_algorithm_simple(sequence_data=None, precursor_width=None, successor_width=None, lagtime_between=None,
                          dependency_evaluator_fn=None, num_rules=None):
    '''
    :param sequence_data: Stream of data (sequence of categorical features)
    :param precursor_width: width of the Rule head
    :param successor_width: width
    :param lagtime_between: time between rule and effect
    :param dependency_evaluator_fn: fundtion to evaluate best rule
    :param num_rules: total Number of expected rules
    :return: rules list as dictionary {}
    '''
    rule_gen_list:List[RuleGenBase] = []

    if sequence_data is not None:
        seq = sequence_data
    else:
        # successor_width = 2
        # precursor_width = 3
        # time_lag_between = 6
        successor_width = 3
        precursor_width = 1
        lagtime_between = 1
        seq = ['S', 'S', 'R', 'I', 'I', 'I', 'R', 'S', 'S', 'S', 'S', 'I', 'R', 'S', 'I', 'I', 'I', 'I', 'R', 'I', 'I',
               'I', 'I', 'R', 'R', 'R', 'R', 'R', 'S', 'I', 'R', 'R', 'R', 'I', 'I', 'I', 'I', 'I', 'R', 'R', 'R', 'R',
               'R', 'R', 'R', 'R', 'R', 'I', ]
        seq = get_sample_seq()
        num_rules = None

    # Rule generation Code Begins
    token_width = lagtime_between if lagtime_between > (
            successor_width + precursor_width) else successor_width + precursor_width

    # Extracting tokens from sequence
    # p_s_token_dict: key is tuple of precursor , successor, width and value is count of the combination
    p_s_token_dict = defaultdict(int)
    p_s_token_list = []
    for i in range(len(seq)):
        if i + token_width >= len(seq):
            break
        p_tok = tuple(seq[i:(i + precursor_width)])
        s_tok = tuple(seq[i + token_width - successor_width: i + token_width])

        # Extra for downstream tasks like writing tokens to file etc
        temp_p_tok = p_tok
        temp_s_tok = s_tok
        if len(p_tok) == 1:
            temp_p_tok = p_tok[0]
        else:
            temp_p_tok = ' '.join(p_tok).strip()
        if len(s_tok) == 1:
            temp_s_tok = s_tok[0]
        else:
            temp_s_tok = ' '.join(s_tok).strip()
        p_s_token_list.append(''.join([temp_p_tok, ' : ', temp_s_tok]))
        p_s_token_dict[(p_tok, s_tok, token_width)] += 1

    total_token_count:int = 0


    for cnt in p_s_token_dict.values():
        total_token_count += cnt
    print(total_token_count)
    write_to_file: bool = False
    if write_to_file:
        pprint(p_s_token_list)
        with open('data/token_data.txt', 'w') as file:
            for item in p_s_token_list:
                file.write(item)
                file.write('\n')


    all_precursor_combination_set = set()
    all_combination_precursors_dict: Dict ={}
    all_combination_p_s_dict: Dict[Set, Tuple] = {}

    for p_s_key, count_val in  p_s_token_dict.items():
        precursor, successor, width = p_s_key
        local_precursor_set = set()
        for i in range(1, len(precursor)+1):
            items = combinations(precursor, i)
            items = list(items)
            for item in items:
                all_precursor_combination_set.add(item)
                local_precursor_set.add(item)

        local_precursor_set =  frozenset(local_precursor_set)
        new_key_p_s = tuple(local_precursor_set, successor)
        all_combination_p_s_dict[new_key_p_s] = count_val
    pprint(all_precursor_combination_set)
    pprint(all_combination_p_s_dict)


    # Extracting only successor tokens from the sequence
    s_token_dict = defaultdict(int)
    for i in range(len(seq)):
        if i + successor_width >= len(seq):
            break
        s_tok = tuple(seq[i:i + successor_width])
        s_token_dict[s_tok] += 1

    # Extracting only precursor tokens from the sequence
    p_token_dict = defaultdict(int)
    for i in range(len(seq)):
        if i + token_width >= len(seq):
            break
        p_tok = tuple(seq[i:i + precursor_width])
        p_token_dict[p_tok] += 1
    p_s_token_prob_dict = defaultdict(int)

    # Probability of Token combination appearing in the data
    for p_s_key in p_s_token_dict.keys():
        p, s, width = p_s_key
        p_s_token_prob_dict[p_s_key] = p_s_token_dict[p_s_key] / p_token_dict[p]

    #  Count that token didnot appear in sequence
    p_s_not_occur_dict = defaultdict(int)
    for tok_key in p_s_token_dict.keys():
        p_token_key, s_token_key, _ = tok_key
        for i in range(len(seq)):
            p_tok_seq = tuple(seq[i:(i + precursor_width)])
            s_tok_seq = tuple(seq[i + token_width - successor_width: i + token_width])
            if p_token_key != p_tok_seq and s_token_key != s_tok_seq:
                p_s_not_occur_dict[tok_key] += 1

    #  Completed using easier way but let it be for now
    g_score_p_s_dict = {}
    for token_key, token_value in p_s_token_dict.items():
        p_token, s_token, _ = token_key
        n1_x_y = token_value
        n2_x_not_y = p_token_dict[p_token] - token_value
        n3_not_x_y = s_token_dict[s_token] - token_value
        n4_not_x_not_y = p_s_not_occur_dict[token_key]
        rule_gen_list.append(RuleGenBase(precursor=p_token, successor=s_token,
                                         only_precursor_count=p_token_dict[p_token],
                                         only_successor_count=s_token_dict[s_token],
                                         both_precursor_successor_count=token_value,
                                         total_data_count=total_token_count))

        g_score_p_s_dict[token_key] = compute_g_score(n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y)
    g_score_p_s_dict_ordered = OrderedDict(sorted(g_score_p_s_dict.items(), key=lambda x: x[1], reverse=True))
    for item in rule_gen_list:
        pprint(item)

    # pprint()
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


def test_msdd_simple():
    rule_dict = msdd_algorithm_simple()
    print(rule_dict)
    # generate_class_from_rule_dict(rule_dict)


def test_msdd_simple2():
    transition = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
                  'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
    SIR = SequenceGenerator(transition_prob=transition)
    sequence_data = SIR.generate_states(current_state='S', nt=2000)
    rule_dict = msdd_algorithm_simple(sequence_data, precursor_width=1, successor_width=1, lagtime_between=1,
                                      dependency_evaluator_fn=None, num_rules=None)
    pprint(rule_dict)
    # generate_class_from_rule_dict_large(rule_dict)


def clean_data(l: List):
    return [x for x in l if x not in ['', ' ', ',', '\n']]


if __name__ == '__main__':
    # test_RuleGenBase()
    test_msdd_simple()
    # test_msdd_simple2()
    # algorithm_multi_sequence_data()
    # test_algorithm_multi_seq_data()
    # create_income_sir_data()
    # read_data_lower_middle_SIR()


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


def generate_class_from_rule_dict_large(rule_dict, class_index=1, has_attr='flu', set_attr='flu'):
    import_str = "from pram.data   import GroupSizeProbe, ProbeMsgMode\nfrom pram.entity import Group, GroupQry, GroupSplitSpec, Site\nfrom pram.rule   import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule\nfrom pram.sim    import Simulation\n\n\n"
    rule_string = "class Autogenerated{}(Rule):\n\tdef apply(self, pop, group, iter, t):\n\t\t".format(class_index)
    rule_string = import_str + rule_string
    condition_string = "\t\tif group.has_attr({'" + has_attr + "': "
    condition_string_2 = "}): return ["
    g_spec_str_1 = "GroupSplitSpec(p= "
    g_spec_str_2 = ", attr_set={ '" + set_attr + "': "
    g_spec_str_3 = " }),"
    condition_string_3 = "]\n"

    for upper_index, (key, val) in enumerate(rule_dict.items()):
        if len(key) < 2:
            rs = condition_string + "'" + key[0] + "'" + condition_string_2
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
                rss = g_spec_str_1 + str(p) + g_spec_str_2 + "'" + child_key[0] + "'" + g_spec_str_3
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


def get_pram_attribute_from_rule(rules_dict: Dict):
    pram_rule_attr = PramRulesArributes()



    # # Extracting only successor tokens from the sequence
    # s_token_dict = defaultdict(int)
    # for i in range(len(seq)):
    #     if i + successor_width >= len(seq):
    #         break
    #     s_tok = tuple(seq[i:i + successor_width])
    #     s_token_dict[s_tok] += 1
    #
    # # Extracting only precursor tokens from the sequence
    # p_token_dict = defaultdict(int)
    # for i in range(len(seq)):
    #     if i + token_width >= len(seq):
    #         break
    #     p_tok = tuple(seq[i:i + precursor_width])
    #     p_token_dict[p_tok] += 1
    # p_s_token_prob_dict = defaultdict(int)
    #
    # # Probability of Token combination appearing in the data
    # for p_s_key in p_s_token_dict.keys():
    #     p, s, width = p_s_key
    #     p_s_token_prob_dict[p_s_key] = p_s_token_dict[p_s_key] / p_token_dict[p]
    #
    # #  Count that token didnot appear in sequence
    # p_s_not_occur_dict = defaultdict(int)
    # for tok_key in p_s_token_dict.keys():
    #     p_token_key, s_token_key, _ = tok_key
    #     for i in range(len(seq)):
    #         p_tok_seq = tuple(seq[i:(i + precursor_width)])
    #         s_tok_seq = tuple(seq[i + token_width - successor_width: i + token_width])
    #         if p_token_key != p_tok_seq and s_token_key != s_tok_seq:
    #             p_s_not_occur_dict[tok_key] += 1
    #
    # #  Completed using easier way but let it be for now
    # g_score_p_s_dict = {}
    # for token_key, token_value in p_s_token_dict.items():
    #     p_token, s_token, _ = token_key
    #     n1_x_y = token_value
    #     n2_x_not_y = p_token_dict[p_token] - token_value
    #     n3_not_x_y = s_token_dict[s_token] - token_value
    #     n4_not_x_not_y = p_s_not_occur_dict[token_key]
    #     rule_gen_list.append(RuleGenBase(precursor=p_token, successor=s_token,
    #                                      only_precursor_count=p_token_dict[p_token],
    #                                      only_successor_count=s_token_dict[s_token],
    #                                      both_precursor_successor_count=token_value,
    #                                      total_data_count=total_token_count))
    #
    #     g_score_p_s_dict[token_key] = compute_g_score(n1_x_y, n2_x_not_y, n3_not_x_y, n4_not_x_not_y)
    # g_score_p_s_dict_ordered = OrderedDict(sorted(g_score_p_s_dict.items(), key=lambda x: x[1], reverse=True))
    # for item in rule_gen_list:
    #     pprint(item)
    #
    # # pprint()
    # if num_rules is None:
    #     num_rules = len(p_s_token_prob_dict)
    # unstructured_rule_dict = {}
    # for i, token_key in enumerate(g_score_p_s_dict_ordered.keys()):
    #     if i > num_rules:
    #         break
    #     else:
    #         rule_left, rule_right, _ = token_key
    #         rule_prob = p_s_token_prob_dict[token_key]
    #         unstructured_rule_dict[(rule_left, rule_right)] = rule_prob
    # rule_dict = {}
    # for token_key, prob in unstructured_rule_dict.items():
    #     rule_left, rule_right = token_key
    #     if rule_dict.get(rule_left) is None:
    #         rule_dict[rule_left] = {}
    #     rule_dict[rule_left][rule_right] = prob
    # # print(rule_dict)
    # return rule_dict
    # return None

    # Upper part
    # if sequence_data is not None:
    #     seq = sequence_data
    # else:
    #     # successor_width = 2
    #     # precursor_width = 3
    #     # time_lag_between = 6
    #     successor_width = 2
    #     precursor_width = 2
    #     lagtime_between = 1
    #     seq = ['S', 'S', 'R', 'I', 'I', 'I', 'R', 'S', 'S', 'S', 'S', 'I', 'R', 'S', 'I', 'I', 'I', 'I', 'R', 'I', 'I',
    #            'I', 'I', 'R', 'R', 'R', 'R', 'R', 'S', 'I', 'R', 'R', 'R', 'I', 'I', 'I', 'I', 'I', 'R', 'R', 'R', 'R',
    #            'R', 'R', 'R', 'R', 'R', 'I', ]
    #     seq = get_sample_seq()
    #     num_rules = None
    #
    # # Rule generation Code Begins
    # token_width = lagtime_between if lagtime_between > (
    #         successor_width + precursor_width) else successor_width + precursor_width
    #
    # # Extracting tokens from sequence
    # # p_s_token_dict: key is tuple of precursor , successor, width and value is count of the combination
    # p_s_token_dict = defaultdict(int)
    # p_s_token_list = []
    # for i in range(len(seq)):
    #     if i + token_width >= len(seq):
    #         break
    #     p_tok = tuple(seq[i:(i + precursor_width)])
    #     s_tok = tuple(seq[i + token_width - successor_width: i + token_width])
    #
    #     # Extra for downstream tasks like writing tokens to file etc
    #     temp_p_tok = p_tok
    #     temp_s_tok = s_tok
    #     if len(p_tok) == 1:
    #         temp_p_tok = p_tok[0]
    #     else:
    #         temp_p_tok = ','.join(p_tok).strip()
    #     if len(s_tok) == 1:
    #         temp_s_tok = s_tok[0]
    #     else:
    #         temp_s_tok = ','.join(s_tok).strip()
    #     p_s_token_list.append(''.join([temp_p_tok, ' : ', temp_s_tok]))
    #     p_s_token_dict[(p_tok, s_tok, token_width)] += 1
    # for cnt in p_s_token_dict.values():
    #     total_token_count += cnt
    # print(total_token_count)
    # write_to_file: bool = True
    # if write_to_file:
    #     pprint(p_s_token_list)
    #     with open('data/token_data.txt', 'w') as file:
    #         for item in p_s_token_list:
    #             file.write(item)
    #             file.write('\n')
    #
