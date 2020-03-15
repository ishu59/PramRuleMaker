import os, sys
from itertools import combinations

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # 'rules' module
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from typing import List, Dict, Union, Tuple
from collections import defaultdict, Set
import numpy as np
from SequenceGenerator import SequenceGenerator
from pprint import pprint
from PramRulesArributes import PramRulesArributes, rules_attr_to_json
from ClassGenerator import PRAMClassGenerator


class RuleGenBase:

    def __init__(self, precursor:Union[List, Tuple, str], successor:Union[List, Tuple, str],
                 only_precursor_count: int, both_precursor_successor_count: int, only_successor_count: int,
                 total_data_count: int, precursor_attribute = None, successor_attribute = None):

        # precursor, successor, only_precursor_count, both_precursor_successor_count, only_successor_count, total_data_count, precursor_attribute

        self.total_data_count = total_data_count
        self.only_successor_count = only_successor_count
        self.both_precursor_successor_count = both_precursor_successor_count
        self.only_precursor_count = only_precursor_count
        self.precursor = precursor
        self.successor = successor
        self.successor_attribute = successor_attribute
        if self.successor_attribute is None:
            self.successor_attribute = ['attribute_'+str(i) for i in range(len(precursor))]
        self.precursor_attribute = precursor_attribute
        if self.precursor_attribute is None:
            self.precursor_attribute = ['attribute_'+str(i) for i in range(len(precursor))]
        self.precursor_str = ' '.join(precursor).strip()
        self.successor_str = ' '.join(successor).strip()
        self.probability = self.both_precursor_successor_count/ self.only_precursor_count
        self.g_score = self._get_g_score()

    def _get_g_score(self):
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

def read_clean_data(file_path):
    if file_path is None:
        raise Exception('File path has not be provided')
        file_path = 'data/sir_sex.txt'
    with open(file_path, 'r') as file:
        data = file.readlines()
    clean_data = []
    for i in range(len(data)):
        data[i] = data[i].rstrip()
        if len((data[i])) < 3:
            continue
        item = data[i]
        prec, succ = item.split(':')
        prec = tuple(prec.split(','))
        succ = tuple(succ.split(','))
        clean_data.append([prec, succ])
    # pprint(clean_data)
    p_s_token_dict = defaultdict(int)
    for prec, succ in clean_data:
        p_s_token_dict[(prec, succ, 1)] += 1
    print(p_s_token_dict)
    return clean_data, p_s_token_dict

def get_clean_data(sequence_data=None, precursor_width=None, successor_width=None, lagtime_between=None):
    '''
    :param sequence_data: path of file for data to be extracted
    :param precursor_width: Width of precursor
    :param successor_width: Width of Successor
    :param lagtime_between: time in between observations
    :return: List of List  [[precursor comma separated][successor comma separated],.....]
    '''
    rule_gen_list: List[RuleGenBase] = []
    if isinstance(sequence_data, str):
        with open(sequence_data, 'r') as f:
            seq = f.read()
        seq = seq.split(',')
    else:
        # successor_width = 2
        # precursor_width = 3
        # time_lag_between = 6
        successor_width = 1
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
    p_s_token_string_list = []
    for i in range(0,len(seq),token_width//2):
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
            temp_p_tok = ','.join(p_tok).strip()
        if len(s_tok) == 1:
            temp_s_tok = s_tok[0]
        else:
            temp_s_tok = ','.join(s_tok).strip()
        p_s_token_string_list.append(''.join([temp_p_tok, ':', temp_s_tok]))
        p_s_token_list.append([p_tok,s_tok] )
        p_s_token_dict[(p_tok, s_tok, token_width)] += 1
    total_token_count: int = 0
    for cnt in p_s_token_dict.values():
        total_token_count += cnt
    print(total_token_count)
    write_to_file: bool = True
    if write_to_file:
        pprint(p_s_token_string_list)
        with open('data/token_data.txt', 'w') as file:
            for item in p_s_token_string_list:
                file.write(item)
                file.write('\n')
    return p_s_token_list, p_s_token_dict

def msdd_algorithm_simple(sequence_data=None, precursor_width=None, successor_width=None, lagtime_between=None,
                          dependency_evaluator_fn=None, num_rules=None, data_type_vertical = False):
    '''
    :param sequence_data:  Path of Stream of data (sequence of categorical features) of format A,B,...(precursor):B,D...(successor)
    :param precursor_width: width of the Rule head
    :param successor_width: width
    :param lagtime_between: time between rule and effect
    :param dependency_evaluator_fn: fundtion to evaluate best rule
    :param num_rules: total Number of expected rules
    :return: rules list as dictionary {}
    '''

    rule_gen_list: List[RuleGenBase] = []
    sequence_data_path = None
    if isinstance(sequence_data, str) :
        sequence_data_path = sequence_data
    if data_type_vertical:
        p_s_token_list, p_s_token_dict = read_clean_data(file_path=sequence_data_path)
    else:
        p_s_token_list, p_s_token_dict = get_clean_data(sequence_data, precursor_width, successor_width, lagtime_between)

    total_token_count:int = 0
    # Contains all the key from all the combination in precursor
    all_precursor_combination_set = set()

    # Contains keys as  tuple of 1 : frozen set all combination of precursor for precursors, 2 successor and value as count of the combination
    all_combination_p_s_dict: Dict[Set, Tuple] = {}

    # contains keys as tuple of 1: individual keys from the set of all the keys 2: Successor and value as count
    ind_combination_precursors_successor_count_dict: Dict =defaultdict(int)

    # contains key : Precursor tuple, value as count of tuple
    ind_precursor_count_dict:Dict = defaultdict(int)

    # contains key : Successor tuple, value as count of tuple
    ind_successor_count_dict: Dict = defaultdict(int)

    for p_s_key, count_val in  p_s_token_dict.items():
        precursor_l, successor_l, width = p_s_key
        local_precursor_set = set()
        for i in range(1, len(precursor_l)+1):
            items = combinations(precursor_l, i)
            items = list(items)
            for item in items:
                all_precursor_combination_set.add(item)
                local_precursor_set.add(item)
        local_precursor_set =  frozenset(local_precursor_set)
        new_key_p_s = tuple([local_precursor_set, successor_l])
        all_combination_p_s_dict[new_key_p_s] = count_val

    token_combination_total_count = 0
    total_successor_count = 0
    total_precursor_count = 0

    for unique_prec_key in all_precursor_combination_set:
        for combined_p_s_key, count_val in all_combination_p_s_dict.items():
            precursor_m, successor_m = combined_p_s_key
            if unique_prec_key in precursor_m:
                new_ind_p_s_key = tuple([unique_prec_key,successor_m])
                ind_combination_precursors_successor_count_dict[new_ind_p_s_key] += count_val
                ind_precursor_count_dict[unique_prec_key] += count_val
                token_combination_total_count += count_val
                ind_successor_count_dict[successor_m] += count_val

    print('token_combination_total_count = ',token_combination_total_count)

    for _, count_val in ind_successor_count_dict.items():
        total_successor_count += count_val
    print('total_successor_count : ', total_successor_count)

    for _, count_val in ind_precursor_count_dict.items():
        total_precursor_count += count_val
    print('total_precursor_count :  ', total_precursor_count)

    for p_s_key, combined_count_val in ind_combination_precursors_successor_count_dict.items():
        prec, succ = p_s_key
        total_prec_count = ind_precursor_count_dict[prec]
        total_succ_count = ind_successor_count_dict[succ]
        rule_gen_list.append(RuleGenBase(precursor = prec, successor = succ, only_precursor_count = total_prec_count,
                                         both_precursor_successor_count = combined_count_val,
                                         only_successor_count =  total_succ_count,
                                         total_data_count = token_combination_total_count))

    rule_gen_list.sort(key=lambda x: x.g_score, reverse=True)
    if num_rules is None:
        num_rules = len(rule_gen_list)
    num_rules = min(num_rules, len(rule_gen_list))
    rule_gen_list = rule_gen_list[:num_rules]
    # for r in rule_gen_list:
    #     pprint(r)

    final_rule_gen_list:Dict[Tuple[str],List[RuleGenBase]] = defaultdict(list)
    pram_rule_attr_list:List[PramRulesArributes] = []

    for r in rule_gen_list:
        final_rule_gen_list[r.precursor].append(r)
    for r_prec, base_rule_list in final_rule_gen_list.items():
        pram_rule_attr:PramRulesArributes = PramRulesArributes(attribute=base_rule_list[0].precursor_attribute,
                                                               precursor=base_rule_list[0].precursor)
        for base_rule in base_rule_list:
            pram_rule_attr.add_new_successor(value=base_rule.successor,
                                             probability=base_rule.probability,
                                             successorAttribute=base_rule.successor_attribute)
        pram_rule_attr_list.append(pram_rule_attr)
    rule_json_value = rules_attr_to_json(pram_rule_attr_list)
    print(rule_json_value)
    return rule_json_value

def test_msdd_with_vertical_data():
    data_path = 'data2/sir_sex.txt'
    json_save_path = 'jsonData/new_rule_from_vert.json'
    output_file_path = './sim2'
    output_file_name = 'SIRGenderVertical'
    num_rules = 9
    rule_dict = msdd_algorithm_simple(sequence_data=data_path, data_type_vertical=True, num_rules=num_rules)
    print(rule_dict)
    with open(json_save_path, 'w') as file:
        file.write(rule_dict)
    pcg = PRAMClassGenerator(json_file_path=json_save_path)
    pcg.generate_class_from_json(data=None, output_file=output_file_name, output_file_path=output_file_path)


def test_msdd_simple():
    # data_path, json_save_path, output_file_path, output_file_name
    data_path = 'data2/sir_gender_flat.txt'
    json_save_path = 'jsonData/new_rule_from_flat.json'
    output_file_path = './sim2'
    output_file_name = 'SIRGenderFlat'
    num_rules = 9
    rule_dict = msdd_algorithm_simple(sequence_data=data_path,
                                      precursor_width=2, successor_width=2,lagtime_between=2,
                                      data_type_vertical=False, num_rules = num_rules)
    with open(json_save_path,'w') as file:
        file.write(rule_dict)
    pcg = PRAMClassGenerator(json_file_path=json_save_path)
    pcg.generate_class_from_json(data=None, output_file=output_file_name, output_file_path=output_file_path)




if __name__ == '__main__':
    print('Starting Task')
    test_msdd_with_vertical_data()
    # test_msdd_simple()
    print('Task completed successfully')


##############################################################################################################
# End of Program
##############################################################################################################






















#
# def test_msdd_simple2():
#     transition = {'S': {'S': 0.7, 'I': 0.3, 'R': 0}, 'I': {'S': 0, 'I': 0.5, 'R': 0.5},
#                   'R': {'S': 0.3, 'I': 0.0, 'R': 0.7}}
#     SIR = SequenceGenerator(transition_prob=transition)
#     sequence_data = SIR.generate_states(current_state='S', nt=2000)
#     rule_dict = msdd_algorithm_simple(sequence_data, precursor_width=1, successor_width=1, lagtime_between=1,
#                                       dependency_evaluator_fn=None, num_rules=None)
#     pprint(rule_dict)
#     # generate_class_from_rule_dict_large(rule_dict)
#
#
# def test_read_clean_data():
#     clean_data, p_s_token_dict = read_clean_data()
#     pprint(clean_data)
#     pprint(p_s_token_dict)




# #
# # if __name__ == '__main__':
# #     print('Starting Task')
# #     # test_msdd_with_vertical_data()
# #     test_msdd_simple()
# #     print('Task completed successfully')
#
#
#
#     # test_read_clean_data
#     # test_msdd_simple2()
#     # algorithm_multi_sequence_data()
#     # test_algorithm_multi_seq_data()
#     # create_income_sir_data()
#     # read_data_lower_middle_SIR()



# def clean_data(l: List):
#     return [x for x in l if x not in ['', ' ', ',', '\n']]
