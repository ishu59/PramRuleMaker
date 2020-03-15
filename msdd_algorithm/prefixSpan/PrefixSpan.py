import os
import json
from pprint import pprint
from collections import defaultdict
from typing import List, Dict, Union

from msdd_algorithm.prefixSpan.Reader import Reader

class PrefixSpan:

    def __init__(self, sequences, minSupport=0.1, maxPatternLength=10):

        minSupport = minSupport * len(sequences)
        self.PLACE_HOLDER = '_'

        freqSequences = self._prefixSpan(
            self.SequencePattern([], None, maxPatternLength, self.PLACE_HOLDER),
            sequences, minSupport, maxPatternLength)

        self.freqSeqs = PrefixSpan.FreqSequences(freqSequences)

    @staticmethod
    def train(sequences, minSupport=0.1, maxPatternLength=10):
        return PrefixSpan(sequences, minSupport, maxPatternLength)

    def freqSequences(self):
        return self.freqSeqs

    class FreqSequences:
        def __init__(self, fs):
            self.fs = fs

        def collect(self):
            return self.fs

    class SequencePattern:
        def __init__(self, sequence, support, maxPatternLength, place_holder):
            self.place_holder = place_holder
            self.sequence = []
            for s in sequence:
                self.sequence.append(list(s))
            self.freq = support

        def append(self, p):
            if p.sequence[0][0] == self.place_holder:
                first_e = p.sequence[0]
                first_e.remove(self.place_holder)
                self.sequence[-1].extend(first_e)
                self.sequence.extend(p.sequence[1:])
            else:
                self.sequence.extend(p.sequence)
                if self.freq is None:
                    self.freq = p.freq
            self.freq = min(self.freq, p.freq)

    def _checkPatternLengths(self, pattern, maxPatternLength):
        for s in pattern.sequence:
            if len(s) > maxPatternLength:
                return False
        return True

    def _prefixSpan(self, pattern, S, threshold, maxPatternLength):
        patterns = []

        if self._checkPatternLengths(pattern, maxPatternLength):
            f_list = self._frequent_items(S, pattern, threshold, maxPatternLength)

            for i in f_list:
                p = self.SequencePattern(pattern.sequence, pattern.freq, maxPatternLength, self.PLACE_HOLDER)
                p.append(i)
                if self._checkPatternLengths(pattern, maxPatternLength):
                    patterns.append(p)

                p_S = self._build_projected_database(S, p)
                p_patterns = self._prefixSpan(p, p_S, threshold, maxPatternLength)
                patterns.extend(p_patterns)

        return patterns

    def _frequent_items(self, S, pattern, threshold, maxPatternLength):
        items = {}
        _items = {}
        f_list = []
        if S is None or len(S) == 0:
            return []

        if len(pattern.sequence) != 0:
            last_e = pattern.sequence[-1]
        else:
            last_e = []
        for s in S:

            # class 1
            is_prefix = True
            for item in last_e:
                if item not in s[0]:
                    is_prefix = False
                    break
            if is_prefix and len(last_e) > 0:
                index = s[0].index(last_e[-1])
                if index < len(s[0]) - 1:
                    for item in s[0][index + 1:]:
                        if item in _items:
                            _items[item] += 1
                        else:
                            _items[item] = 1

            # class 2
            if self.PLACE_HOLDER in s[0]:
                for item in s[0][1:]:
                    if item in _items:
                        _items[item] += 1
                    else:
                        _items[item] = 1
                s = s[1:]

            # class 3
            counted = []
            for element in s:
                for item in element:
                    if item not in counted:
                        counted.append(item)
                        if item in items:
                            items[item] += 1
                        else:
                            items[item] = 1

        f_list.extend([self.SequencePattern([[self.PLACE_HOLDER, k]], v, maxPatternLength, self.PLACE_HOLDER)
                       for k, v in _items.items()
                       if v >= threshold])
        f_list.extend([self.SequencePattern([[k]], v, maxPatternLength, self.PLACE_HOLDER)
                       for k, v in items.items()
                       if v >= threshold])

        # todo: can be optimised by including the following line in the 2 previous lines
        f_list = [i for i in f_list if self._checkPatternLengths(i, maxPatternLength)]

        sorted_list = sorted(f_list, key=lambda p: p.freq)
        return sorted_list

    def _build_projected_database(self, S, pattern):
        """
        suppose S is projected database base on pattern's prefix,
        so we only need to use the last element in pattern to
        build projected database
        """
        p_S = []
        last_e = pattern.sequence[-1]
        last_item = last_e[-1]
        for s in S:
            p_s = []
            for element in s:
                is_prefix = False
                if self.PLACE_HOLDER in element:
                    if last_item in element and len(pattern.sequence[-1]) > 1:
                        is_prefix = True
                else:
                    is_prefix = True
                    for item in last_e:
                        if item not in element:
                            is_prefix = False
                            break

                if is_prefix:
                    e_index = s.index(element)
                    i_index = element.index(last_item)
                    if i_index == len(element) - 1:
                        p_s = s[e_index + 1:]
                    else:
                        p_s = s[e_index:]
                        index = element.index(last_item)
                        e = element[i_index:]
                        e[0] = self.PLACE_HOLDER
                        p_s[0] = e
                    break
            if len(p_s) != 0:
                p_S.append(p_s)

        return p_S

class BasicRule:
    def __init__(self, precurssor , successor, count = None, probability = None):
        self.precurssor = precurssor
        self.count = count
        self.successor = successor
        self.probability = probability

def test_prefix_span():
    r = Reader(data_path='data/file.csv')
    sequences = r.get_seq()
    print(len(sequences))
    model = PrefixSpan.train(sequences, minSupport=0.3, maxPatternLength=2)
    result = model.freqSequences().collect()
    for fs in result:
        print('({} : {})'.format(fs.sequence, fs.freq))

def test_prefix_span_variable():
    r = Reader(data_path='data/variableStates.csv')
    sequences = r.get_seq()
    print(len(sequences))
    model = PrefixSpan.train(sequences, minSupport=0.3, maxPatternLength=2)
    result = model.freqSequences().collect()
    for fs in result:
        if len(fs.sequence) == 1:
            continue
        print('({} : {})'.format(fs.sequence, fs.freq))

def generate_sir():
    r = Reader(data_path='data/file.csv')
    attribute = 'flu'
    sequences = r.get_seq()
    print(len(sequences))
    model = PrefixSpan.train(sequences, minSupport=0.3, maxPatternLength=2)
    result = model.freqSequences().collect()
    result_dict = defaultdict(int)
    result_prob_dict = {}
    pre_count_dict = defaultdict(int)
    post_count_dict = defaultdict(int)
    rule_list: List[BasicRule] = []

    for fs in result:
        # print('{}, {}'.format(fs.sequence, fs.freq))
        pre = tuple(fs.sequence[0])
        if len(fs.sequence) > 1:
            post = tuple(fs.sequence[-1])
            tup_key = (pre, post)
            pre_count_dict[pre]
            rule = BasicRule(precurssor=pre, successor=post, count=fs.freq, probability=None)
            rule_list.append(rule)
        else:
            tup_key = (pre)
        result_dict[tup_key] = fs.freq

    for rule in rule_list:
        numerator = rule.count
        denom_discount = 0
        for other_rule in rule_list:
            if rule == other_rule or rule.precurssor == other_rule.precurssor:
                continue
            if rule.precurssor == other_rule.successor:
                denom_discount += other_rule.count
                print(result_dict[rule.precurssor])
        denominator = result_dict[rule.precurssor] - denom_discount
        rule.probability = numerator / denominator

    rules_dict_list: List[Dict[Union[int, str], Union[int, str]]] = []
    filename = 'rule2'
    for rule in rule_list:
        r_dict = {}
        # print(rule)
        # print(rule[0], list(rule[0])[0])
        precursor_list: List = rule.precurssor
        successor_list: List = rule.successor
        probability: float = rule.probability
        p_attribute: List = []
        s_attribute: List = []
        precursor: List = []
        successor: List = []
        s_dict = {}
        for item in precursor_list:
            # print(item)
            # a, b = item.split('_')
            precursor.append(item)
            p_attribute.append(attribute)
        for item in successor_list:
            # a, b = item.split('_')
            successor.append(item)
            s_attribute.append(attribute)
        s_dict['successorAttribute'] = s_attribute
        s_dict['value'] = successor
        s_dict['probability'] = probability
        s_list = []
        s_list.append(s_dict)
        r_dict['precursor_attribute'] = p_attribute
        r_dict['precursor'] = precursor
        r_dict['successor'] = s_list
        print(r_dict)
        rules_dict_list.append(r_dict)

    pprint(rules_dict_list)
    with open(f'jsonData/{filename}.json', 'w') as file:
        json.dump(rules_dict_list, file)
    pprint(rules_dict_list)
    print(rf'Rules written to {os.getcwd()}/jsonData/{filename}.json')

if __name__ == "__main__":
    print(' Starting Algorithm..')
    # test_prefix_span()
    test_prefix_span_variable()
    # generate_sir()
    # pprint(result_dict)




    # for fs in result:
    #     # print('{}, {}'.format(fs.sequence, fs.freq))
    #     if len(fs.sequence)
    #     pre = tuple(fs.sequence[0])
    #     if len(fs.sequence)>1:
    #         post = tuple(fs.sequence[-1])
    #         tup_key = (pre,post)
    #         pre_count_dict[pre]
    #     else:
    #         tup_key = (pre)
    #     result_dict[tup_key] = fs.freq
    # for key, val in result_dict.items():
    #     if len(key) > 1:
    #         continue

# https://bitbucket.org/lutianming/dm/src/master/PrefixSpan.py
#https://github.com/rangeonnicolas/PrefixSpan/blob/master/PrefixSpan.py

