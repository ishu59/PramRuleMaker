from typing import Dict, List
from pprint import pprint
import json


with open('sampleRule.json') as file:
    data = json.load(file)

class_stmt = '''
from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation


class Autogenerated(Rule):
    def apply(self, pop, group, iter, t):

'''

# why precursor_attribute is list?
# bcoz we can multiple precursor_attribute like flue location other in a single rule
# why precursor_attribute is there in successor list?
# bcoz precursor_attribute A can have rules affecting precursor_attribute B population, can also have multiple attributes
conditional_start = "if "
conditional_head = "group.has_attr("
conditional_stmt = None  # use dict string str(dict) with key from precursor_attribute value from precursor
conditional_and = ' and '
conditional_tail = "): "
consequent_head = "return["
consequent_body_start = "GroupSplitSpec("
consequent_body_end = "),"
consequent_body_prob_text = "p = "
consequent_body_prob_val = None  # just the prob value as float
consequent_body_attr_text = ", attr_set="
consequent_body_attr_val = None  # use successor dict string str(dict) with key from successor precursor_attribute value from successor value
consequent_tail = "]"
comma_str = ', '
rule_buffer = []
for rule in data:
    rule: Dict
    precursor_buffer = ''
    precursor_attribute: List = rule.get("precursor_attribute")
    precursor: List = rule.get("precursor")
    successor: List = rule.get("successor")
    conditions = dict(zip(precursor_attribute, precursor))
    # print(conditions)
    precursor_buffer += conditional_start + conditional_head + str(conditions) + conditional_tail + consequent_head
    successor_buffer = ''
    probability = 0
    for subsequent in successor:
        successor_attribute: List = subsequent.get("precursor_attribute")
        value: List = subsequent.get("value")
        probability = subsequent.get("probability")
        final_successor = dict(zip(successor_attribute, value))
        conseq = consequent_body_start + consequent_body_prob_text + str(probability) + consequent_body_attr_text + str(
            final_successor) + consequent_body_end
        successor_buffer += conseq
        # print(final_successor, probability)
        # print(conseq)
    rule_stmt = precursor_buffer + successor_buffer + consequent_tail
    rule_buffer.append(rule_stmt)

# pprint(rule_buffer)
for item in rule_buffer:
    class_stmt = class_stmt + '\n\t\t' + item

print(class_stmt)

