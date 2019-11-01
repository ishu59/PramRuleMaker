
from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation

class HeartDisease(Rule):

    def apply(self, pop, group, iter, t):
        
        if group.has_attr({'cp': '1'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '0'}),]
        if group.has_attr({'cp': '1'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'hd': 'no'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.7272727272727273, attr_set={'hd': 'no'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.7272727272727273, attr_set={'oldpeak': '0'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.5454545454545454, attr_set={'oldpeak': '0'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.5454545454545454, attr_set={'hd': 'no'}),]
        if group.has_attr({'cp': '4'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'hd': 'yes'}),]
        if group.has_attr({'cp': '3'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '0'}),]
        if group.has_attr({'cp': '3'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'hd': 'no'}),]
        if group.has_attr({'cp': '3'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'hd': 'yes'}),]
        if group.has_attr({'cp': '3'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '4'}),]
        if group.has_attr({'oldpeak': '1'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'hd': 'yes'}),]
        if group.has_attr({'oldpeak': '1'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '4'}),]
        if group.has_attr({'hd': 'yes'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '1'}),]
        if group.has_attr({'cp': '4'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '1'}),]
        if group.has_attr({'hd': 'yes'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '4'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'cp': '3'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'hd': 'yes'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'cp': '4'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'cp': '3'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'hd': 'yes'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.36363636363636365, attr_set={'cp': '4'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '3'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '1'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'hd': 'yes'}),]
        if group.has_attr({'cp': '2'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '4'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '2'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'cp': '2'}),]
        if group.has_attr({'oldpeak': '0'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '1'}),]
        if group.has_attr({'hd': 'no'}): return[GroupSplitSpec(p = 0.2727272727272727, attr_set={'oldpeak': '1'}),]