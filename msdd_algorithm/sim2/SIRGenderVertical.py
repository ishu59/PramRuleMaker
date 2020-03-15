
# Autogenerated file for the pram rule
from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation

class SIRGenderVertical(Rule):

    def apply(self, pop, group, iter, t):
        
        if group.has_attr({'attribute_0': 'S', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.7272727272727273, attr_set={'attribute_0': 'S', 'attribute_1': 'M'}),]
        if group.has_attr({'attribute_0': 'S'}): return [GroupSplitSpec(p = 0.7272727272727273, attr_set={'attribute_0': 'S'}),]
        if group.has_attr({'attribute_0': 'R', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.6903137789904502, attr_set={'attribute_0': 'R', 'attribute_1': 'M'}),GroupSplitSpec(p = 0.3096862210095498, attr_set={'attribute_0': 'S', 'attribute_1': 'M'}),]
        if group.has_attr({'attribute_0': 'R'}): return [GroupSplitSpec(p = 0.6903137789904502, attr_set={'attribute_0': 'R'}),]
        if group.has_attr({'attribute_0': 'I'}): return [GroupSplitSpec(p = 0.4697674418604651, attr_set={'attribute_0': 'I'}),GroupSplitSpec(p = 0.5302325581395348, attr_set={'attribute_0': 'R'}),]
        if group.has_attr({'attribute_0': 'I', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.4697674418604651, attr_set={'attribute_0': 'I', 'attribute_1': 'M'}),GroupSplitSpec(p = 0.5302325581395348, attr_set={'attribute_0': 'R', 'attribute_1': 'M'}),]