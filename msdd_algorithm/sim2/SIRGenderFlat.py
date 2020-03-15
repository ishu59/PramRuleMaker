
# Autogenerated file for the pram rule
from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation

class SIRGenderFlat(Rule):

    def apply(self, pop, group, iter, t):
        
        if group.has_attr({'attribute_0': 'S', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.7272727272727273, attr_set={'attribute_0': 'S', 'attribute_1': 'M'}),]
        if group.has_attr({'attribute_0': 'S'}): return [GroupSplitSpec(p = 0.7272727272727273, attr_set={'attribute_0': 'S'}),]
        if group.has_attr({'attribute_0': 'R'}): return [GroupSplitSpec(p = 0.6898907103825137, attr_set={'attribute_0': 'R'}),GroupSplitSpec(p = 0.31010928961748635, attr_set={'attribute_0': 'S'}),]
        if group.has_attr({'attribute_0': 'R', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.6898907103825137, attr_set={'attribute_0': 'R', 'attribute_1': 'M'}),]
        if group.has_attr({'attribute_0': 'I', 'attribute_1': 'M'}): return [GroupSplitSpec(p = 0.47086247086247085, attr_set={'attribute_0': 'I', 'attribute_1': 'M'}),GroupSplitSpec(p = 0.5291375291375291, attr_set={'attribute_0': 'R', 'attribute_1': 'M'}),]
        if group.has_attr({'attribute_0': 'I'}): return [GroupSplitSpec(p = 0.47086247086247085, attr_set={'attribute_0': 'I'}),GroupSplitSpec(p = 0.5291375291375291, attr_set={'attribute_0': 'R'}),]