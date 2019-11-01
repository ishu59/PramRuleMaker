from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation


class abcd(Rule):

    def apply(self, pop, group, iter, t):

        if group.has_attr({'flu': 'S'}): return [GroupSplitSpec(p=0.7, attr_set={'flu': 'S'}),
                                                 GroupSplitSpec(p=0.3, attr_set={'flu': 'I'}),
                                                 GroupSplitSpec(p=0, attr_set={'flu': 'R'}), ]
        if group.has_attr({'flu': 'I'}): return [GroupSplitSpec(p=0, attr_set={'flu': 'S'}),
                                                 GroupSplitSpec(p=0.5, attr_set={'flu': 'I'}),
                                                 GroupSplitSpec(p=0.5, attr_set={'flu': 'R'}), ]
        if group.has_attr({'flu': 'R'}): return [GroupSplitSpec(p=0.3, attr_set={'flu': 'S'}),
                                                 GroupSplitSpec(p=0, attr_set={'flu': 'I'}),
                                                 GroupSplitSpec(p=0.7, attr_set={'flu': 'R'}), ]