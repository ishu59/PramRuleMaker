from pram.data import GroupSizeProbe, ProbeMsgMode
from pram.entity import Group, GroupQry, GroupSplitSpec, Site
from pram.rule import GoToRule, DiscreteInvMarkovChain, TimeInt, Rule
from pram.sim import Simulation
from pram.data import GroupSizeProbe, ProbeMsgMode, ProbePersistanceDB


class SIRModel(Rule):
    def apply(self, pop, group, iter, t):

        if group.has_attr({'flu': 'S'}): return [GroupSplitSpec(p=0.70, attr_set={'flu': 'S'}),
                                                 GroupSplitSpec(p=0.0, attr_set={'flu': 'I'}), ]

        if group.has_attr({'flu': 'R'}): return [GroupSplitSpec(p=0.70, attr_set={'flu': 'R'}),
                                                 GroupSplitSpec(p=0.30, attr_set={'flu': 'S'}), ]

        if group.has_attr({'flu': 'I'}): return [GroupSplitSpec(p=0.50, attr_set={'flu': 'I'}),
                                                 GroupSplitSpec(p=0.50, attr_set={'flu': 'R'}), ]


probe_grp_size_flu = GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP,
                                            persistance=ProbePersistanceDB(),
                                            memo='Mass distribution across flu status')

p = GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], persistance=ProbePersistanceDB())


(Simulation().
    add_probe(GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP)).
    add_rule(SIRModel()).
    add_group(Group(m=1000, attr={ 'flu': 'S' })).
    run(30)
)
# =======================================================================================================================================
# Sim with plot
# =======================================================================================================================================
#
# (Simulation().
#     add_probe(GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP)).
#     add_rule(SIRModel()).
#     add_probe(p).
#     run(26)
# )
#
#
# series = [
#     { 'var': 'p0', 'lw': 0.75, 'linestyle': '-',  'marker': 'o', 'color': 'red',   'markersize': 0, 'lbl': 'S' },
#     { 'var': 'p1', 'lw': 0.75, 'linestyle': '--', 'marker': '+', 'color': 'blue',  'markersize': 0, 'lbl': 'I' },
#     { 'var': 'p2', 'lw': 0.75, 'linestyle': ':',  'marker': 'x', 'color': 'green', 'markersize': 0, 'lbl': 'R' }
# ]
# p.plot(series, figsize=(16,3))