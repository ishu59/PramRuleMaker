from pram.data import GroupSizeProbe, ProbeMsgMode, ProbePersistanceDB
from pram.entity import Group
from pram.sim    import Simulation
# from Autogenerated1 import Autogenerated1
from abcd import abcd
probe_grp_size_flu = GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP,
                                            persistance=ProbePersistanceDB(), memo='Mass distribution across flu status')

p = GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], persistance=ProbePersistanceDB())

(Simulation().
    add_probes([GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP),p]).

    add_rule(abcd()).
    add_group(Group(m=1000, attr={ 'flu': 'S' })).
    run(26)
)

# (Simulation().
#     add_probe(GroupSizeProbe.by_attr('flu', 'flu', ['S', 'I', 'R'], msg_mode=ProbeMsgMode.DISP)).
#     add_rule(Autogenerated()).
#     add_probe(p).
#     run(26)
# )
#
#

series = [
    { 'var': 'p0', 'lw': 0.75, 'linestyle': '-',  'marker': 'o', 'color': 'red',   'markersize': 0, 'lbl': 'S' },
    { 'var': 'p1', 'lw': 0.75, 'linestyle': '--', 'marker': '+', 'color': 'blue',  'markersize': 0, 'lbl': 'I' },
    { 'var': 'p2', 'lw': 0.75, 'linestyle': ':',  'marker': 'x', 'color': 'green', 'markersize': 0, 'lbl': 'R' }
]
p.plot(series, figsize=(16,3))