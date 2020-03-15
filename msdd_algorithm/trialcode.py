import itertools
from pprint import pprint

d = [i for i in 'ABCD']
d = [i for i in 'AAB']
# d = ['A']
s = set()
for i in range(1,len(d)+1):
    items = itertools.combinations_with_replacement(d,i)
    items = list(items)
    # print(items)
    for item in items:
        s.add(item)
pprint(s)

s = set()
for i in range(1,len(d)+1):
    items = itertools.combinations(d,i)
    items = list(items)
    # print(items)
    for item in items:
        s.add(item)

pprint(s)


d = (i for i in 'AAB')
d = tuple(d)
print(d)
print(list(d))
d = 'ABC'
d = [d]
print(list(d))