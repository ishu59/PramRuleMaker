from msdd_algorithm.prefixSpan.Reader import Reader
from pprint import pprint
r = Reader(data_path='../prefixSpan/data/file.csv')
# print(r.get_seq())
seq = r.get_seq()
# seq = [[['S'], ['SI'], ['RS']], [['I'], ['R'], ['I']], [['S'], ['I'], ['I']], [['R'], ['I'], ['S']], [['S'], ['I'], ['S']], [['R'], ['S'], ['I']], [['S'], ['S'], ['S']]]
# pprint(seq)
max_row = len(seq)
max_col = 0
for row in seq:
    if len(row) > max_col:
        max_col = len(row)

# print(max_row, max_col)
with open('data/spadeData.csv','w') as file:
    for i in range(len(seq)):
        row = seq[i]
        for j in range(len(row)):
            items = row[j]
            # print(i, j, *items)
            data = ''.join(items)
            s = str(i) + ','+ str(j) +','+ data + '\n'
            # print(s)
            file.write(s)