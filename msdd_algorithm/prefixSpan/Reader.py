import csv
from pprint import pprint
import numpy as np

class Reader:
    def __init__(self, data_path):
        self._data_path = data_path

    def get_seq(self):
        seq_set = []
        with open(self._data_path, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                # print ("\n"+str(row))
                seq = []
                for elem in row:
                    seq_elem = list(elem)
                    seq_elem.sort()
                    seq.append(seq_elem)
                seq_set.append(seq)
        return seq_set

if __name__ == '__main__':
    r = Reader(data_path='data/file.csv')
    data = np.array(r.get_seq())
    pprint(data)