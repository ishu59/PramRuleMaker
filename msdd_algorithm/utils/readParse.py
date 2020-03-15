import numpy as np
from typing import Union, Optional, List, Dict
from pprint import pprint

def convert_from_txt(file_path: str, tuple_delimiter: str = ',', element_delimiter: str = ' '
                     ) -> List[List[Union[str, int, float]]]:
    """
    Read file where every line contains agent data where every tuple is separated by tuple_delimiter and
     items are separated by element delimiter
    :param file_path: path of the data file to parsed
    :param tuple_delimiter: separate events by tuple delimiter
    :param element_delimiter: separate elements by elements delimiter
    :return: List of List of elements
    """
    with open(file_path, 'r') as f:
        data = [x.split(tuple_delimiter) for x in f.readlines()]
        data = [([np.array(j.strip().split(element_delimiter)) for j in i]) for i in data]
    return data


def test_convert_from_txt():
    data = convert_from_txt('../data/rule_gen_data.txt', tuple_delimiter='-1', element_delimiter=' ')
    pprint(data)

if __name__ == '__main__':
    test_convert_from_txt()