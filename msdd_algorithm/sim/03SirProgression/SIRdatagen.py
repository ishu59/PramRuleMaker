from algorithm import clean_data
from pprint import pprint

def read_data_lower_middle_SIR():
    with open('data/sir_data.txt', 'r') as file:
        l = file.readlines()
    # print(l[0].split(','))
    l = l[0].split(',')
    l = clean_data(l)
    data = []
    for i in range(len(l)-1):
        data.append([l[i], l[i+1]])
    with open('data/sir_col_data.txt', 'w') as file:
        data_line = []
        for line in data:
            data_line = []
            a = 'flu_' + line[0]
            b = 'flu_' + line[1]
            # c = 'flu_' + line[2]
            file.write(f'{a},{b} \n')

    return l


if __name__ == '__main__':
    read_data_lower_middle_SIR()