import random

class Generator:

    def __init__(self, item_set = None, max_rows = 100, max_cols = 15, max_items = 5):
        self._max_items = max_items
        self._max_cols = max_cols
        self._max_rows = max_rows
        self._item_set = item_set
        if item_set is None:
            self._item_set =  ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

    def generate_seq(self, out_file = 'data/file.csv'):
        with open(out_file, 'w') as file:
            for i in range(0, self._max_rows):
                # elem_limit = random.randrange(1, self._max_cols + 1, 1)
                elem_limit = self._max_cols
                out_string = ""
                for j in range(0, elem_limit):
                    # item_limit = random.randrange(1, self._max_items + 1, 1)
                    item_limit = self._max_items
                    out_elem = ""
                    for k in range(0, item_limit):
                        while (True):
                            item = random.choice(self._item_set)
                            if item not in out_elem:
                                out_elem = out_elem + item
                                break
                    out_string = out_string + out_elem
                    if j != elem_limit - 1:
                        out_string = out_string + ","
                if i != self._max_rows - 1:
                    out_string = out_string + "\n"
                file.write(out_string)
            file.close()

    def create_from_text(self, txt_path = 'data/sir_data.txt'):

        with open('data/sir_data.txt', 'r') as file:
            l = file.readlines()
        # print(l[0].split(','))
        l = l[0].split(',')
        print(l)
        data = []
        with open('data/sir.csv' , 'w') as file:
            for i in range(len(l)-1):
                pre = l[i]
                post = l[i+1]
                data.append([pre,post])
                item = str(pre) +','+ str(post) + '\n'
                file.write(item)

        return data

if __name__ == '__main__':
    # g = Generator(item_set=['S','I','R'], max_items=2, max_cols=2)
    g = Generator(item_set=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'], max_items=4, max_cols=3)
    g.generate_seq(out_file='data/variableStates.csv')
    # g.create_from_text()