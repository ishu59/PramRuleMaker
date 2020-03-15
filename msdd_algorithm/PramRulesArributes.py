from typing import Dict, List, Union
import json


class PramRulesArributes:

    def __init__(self, attribute:Union[List, str], precursor:Union[List, str]):
        # self.attribute = []
        # self.precursor = []
        self._add_precursor(attribute, precursor)
        self.successor = []
        self.successor_child = dict(value=[], successorAttribute=[], probability=0)
        self.rule = []
        self.prob = None

    def _add_precursor(self, attribute:Union[List, str], precursor:Union[List, str]):
        if type(attribute) == str:
            attribute = [attribute]
        if type(precursor) == str:
            precursor = [precursor]
        self.attribute = list(attribute)
        self.precursor = list(precursor)


    def add_new_successor(self, value: List,  probability: float, successorAttribute: List = None):
        if type(value) == str:
            value = [value]
        else:
            value = list(value)
        if successorAttribute is not None:
            if type(successorAttribute) == 'str':
                successorAttribute = [successorAttribute]
            else:
                successorAttribute = list(successorAttribute)
        else:
            successorAttribute = self.attribute

        self.successor.append({'successorAttribute': successorAttribute,
                               'value': value, 'probability': probability})

    def _get_rule(self):
        data = {'attribute': self.attribute, 'precursor': self.precursor, 'successor': self.successor}
        return data

    def get_rule_json(self):
        rule = self._get_rule()
        return json.dumps(rule)

    def get_rule_dict(self):
        data = self._get_rule()
        return data

    def __str__(self):
        data = {'attribute': self.attribute, 'precursor': self.precursor, 'successor': self.successor}
        data = str(data)
        # print(data)
        return data

    def __repr__(self):
        return self.__str__()

def rules_attr_to_json(pram_rules_list: List[PramRulesArributes]):
    data = []
    for item in pram_rules_list:
        data.append(item.get_rule_dict())
    return json.dumps(data)


def test_pram_rule_attr():
    att = ['flu','sex']
    prec = ['S', 'M']
    suc1 =['I', 'M']
    prob_suc1 = 0.7
    suc2 = ['S', 'M']
    prob_suc2 = 0.3
    pra = PramRulesArributes(attribute=att, precursor=prec)
    pra.add_new_successor(value=suc1, probability=prob_suc1)
    pra.add_new_successor(value=suc2, probability=prob_suc2)
    print(pra)
    print(pra.get_rule_dict())
    print(pra.get_rule_dict())
    l:PramRulesArributes = []
    l.append(pra)
    l.append(pra)
    r = rules_attr_to_json(l)
    print(r)


if __name__ == '__main__':
    print(" Testing Pram rule attribute")
    test_pram_rule_attr()

