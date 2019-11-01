from typing import Dict, List
import json


class PramRulesArributes:
    def __init__(self):
        self.attribute = []
        self.precursor = []
        self.successor = []
        self.successor_child = dict(value=[], successorAttribute=[], probability=0)
        self.rule = []
        # self.successorAttribute = []
        # self.value = []
        # self.probability = 0

    def add_new_successor(self, value: List, successorAttribute: List, probability: float):
        self.successor.append({'successorAttribute': successorAttribute,
                               'value': value, 'probability': probability})

    def add_new_rule(self):
        pass

    def __str__(self):
        data = {'attribute': self.attribute, 'precursor': self.precursor, 'successor': self.successor}
        data = str(data)
        # print(data)
        return data

    def __repr__(self):
        return self.__str__()
