from math import ceil
import json
import os
import warnings
from tqdm import tqdm
from utils.readParse import convert_from_txt
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union, Optional
from pprint import pprint
pd.set_option('mode.chained_Assignment', None)

class RuleGenerator:

    def __init__(self):
        self.database: List[List[Optional[int, str, float]]] = None
        self.rules = []  # List containing tuple of valid rules (antecedentSet,consequentSet,support,confidence)
        self.minsup = None
        self.minsupRelative = None
        self.minconf = None
        self.itemlocmap = {}  # Dictionary of item's first and last location in each sequence. {Item: {seqID:[first,last]}}
        self.above_threshold = []  # List of items that are >= the min support
        self.maxAntecedent = None
        self.maxConsequent = None
        self.bin = None  # Interval to bin rules for sorting
        self.fitted = False

    def fit(self, database, minsup, minconf, maxAntecedent=0, maxConsequent=0, bin_interval=0):
        '''
        Generate sequential rules
        
        Arguments:
            database  - List of purchase histories
            
            minsup - Minimum fraction of I => J occurances in all sequences
            
            minconf - Minimum fraction of I => J occurances in sequences with I
            
            maxAntecedent - Maximum number of items in antecedent set (default: 0)
            
            maxConsequent - Maximum number of items in consequent set (default: 0)
            
            bin_interval - Interval to bin rules for sorting (NOT IN USE)
        '''

        # Set all attributes
        self.database: List[List[Optional[int, str, float]]] = database
        self.rules: List = []
        self.minsup = minsup
        self.minsupRelative = ceil(minsup * len(database))
        self.minconf = minconf
        self.itemlocmap = {}
        self.above_threshold = []
        self.maxAntecedent = maxAntecedent
        self.maxConsequent = maxConsequent
        self.bin = bin_interval

        def validate():
            # Check arguments
            assert isinstance(self.minsup,
                              float) and 1 >= self.minsup > 0, 'Support has to be in range (0,1]'
            assert isinstance(self.minconf,
                              float) and 1 >= self.minconf > 0, 'Confidence has to be in range (0,1]'
            assert isinstance(self.maxAntecedent,
                              int) and self.maxAntecedent >= 0, 'maxAntecedent has to be a positive integer'
            assert isinstance(self.maxConsequent,
                              int) and self.maxConsequent >= 0, 'maxConsequent has to be a positive integer'
            assert isinstance(self.bin, (float, int)) and self.bin == 0 or int(
                1 / self.bin) == 1 / self.bin, 'inverse of bin_intervalhas to be an integer. ie. 1has to be divisible by bin_interval'
            if self.bin:
                warnings.warn('bin_interval is not in use and will be ignored.', UserWarning)

        validate()

        # Generate itemlocmap & above_threshold
        def getMinSupItems():
            '''
            Function to generate itemlocmap
            Record items that are above support threshold
            Remove items below threshold of support within database
            '''
            # For each sequence in database
            for seqID, sequence in enumerate(self.database):
                # For each itemset in sequence
                for idx, itemset in enumerate(sequence):
                    # For each item in itemset
                    for item in itemset:
                        # Item not found in itemlocmap yet. Add item into itemlocmap
                        if item not in self.itemlocmap:
                            self.itemlocmap[item] = {seqID: [idx, idx]}
                        # First time item is found in sequence. Add sequence into itemlocmap[item]
                        elif seqID not in self.itemlocmap[item]:
                            self.itemlocmap[item][seqID] = [idx, idx]
                        # Another occurance of item is found. Update last known location of item in sequence.
                        else:
                            self.itemlocmap[item][seqID][1] = idx

            # Generate list of frequentitems
            below_threshold = []
            for item, value in self.itemlocmap.items():
                if len(value) < self.minsupRelative:
                    below_threshold.append(item)
                else:
                    self.above_threshold.append(item)
            # Remove items below min support from database, reduce number of items to loop over
            # self.database = [[np.setdiff1d(itemset,below_threshold) for itemset in seq] for seq in self.database]

        getMinSupItems()

        # Generate Rules
        def genRules():
            '''
            This function first generates valid rules with both antecedent and consequent of size 1.
            Then it will recursively expand the rules through expandLeft and expandRight
            '''
            # Double for loop to run through all possible I => J combination with no repeats
            for i in range(len(self.above_threshold)):
                for j in range(i + 1, len(self.above_threshold)):
                    # items
                    itemI, itemJ = self.above_threshold[i], self.above_threshold[j]
                    # Dictionary of {seqID:[first,last]} for items
                    occurancesI, occurancesJ = self.itemlocmap[itemI], self.itemlocmap[itemJ]
                    # Sequences containing items I & J
                    allseqI, allseqJ = set(occurancesI.keys()), set(occurancesJ.keys())
                    # Sequences containing I => J or J => I
                    allseqIJ, allseqJI = set(), set()

                    # Sequences that contains both I & J
                    allseqboth = set.intersection(allseqI, allseqJ)

                    # Sequences that have I => J or J => I
                    for seqID in allseqboth:
                        if self.itemlocmap[itemI][seqID][0] < self.itemlocmap[itemJ][seqID][1]:
                            allseqIJ.add(seqID)
                        if self.itemlocmap[itemJ][seqID][0] < self.itemlocmap[itemI][seqID][1]:
                            allseqJI.add(seqID)

                    # Check IJ
                    if len(allseqIJ) >= self.minsupRelative:
                        confIJ = len(allseqIJ) / len(occurancesI)
                        antecedentSet = set([itemI])
                        consequentSet = set([itemJ])
                        # Add those with valid support and confidence
                        if confIJ >= self.minconf:
                            self.rules.append(
                                (antecedentSet, consequentSet, len(allseqIJ) / len(self.database), confIJ))
                        # Expand left if possible
                        if not self.maxAntecedent or len(antecedentSet) < self.maxAntecedent:
                            expandLeft(antecedentSet, consequentSet, allseqI, allseqIJ, occurancesJ)
                        # Expand right if possible
                        if not self.maxConsequent or len(consequentSet) < self.maxConsequent:
                            expandRight(antecedentSet, consequentSet, allseqI, allseqJ, allseqIJ, occurancesI,
                                        occurancesJ)

                    # Check JI
                    if len(allseqJI) >= self.minsupRelative:
                        confJI = len(allseqJI) / len(occurancesJ)
                        antecedentSet = set([itemJ])
                        consequentSet = set([itemI])
                        # Add those with valid support and confidence
                        if confJI >= self.minconf:
                            self.rules.append(
                                (antecedentSet, consequentSet, len(allseqJI) / len(self.database), confJI))
                        # Expand left if possible
                        if not self.maxAntecedent or len(antecedentSet) < self.maxAntecedent:
                            expandLeft(antecedentSet, consequentSet, allseqJ, allseqJI, occurancesI)
                        # Expand right if possible
                        if not self.maxConsequent or len(consequentSet) < self.maxConsequent:
                            expandRight(antecedentSet, consequentSet, allseqJ, allseqI, allseqJI, occurancesJ,
                                        occurancesI)

        def expandLeft(antecedentSet, consequentSet, allseqI, allseqIJ, occurancesJ):
            '''
            This function builds on an existing rulebyadding an item to the antecedentset
            '''
            # A dictionary of items C and sequenceIDs where IuC => J
            possibleC = dict()
            # Total number of possible sequences
            seqsLeft = len(allseqIJ)

            # Loop to populate possibleC
            # For each sequenceID where I => J
            for seqID in allseqIJ:
                sequence = self.database[seqID]  # Get sequence
                firstJ, lastJ = occurancesJ[seqID]  # Get last occurance of J in sequene

                # For each itemsetID before the itemset containing the last occurance of J in the sequence
                for itemsetID in range(lastJ):
                    itemset = sequence[itemsetID]  # Get itemset
                    # For each item in itemset
                    for item in itemset:
                        if any([i >= item for i in
                                antecedentSet]) or item in consequentSet or item not in self.above_threshold:
                            # Ensure that the item is not already present in either
                            # antecedent or consequent set.
                            # To prevent repeated rules, only items greater in value
                            # than all items inside antecedent set will be considered
                            continue
                        if item not in possibleC:
                            if seqsLeft >= self.minsupRelative:
                                # items that are not able to meet the
                                # minimum requirement are ignored
                                possibleC[item] = set([seqID])

                        elif len(possibleC[item]) + seqsLeft < self.minsupRelative:
                            # Remove items from possibleC when it is no
                            # longer possible to meet support requirement
                            del possibleC[item]
                            continue

                        else:
                            # Add sequenceID
                            possibleC[item].add(seqID)

                # Decrease max possible sequence left
                seqsLeft -= 1

            # Loop through possibleC to generate valid rules
            for itemC, seqIDs in possibleC.items():
                # Check if minimum support requirement is met
                if len(seqIDs) >= self.minsupRelative:
                    # SeqIDs of IuC 
                    allseqIC = set.intersection(set(self.itemlocmap[itemC].keys()), allseqI)

                    # Confidence of IuC => J
                    # support(IuC => J) / support(IuC)
                    confIC_J = len(seqIDs) / len(allseqIC)

                    # New antecedent set
                    itemsIC = antecedentSet.copy()
                    itemsIC.add(itemC)

                    # Add rule
                    if confIC_J >= self.minconf:
                        self.rules.append((itemsIC, consequentSet, len(seqIDs) / len(self.database), confIC_J))

                    # Expand left if possible
                    if not self.maxAntecedent or len(itemsIC) < self.maxAntecedent:
                        expandLeft(itemsIC, consequentSet, allseqIC, seqIDs, occurancesJ)

        def expandRight(antecedentSet, consequentSet, allseqI, allseqJ, allseqIJ, occurancesI, occurancesJ):
            '''
            This function builds on an existing rule by adding an item to the consequent set
            '''

            # A dictionary of item C and sequenceIDs where I => JuC
            possibleC = dict()
            # Total number of possible sequences
            seqsLeft = len(allseqIJ)

            # Loop to populate possibleC 
            # For each sequenceID where I => J
            for seqID in allseqIJ:
                sequence = self.database[seqID]  # Get sequence
                firstI, lastI = occurancesI[seqID]  # Get first occurance of I in sequence

                # For each itemsetID after the itemset containing the first occurance of I in the sequence
                for itemsetID in range(firstI + 1, len(sequence)):
                    itemset = sequence[itemsetID]  # Get itemset
                    # For each item in itemset
                    for item in itemset:
                        if any([i >= item for i in
                                consequentSet]) or item in antecedentSet or item not in self.above_threshold:
                            # Ensure that the item is not already present in either
                            # antecedent or consequent set
                            # To prevent repeated rules, only item greater in value
                            # than all items inside the consequent set will be considered
                            continue
                        if item not in possibleC:
                            if seqsLeft >= self.minsupRelative:
                                # items that are not able to meet the 
                                # minimum requirement are ignored
                                possibleC[item] = set([seqID])

                        elif len(possibleC[item]) + seqsLeft < self.minsupRelative:
                            # Remove items from possibleC when it is no 
                            # longer possible to meet support requirement
                            del possibleC[item]
                            continue

                        else:
                            # Add sequenceID
                            possibleC[item].add(seqID)

                # Decrease max possible sequence left
                seqsLeft -= 1

            # Loop through possibleC to generate valid rules
            for itemC, seqIDs in possibleC.items():
                # Check if minimum support requirement is met
                if len(seqIDs) >= self.minsupRelative:
                    # SeqIDs of JuC
                    allseqJC = set()
                    # New consequent occurance map
                    occurancesJC = dict()

                    # Loop through the consequent set to find intersection with item C
                    # Update the occurance of consequent set to make sure the last 
                    # occurance is the earliest last occurance among all items
                    for seqID_J in allseqJ:
                        occurancesC = self.itemlocmap[itemC].get(seqID_J, False)
                        if occurancesC:
                            allseqJC.add(seqID_J)
                            firstJ, lastJ = occurancesJ[seqID_J]
                            if occurancesC[1] < lastJ:
                                occurancesJC[seqID_J] = occurancesC
                            else:
                                occurancesJC[seqID_J] = [firstJ, lastJ]

                    # Confidence of I => JuC
                    # support(I => JuC) / support(I)
                    confI_JC = len(seqIDs) / len(allseqI)

                    # New consequent set
                    itemsJC = consequentSet.copy()
                    itemsJC.add(itemC)

                    # Add rule
                    if confI_JC >= self.minconf:
                        self.rules.append((antecedentSet, itemsJC, len(seqIDs) / len(self.database), confI_JC))

                    # Expand left if possible
                    if not self.maxAntecedent or len(antecedentSet) < self.maxAntecedent:
                        expandLeft(antecedentSet, itemsJC, allseqI, seqIDs, occurancesJC)

                    # Expand right if possible
                    if not self.maxConsequent or len(itemsJC) < self.maxConsequent:
                        expandRight(antecedentSet, itemsJC, allseqI, allseqJC, seqIDs, occurancesI, occurancesJC)

        genRules()

        # Sort rules
        def sort_rules():
            '''
            Sort the rules in the following order of importance
            1. Binned confidence
            2. length of antecedent
            3. Confidence
            '''

            def binning(x):
                '''
                Bin confidence by rounding down to nearest 0.05
                binning(0.44) => 0.4
                binning(0.045) => 0.45
                binning(0.51) => 0.5
                '''
                if self.bin == 0:
                    return x
                intervals = 1 / self.bin
                return int(x * intervals) / intervals

            # Sort by confidence
            temp = sorted(self.rules, key=lambda x: x[3], reverse=True)
            # Sort by length of antecedent 
            # temp = sorted(temp,key=lambda x:len(x[0]),reverse=True)
            # Sort by binned confidence
            # temp = sorted(temp,key=lambda x:binning(x[3]),reverse=True)
            # Return sorted rules
            self.rules = temp

        sort_rules()

        # Clear memory
        self.database = None
        self.itemlocmap = {}
        self.above_threshold = []

        # Number of unique consequent items

        # Fitted 
        self.fitted = True

    def write_rules(self, filename, filetype='json'):
        if not self.rules:
            return 'There are no rules to write'
        if filetype == 'txt':
            with open(f'{filename}.txt', 'w') as text_file:
                for rule in self.rules:
                    writable_rule = rule[:-1]
                    text_file.write('{} ==> {} #Probability: {}\n'.format(*rule))
            print(rf'Rules written to {os.getcwd()}\{filename}.txt')
        elif filetype == 'json':
            rules_dict_list: List[Dict[Union[int, str], Union[int, str]]] = []
            for rule in self.rules:
                r_dict = {}
                # print(rule)
                # print(rule[0], list(rule[0])[0])
                precursor_list:List = list(rule[0])
                successor_list:List = list(rule[1])
                probability:float = rule[2]
                p_attribute:List = []
                s_attribute:List = []
                precursor:List = []
                successor:List = []
                s_dict = {}
                for item in  precursor_list:
                    # print(item)
                    a,b = item.split('_')
                    precursor.append(b)
                    p_attribute.append(a)
                for item in successor_list:
                    # print(item)
                    a,b = item.split('_')
                    successor.append(b)
                    s_attribute.append(a)
                s_dict['successorAttribute'] = s_attribute
                s_dict['value'] = successor
                s_dict['probability'] = probability
                s_list = []
                s_list.append(s_dict)
                r_dict['precursor_attribute'] = p_attribute
                r_dict['precursor'] = precursor
                r_dict['successor'] = s_list
                # print(r_dict)
                rules_dict_list.append(r_dict)

            pprint(rules_dict_list)
            with open(f'jsonData/{filename}.json', 'w') as file:
                json.dump(rules_dict_list, file)
            pprint(rules_dict_list)
            print(rf'Rules written to {os.getcwd()}/jsonData/{filename}.json')

        elif filetype == 'csv':
            df = pd.DataFrame(self.rules, columns=['antecedent', 'consequent', 'support', 'confidence'])
            df.to_csv(f'{filename}.csv', index=False)
            print(rf'Rules written to {os.getcwd()}\{filename}.csv')

        else:
            print(f'{filetype} is not a supported file format')


def test_rule_gen_simple():
    data = convert_from_txt('data/rule_gen_data_sir.txt', tuple_delimiter=',')
    rg = RuleGenerator()
    rg.fit(data, minconf=0.2, minsup=0.2, maxAntecedent=1, maxConsequent=1)
    rg.write_rules('file', filetype='txt')

def test_rule_gen_hd_data():
    data = convert_from_txt('data/hdata2.txt', tuple_delimiter=',', element_delimiter='\t')
    rg = RuleGenerator()
    rg.fit(data, minconf=0.2, minsup=0.2, maxAntecedent=1, maxConsequent=1)
    rg.write_rules('hdfile', filetype='json')

def test_rule_gen_sir_data():
    data = convert_from_txt('data/sir_rgen_data.txt', tuple_delimiter=',')
    rg = RuleGenerator()
    rg.fit(data, minconf=0.2, minsup=0.2, maxAntecedent=1, maxConsequent=1)
    rg.write_rules('SIRNew1', filetype='json')



if __name__ == '__main__':
    print('======================== Running Program =====================================')
    # test_rule_gen_simple()
    test_rule_gen_hd_data()
    # test_rule_gen_sir_data()









def other():
    # data = convert_from_txt('data/rule_gen_data_sir.txt', tuple_delimiter=',')
    data = convert_from_txt('data/hdata2.txt', tuple_delimiter=',', element_delimiter='\t')
    rg = RuleGenerator()
    rg.fit(data, minconf=0.2, minsup=0.2, maxAntecedent=1, maxConsequent=1)
    # rg.write_rules('file', filetype='txt')
    rg.write_rules('hdfile', filetype='json')
