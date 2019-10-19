from utils.helpers import *
from  RuleGrowth import  RuleGrowth

database = convert_from_txt('sample/word.txt')
rg =  RuleGrowth()
rg.fit(database,minconf=0.6, minsup=0.6)
# print(rg.save_rules('out.txt'))
rg.write_rules('file.txt', filetype='txt')
print(rg.predict(database))