
from indra.sources import eidos
from indra.assemblers.cag import CAGAssembler
from indra.assemblers.pysb import PysbAssembler
import pickle




# pa = PysbAssembler()
# pa.add_statements(ep.statements)
# pa.make_model(policies='two_step')
# for monomer in pa.model.monomers:
#     print(monomer)
# for rule in pa.model.rules:
#     print(rule)
# for parameter in pa.model.parameters:
#     print(parameter)

def dashes(keyword: str = None, pattern: str = '=', num: int = 10):
    if keyword is None:
        keyword = ''
    start_line = ' '.join([num*pattern,keyword,num*pattern])
    end_line = ''.join(len(start_line)*pattern)
    return start_line, end_line

def eidos_process_statements(sentence=None,webservice=None,use_webService = False):
    print('Running...rul    ')
    ep = None
    if use_webService:
        if webservice is None:
            webservice = 'http://localhost:9000'
        if sentence is not None:
            ep = eidos.process_text(sentence, webservice=webservice)
    else:
        ep = eidos.process_text(sentence)
    return ep

def indra_processing(statements):
    pa = PysbAssembler()
    pa.add_statements(statements)
    # pa.make_model(policies='two_step')
    pa.make_model()
    print(10*'=', 'Monomer', 10*'=')
    for monomer in pa.model.monomers:
        print(monomer)
    for rule in pa.model.rules:
        print(rule)
    for parameter in pa.model.parameters:
        print(parameter)
    return pa

def test():
    sent = 'Conflict causes displacement, which leads to hunger'

    ep = eidos_process_statements(sent)
    pa = indra_processing(ep.statements)
    print(pa)

if __name__ == '__main__':
    print(' Processing.. ')
    test()














# saved_file = 'savedFiles/sent1.pkl'
    # try:
    #     ep = pickle.load(open(saved_file, 'rb'))
    # except (OSError, IOError) as e:
    #     ep = eidos_process_statements(sent)
    #     pickle.dump(ep, open(saved_file, 'wb'))