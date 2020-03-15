import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import f, chisquare
def KL(a, b):
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)

    return np.sum(np.where(a != 0, a * np.log(a / b), 0))


gen = np.array([0.37,0.25,0.38])
actual = np.array([0.38, 0.24, 0.38])

print(KL(actual, gen ))

def kl_divergence(p, q):

    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

print(kl_divergence(actual,gen))

actual = [0.7,0.3,0.5,0.5,0.7,0.3]
predicted = [0.68,0.32,0.52,0.48,0.69,0.31]
print(r2_score(actual,predicted))

actual = [384.6, 230.8, 384.6]
generated = [369.9,247.1,383]


print(chisquare(actual, generated))
