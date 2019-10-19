# from sklearn.datasets import load_iris
# from sklearn import tree
# import matplotlib.pyplot as plt
# iris = load_iris()
# clf = tree.DecisionTreeClassifier()
# clf = clf.fit(iris.data, iris.target)
#
# tree.plot_tree(clf.fit(iris.data, iris.target))
# plt.show()

from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree.export import export_text
from sklearn import tree
iris = load_iris()
X = iris['data']
y = iris['target']
decision_tree = DecisionTreeClassifier(random_state=0, max_depth=2)
decision_tree = decision_tree.fit(X, y)
tree.plot_tree(decision_tree)
r = export_text(decision_tree, feature_names=iris['feature_names'])
print(r)