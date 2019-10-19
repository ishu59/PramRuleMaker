# Load libraries
import pandas as pd
from sklearn.tree import DecisionTreeClassifier # Import Decision Tree Classifier
from sklearn.model_selection import train_test_split # Import train_test_split function
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn.tree import export_graphviz
from sklearn.externals.six import StringIO
from IPython.display import Image
import pydotplus
from sklearn.tree import _tree


def dt_default(X_train, X_test, y_train, y_test ):
    # Create Decision Tree classifer object
    clf = DecisionTreeClassifier()
    # Train Decision Tree Classifer
    clf = clf.fit(X_train, y_train)
    # Predict the response for test dataset
    y_pred = clf.predict(X_test)
    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
    return clf

def extract_rules(clf):
    pass

def rules_to_pycode():
    pass

def dt_entropy_depth(X_train, X_test, y_train, y_test, criteria = "entropy", depth = 1):
    # Create Decision Tree classifer object
    clf = DecisionTreeClassifier(criterion=criteria, max_depth= depth)
    # Train Decision Tree Classifer
    clf = clf.fit(X_train, y_train)
    # Predict the response for test dataset
    y_pred = clf.predict(X_test)
    # Model Accuracy, how often is the classifier correct?
    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
    return clf


def visualize_dt(clf, feature_cols):
    dot_data = StringIO()
    export_graphviz(clf, out_file=dot_data,
                    filled=True, rounded=True,
                    special_characters=True, feature_names=feature_cols, class_names=['0', '1'])
    graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
    graph.write_png('diabetes.png')
    Image(graph.create_png())

def program():
    col_names = ['pregnant', 'glucose', 'bp', 'skin', 'insulin', 'bmi', 'pedigree', 'age', 'label']
    # pima = pd.read_csv("data/diabetes.csv", header=None, names=col_names)
    pima = pd.read_csv("data/diabetes.csv", header= 1, names=col_names)
    print(pima.head())
    feature_cols = ['pregnant', 'insulin', 'bmi', 'age', 'glucose', 'bp', 'pedigree']
    X = pima[feature_cols]  # Features
    y = pima.label  # Target variable
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1)  # 70% training and 30% test
    c1 = dt_default(X_train, X_test, y_train, y_test)
    c2 = dt_entropy_depth(X_train, X_test, y_train, y_test)
    print(40*'=', 'auto generated code starts', 40*'=')
    tree_to_code(c2, feature_cols)
    print(40*'=', 'auto generated code ends', 40*'=')


def tree_to_code(tree, feature_names):
    '''
    Outputs a decision tree model as a Python function
    Parameters:
    -----------
    tree: decision tree model
        The decision tree to represent as a function
    feature_names: list
        The feature names of the dataset used for building the decision tree
    '''
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]
    print("def tree({}):".format(", ".join(feature_names)))

    def recurse(node, depth):
        indent = "  " * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            print("{} if {} <= {}:".format(indent, name, threshold))
            recurse(tree_.children_left[node], depth + 1)
            print("{} else:  # if {} > {}".format(indent, name, threshold))
            recurse(tree_.children_right[node], depth + 1)
        else:
            print("{}return {}".format(indent, tree_.value[node]))

    recurse(0, 1)

if __name__ == '__main__':
    program()