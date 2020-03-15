from collections import defaultdict

class Node:
    def __init__(self, character, is_leaf = False):
        self.children = dict()
        self.character = character
        self.is_leaf = is_leaf


class TrieStructure:
    def __init__(self):
        self.root = Node("")

    def add(self, word):
        cur = self.root
        for character in word:
            if character not in cur.children:
                cur.children[character] = Node(character)
            cur = cur.children[character]
        cur.is_leaf = True

    def search(self, word, node=None):
        cur = node
        if not cur:
            cur = self.root
        for i, character in enumerate(word):
            if character == "*":
                if i == len(word) - 1:
                    for child in cur.children.values():
                        if child.is_terminal:
                            return True
                    return False
                for child in cur.children.values():
                    if self.search(word[i+1:], child) == True:
                        return True
                return False
            if character not in cur.children:
                return False
            cur = cur.children[character]
        return cur.is_leaf


def test_Trie():
    t =TrieStructure()
    t.add('hello')
    t.add('help')
    t.add('help')
    t.add('herman')
    t.add('german')
    t.add('violin')
    t.add('umbrella')

    z = t.search('help')
    print(z)
    z = t.search('hell')
    print(z)
    z = t.search('***lin')
    print(z)
    z = t.search('world')
    print(z)

if __name__ == '__main__':
    test_Trie()