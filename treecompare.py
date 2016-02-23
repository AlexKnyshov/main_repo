from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
import sys
# store and return all _BitStrings
def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]  
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs
    
def compare(tree1, tree2):
    term_names1 =  [term.name for term in tree1.get_terminals()]
    term_names2 =  [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same 
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = '%d_%s' % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

infile1 = sys.argv[1]
infile2 = sys.argv[2]
tree1 = Phylo.read(infile1, "newick")
tree2 = Phylo.read(infile2, "newick")
print compare(tree1, tree2)

query = "Chinannus"#"Chinannus_sp2_"
print "query:", query

print "searching tree1..."
cladelist = []
for cladename in tree1.find_clades():
    for term in cladename.get_terminals():
        if query in term.name:
            cladelist.append(term)
cladeset = list(set(cladelist))

#tree2
print "searching tree2..."
cllst2 = []
for cladename in tree2.find_clades():
    for term in cladename.get_terminals():
        if query in term.name:
            cllst2.append(term)
cladeset2 = list(set(cllst2))


clade1 = tree1.common_ancestor(cladeset)
clade2 = tree2.common_ancestor(cladeset2)

#Phylo.draw_ascii(clade1)
#Phylo.draw_ascii(clade2)

print "before collapse:", compare(clade1, clade2)

print "collapsing..."
clade1.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)
clade2.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)


#Phylo.draw_ascii(clade1)
#Phylo.draw_ascii(clade2)

#compare clades
print "after collapse:", compare(clade1, clade2)
