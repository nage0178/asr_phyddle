import dendropy
import sys
import pandas as pd

treeFile = sys.argv[1]
outTreeFile = sys.argv[2]
outASRFile = sys.argv[3]
outTreeFinal = sys.argv[4]

# Example Newick string for a tree with unlabeled internal nodes
# (A:0.1, B:0.2):0.3, (C:0.4, D:0.5):0.6);
# Or load from a file:
tree = dendropy.Tree.get(
    path=treeFile,
    schema="nexus")

# A counter for generic labels (optional, you could use specific names)
internal_node_counter = 1

# Iterate through all nodes in the tree
for node in tree.postorder_node_iter():
    # Check if the node is an internal node (not a leaf node)
    if not node.is_leaf():
        # Assign a label to the internal node
        node.label = f"Node_{internal_node_counter}"
        internal_node_counter += 1

## Print the tree in Newick format to verify labels
#print(tree.as_string(schema="newick"))

# Write the tree to a file
tree.write(
    path= outTreeFile,
    schema="nexus")

name = []
loc = []
for node in tree.postorder_node_iter():
    if not node.is_leaf():
        name.append(node.label)
        for annotation in node.annotations:
            if annotation.name == 'location':
                loc.append(annotation.value[0])
df = pd.DataFrame(loc, name)
   # if node.is_leaf():
   #     print(node.taxon)
   #     for annotation in node.annotations:
   #         if annotation.name == 'location':
   #             print(annotation.value)
df.to_csv(outASRFile, header = False )


reroot=None
for node in tree.postorder_node_iter():
    if len(node.child_nodes()) == 1 and node.is_internal():
        child = node.child_nodes()[0]
        parent = node.parent_node

        if parent: 
            parent.remove_child(node)
            parent.add_child(child)
            child.edge_length += node.edge_length
        else: 
            reroot=node
if reroot: 
    tree.reroot_at_node(reroot, update_bipartitions=False)

    for node in tree.postorder_node_iter():
        if len(node.child_nodes()) == 1 and node.is_internal():
            child = node.child_nodes()[0]
            parent = node.parent_node
    
            if parent: 
                parent.remove_child(node)
                parent.add_child(child)
                child.edge_length += node.edge_length

tree.write(path = outTreeFinal, schema= "nexus")
