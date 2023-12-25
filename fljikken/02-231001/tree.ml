type 'a tree =
        | Leaf
        | Node of 'a * 'a tree * 'a tree

let node3=Node (3,Leaf,Leaf)
let node4=Node (4,Leaf,Leaf)
let node5=Node (5,Leaf,Leaf)
let node6=Node (6,Leaf,Leaf)
let node1=Node (1,node3,node4)
let node2=Node (2,node5,node6)
let node0=Node (0,node1,node2)
