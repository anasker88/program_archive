type complex = { re : float; im : float; }
let prod x y = {re = x.re *. y.re -. x.im *. y.im ; im =x.re *. y.im +. x.im *. y.re}
type str_tree =
        | Leaf
        | Node of string * str_tree * str_tree
let leaf1=Leaf
let node1=Node ("hello",Leaf,Leaf)
let node2=Node("konnichiha",Leaf,leaf1)
type ib_list = INil| ICons of int * bi_list
 and bi_list = BNil| BCons of bool * ib_list;;
let ib1=INil
let ib2=ICons (2,BNil)
let ib3=ICons (3,BCons(true, INil))
