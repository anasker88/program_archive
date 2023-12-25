type 'a tree =
        | Leaf
        | Node of 'a * 'a tree * 'a tree
(*行きがけ順の関数。自分、左の子、右の子の順に探索*)
let rec pre_order node=match node with
                        |Leaf -> []
                        |Node (self, ch1, ch2) -> self :: (pre_order ch1 @ pre_order ch2)
(*通りがけ順の関数。左の子、自分、右の子の順に探索*)
let rec in_order node=match node with
                        |Leaf -> []
                        |Node (self, ch1, ch2) -> in_order ch1 @ (self :: pre_order ch2)
(*帰りがけ順の関数。左の子、右の子、自分の順に探索*)
let rec post_order node=match node with
                        |Leaf -> []
                        |Node (self, ch1, ch2) -> (post_order ch1 @ pre_order ch2) @ [self]

