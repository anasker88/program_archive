type 'a tree =
        | Leaf
        | Node of 'a * 'a tree * 'a tree
(*末尾再帰で定義可能*)
let level_order node =
        (*loopは探索済みノードのリスト,探索待ちのノードキューのリストをタプルで管理*)
        (*現在地のノードを探索済みにし、子を探索待ちキューに追加*)
        let rec loop result_queue = match result_queue with
                                        |(rsl,[]) -> rsl
                                        |(rsl,(Node (self,ch1,ch2))::que) -> loop (rsl@[self] , if ch1=Leaf &&ch2=Leaf then que
                                                else if ch1=Leaf then que@[ch2]
                                                else if ch2=Leaf then que@[ch1]
                                                else que@(ch1 ::[ch2]))
                                        |_ -> []
        in loop ([],[node])
