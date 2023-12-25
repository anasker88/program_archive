type order = LT | EQ | GT

module type ORDERED_TYPE =
sig
  type t
  val compare : t -> t -> order
end

exception GetMinError

module MakeMap =
  functor (T : ORDERED_TYPE) -> struct
    type t =(*keyに関する二分探索木。葉もしくは左右の子を持つノード*)
      |Leaf
      |Node of T.t * string * t * t 
    let empty = Leaf

    let rec add mykey ele tree = match tree with
      |Leaf -> Node (mykey,ele,Leaf,Leaf)(*末端ならノードを追加*)
      |Node (key,cur,left,right) -> (match (T.compare key mykey) with(*それ以外なら比較し、適当な方に挿入*)
            |GT ->Node (key,cur, (add mykey ele left),right)
            |LT ->Node (key,cur, left,(add mykey ele right))
            |EQ ->Node (key,ele, left, right))(*すでにあるなら上書き*)

    let rec dltmin tree = match tree with(*最小ノードを消去*)
            | Leaf -> Leaf(*葉を引数に取るのは想定外、とりあえずLeafを返す*)
            | Node (key,cur, Leaf, right) -> right(*左が葉ならそこで最小*)
            | Node (key,cur, left, right) -> Node (key,cur, dltmin left, right)(*それ以外なら左に潜る*)
          
    let rec getmin tree = match tree with(*上と同様にして最小ノードを取得*)
            | Node (key,cur, Leaf, right) -> key
            | Node (key,cur, left, right) -> getmin left
            | Leaf -> raise GetMinError
            
    let rec remove ele tree = match tree with(*treeからeleを１つ削除*)
    | Leaf -> Leaf(*見つからないならそのまま*)
    |Node (key,cur, left, right) -> (match (T.compare key ele) with
    |GT ->Node (key,cur, (remove ele left),right)(*左へ探索*)
    |LT ->Node (key,cur, left,(remove ele right))(*右へ探索*)
    |EQ -> if left=Leaf then right 
           else if right=Leaf then left(*子のいずれかが葉ならもう一方で置換して終わり*)
           else Node (getmin left,cur, dltmin left, right))(*それ以外なら右から最小ノードを除去し、その最小ノードで置換*)


exception Not_found
    let rec lookup mykey tree = match tree with
    | Leaf -> raise Not_found(*末端まで来たらNot_found*)
    |Node (key,cur, left, right)-> (match (T.compare key mykey) with
        |GT -> lookup mykey left
        |LT -> lookup mykey right
        |EQ -> cur)
  end
