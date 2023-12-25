type order = LT | EQ | GT

module type ORDERED_TYPE =
sig
  type t
  val compare : t -> t -> order
end

module type MULTISET2 =
  functor (T : ORDERED_TYPE) ->
    sig
      type t
      val empty : t
      val add    : T.t -> t -> t
      val remove : T.t -> t -> t
      val count  : T.t -> t -> int
    end

exception GetMinError(*getminで対象がLeafならerror*)

module AbstMultiset2 : MULTISET2 =
  functor (T : ORDERED_TYPE) -> struct

    type t =(*二分探索木。葉もしくは左右の子を持つノード*)
      |Leaf
      |Node of T.t *t * t 
    let empty = Leaf

    let rec add ele tree = match tree with
      |Leaf -> Node (ele,Leaf,Leaf)(*末端ならノードを追加*)
      |Node (cur,left,right) -> (match (T.compare cur ele) with(*それ以外なら比較し、適当な方に挿入*)
            |GT ->Node (cur, (add ele left),right)
            |_  ->Node (cur, left,(add ele right)))

    let rec dltmin tree = match tree with(*最小ノードを消去*)
            | Leaf -> Leaf(*葉を引数に取るのは想定外、とりあえずLeafを返す*)
            | Node (cur, Leaf, right) -> right(*左が葉ならそこで最小*)
            | Node (cur, left, right) -> Node (cur, dltmin left, right)(*それ以外なら左に潜る*)
          
    let rec getmin tree = match tree with(*上と同様にして最小ノードを取得*)
            | Node (cur, Leaf, right) -> cur
            | Node (cur, left, right) -> getmin left
            | Leaf -> raise GetMinError
            
    let rec remove ele tree = match tree with(*treeからeleを１つ削除*)
    | Leaf -> Leaf(*見つからないならそのまま*)
    |Node (cur, left, right) -> (match (T.compare cur ele) with
    |GT ->Node (cur, (remove ele left),right)(*左へ探索*)
    |LT ->Node (cur, left,(remove ele right))(*右へ探索*)
    |EQ -> if left=Leaf then right 
           else if right=Leaf then left(*子のいずれかが葉ならもう一方で置換して終わり*)
           else Node (getmin left, dltmin left, right))(*それ以外なら右から最小ノードを除去し、その最小ノードで置換*)

    let rec count ele tree = match tree with
    | Leaf -> 0(*末端まで来たら0*)
    |Node (cur, left, right)-> (if cur=ele then 1 else 0) + count ele left + count ele right(*自分を判定しながら順に潜っていく*)
  end
