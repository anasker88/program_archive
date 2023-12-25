module type SEMIRING = 
sig
  type t
  val add : t -> t -> t
  val mul : t -> t -> t
  val unit : t
  val zero : t
end

exception DifferentSize
module Matrix (T: SEMIRING)=
struct
  type vec = T.t list(*ベクトルはリストで表現*)
  type mat = T.t list list(*行列はリストのリストで表現*)
  let rec add_vec (v1 : vec) (v2 : vec) = match v1,v2 with(*ベクトルの足し算*)
    | [],[] -> [] 
    | x :: xs , y :: ys -> (T.add x y) :: (add_vec xs ys)
    | _ , _ -> raise DifferentSize
  let rec mul_vec (v1 : vec) (v2 : vec) = match v1,v2 with(*ベクトルの掛け算*)
    | [],[] -> T.zero    | x :: xs , y :: ys -> T.add (T.mul x y)  (mul_vec xs ys)
    | _ , _ -> raise DifferentSize
  let rec add_mat (m1 : mat) (m2 : mat) = match m1,m2 with
    | [],[] -> [] 
    | vx :: vxs , vy :: vys -> (add_vec vx vy) :: (add_mat vxs vys)
    | _ , _ -> raise DifferentSize
  
  let rec trans (ma : mat) = match ma with(*行列を転置*)
      | [] -> []
      | [] :: _ -> []
      | m -> (List.map List.hd m) :: trans (List.map List.tl m) (*一列目の転置を残りの転置の上につける*)
  let rec mul_vecmat (ve : vec) (ma : mat) = match (trans ma) with(*ベクトルと行列の掛け算*)
| [] -> []
|  vx :: vxs -> (mul_vec ve vx) :: (mul_vecmat ve (trans vxs))
  let rec mul_mat (m1 : mat) (m2 : mat) = match m1 with (*m1を1行ずつm2とかけて行列の積を計算*)
  | [] -> []
  | vx :: vxs -> mul_vecmat vx m2 :: (mul_mat vxs m2)
end

module BoolSemiring =(*Boolの半環*)
struct
  type t =bool
  let add  x y= x || y (*和はor*)
  let mul x y = x && y (*積はand*)
  let zero = false
  let unit = true
end

module BoolMatrix = Matrix(BoolSemiring) (*Bool行列*)

module TropSemiring =(*Tropの半環*)
struct
  type t = INF| NUM of int
  let add  x y = match x,y with(*和はmin*)
    | INF , q -> q
    | p , INF -> p
    | NUM p, NUM q ->NUM (min p q) 
  let mul x y = match x,y with
  | INF , _ -> INF
  | _ , INF -> INF
  | NUM p , NUM q -> NUM ( p + q ) (*積は和*)
  let zero = INF
  let unit = NUM 0
end

module TropMatrix = Matrix(TropSemiring) (*Trop行列*)