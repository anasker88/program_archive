(*Zは0を表し、Sは+1を表す*)
type nat = Z | S of nat
(*x+y=((x-1)+y)+1*)
let rec add x y = match x with
                | Z -> y
                | S z -> S (add z y)
(*x-1を求める関数des*)
let des x = match x with
                | Z -> Z
                | S y -> y
(*x-y=(x-1)-(y-1)*)
let rec sub x y = match y with
                | Z -> x
                | S z -> sub (des x) z
(* x*y=y+(x-1)*y *)
let rec mul x y = match x with
                | Z -> Z
                | S z -> add y (mul z y)
(* x^y=x*x^(y-1) *)
let rec pow x y = match y with
                | Z -> S(Z)
                | S z -> mul x (pow x z)
(*nat->intの変換。Zなら0に、そうでないならSを取り除き、代わりに結果に１を足す*)
let rec n2i x = match x with
                | Z -> 0
                | S y -> 1 + n2i (y)
(*int->natの変換。負ならZにし、それ以外なら1を引いて代わりに結果にSを作用させる*)
let rec i2n x = if x <= 0 then Z else S (i2n (x-1))

