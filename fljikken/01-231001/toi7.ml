(*リストの各要素に関数ｆを適用する関数map*)
let rec map f xs =
        match xs with
        | [] -> []
        |y :: ys -> f y :: map f ys;;
(*リストのリストを受け取り、それらを結合して一つのリストにする関数concat*)
let rec concat lst=
        match lst with
        | [[]] -> []
        | x :: y -> x @ (concat y)
        | _ -> []
(*リストと要素を受け取り、リストのどこかに要素ｘを挿入したリストすべてのリストを返す関数insert_all*)
let rec insert_all x lst =
  match lst with
  | [] -> [[x]]
  | y :: ys ->
      (x :: lst) :: map (fun l -> y :: l) (insert_all x ys)

(*左端の要素xを取り出し、それ以外のリストxsの順列にxを挿入した全パターンを求める*)
let rec perm lst =
  match lst with
  | [] -> [[]]
  | x :: xs ->concat (map (insert_all x) (perm xs))
