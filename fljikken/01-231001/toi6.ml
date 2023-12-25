(*以下２つの関数は問５のもの*)
let rec fold_right f xs e = match xs with
                                |[] -> e
                                |y :: ys -> f y (fold_right f ys e)
let rec fold_left f e xs = match xs with
                                |[] -> e
                                |y :: ys -> fold_left f (f e y) ys
(*リストxsと要素xを受け取り、xが関数fを満たすならリストに右から加える*)
let rec add_end_filter  f xs x = match xs with
                        |[] -> if (f x)  then [x] else [] 
                        |y ::ys -> y :: ( add_end_filter f  ys x )
(*fold_leftを用いてysの要素を左から順にxsの右に加える*)
let append_left xs ys = fold_left (add_end_filter (fun x-> true)) xs ys
(*fold_leftを用いてxsの要素を左からとりだし、fを満たすものを空リストに右から加えていく*)
let filter_left f xs =fold_left  (add_end_filter f)  [] xs

(*fold_rightを用いてxsの要素を右から順にysの左に加える*)
let append_right xs ys = fold_right (fun p qs -> p::qs) xs ys
(*fold_rightを用いてxsの要素を右からとりだし、fを満たすものを空リストに左から加えていく*)
let filter_right f xs =fold_right (fun  p qs -> if f p then p::qs else qs) xs []
