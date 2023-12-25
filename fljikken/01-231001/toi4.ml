(*空リストなら初期値eを返し、そうでないなら左端を除いたリストを与えた時の出力をfに右から与える*)
let rec fold_right f xs e = match xs with 
                                |[] -> e
                                |y :: ys -> f y (fold_right f ys e);;
(*空リストなら初期値eを返し、そうでないなら初期値をf e (左端の要素)に書き換え、左端を除いて改めて計算*)
let rec fold_left f e xs = match xs with
                                |[] -> e
                                |y :: ys -> fold_left f (f e y) ys;;
