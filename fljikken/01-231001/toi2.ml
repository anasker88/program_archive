(*関数適用の優先度に注意*)
let twice f x = f (f x);;
(*再帰により繰り返しfを適用*)
let rec repeat f n x =
        if n=0 then x
        else f (repeat f (n-1) x)
