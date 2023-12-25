(*xsの左端を取り出し、それ以外とysの結合の左端に加える*)
let rec append xs ys = match xs with
        |[] -> ys
        |z :: zs -> z :: (append zs ys);;
(*xsの左端を取り出し、それ以外にfilterを適用したものに状況に応じて加える*)
let rec filter f xs=match xs with 
        |[] -> []
        |y :: ys -> if f y then y :: (filter f ys) else (filter f ys);;
