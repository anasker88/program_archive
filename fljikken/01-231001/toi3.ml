let rec fix f x = f (fix f) x;;
(*0以外なら前までの和とxの和、0なら0で確定*)
let sum_to_fix x = fix (fun f x -> if x=0 then 0 else x + f (x-1) ) x;;
(*xが２以上y以下の整数で割り切れるか否かの関数。yが１ならfalseで確定、それ以外なら(y-1まで見た時の真偽)||(x|y)*)
let dividable_fix x y = fix (fun f y -> if y=1 then false else x mod y = 0 || f (y-1)) y;;
(*x-1まででxが割り切れないなら素数*)
let is_prime_fix x = if x= 1 then false else not ( dividable_fix x (x-1));;
(*ユークリッドの互除法。割り切れるなら確定、そうでないなら余りを計算して再度適用*)
let gcd_fix x y  = fix (fun f x y ->
        if x>y then (if x mod y = 0 then y else f (x mod y)  y)
        else (if y mod x = 0 then x else f x ( y mod x))) x y;;
