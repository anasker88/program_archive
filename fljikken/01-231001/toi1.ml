(*帰納的に加算を定義*)
let rec sum_to x =
        if x=0 then 0
        else x + sum_to (x - 1) ;;
(*2以上y以下の整数でｘが割り切れるか*)
let rec dividable x y=
        if y=1 then false
        else x mod y = 0 || dividable x (y-1);;
(*x-1までで割り切れなければ素数*)
let is_prime x = 
        if x=1 then false
        else not (dividable x (x-1));;
(*ユークリッドの互除法*)
let rec gcd x y =
        if x > y then
                if x mod y = 0 then y
                else gcd (x mod y) y
        else 
                if y mod x = 0 then x
                else gcd x (y mod x);;

