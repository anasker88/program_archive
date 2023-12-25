let rec map f xs =
        match xs with 
        | [] -> []
        |y :: ys -> f y :: map f ys;;
