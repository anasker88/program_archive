let rec sigma f x =
        if x=0 then f x
        else
                sigma  f ( x- 1 )   + f x;;
