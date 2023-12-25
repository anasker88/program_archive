type tyvar = int

and ty =
  | TyInt
  | TyBool
  | TyFun of ty * ty
  | TyVar of tyvar

  let new_tyvar = let x= ref 0 in fun () -> (x := !x+1;!x)

let rec print_type ty = match ty with
| TyInt -> print_string "Int "
| TyBool -> print_string "Bool "
| TyFun (ty1,ty2)-> print_string "( ";print_type ty1;print_string "-> "; print_type ty2;print_string(") ")
| TyVar num -> print_string("a");print_int(num)