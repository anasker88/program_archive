type tyvar = int

and ty =
  | TyInt
  | TyBool
  | TyFun of ty * ty
  | TyVar of tyvar


  type subst = (tyvar*ty) list 
  let new_tyvar = let x= ref 0 in fun () -> (x := !x+1;!x)

let rec print_type ty = match ty with
| TyInt -> print_string "Int "
| TyBool -> print_string "Bool "
| TyFun (ty1,ty2)-> print_string "( ";print_type ty1;print_string "-> "; print_type ty2;print_string(") ")
| TyVar num -> print_string("a");print_int(num)

let rec ty_subst (sbt:subst) t=
  match t with
  | TyInt -> TyInt
  | TyBool ->TyBool
  | TyFun (a,b) -> TyFun (ty_subst sbt a,ty_subst sbt b)
  | TyVar num -> try List.assoc num sbt 
      with Not_found -> TyVar num


