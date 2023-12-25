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

let rec compose sbt1 sbt2 =
  match sbt2 with
  | [] -> sbt1
  | (tyv,ty) :: rest-> (tyv,ty_subst sbt1 ty) :: List.filter (fun (tyvar,_)-> tyvar!=tyv) (compose sbt1 rest)


let rec replace cons (a,t)=match cons with 
| [] -> []
| (ty1,ty2)::rest -> let fin_ty1 = ty_subst [a,t] ty1 in
  let fin_ty2 = ty_subst [a,t] ty2 in
  (fin_ty1,fin_ty2) ::replace rest (a,t)

exception ConstraintConflict


(*tyvarがあるtyの中に現れているか*)
let rec exist_in_expr_type ty tyvar =match ty with
    |TyInt -> false
    |TyBool -> false
    |TyFun (a,b) -> (exist_in_expr_type a tyvar) || (exist_in_expr_type b tyvar)
    |TyVar a -> a=tyvar

let rec unify cons = match cons with 
|[] -> []
| (ty1,ty2) :: rest -> 
  if ty1=ty2 then unify rest
  else match (ty1,ty2) with
  | (TyFun (s1,t1),TyFun (s2,t2)) -> unify ((s1,s2)::(t1,t2)::rest)
  | (TyVar a,t) -> if not (exist_in_expr_type t a) then compose (unify (replace rest (a,t))) [(a,t)] else raise ConstraintConflict
  | (t,TyVar a) -> if not (exist_in_expr_type t a) then compose (unify (replace rest (a,t))) [(a,t)] else raise ConstraintConflict
  | _ -> raise ConstraintConflict