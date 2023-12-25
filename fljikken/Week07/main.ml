open TySyntax
open ConstraintSolver

let alphax = new_tyvar ()
let alphaf = new_tyvar ()
let alphag = new_tyvar ()
let alphagx = new_tyvar ()
let alphafgx = new_tyvar ()
let myconst =[(TyVar alphag,TyFun(TyVar alphax,TyVar alphagx));(TyVar alphaf,TyFun(TyVar alphagx,TyVar alphafgx))]

let sbst= unify myconst
let ty = (ty_subst sbst (TyFun(TyVar alphaf,TyFun(TyVar alphag,TyFun(TyVar alphax,TyVar alphafgx)))))
let () =
  print_type ty;
  print_endline ""
