open Syntax
open TySyntax
open ConstraintSolver

exception Unbound

exception TypeErr

type constraints = (ty * ty) list

type tyenv = (name * ty) list


let empty_tyenv = []
let extend x v env = (x, v) :: env

let lookup x tyenv =
  try List.assoc x tyenv with Not_found -> raise Unbound


let rec infer_expr tyenv e = match e with
  | EConstInt _ -> 
      (TyInt,[])
  | EConstBool _->
      (TyBool,[])
  | EVar      x ->
      ((try lookup x tyenv
      with
      | Unbound -> raise TypeErr),[])
  | EAdd  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyInt,(ty1,TyInt)::(ty2,TyInt)::con1@con2)
  | ESub  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyInt,(ty1,TyInt)::(ty2,TyInt)::con1@con2)
  | EMul  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyInt,(ty1,TyInt)::(ty2,TyInt)::con1@con2)
  | EDiv  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyInt,(ty1,TyInt)::(ty2,TyInt)::con1@con2)
  | EEq   (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyBool,(ty1,ty2)::con1@con2)
  | ELt   (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      (TyBool,(ty1,TyInt)::(ty2,TyInt)::con1@con2)
  | EIf   (e1,e2,e3) ->
    let (ty1,con1)= infer_expr tyenv e1 in 
    let (ty2,con2)= infer_expr tyenv e2 in
    let (ty3,con3)= infer_expr tyenv e3 in
      (ty2,(ty1,TyBool)::(ty2,ty3)::con1@con2@con3)
  | ELet  (x,e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr (extend x ty1 tyenv) e2 in
       (ty2,con1@con2)
  | EFun   (x,e) ->
      let a = new_tyvar () in
      let (ty,con)= infer_expr (extend x (TyVar a) tyenv) e in
       (TyFun(TyVar a,ty),con) 
  | EApp  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      let a = new_tyvar () in
       (TyVar a,(ty1,TyFun(ty2,TyVar a))::con1@con2)    
  | ELetRec  (f,x,e1,e2) ->
      let a = new_tyvar () in
      let b = new_tyvar () in
      let (ty1,con1)= infer_expr (extend x (TyVar a) (extend f (TyFun (TyVar a,TyVar b)) tyenv)) e1 in
      let (ty2,con2)= infer_expr (extend f (TyFun(TyVar a,TyVar b)) tyenv) e2 in
      (ty2,(ty1,(TyVar b))::con1@con2)

      
let infer_cmd tyenv c =
  match c with
  | CExp e -> 
    let (ty,con)=infer_expr tyenv e in 
    let sbst=unify con in
    let fin_ty = ty_subst sbst ty in
    (fin_ty,tyenv)
  | CDecl (x,e) ->
    let (ty,con)=infer_expr tyenv e in 
    let sbst=unify con in
    let fin_ty = ty_subst sbst ty in
    (fin_ty,extend x fin_ty tyenv)
  | CRecDecl (id,x,e) ->
    let a = new_tyvar () in
    let b = new_tyvar () in
    let (ty,con)= infer_expr (extend x (TyVar a) (extend id (TyFun(TyVar a,TyVar b)) tyenv)) e in
    let sbst= unify ((TyVar b,ty)::con) in
    let fin_ty= ty_subst sbst (TyFun(TyVar a,TyVar b)) in
    ( fin_ty, extend id fin_ty tyenv)