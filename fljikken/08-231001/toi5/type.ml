open Syntax
open Const_solver

exception Unbound

exception TypeErr

type constraints = (ty * ty) list

type type_schema = (tyvar list*ty)(*内部の自由変数のリストと自分の型*)

type type_env = (name * type_schema) list


let empty_tyenv = []
let extend x v env = (x, v) :: env

let lookup x tyenv =
  try List.assoc x tyenv with Not_found -> raise Unbound

let rec connect_as_set li1 li2 = match li1 with
    | [] -> li2
    | ele :: rest -> if List.mem ele li2 then connect_as_set rest li2 else connect_as_set rest (ele :: li2)



(*tyvarがあるtyの中に現れているか*)
let rec exist_in_expr_type ty tyvar =match ty with
    |TyInt -> false
    |TyBool -> false
    |TyFun (a,b) -> (exist_in_expr_type a tyvar) || (exist_in_expr_type b tyvar)
    |TyVar a -> a=tyvar

(*tyvarがtyenvの中に現れているか*)
let rec exist_in_tyenv (tyenv: type_env) tyvar = match tyenv with
| [] -> false
| (_,(_,a)):: rest -> (exist_in_expr_type a tyvar) || exist_in_tyenv rest tyvar

(*型環境をもとに型を型スキームに*)
let rec generize tyenv ty = match ty with
    | TyInt -> (([],TyInt):type_schema)
    | TyBool -> ([],TyBool)
    | TyFun (a,b)-> let (tvli1,ty1)=generize tyenv a in
                    let (tvli2,ty2)=generize tyenv b in
                    ((connect_as_set tvli1 tvli2),TyFun (ty1,ty2))
    | TyVar a -> if (exist_in_tyenv tyenv a) then ([],TyVar a ) else ([a],TyVar a) 

(*型スキームの各自由変数を新たな型変数に置き換え、型を生成*)
let rec instantiate (tysche:type_schema) = match tysche with
|([],ty)-> ty
|(a::rest,ty)->instantiate (rest,ty_subst [a,TyVar (new_tyvar())] ty)

let rec get_type_vars ty = match ty with
    | TyInt -> []
    | TyBool -> []
    | TyFun (a,b) -> connect_as_set (get_type_vars a) (get_type_vars b)
    | TyVar a -> [a]

let rec tyenv_subst subst (tyenv:type_env)= match tyenv with
| [] -> []
| (_,(_,ty))::rest -> (ty_subst subst ty):: tyenv_subst subst rest

let rec infer_expr (tyenv:type_env) e = match e with
  | EConstInt _ -> 
      ((TyInt),([]:constraints))
  | EConstBool _->
      (TyBool,[])
  | EVar      x ->
      (instantiate (try lookup x tyenv
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
      let sigma = unify con1 in
      let s1 = ty_subst sigma ty1 in
      let gamma= tyenv_subst sigma tyenv in
      let p = List.filter (fun x -> not (List.mem x (List.concat (List.map get_type_vars gamma)))) (get_type_vars s1) in
      let (ty2,con2)= infer_expr (extend x (p,s1) tyenv) e2 in
       (ty2,con1@con2)
  | EFun   (x,e) ->
      let a = new_tyvar () in
      let (ty,con)= infer_expr (extend x ([],(TyVar a)) tyenv) e in
       (TyFun(TyVar a,ty),con) 
  | EApp  (e1,e2) ->
      let (ty1,con1)= infer_expr tyenv e1 in 
      let (ty2,con2)= infer_expr tyenv e2 in
      let a = new_tyvar () in
       (TyVar a,(ty1,TyFun(ty2,TyVar a))::con1@con2)    
  | ELetRec  (f,x,e1,e2) ->
      let a = new_tyvar () in
      let b = new_tyvar () in
      let (ty1,con1)= infer_expr (extend x ([], (TyVar a)) (extend f ([],(TyFun (TyVar a,TyVar b))) tyenv)) e1 in
      let sigma = unify ((ty1,(TyVar b))::con1) in
      let s1 = ty_subst sigma (TyFun (TyVar a,TyVar b)) in
      let gamma= tyenv_subst sigma tyenv in
      let p = List.filter (fun x -> (not (List.mem x (List.concat (List.map get_type_vars gamma))))) (get_type_vars s1) in
      let (ty2,con2)= infer_expr (extend f (p,s1) tyenv) e2 in
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
    (fin_ty,extend x (generize tyenv fin_ty) tyenv)
  | CRecDecl (id,x,e) ->
    let a = new_tyvar () in
    let b = new_tyvar () in
    let (ty,con)= infer_expr (extend x ([],(TyVar a)) (extend id ([],(TyFun(TyVar a,TyVar b))) tyenv)) e in
    let sbst= unify ((TyVar b,ty)::con) in
    let fin_ty= ty_subst sbst (TyFun(TyVar a,TyVar b)) in
    ( fin_ty, extend id (generize tyenv fin_ty) tyenv)