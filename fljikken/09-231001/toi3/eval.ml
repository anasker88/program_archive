open Syntax

exception Unbound


let empty_env = []
let extend x v env = (x, v) :: env

let lookup x env =
  try List.assoc x env with Not_found -> raise Unbound

exception EvalErr


let rec find_match pat va = match pat with
| PInt i -> (match va with
  | VInt j -> if i=j then Some [] else None
  | _ -> None)
| PBool i -> (match va with
  | VBool j -> if i=j then Some [] else None
  | _ -> None)
|PVar s -> Some [s,va]
|PPair (p1,p2)-> (match va with
  |VPair  (v1,v2) -> (match (find_match p1 v1),(find_match p2 v2)  with
    | None,_-> None
    | _,None -> None
    | Some a ,Some b -> Some (a@b))
  | _ -> None)
|PNil->(match va with
  | Vlist [] -> Some []
  | _ -> None)
|PCons (c1,c2)-> (match va with
  |Vlist  (v1::v2) -> (match (find_match c1 v1),(find_match c2 (Vlist v2))  with
    | None,_-> None
    | _,None -> None
    | Some a ,Some b -> Some (a@b))
  | _ -> None)

let rec eval_expr env e =
  match e with
  | EConstInt i -> 
      VInt i
  | EConstBool b ->
      VBool b
  | EVar x ->
      (try
         lookup x env
       with
       | Unbound -> raise EvalErr)
  | EAdd (e1,e2) ->
      let v1 = eval_expr env e1 in
      let v2 = eval_expr env e2 in
      (match v1, v2 with
       | VInt i1, VInt i2 -> VInt (i1 + i2)
       | _ -> raise EvalErr)
  | ESub (e1,e2) ->
      let v1 = eval_expr env e1 in
      let v2 = eval_expr env e2 in
      (match v1, v2 with
       | VInt i1, VInt i2 -> VInt (i1 - i2)
       | _ -> raise EvalErr)
   | EMul (e1,e2) ->
      let v1 = eval_expr env e1 in
      let v2 = eval_expr env e2 in
      (match v1, v2 with
       | VInt i1, VInt i2 -> VInt (i1 * i2)
       | _ -> raise EvalErr)
  | EDiv (e1,e2) ->
        let v1 = eval_expr env e1 in
        let v2 = eval_expr env e2 in
        (match v1, v2 with
         | _,VInt 0 -> raise EvalErr
         | VInt i1, VInt i2 -> VInt (i1 / i2)
         | _ -> raise EvalErr)
  | EEq (e1,e2) ->
      let v1 = eval_expr env e1 in
      let v2 = eval_expr env e2 in
      (match v1, v2 with
       | VInt i1,  VInt i2  -> VBool (i1 = i2)
       | VBool b1, VBool b2 -> VBool (b1=b2)
       | _ -> raise EvalErr)
  | ELt (e1,e2) ->
      let v1 = eval_expr env e1 in
      let v2 = eval_expr env e2 in
      (match v1, v2 with
       | VInt i1,  VInt i2  -> VBool (i1 < i2)
       | _ -> raise EvalErr)
  | EIf (e1,e2,e3) ->
      let v1 = eval_expr env e1 in
      (match v1 with
       | VBool b ->
           if b then eval_expr env e2 else eval_expr env e3
       | _ -> raise EvalErr)
  | ELet (x,e1,e2) ->
        let v1 = eval_expr env e1 in
        eval_expr (env |> extend x v1) e2
  | EFun (x,e) ->
        VFun (x,e,env)
  | ELetRec (f,x,e1,e2) ->
        let env' = extend f (VRecFun (f,x,e1,env)) env 
            in eval_expr env' e2
  | EApp (e1,e2) ->
        let v1 = eval_expr env e1 in
        let v2 = eval_expr env e2 in
        (match v1 with
          | VFun (x,e,oenv) ->eval_expr (extend x v2 oenv) e
          | VRecFun (f,x,e,oenv) ->
                let env'=extend x v2 (extend f (VRecFun (f,x,e,oenv)) oenv) 
                    in eval_expr env' e
          | _ -> raise EvalErr)
  | EPair (e1,e2) ->
        let v1 = eval_expr env e1 in
        let v2 = eval_expr env e2 in
        VPair (v1,v2)
  | ENil -> Vlist []
  | ECons (e1,eli) -> let v = eval_expr env e1 in
    (match eval_expr env eli with
    Vlist vli ->Vlist (v::vli)
    | _ -> raise EvalErr)
  |EMatch (e1,patli)->
    let v1=eval_expr env e1 in
    (match patli with
    | [] -> raise EvalErr
    | (pat,e2)::rest ->(match find_match pat v1 with
        |None ->eval_expr env (EMatch (e1,rest))
        |Some addenv -> eval_expr (addenv@env) e2)
  )

    

let eval_command env c =
  match c with
  | CExp e -> ("-", env, eval_expr env e)
  | CDecl (x,e) -> let v=eval_expr env e in("val "^x,env |> extend x v, eval_expr env e)
  | CRecDecl (id,x,e) ->let v= VRecFun(id,x,e,env) in("val "^id,env |> extend id v, VRecFun(id,x,e,env))