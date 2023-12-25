open Syntax

exception Unbound


let empty_env = []
let extend x v env = (x, v) :: env

let lookup x env =
  try List.assoc x env with Not_found -> raise Unbound

exception EvalErr

let rec eval_expr (env:env) e =
  match e with
  | EConstInt i -> 
      VInt i
  | EConstBool b ->
      VBool b
  | EVar x ->
      let Thunk(e,oenv)=(try
         lookup x env
       with
       | Unbound -> raise EvalErr) in (eval_expr !(oenv) e)
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
        let oenv=ref[] in
        oenv:=env;
        eval_expr (env |> extend x (Thunk(e1,oenv))) e2
  | EFun (x,e) ->
        VFun (x,e,env)
  | ELetRec (f,x,e1,e2) ->
        let oenv = ref[] in
        oenv:=extend f (Thunk(EFun (x,e1),oenv)) env;
        eval_expr (extend f (Thunk(EFun (x,e1),oenv)) env) e2 
        
  | EApp (e1,e2) ->
        let v1 = eval_expr env e1 in
        let env'=ref[] in
        env':=env;
        (match v1 with
          | VFun (x,e,oenv) ->eval_expr (extend x (Thunk(e2,env')) oenv) e
          | _ -> raise EvalErr)
  | EPair (e1,e2) ->
        let v1 = eval_expr env e1 in
        let v2 = eval_expr env e2 in
        VPair (v1,v2)
  | ENil -> Vlist []
  | ECons (e1,eli) -> let v = eval_expr env e1 in
    match eval_expr env eli with
    Vlist vli ->Vlist (v::vli)
    | _ -> raise EvalErr

    

let eval_command env c =
  match c with
  | CExp e -> ("-", env, eval_expr env e)
  | CDecl (x,e) -> let env'=ref[] in env':=env;("val "^x,env |> extend x (Thunk(e,env')), eval_expr env e)
  | CRecDecl (id,x,e) ->let oenv=ref[] in oenv:=extend id (Thunk(EFun (x,e),oenv)) env; ("val "^id,!(oenv), eval_expr (!(oenv)) (EFun (x,e)))