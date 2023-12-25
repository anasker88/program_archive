open Syntax

exception Unbound


let empty_env = []
let extend x v env = (x, v) :: env

let lookup x env =
  try List.assoc x env with Not_found -> raise Unbound

exception EvalErr

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
        VFun (x,e,ref env)
  | ELetRec(letlst,e2) ->
        let oenv=ref [] in
        let rec defenv lst curenv = match lst with
        | [] -> curenv
        | (fi,xi,ei) :: rest -> let vi = VFun(xi,ei,oenv) in (defenv rest (extend fi vi curenv)) in
            (oenv := defenv letlst env;
            eval_expr (defenv letlst env) e2)
  | EApp (e1,e2) ->
        match eval_expr env e1 with
        |VFun (x,e,oenv)  ->
            let v2 = eval_expr env e2 in
            eval_expr (extend x v2 (!oenv)) e
        |_ -> raise EvalErr

let eval_command env c =
  match c with
  | CExp e -> (["-",eval_expr env e],env)
  | CDecl (x,e) -> let v=eval_expr env e in(["val "^x, eval_expr env e],env |> extend x v)
  | CRecDecl letlst -> let oenv = ref [] in
   let rec defret lst = match lst with
  | [] -> ([],env)
  | (fi,xi,ei) :: rest -> let vi = VFun(xi,ei,oenv) in 
  match defret rest with
  | (setlst,newenv) -> (("val "^fi,vi)::setlst,extend fi vi newenv) in
    match defret letlst with 
    |(setlst,newenv) -> oenv :=newenv;(setlst,newenv)