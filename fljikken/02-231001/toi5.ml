type value = VInt of int | VBool of bool

exception Eval_error;;

type expr =
  | EConstInt of int(*数値*)
  | EAdd of expr * expr(* E + E *)
  | ESub of expr * expr(* E - E *)
  | EMul of expr * expr(* E * E *) 
  | EDiv of expr * expr(* E / E *)
  | EConstBool of bool(* true , false *)
  | EEqual of expr * expr(* E = E*)
  | ELs of expr * expr(* E < E *)
  | EIf of expr * expr*expr(* if E then E else E *)

let add e1 e2 = match e1,e2 with
              | VInt e1 , VInt e2 -> VInt (e1 + e2)
              | _ , _ -> raise Eval_error

let sub e1 e2 = match e1,e2 with
              | VInt e1 , VInt e2 -> VInt (e1 - e2)
              | _ , _ -> raise Eval_error

let mul e1 e2 = match e1,e2 with
              | VInt e1 , VInt e2 -> VInt (e1 * e2)
              | _ , _ -> raise Eval_error

let div e1 e2 = match e1,e2 with
              | VInt e1 , VInt e2 -> VInt (e1 / e2)
              | _ , _ -> raise Eval_error

let equal e1 e2 = match e1,e2 with
              | VInt e1, VInt e2 -> VBool (e1 = e2)
              | VBool e1 , VBool e2 -> VBool (e1 = e2)
              | _ , _ -> raise Eval_error

let ls e1 e2 = match e1,e2 with
              | VInt e1 , VInt e2 -> VBool (e1 < e2)
              | _ , _ -> raise Eval_error

let myif e1 e2 e3 = match e1 with
              | VBool i -> if i then e2 else e3 
              | _ -> raise Eval_error

let rec eval e  = match e with
                |EConstInt i -> VInt i
                |EAdd (e1,e2) -> add (eval e1)  (eval e2)
                |ESub (e1,e2) -> sub (eval e1)  (eval e2)
                |EMul (e1,e2) -> mul (eval e1)  (eval e2)
                |EDiv (e1,e2) -> div (eval e1)  (eval e2)
                |EConstBool i -> VBool i
                |EEqual (e1,e2) -> equal (eval e1)  (eval e2)
                |ELs (e1,e2) -> ls (eval e1)  (eval e2)
                |EIf (e1,e2,e3) -> myif (eval e1)  (eval e2) (eval e3)
