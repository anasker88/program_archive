type name = string
type pattern =
  | PInt of int 
  | PBool of bool
  | PVar of name
  | PPair of pattern * pattern
  | PNil 
  | PCons of pattern * pattern

type value =
  | VInt of int
  | VBool of bool
  | VFun of name * expr * env
  | VRecFun of name * name * expr * env
  | VPair of value * value
  | VNil
  | VCons of value * value
and env = (name * value) list
and expr =
  | EConstInt  of int
  | EConstBool of bool
  | EVar       of name
  | EAdd       of expr * expr
  | ESub       of expr * expr
  | EMul       of expr * expr
  | EDiv       of expr * expr
  | EEq        of expr * expr
  | ELt        of expr * expr
  | EIf        of expr * expr * expr
  | ELet       of name * expr * expr
  | EFun       of name * expr
  | EApp       of expr * expr 
  | ELetRec    of name * name * expr * expr

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
  | VNil -> Some []
  | _ -> None)
|PCons (c1,c2)-> (match va with
  |VCons  (v1,v2) -> (match (find_match c1 v1),(find_match c2 v2)  with
    | None,_-> None
    | _,None -> None
    | Some a ,Some b -> Some (a@b))
  | _ -> None)

