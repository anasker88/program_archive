
(* The type of tokens. *)

type token = 
  | THEN
  | SEMISEMI
  | RPAR
  | PLUS
  | LT
  | LPAR
  | INT of (int)
  | IF
  | ID of (string)
  | EQ
  | ELSE
  | BOOL of (bool)

(* This exception is raised by the monolithic API functions. *)

exception Error

(* The monolithic API. *)

val toplevel: (Lexing.lexbuf -> token) -> Lexing.lexbuf -> (Syntax.command)
