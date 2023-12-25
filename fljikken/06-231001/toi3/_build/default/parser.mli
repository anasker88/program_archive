
(* The type of tokens. *)

type token = 
  | TIMES
  | THEN
  | SEMISEMI
  | RPAR
  | REC
  | PLUS
  | MINUS
  | LT
  | LPAR
  | LET
  | INT of (int)
  | IN
  | IF
  | ID of (string)
  | FUN
  | EQ
  | ELSE
  | DIV
  | BOOL of (bool)
  | ARROW
  | AND

(* This exception is raised by the monolithic API functions. *)

exception Error

(* The monolithic API. *)

val toplevel: (Lexing.lexbuf -> token) -> Lexing.lexbuf -> (Syntax.command)
