open Syntax
open Eval

let rec read_eval_print env =
  print_string "# ";
  flush stdout;
  let cmd = Parser.toplevel Lexer.main (Lexing.from_channel stdin) in
  let ret = eval_command env cmd in
  let rec loop  x=
  match x with
  | ([],newenv) -> read_eval_print newenv
  | ((id,v)::rest,newenv) -> (Printf.printf "%s = " id;
  print_value v;
  print_newline ();
  loop (rest,newenv))
in loop ret
  

let initial_env =
  empty_env
  |> extend "i" (VInt 1)
  |> extend "v" (VInt 5)
  |> extend "x" (VInt 10)

let _ = read_eval_print initial_env
