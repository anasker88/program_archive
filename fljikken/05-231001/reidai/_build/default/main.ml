open Syntax
open Eval

let rec read_eval_print env =
  print_string "# ";
  flush stdout;
  let cmd = Parser.toplevel Lexer.main (Lexing.from_channel stdin) in
  let (id, newenv, v) = eval_command env cmd in
  (Printf.printf "%s = " id;
   print_value v;
   print_newline ();
   read_eval_print newenv)

let initial_env =
  empty_env
  |> extend "i" (VInt 1)
  |> extend "ii" (VInt 2)
  |> extend "iii" (VInt 3)
  |> extend "iv" (VInt 4)
  |> extend "v" (VInt 5)
  |> extend "x" (VInt 10)

let _ = read_eval_print initial_env
