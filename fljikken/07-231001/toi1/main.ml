open Syntax
open Eval
open Type
open TySyntax

let rec read_eval_print tyenv env  =
  print_string "# ";
  flush stdout;
  let cmd = Parser.toplevel Lexer.main (Lexing.from_channel stdin) in
  let (ty, newtyenv ) = infer_cmd tyenv  cmd in
  let (id, newenv, v) = eval_command env cmd in
  Printf.printf "%s : " id;
  print_type ty;
  Printf.printf(" = ");
  print_value v;
  print_newline ();
  read_eval_print newtyenv newenv

let initial_env =
  empty_env
  |> extend "i" (VInt 1)
  |> extend "v" (VInt 5)
  |> extend "x" (VInt 10)

  let initial_tyenv= empty_tyenv

let _ = read_eval_print initial_tyenv initial_env
