open Effect
open Effect.Deep

module type MEMO_AND_WRITE = sig
  val memo : (int -> int) -> int -> int
  val write : string -> unit
  val run : (int -> int) -> int -> int * string
end

module M : MEMO_AND_WRITE = struct

  type _ Effect.t += Memo : ((int -> int) *  int) -> int Effect.t
                   | Write : string -> unit Effect.t


  let memo f n = perform (Memo (f, n))

  let write s = perform (Write s)


  let rec with_memo_handler f n =
    match_with f n
        { retc = (fun r -> failwith "not implemented");
          exnc = (fun e -> raise e);
          effc = (fun (type b) (eff: b Effect.t) ->
            match eff with
             Memo (f, n) -> Some (fun (k: (b,_) continuation) ->
                  failwith "not implemented")
            | _ -> None);
        }

  let run_memo f n : int = failwith "not implemented"

  let with_write_handler f n : int * string =
    match_with f n
        { retc = (fun r -> failwith "not implemented");
          exnc = (fun e -> raise e);
          effc = (fun (type b) (eff: b Effect.t) ->
            match eff with
            | Write s -> Some (fun (k: (b,_) continuation) ->
                  failwith "not implemented")
            | _ -> None);
        }

 let run f n = with_write_handler (run_memo f) n
end


let rec fib n =
  M.write ("Fib(" ^ (string_of_int n) ^")\n");
  if n <= 1 then
    n
  else
    let x = M.memo fib (n - 2) in
    let y = M.memo fib (n - 1) in
    x + y

let _ =
  let (x, s) = M.run fib 10 in
  assert (x = 55);
  print_string s
