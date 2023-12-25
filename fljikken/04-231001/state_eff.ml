open Effect
open Effect.Deep

module type STATE = sig
  type t
  val get : unit -> t
  val put : t -> unit
  val run1 : (unit -> unit) -> t -> unit
  val run2 : (unit -> unit) -> t -> unit
end

module State (S : sig type t end) : STATE with type t = S.t = struct

  type t = S.t

  type _ Effect.t += Get : t Effect.t
                   | Put : t -> unit Effect.t

  let get () = perform Get

  let put v  = perform (Put v)

  let run1 f init =
      let g = match_with f ()
        { retc = (fun r -> fun s -> (r, s));
          exnc = (fun e -> raise e);
          effc = (fun (type b) (eff: b Effect.t) ->
            match eff with
            | Get -> Some (fun (k: (b,_) continuation) ->
                    fun (s : t) -> (continue k s) s)
            | Put v -> Some (fun (k: (b,_) continuation) ->
                    fun (s : t) -> (continue k ()) v)
            | _ -> None);
        }
      in ignore (g init)

  let run2 f init =
      let buff : t ref = ref init in
      match_with f ()
        { retc = (fun r -> r);
          exnc = (fun e -> raise e);
          effc = (fun (type b) (eff: b Effect.t) ->
            match eff with
            | Get -> Some (fun (k: (b,_) continuation) ->
                    continue k !buff)
            | Put v -> Some (fun (k: (b,_) continuation) ->
                    buff := v; continue k ())
            | _ -> None);
        }
end

module IS = State (struct type t = int end)
module SS = State (struct type t = string end)

let example () : unit =
  let open Printf in
  printf "%d\n" (IS.get ());
  IS.put 10;
  printf "%d\n" (IS.get ());
  IS.put 20;
  printf "%d\n" (IS.get ());
  printf "%s\n" (SS.get ());
  SS.put "hello world";
  printf "%s\n" (SS.get ())

let _ = IS.run1 (fun () -> SS.run2 example "forty two") 42
