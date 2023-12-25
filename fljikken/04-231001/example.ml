(* Definition of "the" list monad *)
type 'a m = 'a list

let (>>=) (x : 'a m) (f : 'a -> 'b m) : 'b m =
  List.concat (List.map f x)

let return (x : 'a) : 'a m = [x]

let (let*) = (>>=)

let guard (x : bool) : unit m =
  if x then return () else []



(** Examples of boolean functions *)
(** SAT *)
let phi x y z =
  not(x)
  && (x || not(y))
  && (x || y || z)

(** UNSAT *)
let psi x y z =
  not(x)
  && (x || not(y))
  && (x || y || z)
  && (x || y || not(z))

  let sat f =
    [true;false] >>= (fun x ->
    [true;false] >>= (fun y ->
    [true;false] >>= (fun z ->
    (guard (f x y z)) >>= (fun _ ->
    return (x, y, z)))))