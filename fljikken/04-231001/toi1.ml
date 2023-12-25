type 'a m = ('a, string) result

let (>>=) (x : 'a m) (f : 'a -> 'b m) : 'b m =
  match x with
  | Ok(v) -> f v
  | Error msg -> Error msg

let return (x : 'a) : 'a m = Ok x

let (let*) = (>>=)

let err msg = Error msg

let myDiv x y : int m =
  if y = 0 then
    Error "Division by Zero"
  else
    return (x / y)

let rec eLookup (key : 'a)  (t : ('a * 'b) list) : 'b m = match t with
    |[] -> Error "Not Found" (*最後まで見つからなければNotFound*)
    |(v1 , v2) :: rest -> if v1=key then return v2 (*keyが見つかったら相方を返す*) 
                          else eLookup key rest 


let lookupDiv (kx : 'a) (ky : 'a) (t : ('a * int) list) : int m = 
  let* x = eLookup kx t in
  let* y = eLookup ky t in
  myDiv x y

(** Tests *)
let table = [("x", 6); ("y", 0); ("z", 2)]

(* same as is_error *)
let isErr x =
  match x with
  | Ok(_) -> false
  | Error(_) -> true

let () =
  let b1 = isErr (lookupDiv "x" "y" table) in
  let b2 = (Ok 3 = lookupDiv "x" "z" table) in
  let b3 = isErr (lookupDiv "x" "a" table) in
  let b4 = isErr (lookupDiv "a" "z" table) in
  if b1 && b2 && b3 && b4 then
    print_string "ok\n"
  else
    print_string "wrong\n"
