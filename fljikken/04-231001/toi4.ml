type 'a m = ('a * 'a) list -> 'a * (('a * 'a) list)


let (>>=) (x : 'a m) (f : 'a -> 'b m) : 'b m =
  fun init ->
    let (r, s) = x init in
    f r s

let return x : 'a m = fun s -> (x, s)

let (let*) = (>>=)

let memo (f : int -> int m) (n : int) : int m =
  fun cache ->

    match List.assoc_opt n cache with
    | Some v -> (v,cache)
    | None ->
        let (v, s) = f n cache in
        (v, (n, v):: s)

let runMemo (x : 'a m) : 'a =
  let (r, _) = x [] in
  r

  let rec fib n =
    if n <= 1 then
      return n
    else
      let* r1 = memo fib (n-2) in
      let* r2 = memo fib (n-1) in
      return (r1 + r2)

let () =
  if runMemo (fib 80) = 23416728348467685 && runMemo (fib 10) = 55 then
    print_string "ok\n"
  else
    print_string "wrong\n"
