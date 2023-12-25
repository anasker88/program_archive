type 'a m = 'a list

let (>>=) (x : 'a m) (f : 'a -> 'b m) : 'b m =
  List.concat (List.map f x)

let return (x : 'a) : 'a m = [x]

let (let*) = (>>=)

let guard (x : bool) : unit m =
  if x then return () else []



(** check if "banana + banana = sinamon" *)
let test_banana ba na si mo n =
  (100 * ba + 10 * na + na
   + 100 * ba + 10 * na + na
   = 1000 * si + 100 * na + 10 * mo + n)

(** check if "send + more = money" *)
let test_money s e n d m o r y =
  (1000 * s + 100 * e + 10 * n + d
   + 1000 * m + 100 * o + 10 * r + e
   = 10000 * m + 1000 * o + 100 * n + 10 * e + y)

let calc_banana f =
  let* ba = [0;1;2;3;4;5;6;7;8;9] in
  let* na = [0;1;2;3;4;5;6;7;8;9] in
  let* si = [0;1;2;3;4;5;6;7;8;9] in
  let* mo = [0;1;2;3;4;5;6;7;8;9] in
  let* n = [0;1;2;3;4;5;6;7;8;9] in
  let* _=guard (f ba na si mo n) in
  return (ba,na,si,mo,n)

let ans_banana=calc_banana test_banana



let calc_money f =
  let e_li=[0;1;2;3;4;5;6;7;8;9] in
  let* e = e_li in
  let n_li= List.filter (fun x -> x <> e) [0;1;2;3;4;5;6;7;8;9] in
  let* n = n_li in
  let d_li= List.filter (fun x -> x <> n) n_li in
  let* d = d_li in
  let o_li= List.filter (fun x -> x <> d) d_li in
  let* o = o_li in
  let r_li= List.filter (fun x -> x <> o) o_li in
  let* r = r_li in
  let y_li= List.filter (fun x -> x <> r) r_li in
  let* y = y_li in
  let m_li= List.filter (fun x -> x <> y && x <> 0) y_li in
  let* m = m_li in
  let s_li= List.filter (fun x -> x <> m) m_li in
  let* s = s_li in
  let* _=guard (f s e n d m o r y) in
  return (s,e,n,d,m,o,r,y)

let ans_money=calc_money test_money