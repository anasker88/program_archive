exception PopError(*空スタックに対するpopはエラー*)

module type ABSTSTACK = sig
  type 'a t
  val pop: 'a t ->('a * 'a t)
  val push: 'a -> 'a t -> 'a t
  val empty: 'a t
  val size: 'a t -> int  
end

module AbstStack :ABSTSTACK = struct
  type 'a t ='a list
  let empty = [](*空スタック*)
  let pop xs= match xs with
    | [] -> raise PopError(*空スタックならエラー*)
    | x :: ys -> (x, ys) (*先頭を取り出し*)
  let push x lst = x :: lst(*先頭に挿入*)
  let rec size xs = match xs with
  | [] -> 0
  | _ :: ys -> 1 + size ys(*一要素ずつ取り出して長さ+1*)
end

