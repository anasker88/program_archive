let r = ref 0
let change x =let y = (!r) in (r:= x);y
let f x = change x
