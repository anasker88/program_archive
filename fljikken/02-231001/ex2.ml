type iexpr =
        | EConstInt of int
        | EAdd of iexpr * iexpr
        | ESub of iexpr * iexpr
        | EMul of iexpr * iexpr
let rec  eval x = match x with
        |EConstInt y-> y
        |EAdd(y,z)->eval y+eval z
        |ESub(y,z)->eval y-eval z
        |EMul(y,z)->eval y*eval z
