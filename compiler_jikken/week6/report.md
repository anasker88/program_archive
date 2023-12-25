# コンパイラ実験第六回課題レポート　阿部 桃大

## 課題1

以下のocamlプログラムをrisc-vアセンブリに変換する。なお、risc-vの各命令については、[risc-vの仕様書](https://riscv.org/specifications/)を参照した。

```ocaml
let rec gcd m n =
if m <= 0 then n else
if m <= n then gcd m (n – m)
else gcd n (m – n) in
print_int (gcd 21600 337500)
```

これを、末尾最適化を適用してrisc-vアセンブリに変換すると、以下のようになる。

```riscv
main:
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    addi s0, sp, 16
    li a0, 21600
    li a1, 337500
    jal gcd
    jal print_int
    lw ra, 12(sp)
    lw s0, 8(sp)
    addi sp, sp, 16
    ret
gcd:
    bgtz a0, L1
    mv a0, a1
    ret
L1:
    bgt a0, a1, L3
    sub a1, a1, a0
    j gcd
L3:
    sub a0, a0, a1
    j gcd
```

`gcd`内で再帰的に`gcd`を呼ぶ際に`j`命令を使っているため、raレジスタの値を変化させることなく再帰的に呼び出すことができる。これにより、`gcd`内で再帰的に`gcd`を呼ぶ際に、raレジスタの値を保存する必要がなくなり、末尾最適化を行うことができる。

こうした最適化を行わない場合、以下のようなアセンブリが生成されると考えられる。

```riscv
main:
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    addi s0, sp, 16
    li a0, 21600
    li a1, 337500
    jal gcd
    jal print_int
    lw ra, 12(sp)
    lw s0, 8(sp)
    addi sp, sp, 16
    ret
gcd:
    bgtz a0, L1
    mv a0, a1
    j L2
L1:
    bgt a0, a1, L3
    sub a1, a1, a0
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    jal gcd
    lw ra, 12(sp)
    lw s0, 8(sp)
    addi sp, sp, 16
    j L2
L3:
    sub a0, a0, a1
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    jal gcd
    lw ra, 12(sp)
    lw s0, 8(sp)
    addi sp, sp, 16
L2:
    ret
```

ここでは、`gcd`内で再帰的に`gcd`を呼ぶ際に、スタックにraの値を積んでから`jal`命令で呼び出しを行い、呼び出し処理が終わるとスタックに積んだraの値を読みだして`ret`命令で戻るという処理を行っている。これにより、`gcd`内で再帰的に`gcd`を呼ぶ際に、raレジスタの値を保存する必要が生じており、処理数が増加しているだけでなく無駄にメモリを消費している。

## 課題2

```ocaml
let rec ack x y =
if x <= 0 then y + 1 else
if y <= 0 then ack (x-1) 1
else ack (x-1) (ack x (y-1))
in
print_int (ack 3 10)
```

このコードをCPS変換すると以下のようになる。

```ocaml
let rec ack' x y k =
if x <= 0 then k (y + 1) else
if y <= 0 then ack' (x-1) 1 k
else ack' x (y-1) (fun v -> ack' (x-1) v k)
in
ack' 3 10 (fun r -> print_int r)
```

2行目では元の`ack`関数では値`y+1`を返すので、継続`k`に`y+1`を渡している。3行目では元の`ack`関数では`ack (x-1) 1`という再帰呼び出しを行うので、`ack'`に継続`k`を渡して呼び出しを行っている。4行目では元の`ack`関数では`ack`を2回再帰呼び出しするが、先に内側の`ack`関数が計算されてその結果が外側の`ack`関数に渡されるので、外側の`ack`関数を継続`k`を持ち、引数から`ack'`を計算する継続として内側の`ack`関数に対応する`ack'`関数に渡して呼び出すことで、CPS変換を行っている。

これにより、すべての再帰呼び出しが末尾再帰呼び出しになっているため、末尾最適化を行うことができる。そのため、実行時のスタックの深さや命令数を削減することができると期待される。

## 課題4

まずは、λ抽象について説明する。

mylexer.mllに新たにトークンFUN,ARROWを追加した。

mylexer.mll(33~36行目)

```ocaml
| "fun"
    { FUN }
| "->"
    { ARROW }
```

これにより、`fun`と`->`をトークンとして認識するようになった。こうして認識したトークンは、myparser.mlyにおいて,`let rec`文と同様に解釈するようにした。

myparser.mly(193~198行目)

```ocaml
lambda:
| FUN formal_args ARROW exp
    %prec prec_fun
    { let f = Id.genid "f" in
        fun_args := (f, List.length $2) :: !fun_args;
        LetRec({ name = addtyp f; args = $2; body = $4 }, Var(f)) }
```

これにより、`fun`と`->`を用いたλ抽象を解釈することができるようになった。

次に、部分適用の実装について説明する。該当する変更点は以下の通りである。

myparser.mly(193~211行目)

```ocaml
lambda:
| FUN formal_args ARROW exp
    %prec prec_fun
    { let f = Id.genid "f" in
        fun_args := (f, List.length $2) :: !fun_args;
        LetRec({ name = addtyp f; args = $2; body = $4 }, Var(f)) }

let_partial:
| LET IDENT EQUAL exp
    %prec prec_let
    { try let x = is_fun $4 in
        fun_args := ($2, List.assoc x !fun_args) :: !fun_args;
        ($2, $4)
     with Notfunction -> ($2, $4) }

fundef_f:
| IDENT formal_args
    { fun_args := ($1, List.length $2) :: !fun_args;
        (addtyp $1, $2)}
```

関数定義の際に、その関数の引数の数を`ref`の`fun_args`に保存するようにした。そして、関数適用の際に、とっている引数の数が関数の引数の数よりも少ない場合は、部分適用を行い、`let rec`文として関数を定義するようにした。これにより、部分適用を行うことができるようになった。

ただし、今回の実装ではクロージャ―変換と部分適用が同時に起こる場合に対応していない。また、最初の関数定義の際にすべての引数が明示されていないと部分適用を行うことができない。例えば以下のような場合である。

```ocaml
let rec make_adder x =
  let rec adder y = x + y in
  adder in
print_int (make_adder 3 7)
```
