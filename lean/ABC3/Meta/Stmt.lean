/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Lean

/-!
# `stmt%` —— 定数の**型そのもの**を項として取る

## なぜ要るか

★「証明を書く前に依存を型で出す」ためである。

Skeleton の主張どうしは、statement のうえでほぼ繋がっていない
(2026-09-09 実測: 247 件中 **2 件**しか他の主張に触れない。`tools/stmt-edges.py`)。
★これは原典に忠実な結果である——原典でも「Lemma 3.1」は Proposition 3.4 の
**証明**に出て**主張**には出ない。

辺を証明の前に出すには、還元を**含意として**書けばよい:

```lean
theorem Y_of_X (h : <X の型>) : <Y の型> := sorry
theorem chain_Y : <Y の型> := Y_of_X X   -- ★ここが型検査され、定数 X を消費する
```

★問題は `<X の型>` を**書き写さねばならない**ことである。長い statement では
数十行になり、写し間違いも起きる(2026-09-09 に本体が 2 回落ちた)。

`stmt% X` は定数 `X` の型を**そのまま**項として返すので、書き写しが要らない。

```lean
theorem Y_of_X (h : stmt% X) : stmt% Y := sorry
theorem chain_Y : stmt% Y := Y_of_X @X
```

## ★これが作る辺と、作らない辺

- `stmt% X` 自体は `X` への依存を**作らない**(型を複製するだけ)。
- 辺を作るのは `:= Y_of_X @X` の側である。★意図した依存だけが辺になる。
-/

namespace ABC3.Meta

open Lean Elab Term Meta

/-- `stmt% X` —— 定数 `X` の型を項として返す。`X` が命題なら結果は `Prop`。 -/
elab "stmt% " n:ident : term => do
  let c ← realizeGlobalConstNoOverload n
  let info ← getConstInfo c
  let us ← info.levelParams.mapM fun _ => mkFreshLevelMVar
  return info.instantiateTypeLevelParams us

end ABC3.Meta
