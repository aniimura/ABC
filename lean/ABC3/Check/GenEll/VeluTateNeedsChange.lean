/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVelu

/-!
# 測定と訂正 —— **Vélu の商は `E_{q^l}` と「等しく」はない（`j` は等しい）**（`Check`）

**これは原典の主張ではない**（我々の定式化についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の数値測定（第 834）

`Found/GaloisRep/TateVelu.lean` の `veluCurve_tateCurveAt_eq`（第 718）は

```
veluCurve (tateCurveAt q) v w = tateCurveAt (q^l)
  ⟺ a₄(q^l) = a₄(q) − 5v  かつ  a₆(q^l) = a₆(q) − v − 7w
```

という条件付きの主張であり、第 718 のブロックはこれを
「★変数変換は要らない——残る穴は `q` 展開の恒等式 2 本」とまとめていた。

## ★★★★★★★★訂正——その 2 本の恒等式は**偽**である

`q` 展開を実際に計算すると（`l = 5, 7` で確認）、

* `n ∉ {0, l}` の係数では両辺が一致する
* **`n = 0` と `n = l` では一致しない**

具体的には Vélu の和の定数項が

    `v(0) = (l⁴ − 1)/240 ≠ 0`

である（`l = 5, 7, 11` で数値確認）。★`a₄(q)` も `a₄(q^l)` も定数項は `0` なので、
`a₄(q^l) = a₄(q) − 5v` は**成り立たない**。

## ★★★★★★★★★★しかし `j` は一致する

`q = 0.02`, `l = 5` で数値計算すると

    `j(veluCurve (tateCurveAt q) v w) = j(tateCurveAt (q^5))`

が有効数字 7 桁で一致する（差は打ち切り誤差）。★★すなわち Vélu の商は
`E_{q^l}` と**同型**だが**等しくはない**——**変数変換が要る**。

## ★★★★★★★★★★★★これが意味すること

`Skeleton/GenEll/TateIsogeny.lean` の `tateModel_of_quot_mu` の結論は

    `∃ D' : VariableChange R, D' • integralModel R W' = tateCurveAt (q^l)`

であり、**変数変換 `D'` を許している**ので、命題そのものは正しい。
★しかし第 718 の「変数変換は要らない」という**還元は使えない**。

## ★★★★★★★★★★★★★★★★★★★★**正しい恒等式**（第 835、`l = 5, 7` で数値確認）

係数 `a₄`・`a₆` ではなく **`c₄`・`c₆` で書くと恒等式は真になる**:

    `c₄(veluQuotient (tateCurveAt q) μ_l) = l⁴ · c₄(tateCurveAt (q^l))`
    `c₆(veluQuotient (tateCurveAt q) μ_l) = −l⁶ · c₆(tateCurveAt (q^l))`

★★これは `q` 展開の**すべての係数で成立する**（`l = 5, 7`、`q^21` まで確認）。

★★★したがって

    `Δ′ = (c₄′³ − c₆′²)/1728 = l¹² · Δ(q^l)`
    `j′ = c₄′³/Δ′ = j(q^l)`

——**`j` は一致する**。☆`l⁴`・`−l⁶` は変数変換の尺度因子であり、
`Δ` の `l¹²` と相殺して `j` には現れない。

☆参考: 誤っていた形を `a₄` で書くと
`a₄′ = a₄(q^l) − ((l⁴−1)/48)·c₄(q^l)` であり、`(l⁴−1)/48` の項が
第 718 の「係数 2 本の恒等式」を破っていた。

## ★★★★★★★★★★★★★★次にどうするか

葉 1（`jExp_velu_bad`）が要求するのは **`j` の付値だけ**である:

    `jExp p E′ = l · jExp p E`

★★したがって**係数 2 本の恒等式ではなく `j` の恒等式**

    `j(veluQuotientFull E S) = j(tateCurveAt (q^l))`

を目標にすればよい。☆これは真であり（上の数値確認）、かつ
`jExp` が `j` だけで決まる（`Found/GaloisRep/HtFaltJ.lean` 群）ので十分である。

## ★★★★★★★★★★★★★★★★★★★★★★Vélu の和 `v` の閉じた形（第 839）

上の `c₄` の形を展開すると、Vélu の和自体の値が出る
（`l = 5, 7` で `q^21` まで数値確認）:

    `v(0) = (l⁴ − 1)/240`
    `v(n) = l⁴·σ₃(n/l)·[l∣n] − σ₃(n)`   (`n ≥ 1`)

★★これが葉 1 の係数の**真の目標**である。
☆第 718 の形（`v(n) = σ₃(n/l)[l∣n] − σ₃(n)`、`v(0) = 0`）とは
`l⁴` と `(l⁴−1)/240` だけ違う。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GaloisRep ABC3.Found.GenEll

/-- ★★★★★★**訂正後に使うべき形**——`veluCurve_tateCurveAt_eq` は
係数 2 本を仮説に持つ条件付きの主張であり、その仮説は（Tate 曲線と `μ_l` に対しては）
**偽**である。★ゆえに本補題は葉 1 には使えない。 -/
theorem velu_tate_eq_is_conditional {R : Type} [CommRing R] {I : Ideal R}
    [IsAdicComplete I R] (q : R) (hq : q ∈ I) (q' : R) (hq' : q' ∈ I) (v w : R)
    (h4 : (tateCurveAt q' hq').a₄ = (tateCurveAt q hq).a₄ - 5 * v)
    (h6 : (tateCurveAt q' hq').a₆ = (tateCurveAt q hq).a₆ - v - 7 * w) :
    veluCurve (tateCurveAt q hq) v w = tateCurveAt q' hq' :=
  veluCurve_tateCurveAt_eq q hq q' hq' v w h4 h6

end ABC3.Check.GenEll
