/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalRepClosed

/-!
# 界面の測定と訂正 —— **`α` の仮説は `∃ e` でよい（`∀ e` では循環する）**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の測定と訂正（第 824-826）

`Found/GenEll/EllModuliGalois.lean` の `imageContainsSL2J_of_alpha` は、以前

```
halpha : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
  α ∈ (galRep E.W l e).range.map (glRedPadic l)
```

という仮説を持っていた。★★**これは循環する。**

## ★★★なぜか

`galRep E.W l e` は基底 `e` に依存し、基底を取り替えると**共役**になる
（`Found/GaloisRep/GalRepBasis.lean` の `galRep_basisChange`、第 825）。
★したがって `∀ e` は「`α` が像の**すべての共役**に入る」ことを意味し、
それは `SL₂ ⊆ 像` を既に知っていないと言えない。

一方、原文の局所理論が与えるのは **Tate 一意化に適合した基底
`(ζ_l, q^{1/l})` について**の主張である。★★ゆえに正しい形は `∃ e₀` である。

## ★★★★訂正（第 826）

仮説を `∃ e₀` に弱め、`sl2_range_basis_transfer`（第 825）で
すべての基底へ移すようにした。★★**結論 `ImageContainsSL2J E l` は変わらない**
——`SL₂` は `GL₂` の**正規部分群**なので共役で保たれるからである。

## ☆同じ形の測定

* `Check/GenEll/EllModuliDegInfPos.lean`（第 745）——界面は `deg∞ > 0` を強制する
* `Check/GenEll/LcyclicExcTooStrong.lean`（第 754-755）——`mem_lcyclicExc` は `l` の下界を落としていた
* `Check/GenEll/ImageSL2NeedsL5.lean`（第 776）——`5 ≤ l` を落としていた
* 本ファイル（第 826）——`α` の仮説が `∀ e` で循環していた

★どれも **witness を実際に作ろうとして初めて見えた**ものである。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Interface.GaloisRep

open scoped Classical

/-- ★★★★★★**訂正後に界面が主張していること**——仮説は `∃ e₀` である。 -/
theorem alpha_needs_one_basis (E : SSCurve) (l : ℕ) [Fact l.Prime] (hl5 : 5 ≤ l)
    (halpha : ∃ e₀ : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (Matrix.SpecialLinearGroup.toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e₀).range).map (glRedPadic l))
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l :=
  imageContainsSL2J_of_alpha' E l hl5 halpha hno

end ABC3.Check.GenEll
