/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.RingTheory.PowerSeries.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Weierstrass
import ABC3.Found.GaloisRep.TateSeries.Definition33

/-!
# TateSeries —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep



/-! ## ★★★★★★★★★★★★`c₄ = E₄`・`c₆ = −E₆` -/

/-- ★★★★★★★★★★★★**Tate 曲線の `c₄` は `E₄`**——★**無条件**（第 1254）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`a₁ = 1`・`a₂ = a₃ = 0`・`a₄ = −5 σ₃` なので
`c₄ = b₂² − 24 b₄ = 1 − 48·(−5 σ₃) = 1 + 240 σ₃`。 -/
theorem c₄_tateCurve :
    tateCurve.c₄ = 1 + 240 * sigmaSeries 3 := by
  simp only [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    tateCurve, tateA4, map_neg, map_ofNat]
  ring

/-- ★★★★★★★★★★★★**Tate 曲線の `c₆` は `−E₆`**——★**無条件**（第 1254）。

☆`c₆ = −b₂³ + 36 b₂ b₄ − 216 b₆ = −1 + 72 a₄ − 864 a₆`
であり、`12 a₆ = −(5σ₃ + 7σ₅)` から `−1 + 504 σ₅` になる。 -/
theorem c₆_tateCurve :
    tateCurve.c₆ = -1 + 504 * sigmaSeries 5 := by
  have h12 : (12 : PowerSeries ℤ) * tateA6
      = -(5 * sigmaSeries 3 + 7 * sigmaSeries 5) := by
    have h := twelve_mul_tateA6
    simpa [map_ofNat] using h
  have hc : tateCurve.c₆ = -1 + 72 * tateA4 - 72 * ((12 : PowerSeries ℤ) * tateA6) := by
    simp only [WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
      WeierstrassCurve.b₆, tateCurve]
    ring
  rw [hc, h12, tateA4]
  simp only [map_neg, map_ofNat]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def c₄_tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の c₄ は E₄ = 1 + 240σ₃。★無条件)",
    sectionId := "genell-lemma-3-2" }

def c₆_tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の c₆ は −E₆ = −1 + 504σ₅。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
