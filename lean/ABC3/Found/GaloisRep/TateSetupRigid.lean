/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TatePhi
import ABC3.Meta.Claim

/-!
# 第 1275 ブロック —— **`TateSetup` の底環 `R` は `K` を決めてしまう**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★☆★★☆測定——`σA : K →ₐ[R] K` は**恒等射しかない**

`TateSetup` の `hmem0` は「`v ≥ 0` の元は `R` から来る」と言う。
★どの `x : Kˣ` についても `v(x) ≥ 0` か `v(x⁻¹) ≥ 0` のどちらかなので、
**`K` は `R` の像の分数体**である。

☆したがって `R`-代数準同型 `K →ₐ[R] K` は `R` の像の上で恒等、
逆元も保つので **`K` 全体で恒等**である。

★★★これは `tatePhi_pointMap`（在庫）の `σA : K →ₐ[R] K` が
**恒等射のときしか内容を持たない**ことを意味する
——同変性は形の上では書けているが、**非自明な Galois 元には当たらない**。

☆同じ理由で `tate_stab_of_pointStab`・`lemma_3_2_i_tate` の
`σ : L ≃ₐ[K] L` も、`IsScalarTower R K L` と `TateSetup R I L` から `K = L` が出るので
**恒等射しか動かない**。

★★正しい形は `tatePhi_map`（在庫、`σR : R →+* R` と `σK : K →+* K` の対を受ける）である
——そちらは `R` を動かす自己同型を許す。☆本ブロックは**測定だけ**を記録する。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★
**`TateSetup` があれば `K` は `R` の像の分数体**——★**無条件**（第 1275）。 -/
theorem tateSetup_exists_of_ne_zero (S : TateSetup R I K) (x : K) (hx : x ≠ 0) :
    (∃ y : R, algebraMap R K y = x) ∨ (∃ y : R, algebraMap R K y = x⁻¹) := by
  rcases le_total 0 (vAdd S.v (Units.mk0 x hx)) with h | h
  · exact Or.inl (S.hmem0 (Units.mk0 x hx) h)
  · refine Or.inr ?_
    have h' : 0 ≤ vAdd S.v (Units.mk0 x hx)⁻¹ := by
      rw [vAdd_inv]
      omega
    obtain ⟨y, hy⟩ := S.hmem0 (Units.mk0 x hx)⁻¹ h'
    exact ⟨y, by simpa using hy⟩

/-- ★★★★★★★★★★★★★★★★
**`R`-代数準同型 `K →ₐ[R] K` は恒等射しかない**——★**無条件**（第 1275）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★これが `tatePhi_pointMap`（在庫）の同変性が
**非自明な Galois 元に当たらない**ことの証拠である。 -/
theorem tateSetup_algHom_eq_id (S : TateSetup R I K) (σA : K →ₐ[R] K) (x : K) :
    σA x = x := by
  by_cases hx : x = 0
  · simp [hx]
  rcases tateSetup_exists_of_ne_zero S x hx with ⟨y, hy⟩ | ⟨y, hy⟩
  · rw [← hy]
    exact σA.commutes y
  · have hxi : x⁻¹ ≠ 0 := inv_ne_zero hx
    have h1 : σA x⁻¹ = x⁻¹ := by
      rw [← hy]
      exact σA.commutes y
    have h2 : σA x * σA x⁻¹ = 1 := by
      rw [← map_mul, mul_inv_cancel₀ hx, map_one]
    rw [h1] at h2
    field_simp at h2
    exact h2

/-! ## ★出典の紐付け(`.src`) -/

def tateSetup_exists_of_ne_zero.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(TateSetup があれば K は R の像の分数体。★無条件)",
    sectionId := "genell-def-3-3" }

def tateSetup_algHom_eq_id.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(R-代数準同型 K →ₐ[R] K は恒等射しかない。★無条件)",
    sectionId := "genell-def-3-3" }

def tateSetup_algHom_eq_id.needs : List ProofObligation :=
  [ .implicitStep
      ("★★☆**2026-09-02（第 1275）**——**測定**である。" ++
       "☆`tatePhi_pointMap`（在庫）の `σA : K →ₐ[R] K` は本定理により恒等射しかないので、" ++
       "同変性は**非自明な Galois 元には当たらない**。" ++
       "★正しい形は `tatePhi_map`（在庫、`σR : R →+* R` と `σK : K →+* K` の対）であり、" ++
       "点の側の同変性もその対から書き直す必要がある。") 3 ]

end ABC3.Found.GaloisRep
