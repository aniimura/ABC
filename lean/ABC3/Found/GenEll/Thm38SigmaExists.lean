/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.QDomain
import Mathlib.FieldTheory.Galois.Infinite
import ABC3.Meta.Claim

/-!
# 第 1280 ブロック —— **`π` を動かす `σ` は存在する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——訂正後の道筋の段 5

第 1279 で `π`（`π^l = Q`）は体を拡げずに作れた。
★次は「`π` を動かす Galois 元がある」ことである。

☆機構は 2 段だけ:

| 段 | 内容 |
|---|---|
| 1 | `l ∤ v(Q)` なら `π` は**基礎体に無い**（あれば `v(Q) = l·v(y)`） |
| 2 | 代数閉包は基礎体上 Galois なので、基礎体に無い元は**どれかの `σ` が動かす** |

★★★段 2 は mathlib の `InfiniteGalois.mem_range_algebraMap_iff_fixed` そのものである。

☆このとき `σ(π)^l = σ(Q) = Q = π^l` なので `σ(π)/π` は `l` 乗根であり、
`σ(π) ≠ π` だから**原始 `l` 乗根**になる——これが `α` の非自明性の源である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

/-- ★★★★★★★★★★★★
**`l ∤ v(Q)` なら `π` は基礎体に無い**——★**無条件**（第 1280）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`π = ι y` なら `y^l = Q` なので `v(Q) = l·v(y)` になってしまう。 -/
theorem not_mem_range_of_not_dvd_vAdd {K₀ K : Type} [Field K₀] [Field K] [Algebra K₀ K]
    (v₀ : K₀ˣ →* Multiplicative ℤ) {l : ℕ} (Q₀ : K₀ˣ)
    (hnd : ¬ ((l : ℤ) ∣ vAdd v₀ Q₀))
    (π : Kˣ) (hπ : π ^ l = Units.map (algebraMap K₀ K : K₀ →* K) Q₀) :
    ((π : K)) ∉ Set.range (algebraMap K₀ K) := by
  rintro ⟨y₀, hy₀⟩
  have hy0ne : y₀ ≠ 0 := by
    intro h
    rw [h, map_zero] at hy₀
    exact (π.ne_zero) hy₀.symm
  set y : K₀ˣ := Units.mk0 y₀ hy0ne with hy
  have hmap : Units.map (algebraMap K₀ K : K₀ →* K) y = π := Units.ext hy₀
  have hinj : Function.Injective (Units.map (algebraMap K₀ K : K₀ →* K)) :=
    Units.map_injective (algebraMap K₀ K).injective
  have hyl : y ^ l = Q₀ := by
    refine hinj ?_
    rw [map_pow, hmap, hπ]
  refine hnd ⟨vAdd v₀ y, ?_⟩
  have : vAdd v₀ (y ^ l) = (l : ℤ) * vAdd v₀ y := by
    rw [← zpow_natCast y l, vAdd_zpow]
  rw [← hyl, this]

/-- ★★★★★★★★★★★★★★★★
**基礎体に無い元はどれかの `σ` が動かす**——★**無条件**（第 1280）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆mathlib の `InfiniteGalois.mem_range_algebraMap_iff_fixed` の対偶である。 -/
theorem exists_algEquiv_move {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]
    [IsGalois K₀ M] (x : M) (hx : x ∉ Set.range (algebraMap K₀ M)) :
    ∃ σ : M ≃ₐ[K₀] M, σ x ≠ x := by
  by_contra hcon
  push_neg at hcon
  exact hx ((InfiniteGalois.mem_range_algebraMap_iff_fixed x).2 hcon)

/-- ★★★★★★★★★★★★★★★★★★★★
**`l ∤ v(Q)` なら `π` を動かす `σ` がある**——★**無条件**（第 1280）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `α` の**非自明性**の源である。 -/
theorem exists_algEquiv_move_pi {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]
    [IsGalois K₀ M] (v₀ : K₀ˣ →* Multiplicative ℤ) {l : ℕ} (Q₀ : K₀ˣ)
    (hnd : ¬ ((l : ℤ) ∣ vAdd v₀ Q₀))
    (π : Mˣ) (hπ : π ^ l = Units.map (algebraMap K₀ M : K₀ →* M) Q₀) :
    ∃ σ : M ≃ₐ[K₀] M, σ (π : M) ≠ (π : M) :=
  exists_algEquiv_move _ (not_mem_range_of_not_dvd_vAdd v₀ Q₀ hnd π hπ)

/-! ## ★出典の紐付け(`.src`) -/

def not_mem_range_of_not_dvd_vAdd.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l ∤ v(Q) なら π は基礎体に無い。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_move.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体に無い元はどれかの σ が動かす。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_move_pi.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l ∤ v(Q) なら π を動かす σ がある。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_move_pi.needs : List ProofObligation :=
  [ .citation "[mathlib]" "InfiniteGalois.mem_range_algebraMap_iff_fixed"
      (.inMathlib "InfiniteGalois.mem_range_algebraMap_iff_fixed") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1280）**——訂正後の道筋（第 1277）の段 5 である。" ++
       "☆`σ(π)^l = σ(Q) = Q = π^l` なので `σ(π)/π` は `l` 乗根であり、" ++
       "`σ(π) ≠ π` だから**原始 `l` 乗根**になる。" ++
       "★★★これが `α` の非自明性の源である。") 2 ]

end ABC3.Found.GenEll
