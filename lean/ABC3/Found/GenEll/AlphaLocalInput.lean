/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.TorsionNotRational
import ABC3.Found.GenEll.UnipFromRational
import ABC3.Found.GenEll.NeFromCoord
import ABC3.Meta.Claim

/-!
# 第 1307 ブロック —— **局所データから `h2`・`h1` を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——II 側の**局所側の完成形**

☆受け取るのは基礎体 `K₀` についての 3 つだけ:

| # | 入力 | Tate 曲線での出どころ |
|---|---|---|
| 1 | `ζ₀ ∈ K₀` が原始 `l` 乗根 | 大域で `ζ_l` を添加してよい |
| 2 | `P₀ ∈ E(K₀)` が位数 `l` | `μ_l ⊂ E(K₀)`（第 1297） |
| 3 | `K₀` の `l`-捩れはたかだか `l` 個 | `l ∤ v(Q)`（第 1304） |

★★★これだけで、**ある `σ` について `h2` と `h1` が同時に出る**
——`alpha_mem_map_of_galTate`（第 1237）の入力である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

/-- ★★★★★★★★★★★★★★★★
**基礎体の `l`-捩れが少なければ、動かされる `l`-捩れ点と `σ` がある**——★**無条件**（第 1307）。 -/
theorem exists_sigma_moved_of_card_le [IsAlgClosed M] [CharZero K₀] [IsGalois K₀ M]
    (W : WeierstrassCurve K₀) [W.IsElliptic]
    [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) (hl : l.Prime)
    (hcard : ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∃ (σ : M ≃ₐ[K₀] M) (Q : (W.baseChange M).toAffine.Point),
      l • Q = 0 ∧ galPoint W σ Q ≠ Q := by
  obtain ⟨Q, hQ, hnr⟩ := exists_torsion_not_rational (M := M) W l hl hcard
  obtain ⟨σ, hσ⟩ := exists_galPoint_ne_of_not_rational W Q hnr
  exact ⟨σ, Q, hQ, hσ⟩

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**局所データから `h2`・`h1` を同時に出す**——★**無条件**（第 1307）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`h2` は第 1300 が**どの `σ` についても**与えるので、
`h1` の `σ`（第 1307 前半）をそのまま使えばよい。 -/
theorem exists_sigma_h2_h1_of_local [IsAlgClosed M] [CharZero K₀] [IsGalois K₀ M]
    (W : WeierstrassCurve K₀) [W.IsElliptic]
    [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange M) l ≃+ (Fin 2 → ℤ_[l]))
    (ζ₀ : K₀) (hζ : IsPrimitiveRoot (algebraMap K₀ M ζ₀) l)
    (P₀ : W.toAffine.Point) (hP₀ : addOrderOf P₀ = l)
    (hcard : ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∃ σ : M ≃ₐ[K₀] M,
      (∀ x : tateModule (W.baseChange M) l, ∃ u : tateModule (W.baseChange M) l,
          galTate W l σ (galTate W l σ x) + x
            = galTate W l σ x + galTate W l σ x + l • u) ∧
        (∃ x : tateModule (W.baseChange M) l, ∀ u : tateModule (W.baseChange M) l,
          galTate W l σ x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  haveI : CharZero M := charZero_of_injective_algebraMap (algebraMap K₀ M).injective
  obtain ⟨σ, Q, hQ, hmoved⟩ := exists_sigma_moved_of_card_le (M := M) W l hl hcard
  refine ⟨σ, ?_, ?_⟩
  · exact galTate_unipotent_of_rational W l e σ ζ₀ hζ P₀ hP₀
  · exact exists_galTate_ne_of_galPoint W l σ Q hQ hmoved

/-! ## ★出典の紐付け(`.src`) -/

def exists_sigma_moved_of_card_le.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体の l-捩れが少なければ動かされる点と σ がある。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_h2_h1_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所データから h2・h1 を同時に出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_h2_h1_of_local.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galTate_unipotent_of_rational(第 1300、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galTate_unipotent_of_rational") 1,
    .citation "[ABC3]" "exists_torsion_not_rational(第 1306、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_torsion_not_rational") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1307）**——II 側の**局所側の完成形**である。" ++
       "☆受け取るのは `ζ₀ ∈ K₀`（原始 `l` 乗根）・`P₀ ∈ E(K₀)`（位数 `l`）・" ++
       "「`K₀` の `l`-捩れはたかだか `l` 個」の 3 つだけ。" ++
       "★★★どれも `K₀` の上の Tate 一意化から出る（第 1297・1304）。") 3 ]

end ABC3.Found.GenEll
