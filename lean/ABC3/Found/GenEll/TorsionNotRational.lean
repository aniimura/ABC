/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.NotRationalMoved
import ABC3.Found.GaloisRep.TorsionCard
import ABC3.Found.GaloisRep.TorsionAll
import ABC3.Meta.Claim

/-!
# 第 1306 ブロック —— **基礎体の `l`-捩れが少なければ、有理でない `l`-捩れ点がある**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`h1` の入力の**個数の勘定**

代数閉体の上では `E[l]` はちょうど `l²` 個（`torsion_card`、在庫）。
★基礎体の上の `l`-捩れがたかだか `l` 個（第 1304）なら、`l ≥ 2` より `l < l²` なので、
**有理でない `l`-捩れ点が必ずある**。

☆それを第 1305 に渡せば「動かす `σ`」が取れ、第 1301 で `h1` になる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**基礎体の `l`-捩れが `l` 個以下なら、有理でない `l`-捩れ点がある**——★**無条件**（第 1306）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆代数閉体の上では `l²` 個あるので、個数が合わない。 -/
theorem exists_torsion_not_rational [IsAlgClosed M] [CharZero K₀]
    (W : WeierstrassCurve K₀) [W.IsElliptic]
    [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) (hl : l.Prime)
    (hcard : ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∃ P : (W.baseChange M).toAffine.Point, l • P = 0 ∧
      ¬ ∃ P₀ : W.toAffine.Point, rhPoint (algebraMap K₀ M) W P₀ = P := by
  classical
  haveI : CharZero M := charZero_of_injective_algebraMap (algebraMap K₀ M).injective
  haveI hEM : (W.map (algebraMap K₀ M)).IsElliptic :=
    inferInstanceAs ((W.baseChange M).IsElliptic)
  by_contra hcon
  push_neg at hcon
  -- すべての `l`-捩れ点は基礎体から来る
  have hchar : ∀ k : ℕ, 1 ≤ k → k ≤ l → ((k : M) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  have hfin : {P : (W.baseChange M).toAffine.Point | l • P = 0}.Finite :=
    finite_torsion (W.baseChange M) l hl.pos (Nat.cast_ne_zero.2 hl.ne_zero)
  have hcard2 : Nat.card {P : (W.baseChange M).toAffine.Point // l • P = 0} = l ^ 2 :=
    torsion_card (W.baseChange M) (W.baseChange M).isUnit_Δ.ne_zero l hl.pos hchar
  -- 選択関数
  set g : (W.baseChange M).toAffine.Point → W.toAffine.Point := fun P =>
    if h : ∃ P₀ : W.toAffine.Point, rhPoint (algebraMap K₀ M) W P₀ = P then h.choose else 0
    with hg
  have hgspec : ∀ P : (W.baseChange M).toAffine.Point, l • P = 0 →
      rhPoint (algebraMap K₀ M) W (g P) = P := by
    intro P hP
    have hex : ∃ P₀ : W.toAffine.Point, rhPoint (algebraMap K₀ M) W P₀ = P := hcon P hP
    rw [hg]
    simp only [dif_pos hex]
    exact hex.choose_spec
  set S := hfin.toFinset with hS
  have hmemS : ∀ P : (W.baseChange M).toAffine.Point, P ∈ S ↔ l • P = 0 := by
    intro P
    rw [hS, Set.Finite.mem_toFinset]
    rfl
  have hinj : Set.InjOn g S := by
    intro a ha b hb hab
    have h1 := hgspec a ((hmemS a).1 ha)
    have h2 := hgspec b ((hmemS b).1 hb)
    rw [← h1, ← h2, hab]
  have hcardS : (S.image g).card = S.card := Finset.card_image_of_injOn hinj
  have htors : ∀ p ∈ S.image g, l • p = 0 := by
    intro p hp
    obtain ⟨P, hPS, rfl⟩ := Finset.mem_image.1 hp
    have hP := (hmemS P).1 hPS
    refine rhPoint_injective (algebraMap K₀ M) W ?_
    rw [rhPoint_nsmul, rhPoint_zero, hgspec P hP]
    exact hP
  have hle := hcard (S.image g) htors
  rw [hcardS] at hle
  -- 個数の矛盾
  have hSc : S.card = l ^ 2 := by
    haveI := hfin.fintype
    rw [hS, Set.Finite.card_toFinset, ← Nat.card_eq_fintype_card]
    exact hcard2
  rw [hSc] at hle
  have h2 : 2 ≤ l := hl.two_le
  nlinarith [hle, h2]

/-! ## ★出典の紐付け(`.src`) -/

def exists_torsion_not_rational.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体の l-捩れが l 個以下なら有理でない l-捩れ点がある。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_torsion_not_rational.needs : List ProofObligation :=
  [ .citation "[ABC3]" "torsion_card(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.torsion_card") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1306）**——`h1` の入力の**個数の勘定**である。" ++
       "☆代数閉体の上では `l²` 個あるので、基礎体の上が `l` 個以下なら合わない。") 2 ]

end ABC3.Found.GenEll
