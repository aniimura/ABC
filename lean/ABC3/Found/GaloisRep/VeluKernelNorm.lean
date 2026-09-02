/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DiscIdentitySemistable
import ABC3.Found.GaloisRep.TorsionIntegralGood
import ABC3.Meta.Claim

/-!
# 第 1386 ブロック —— **Vélu の核のノルム `N` と、その整性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1385 が要る `N`

第 1385（`semistableAt_of_disc_pow_eq`）は

    Δ(E)^l = Δ(E′) · N^4,   `N` は `p` 上整

という形の恒等式を受ける。★本ブロックはその `N` を

    N ≔ ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃)

と定義し、**核の座標が整なら `N` も整**であることを示す。

☆核の座標の整性は第 1073（`pointCoords_mem_of_addOrderOf_prime`、証明済み）が与える
——ただし `p ∤ l` が要る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset ABC3.Meta
open ABC3.Found.GenEll

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★**Vélu の核のノルム**（第 1386）。

☆`2y + a₁x + a₃` は `P ↦ −P` で符号が変わる量であり、
`P` が 2-捩れでないことと `≠ 0` は同値である。 -/
noncomputable def veluKernelNorm (W : WeierstrassCurve L) (S : Finset (L × L)) : L :=
  ∏ q ∈ S, (2 * q.2 + W.a₁ * q.1 + W.a₃)

/-- ★★★★**`p` 上整な曲線の `a₁` は整**（第 1386）。 -/
theorem a₁_mem_primeSubring (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] : W.a₁ ∈ primeSubring p := by
  rw [← WeierstrassCurve.integralModel_a₁_eq (primeSubring p) W]
  exact ((WeierstrassCurve.integralModel (primeSubring p) W).a₁).2

/-- ★★★★**`p` 上整な曲線の `a₃` は整**（第 1386）。 -/
theorem a₃_mem_primeSubring (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] : W.a₃ ∈ primeSubring p := by
  rw [← WeierstrassCurve.integralModel_a₃_eq (primeSubring p) W]
  exact ((WeierstrassCurve.integralModel (primeSubring p) W).a₃).2

/-- ★★★★★★★★★★★★
**核の座標が整なら `N` も整**——★**無条件**（第 1386）。 -/
theorem veluKernelNorm_mem_primeSubring (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (S : Finset (L × L))
    (hS : ∀ q ∈ S, q.1 ∈ primeSubring p ∧ q.2 ∈ primeSubring p) :
    veluKernelNorm W S ∈ primeSubring p := by
  refine prod_mem ?_
  intro q hq
  obtain ⟨hx, hy⟩ := hS q hq
  have h2 : (2 : L) ∈ primeSubring p := by
    have : (2 : L) = 1 + 1 := by norm_num
    rw [this]
    exact add_mem (one_mem _) (one_mem _)
  exact add_mem (add_mem (mul_mem h2 hy) (mul_mem (a₁_mem_primeSubring p W) hx))
    (a₃_mem_primeSubring p W)

/-- ★★★★★★★★★★★★★★★★
**Vélu の核の座標は整**——`p ∤ l` のとき（第 1386）。

☆第 1073（`pointCoords_mem_of_addOrderOf_prime`、証明済み）を
`veluQuotientFull` が受ける形の `Finset` に直しただけである。 -/
theorem veluKernel_coords_mem (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    ∀ q ∈ ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)),
      q.1 ∈ primeSubring p ∧ q.2 ∈ primeSubring p := by
  intro q hq
  obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hq
  rw [mem_erase, mem_range] at hk
  exact pointCoords_mem_of_addOrderOf_prime p E hl hlu Q hQ hk.1 hk.2

/-! ## ★出典の紐付け(`.src`) -/

def veluKernelNorm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の核のノルム N の定義)",
    sectionId := "genell-lemma-3-5" }

def veluKernelNorm_mem_primeSubring.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核の座標が整なら N も整。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluKernel_coords_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の核の座標は整——p ∤ l のとき。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluKernelNorm_mem_primeSubring.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pointCoords_mem_of_addOrderOf_prime(第 1073、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_mem_of_addOrderOf_prime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1386）**——第 1385 が要る `N` を作る段である。" ++
       "☆残るのは恒等式 `Δ(E)^l = Δ(E′)·N^4` 1 本と、`p ∣ l` の場合である。") 17 ]

end ABC3.Found.GaloisRep
