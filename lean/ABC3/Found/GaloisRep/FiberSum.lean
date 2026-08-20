import ABC3.Found.GaloisRep.DvdCount
import ABC3.Found.GaloisRep.TorsionStructure

/-!
# Galois (G5) 第 150 ブロック —— **★★★★★ファイバーの総和**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★D2 の心臓——ファイバーの総和が 0 になる

D2(`J` が単項)の類の計算で使うのは次の 2 つである:

    Σ_{T ∈ E[n]} T = 0,        Σ_{Q ∈ [n]⁻¹(P)} Q = n²·Q₀ = n·P = 0

★本ブロックはその両方を出す。

### ★★★★(G1) が効く

**(G1) は達成済**であり `E[n] ≅ (ℤ/n)²` を与える(`torsion_addEquiv`)。
★第 140 の `sum_univ_eq_zero_of_addEquiv`(有限可換群の直積では全元の和が 0)を
当てるだけで `Σ_{T ∈ E[n]} T = 0` が出る。
★★`#E[n] = n²` も同じ同型から出る。

### ★★★ファイバー

`n·Q₀ = P` なる `Q₀` を 1 つ取ると、ファイバーは `Q₀ + E[n]` であり

    Σ_{T ∈ E[n]} (Q₀ + T) = #E[n]·Q₀ + Σ_T T = n²·Q₀ = n·(n·Q₀) = n·P

★`P` が `n` 等分点なら **これは 0** である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sum_torsion_eq_zero` | ★★★★★**`Σ_{T ∈ E[n]} T = 0`** |
| `card_torsion` | ★★★★**`#E[n] = n²`** |
| `sum_fiber_eq` | ★★★★★**`Σ_{T ∈ E[n]} (Q₀ + T) = n·(n·Q₀)`** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★**`n` 等分点の総和は 0 である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★(G1) の `E[n] ≅ (ℤ/n)²` に第 140 の
`sum_univ_eq_zero_of_addEquiv`(有限可換群の直積では全元の和が 0)を当てる。 -/
theorem sum_torsion_eq_zero [IsAlgClosed F] [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)] :
    ∑ T : torsionPoints W n, (T : W.toAffine.Point) = 0 := by
  haveI : NeZero n := ⟨by omega⟩
  have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  obtain ⟨e⟩ := torsion_addEquiv W hΔ n hn hchar
  have hzero : ∑ T : torsionPoints W n, T = 0 :=
    sum_univ_eq_zero_of_addEquiv (A := torsionPoints W n) e.symm
  have hmap := congrArg (AddSubgroup.subtype (torsionPoints W n)) hzero
  rw [map_sum, map_zero] at hmap
  exact hmap

/-- ★★★★**`#E[n] = n²`**——(G1) の同型から。 -/
theorem card_torsion [IsAlgClosed F] [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    Nat.card (torsionPoints W n) = n * n := by
  haveI : NeZero n := ⟨by omega⟩
  have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  obtain ⟨e⟩ := torsion_addEquiv W hΔ n hn hchar
  rw [torsionPoints_eq_ker, ← Nat.card_congr e.toEquiv, Nat.card_prod, Nat.card_zmod]

/-- ★★★★★**ファイバーの総和**——`Σ_{T ∈ E[n]} (Q₀ + T) = n·(n·Q₀)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`n·Q₀ = P` かつ `n·P = 0` なら **これは 0** である。
★★これが D2(因子の類が自明であること)の計算の核になる。 -/
theorem sum_fiber_eq [IsAlgClosed F] [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)] (Q₀ : W.toAffine.Point) :
    ∑ T : torsionPoints W n, (Q₀ + (T : W.toAffine.Point)) = n • (n • Q₀) := by
  rw [Finset.sum_add_distrib, sum_torsion_eq_zero W n hn hchar, add_zero,
    Finset.sum_const, Finset.card_univ, ← Nat.card_eq_fintype_card,
    card_torsion W n hn hchar, mul_smul]

/-! ## ★出典の紐付け(`.src`) -/

def sum_torsion_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——n 等分点の総和が 0 でありファイバーの総和が n·P であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
