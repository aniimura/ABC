import ABC3.Found.GaloisRep.PointPlace
import ABC3.Found.GaloisRep.FiberSum

/-!
# Galois (G5) 第 173 ブロック —— **★★★★★★台を 2 つに分けて総和を 0 にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★係数が一致する必要はない

D2 の残りは `Σ_v e_v · Q_v = 0` の組み立てである。★台は

* `[n]⁻¹(P)`(場合 A、`count > 0`)
* `E[n] ∖ {O}`(場合 B、`count < 0`)

の 2 つに分かれる(第 171 で判定が `n·Q_v` の言葉になった)。

★★**2 つの係数 α, β が一致する必要はない**:

    S = α·(Σ_{[n]⁻¹(P)} Q) + β·(Σ_{T ∈ E[n]∖{O}} T) = α·0 + β·0 = 0

`Σ_{[n]⁻¹(P)} Q = n·(n·Q₀) = n·P = 0`(第 150 + `n·P = 0`)と
`Σ_{T ∈ E[n]} T = 0`(第 150)がそれぞれ**独立に** 0 だからである。

★★★これは**分岐指数の大域的な一定性を示さずに済む**ということである。
通常の証明は `[n]` の不分岐性(`e = 1`)を経由するが、ここでは要らない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sum_smul_const_eq_zero` | ★★係数が定数で総和が 0 なら 0 |
| `sum_split_eq_zero` | ★★★★★★**台を 2 つに分けて総和を 0 に** |
| `sum_fiber_zero` | ★★★★★ファイバーの総和は 0 |
| `sum_torsion_erase_zero` | ★★★★`E[n] ∖ {O}` の総和は 0 |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep
open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★抽象的な総和の補題 -/

/-- ★★係数が定数で、点の総和が 0 なら、重みつき総和も 0。 -/
theorem sum_smul_const_eq_zero {G : Type} [AddCommGroup G] (s : Finset G) (α : ℤ)
    (ε : G → ℤ) (hconst : ∀ Q ∈ s, ε Q = α) (hsum : ∑ Q ∈ s, Q = 0) :
    ∑ Q ∈ s, ε Q • Q = 0 := by
  rw [Finset.sum_congr rfl (fun Q hQ => by rw [hconst Q hQ]), ← Finset.smul_sum, hsum, smul_zero]

/-- ★★★★★★**台を 2 つに分けて総和を 0 にする**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★2 つの係数 `α`、`β` が一致する必要はない。 -/
theorem sum_split_eq_zero {G : Type} [AddCommGroup G] [DecidableEq G]
    (s s₁ s₂ : Finset G) (ε : G → ℤ) (α β : ℤ)
    (hsub : s ⊆ s₁ ∪ s₂) (hdisj : Disjoint s₁ s₂)
    (hzero : ∀ Q ∈ s₁ ∪ s₂, Q ∉ s → ε Q = 0)
    (h1 : ∀ Q ∈ s₁, ε Q = α) (h2 : ∀ Q ∈ s₂, ε Q = β)
    (hs1 : ∑ Q ∈ s₁, Q = 0) (hs2 : ∑ Q ∈ s₂, Q = 0) :
    ∑ Q ∈ s, ε Q • Q = 0 := by
  have hext : ∑ Q ∈ s, ε Q • Q = ∑ Q ∈ s₁ ∪ s₂, ε Q • Q :=
    Finset.sum_subset hsub (fun Q hQ hQs => by rw [hzero Q hQ hQs, zero_smul])
  rw [hext, Finset.sum_union hdisj,
    sum_smul_const_eq_zero s₁ α ε h1 hs1, sum_smul_const_eq_zero s₂ β ε h2 hs2, add_zero]

/-! ## ★★★★★2 つの台の総和 -/

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve F)
  [DecidableEq W.toAffine.Point]

/-- ★★★★★**ファイバーの総和は 0**(`n·P = 0` のとき)。 -/
theorem sum_fiber_zero [IsAlgClosed F] [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)] {P : W.toAffine.Point} (Q₀ : W.toAffine.Point)
    (hQ₀ : n • Q₀ = P) (hP : n • P = 0) :
    ∑ Q ∈ (Finset.univ : Finset (torsionPoints W n)).image
      (fun T : torsionPoints W n => Q₀ + (T : W.toAffine.Point)), Q = 0 := by
  rw [Finset.sum_image (by
    intro a _ b _ hab
    exact Subtype.ext (add_right_injective Q₀ hab))]
  rw [sum_fiber_eq W n hn hchar Q₀, hQ₀, hP]

/-- ★★★★**`E[n] ∖ {O}` の総和は 0**。 -/
theorem sum_torsion_erase_zero [IsAlgClosed F] [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)] :
    ∑ Q ∈ ((Finset.univ : Finset (torsionPoints W n)).image
      (fun T : torsionPoints W n => (T : W.toAffine.Point))).erase 0, Q = 0 := by
  have herase : ∀ s : Finset W.toAffine.Point, ∑ Q ∈ s.erase 0, Q = ∑ Q ∈ s, Q :=
    fun s => Finset.sum_erase s rfl
  rw [herase, Finset.sum_image (by
    intro a _ b _ hab
    exact Subtype.ext hab)]
  exact sum_torsion_eq_zero W n hn hchar

/-! ## ★出典の紐付け(`.src`) -/

def sum_split_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——台を 2 つに分けて総和を 0 にすること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
