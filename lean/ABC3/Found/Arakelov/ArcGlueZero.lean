import ABC3.Found.Arakelov.ArcGlue

/-!
# Arakelov (C3) 第 263 ブロック —— ★★★★**貼り合わせたノルムの非退化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★1 の分割の「どこかで正」が非退化を与える

`gluedNorm p w = ∑ᶠ i, ρ i p * extNorm i p w` は各項が非負なので、
★**1 つでも正の項があれば和は正**である。

★★1 の分割は各点で `∑ ρ i p = 1` なので `ρ i₀ p > 0` なる `i₀` が在り、
subordinate なので `p ∈ U i₀`、したがって `extNorm i₀ p w = genNorm … w` である。

★★★これで `gluedNorm p w = 0 → w = 0` が出る(`single_le_finsum` で 1 項を取り出す)。

## ★逆向きは仮定が要らない

`genNorm … 0 = 0` は**非消滅の仮定なしで**成り立つ(`genNorm_zero`)——
分子が `nrm(0) = 0` だからである。★これで `w = 0 → gluedNorm = 0` は
各項ごとに閉じる。

| 定理 | 内容 |
|---|---|
| `genNorm_zero` | ★`genNorm … 0 = 0`(仮定なし) |
| `gluedNorm_eq_zero_iff` | ★★★★**非退化** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite
open scoped Classical

variable {X : Scheme.{0}} {ι : Type}

theorem genNorm_zero (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤) :
    genNorm F hF V g p h 0 = 0 := by
  show (arcMetricOf F hF).nrm p 0 / _ = 0
  rw [((arcMetricOf F hF).eq_zero_iff p 0).2 rfl, zero_div]

theorem gluedNorm_eq_zero_iff (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ) (hρ : ∀ i p, 0 ≤ ρ i p)
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (i₀ : ι) (hpos : 0 < ρ i₀ p) (h₀ : p ⁻¹ᵁ (U i₀) = ⊤)
    (hg₀ : arcEvalOnTop F p (U i₀) h₀ (g i₀) ≠ 0)
    (w : ↥(arcFiber p F))
    (hfw : (Function.support (fun i => ρ i p * extNorm F hF U g i p w)).Finite) :
    gluedNorm F hF U g ρ p w = 0 ↔ w = 0 := by
  constructor
  · intro hz
    have hle : ρ i₀ p * extNorm F hF U g i₀ p w
        ≤ ∑ᶠ i, ρ i p * extNorm F hF U g i p w :=
      single_le_finsum i₀ hfw (fun j => mul_nonneg (hρ j p) (extNorm_nonneg F hF U g j p w))
    rw [show (∑ᶠ i, ρ i p * extNorm F hF U g i p w) = gluedNorm F hF U g ρ p w from rfl, hz] at hle
    have hnn : 0 ≤ ρ i₀ p * extNorm F hF U g i₀ p w :=
      mul_nonneg (hρ i₀ p) (extNorm_nonneg F hF U g i₀ p w)
    have hzero : ρ i₀ p * extNorm F hF U g i₀ p w = 0 := le_antisymm hle hnn
    have he : extNorm F hF U g i₀ p w = 0 := by
      rcases mul_eq_zero.1 hzero with h | h
      · exact absurd h (ne_of_gt hpos)
      · exact h
    rw [show extNorm F hF U g i₀ p w = genNorm F hF (U i₀) (g i₀) p h₀ w from
      by simp only [extNorm, dif_pos h₀]] at he
    exact (genNorm_eq_zero_iff F hF (U i₀) (g i₀) p h₀ hg₀ w).1 he
  · intro hw
    show (∑ᶠ i, ρ i p * extNorm F hF U g i p w) = 0
    have : ∀ i, ρ i p * extNorm F hF U g i p w = 0 := by
      intro i
      rw [hw]
      by_cases hi : p ⁻¹ᵁ (U i) = ⊤
      · simp only [extNorm, dif_pos hi]
        rw [genNorm_zero F hF (U i) (g i) p hi, mul_zero]
      · simp only [extNorm, dif_neg hi, mul_zero]
    simp only [this]
    exact finsum_zero


/-! ## ★出典の紐付け(`.src`) -/

def gluedNorm_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——貼り合わせたノルムの非退化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
