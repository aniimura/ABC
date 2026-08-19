import ABC3.Found.Arakelov.ArcXNorm

/-!
# Arakelov (C3) 第 272 ブロック —— ★★★★★**貼り合わせ(`xNorm` 版)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 262–263 の書き直し —— **仮定が 1 つ消えた**

第 262–263 は `genNorm`(迂回路)に対して書いてあった。
★第 270–271 で `xNorm`(橋から直接)が出たので書き直す。

★★★**違いは仮定が 1 つ消えたこと**である:

| 版 | `gluedNorm = 0 → w = 0` に要るもの |
|---|---|
| 第 263(`genNorm`) | `ρ i₀ p > 0` + `p ∈ U i₀` + ★**生成切断の非消滅** |
| ★本ブロック(`xNorm`) | `ρ i₀ p > 0` + `p ∈ U i₀` のみ |

★★迂回路では「比の分母が 0 でない」ことに非消滅が要ったが、
橋があれば `xNorm` は最初から非退化である。

| 定義・定理 | 内容 |
|---|---|
| `extNormX` / `gluedNormX` | ★0 で延長・1 の分割で貼る |
| `extNormX_nonneg` / `_smul` | ★各項の法則 |
| `gluedNormX_nonneg` / `_smul` | ★★和の法則 |
| `gluedNormX_eq_zero_iff` | ★★★★**非退化**(仮定 2 つ) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open scoped Classical

variable {X : Scheme.{0}} {ι : Type}
/-- ★局所ノルムを 0 で延長する。 -/
noncomputable def extNormX (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (i : ι) (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  if h : p ⁻¹ᵁ (U i) = ⊤ then xNorm (U i) F (e i) p h w else 0

/-- ★★貼り合わせたノルム。 -/
noncomputable def gluedNormX (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  ∑ᶠ i, ρ i p * extNormX F U e i p w

theorem extNormX_nonneg (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (i : ι) (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    0 ≤ extNormX F U e i p w := by
  by_cases h : p ⁻¹ᵁ (U i) = ⊤
  · simp only [extNormX, dif_pos h]
    exact xNorm_nonneg (U i) F (e i) p h w
  · simp only [extNormX, dif_neg h]
    exact le_rfl

theorem extNormX_smul (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (i : ι) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    extNormX F U e i p (c • w) = ‖c‖ * extNormX F U e i p w := by
  by_cases h : p ⁻¹ᵁ (U i) = ⊤
  · simp only [extNormX, dif_pos h]
    exact xNorm_smul (U i) F (e i) p h c w
  · simp only [extNormX, dif_neg h]
    rw [mul_zero]

theorem gluedNormX_nonneg (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ) (hρ : ∀ i p, 0 ≤ ρ i p)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    0 ≤ gluedNormX F U e ρ p w :=
  finsum_nonneg (fun i => mul_nonneg (hρ i p) (extNormX_nonneg F U e i p w))

theorem gluedNormX_smul (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ)
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    gluedNormX F U e ρ p (c • w) = ‖c‖ * gluedNormX F U e ρ p w := by
  show (∑ᶠ i, ρ i p * extNormX F U e i p (c • w))
    = ‖c‖ * ∑ᶠ i, ρ i p * extNormX F U e i p w
  rw [mul_finsum]
  refine finsum_congr (fun i => ?_)
  rw [extNormX_smul]
  ring


theorem gluedNormX_eq_zero_iff (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ) (hρ : ∀ i p, 0 ≤ ρ i p)
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (i₀ : ι) (hpos : 0 < ρ i₀ p) (h₀ : p ⁻¹ᵁ (U i₀) = ⊤)
    (w : ↥(arcFiber p F))
    (hfw : (Function.support (fun i => ρ i p * extNormX F U e i p w)).Finite) :
    gluedNormX F U e ρ p w = 0 ↔ w = 0 := by
  constructor
  · intro hz
    have hle : ρ i₀ p * extNormX F U e i₀ p w
        ≤ ∑ᶠ i, ρ i p * extNormX F U e i p w :=
      single_le_finsum i₀ hfw (fun j => mul_nonneg (hρ j p) (extNormX_nonneg F U e j p w))
    rw [show (∑ᶠ i, ρ i p * extNormX F U e i p w) = gluedNormX F U e ρ p w from rfl, hz] at hle
    have hzero : ρ i₀ p * extNormX F U e i₀ p w = 0 :=
      le_antisymm hle (mul_nonneg (hρ i₀ p) (extNormX_nonneg F U e i₀ p w))
    have he : extNormX F U e i₀ p w = 0 := by
      rcases mul_eq_zero.1 hzero with h | h
      · exact absurd h (ne_of_gt hpos)
      · exact h
    rw [show extNormX F U e i₀ p w = xNorm (U i₀) F (e i₀) p h₀ w from
      by simp only [extNormX, dif_pos h₀]] at he
    exact (xNorm_eq_zero_iff (U i₀) F (e i₀) p h₀ w).1 he
  · intro hw
    show (∑ᶠ i, ρ i p * extNormX F U e i p w) = 0
    have hall : ∀ i, ρ i p * extNormX F U e i p w = 0 := by
      intro i
      by_cases hi : p ⁻¹ᵁ (U i) = ⊤
      · simp only [extNormX, dif_pos hi]
        rw [(xNorm_eq_zero_iff (U i) F (e i) p hi w).2 hw, mul_zero]
      · simp only [extNormX, dif_neg hi, mul_zero]
    simp only [hall]
    exact finsum_zero


/-! ## ★出典の紐付け(`.src`) -/

def gluedNormX_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——貼り合わせ、xNorm 版)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
