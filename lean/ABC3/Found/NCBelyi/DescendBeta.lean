/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.NestedInductionBeta
import ABC3.Found.NCBelyi.SeparationStep
import ABC3.Found.NCBelyi.RatSeparation
import ABC3.Found.NCBelyi.DescendData
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`Lemma 2.4` の `β` つき降下段（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 降下段が丸ごと取れた

`Lemma 2.4` の帰納段は、原文では次の 1 続きである:

> applying an automorphism (with rational coefficients!) as in Lemma 2.3 and then
> multiplying by some positive rational number, we may assume that |α| ≤1, ...
> Thus, for a suitable choice of C, it follows that f0(β) /∈S′ ...
> replacing S by S′, β by f0(β), applying the induction hypothesis

★★★**本ファイルはこの 1 続きを 1 本の定理にする**（`exists_descend_data_beta`）:

    `S`（`Gal`-安定・代数的・`m(S) > 0`）と `β ∈ ℚ∖S` から
    `T`（`Gal`-安定・代数的）と `β′ ∈ ℚ∖T` が作れて、**測度が下がる**。

## ★組み立ての順序

| 段 | 使うもの |
|---|---|
| 正規化（`‖α‖ ≤ 1`、`‖β‖ ≥ 4`） | `RatSeparation.lean` の `exists_rat_normalization`（第 409） |
| 正規化で測度が壊れない | `MobiusRedDeg.lean`（`§9-983`）＋ 本ファイルの `dSum_mobius` |
| 正規化で `Gal`-安定性が壊れない | `MobiusConj.lean`（`§9-985`） |
| 最大層の点 `x₀` と `f₀ ≔ minpoly x₀` | `DescendData.lean` の `exists_descend_data`（第 416） |
| **`f₀(β) ∉ S′`** | `SeparationStep.lean`（`§9-984`） |
| `S′` が `Gal`-安定 | ★本ファイルの `conjSet_aeval_subset` |
| 測度が下がる | `NestedInduction.lean` の `measure_lt`（第 408） |

## ★★★★★`conjSet_aeval_subset` —— `§9-985` の多項式版

    **`conjSet (p(x)) ⊆ p(conjSet x)`**   （`p ∈ ℚ[X]`）

★機構は**埋め込みの値域**である——mathlib の

    `Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly :
       (Set.range fun ψ : K →ₐ[F] A => ψ x) = (minpoly F x).rootSet A`

を `K ≔ ℚ⟮x⟯`、`A ≔ ℂ` で使う。★★`z` が `p(x)` の共役なら
`ψ : ℚ⟮x⟯ →ₐ[ℚ] ℂ` が取れて `ψ(p(gen)) = z`、
そして `ψ` は `ℚ`-代数準同型だから `z = p(ψ gen)` かつ `ψ gen ∈ conjSet x` である。

## ★★★これで `Lemma 2.4` に残るのは**配管 1 つだけ**

★☆**合成を `ℙ¹` の有理写像として組み立てること**——原文の `f(x) ∈ ℚ(x)` は
Möbius と多項式の交互合成である。`Separation.lean` の `P1C` が受け皿。
★★数値の核・測度・`Gal`-安定性・`β` を運ぶ帰納・**降下段**はすべて取れた。
-/

namespace ABC3.Found.NCBelyi

open Polynomial IntermediateField

/-! ## ★★★★★★★★★★★★共役は多項式の像で写る（`§9-985` の多項式版） -/

/-- ★★★★★★★★★★★★**`conjSet (p(x)) ⊆ p(conjSet x)`**（`p ∈ ℚ[X]`）。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

★原文が `S′ ≔ f₀(S) ∪ f₀(S₀)` を `Gal`-安定と言い切るところの中身である。
★★機構は**埋め込みの値域**——mathlib の
`Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly` を `ℚ⟮x⟯ →ₐ[ℚ] ℂ` で使う。 -/
theorem conjSet_aeval_subset (p : ℚ[X]) (x : ℂ) (hx : IsIntegral ℚ x) :
    conjSet (aeval x p) ⊆ (conjSet x).image (fun w => aeval w p) := by
  classical
  haveI : FiniteDimensional ℚ ℚ⟮x⟯ := IntermediateField.adjoin.finiteDimensional hx
  haveI : Algebra.IsAlgebraic ℚ ℚ⟮x⟯ := Algebra.IsAlgebraic.of_finite ℚ ℚ⟮x⟯
  set gen : ℚ⟮x⟯ := ⟨x, mem_adjoin_simple_self ℚ x⟩ with hgen
  have hinj : Function.Injective (algebraMap ℚ⟮x⟯ ℂ) := (algebraMap ℚ⟮x⟯ ℂ).injective
  have hgx : (algebraMap ℚ⟮x⟯ ℂ) gen = x := rfl
  have hmap : (algebraMap ℚ⟮x⟯ ℂ) ((aeval gen) p) = aeval x p :=
    (Polynomial.aeval_algHom_apply (IsScalarTower.toAlgHom ℚ ℚ⟮x⟯ ℂ) gen p).symm
  intro z hz
  have hmy : minpoly ℚ (aeval gen p) = minpoly ℚ (aeval x p) := by
    have h : minpoly ℚ (algebraMap ℚ⟮x⟯ ℂ (aeval gen p)) = minpoly ℚ (aeval gen p) :=
      minpoly.algebraMap_eq hinj _
    rw [← h, hmap]
  have hzroot : z ∈ (minpoly ℚ (aeval gen p)).rootSet ℂ := by
    rw [Polynomial.mem_rootSet]
    refine ⟨minpoly.ne_zero (IsIntegral.of_finite ℚ _), ?_⟩
    rw [hmy]
    exact (mem_conjSet_iff (isIntegral_aeval hx p)).1 hz
  rw [← Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly ℂ (aeval gen p)] at hzroot
  obtain ⟨ψ, hψ⟩ := hzroot
  refine Finset.mem_image.2 ⟨ψ gen, ?_, ?_⟩
  · refine (mem_conjSet_iff hx).2 ?_
    have hgr : ψ gen ∈ (minpoly ℚ gen).rootSet ℂ := by
      rw [← Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly ℂ gen]
      exact ⟨ψ, rfl⟩
    have h2 : minpoly ℚ gen = minpoly ℚ x := by
      have h : minpoly ℚ (algebraMap ℚ⟮x⟯ ℂ gen) = minpoly ℚ gen := minpoly.algebraMap_eq hinj _
      rw [← h, hgx]
    rw [h2, Polynomial.mem_rootSet] at hgr
    exact hgr.2
  · rw [← hψ]
    exact Polynomial.aeval_algHom_apply ψ gen p

/-- ★★★★★**多項式の像も `Gal`-安定**。 -/
theorem isConjStable_image_aeval (p : ℚ[X]) (S : Finset ℂ)
    (hint : ∀ x ∈ S, IsIntegral ℚ x) (hstab : IsConjStable S) :
    IsConjStable (S.image (fun x => aeval x p)) := by
  classical
  intro y hy z hz
  obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
  obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 (conjSet_aeval_subset p x (hint x hxS) hz)
  exact Finset.mem_image.2 ⟨w, hstab x hxS hw, rfl⟩

/-- ★★★**`ℚ`-多項式の根の集まりは `Gal`-安定**。 -/
theorem isConjStable_rootSet (g : ℚ[X]) (hg : g ≠ 0) :
    IsConjStable (((g.map (algebraMap ℚ ℂ)).roots).toFinset) := by
  classical
  have hmapne : (g.map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hg
  intro x hx z hz
  rw [Multiset.mem_toFinset] at hx
  have hxr : aeval x g = 0 := by
    have h := (Polynomial.mem_roots hmapne).1 hx
    simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using h
  have hxint : IsIntegral ℚ x := isIntegral_of_aeval_eq_zero hg hxr
  have hdvd : minpoly ℚ x ∣ g := minpoly.dvd ℚ x hxr
  have hzr : aeval z g = 0 := by
    obtain ⟨r, hr⟩ := hdvd
    rw [hr, map_mul, (mem_conjSet_iff hxint).1 hz, zero_mul]
  rw [Multiset.mem_toFinset]
  refine (Polynomial.mem_roots hmapne).2 ?_
  simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using hzr

/-- ★`Gal`-安定性は合併で保たれる。 -/
theorem isConjStable_union {A B : Finset ℂ} (hA : IsConjStable A) (hB : IsConjStable B) :
    IsConjStable (A ∪ B) := by
  intro x hx z hz
  rcases Finset.mem_union.1 hx with h | h
  · exact Finset.mem_union_left _ (hA x h hz)
  · exact Finset.mem_union_right _ (hB x h hz)

/-! ## ★★★★★★`ℚ`-Möbius は `d(S)` も変えない -/

open scoped Classical in
/-- ★★★★★★**`ℚ`-Möbius は `d(S)` を変えない**——単射で `redDeg` を保つから。 -/
theorem dSum_mobius (lam nu mu : ℚ) (hnu : nu ≠ 0) (S : Finset ℂ)
    (hlam : ((lam : ℂ)) ∉ S) :
    dSum (S.image (fun x => (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ))) = dSum S := by
  classical
  have hnuC : ((nu : ℂ)) ≠ 0 := by exact_mod_cast hnu
  set σ : ℂ → ℂ := fun x => (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ) with hσ
  have hmax : maxRedDeg (S.image σ) = maxRedDeg S := maxRedDeg_mobius lam nu mu hnu S hlam
  have hinjOn : ∀ x y : ℂ, σ x = σ y → x = y := by
    intro x y h
    simp only [hσ, add_left_inj] at h
    have h2 : ((x - (lam : ℂ))⁻¹) = ((y - (lam : ℂ))⁻¹) := mul_left_cancel₀ hnuC h
    exact sub_left_inj.1 (inv_inj.1 h2)
  have htop : topLayer (S.image σ) = (topLayer S).image σ := by
    ext y
    simp only [mem_topLayer_iff, Finset.mem_image, hmax]
    constructor
    · rintro ⟨⟨x, hxS, rfl⟩, hdeg⟩
      have hx' : x ≠ (lam : ℂ) := fun hc => hlam (hc ▸ hxS)
      refine ⟨x, ⟨hxS, ?_⟩, rfl⟩
      rwa [redDeg_mobius lam nu mu hnu x hx'] at hdeg
    · rintro ⟨x, ⟨hxS, hdeg⟩, rfl⟩
      have hx' : x ≠ (lam : ℂ) := fun hc => hlam (hc ▸ hxS)
      exact ⟨⟨x, hxS, rfl⟩, by rw [redDeg_mobius lam nu mu hnu x hx']; exact hdeg⟩
  rw [dSum_eq, dSum_eq, hmax, htop, Finset.card_image_of_injOn]
  intro x _ y _ h
  exact hinjOn x y h

/-! ## ★★★★★★★正規化を `ℂ` の言葉に直す -/

/-- ★★★★★★★**正規化**（`Lemma 2.3` ＋ 正の有理数倍）を `ℂ` の言葉で。

原文 (NCBelyi p.5):
> multiplying by some positive rational number, we may assume that |α| ≤1, for all

    `‖c/(α−λ)‖ ≤ 1`（`∀ α ∈ S`）  かつ  `‖c/(β−λ)‖ ≥ 4`

★`RatSeparation.lean` の `exists_rat_normalization` を `C = 4` で引き、
`P1C` の言葉から `ℂ` の言葉へ移しただけである。 -/
theorem exists_normalization_complex (S : Finset ℂ) (β : ℚ) (hβ : ((β : ℂ)) ∉ S) :
    ∃ (lam c : ℚ), 0 < c ∧ ((lam : ℂ)) ∉ S ∧ (lam : ℚ) ≠ β
      ∧ (∀ α ∈ S, ‖(c : ℂ) * (α - (lam : ℂ))⁻¹‖ ≤ 1)
      ∧ 4 ≤ ‖(c : ℂ) * ((β : ℂ) - (lam : ℂ))⁻¹‖ := by
  classical
  have hb : (some ((β : ℂ))) ∉ S.image some := by
    intro h
    obtain ⟨α, hα, heq⟩ := Finset.mem_image.1 h
    exact hβ (by rwa [Option.some_inj.1 heq] at hα)
  obtain ⟨lam, c, hc0, hne, hnotmem, hle1, hge⟩ :=
    exists_rat_normalization (S.image some) 4 (by norm_num) β hb
  have hcC : ‖((c : ℂ))‖ = (c : ℝ) := by
    rw [Complex.norm_ratCast]
    have h : (0:ℝ) < (c:ℝ) := by exact_mod_cast hc0
    exact abs_of_pos h
  refine ⟨lam, c, hc0, ?_, ?_, ?_, ?_⟩
  · intro h
    exact hnotmem (Finset.mem_image.2 ⟨(lam : ℂ), h, rfl⟩)
  · intro h
    exact hne (by rw [h])
  · intro α hα
    have h := hle1 (some (α : ℂ)) (Finset.mem_image.2 ⟨α, hα, rfl⟩)
    rw [absInvShift_some] at h
    rw [norm_mul, hcC, one_div] at *
    exact h
  · rw [absInvShift_some] at hge
    rw [norm_mul, hcC, one_div] at *
    exact hge

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`β` つきの降下段 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`Lemma 2.4` の降下段（`β` つき）**。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

`S`（`Gal`-安定・代数的・`m(S) > 0`）と `β ∈ ℚ`（`β ∉ S`）から、

* `T`（`Gal`-安定・代数的）と `β′ ∈ ℚ` で **`β′ ∉ T`**
* **測度が下がる**（`m(T) < m(S)`、または `m(T) = m(S)` かつ `d(T) < d(S)`）

が作れる。★★★これが原文の
『replacing S by S′, β by f₀(β), applying the induction hypothesis』の**前半**である。

★具体的には `T = f₀(σ(S)) ∪ f₀(σ(S)₀)`、`β′ = f₀(σ(β))`
（`σ` は正規化の Möbius、`f₀` は最大層の点の最小多項式）。 -/
theorem exists_descend_data_beta (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x)
    (hstab : IsConjStable S) (hS : 0 < maxRedDeg S) (β : ℚ) (hβ : ((β : ℂ)) ∉ S) :
    ∃ (T : Finset ℂ) (β' : ℚ),
      (∀ x ∈ T, IsIntegral ℚ x) ∧ IsConjStable T ∧ ((β' : ℂ)) ∉ T
      ∧ (maxRedDeg T < maxRedDeg S ∨ (maxRedDeg T = maxRedDeg S ∧ dSum T < dSum S)) := by
  classical
  obtain ⟨lam, c, hc0, hlamS, hlamβ, hnorm1, hnorm4⟩ := exists_normalization_complex S β hβ
  have hcne : c ≠ 0 := ne_of_gt hc0
  set σ : ℂ → ℂ := fun x => (c : ℂ) * (x - (lam : ℂ))⁻¹ + (((0:ℚ)) : ℂ) with hσ
  set S₁ : Finset ℂ := S.image σ with hS1
  have hlamne : ∀ x ∈ S, x ≠ (lam : ℂ) := fun x hx hc => hlamS (hc ▸ hx)
  have hS1int : ∀ y ∈ S₁, IsIntegral ℚ y := by
    intro y hy
    obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
    exact isIntegral_mobius lam c 0 (hint x hxS)
  have hS1stab : IsConjStable S₁ := isConjStable_mobius lam c 0 hcne S hint hlamne hstab
  have hS1max : maxRedDeg S₁ = maxRedDeg S := maxRedDeg_mobius lam c 0 hcne S hlamS
  have hS1d : dSum S₁ = dSum S := dSum_mobius lam c 0 hcne S hlamS
  have hS1norm : ∀ y ∈ S₁, ‖y‖ ≤ 1 := by
    intro y hy
    obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
    simpa [hσ] using hnorm1 x hxS
  have hS1pos : 0 < maxRedDeg S₁ := by rw [hS1max]; exact hS
  obtain ⟨E, x₀, hx₀S, hx₀top, hno, hE, hx₀drop, hEdef⟩ := exists_descend_data S₁ hS1int hS1pos
  have hx₀int : IsIntegral ℚ x₀ := hS1int x₀ hx₀S
  have hdeg : 0 < (minpoly ℚ x₀).natDegree := by
    rw [natDegree_minpoly_eq_redDeg_succ hx₀int]; omega
  have hd0 : derivative (minpoly ℚ x₀) ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := minpoly ℚ x₀)).1 hc
    omega
  have hmapne : ((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
  set β₁ : ℚ := c / (β - lam) with hβ1
  have hβlam : (β : ℚ) - lam ≠ 0 := sub_ne_zero.mpr (fun h => hlamβ h.symm)
  have hβ1C : ((β₁ : ℂ)) = σ ((β : ℂ)) := by
    rw [hβ1, hσ]
    push_cast
    rw [div_eq_mul_inv, add_zero]
  have hβ1norm : 3 < ‖((β₁ : ℂ))‖ := by
    rw [hβ1C, hσ]
    simp only [Rat.cast_zero, add_zero]
    linarith [hnorm4]
  have hRint : ∀ w ∈ ((((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)).roots).toFinset),
      IsIntegral ℚ w := by
    intro w hw
    rw [Multiset.mem_toFinset] at hw
    have hwr : aeval w (derivative (minpoly ℚ x₀)) = 0 := by
      have h2 := (Polynomial.mem_roots hmapne).1 hw
      simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using h2
    exact isIntegral_of_aeval_eq_zero hd0 hwr
  refine ⟨S₁.image (fun y => aeval y (minpoly ℚ x₀)) ∪ E, (minpoly ℚ x₀).eval β₁,
    ?_, ?_, ?_, ?_⟩
  · -- ★代数的
    intro x hx
    rcases Finset.mem_union.1 hx with h | h
    · obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 h
      exact isIntegral_aeval (hS1int y hy) _
    · rw [hEdef] at h
      obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 h
      exact isIntegral_aeval (hRint w hw) _
  · -- ★★`Gal`-安定
    rw [hEdef]
    exact isConjStable_union (isConjStable_image_aeval _ S₁ hS1int hS1stab)
      (isConjStable_image_aeval _ _ hRint
        (isConjStable_rootSet (derivative (minpoly ℚ x₀)) hd0))
  · -- ★★★`β′ ∉ T` —— 分離
    have hcast : ((((minpoly ℚ x₀).eval β₁ : ℚ)) : ℂ) = aeval ((β₁ : ℂ)) (minpoly ℚ x₀) := by
      have h : ((β₁ : ℂ)) = algebraMap ℚ ℂ β₁ := by simp
      rw [h, Polynomial.aeval_algebraMap_apply]
      simp
    rw [hcast, hEdef]
    exact separation_step_not_mem S₁ hS1stab hS1norm hx₀S hx₀int hdeg hβ1norm
  · -- ★★★★測度が下がる
    have hmeas := measure_lt S₁ E (fun y => aeval y (minpoly ℚ x₀)) hno hE hx₀S hx₀top
      hx₀drop hS1pos
    rw [hS1max, hS1d] at hmeas
    exact hmeas

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def conjSet_aeval_subset.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(conjSet (p(x)) ⊆ p(conjSet x))",
    sectionId := "ncbelyi-lemma-2-4" }

def exists_normalization_complex.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(正規化を ℂ の言葉で——‖α‖ ≤ 1 かつ ‖β‖ ≥ 4)",
    sectionId := "ncbelyi-lemma-2-4" }

def dSum_mobius.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ-Möbius は d(S) を変えない)",
    sectionId := "ncbelyi-lemma-2-4" }

def exists_descend_data_beta.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(β つきの降下段)",
    sectionId := "ncbelyi-lemma-2-4" }

def exists_descend_data_beta.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_rat_normalization(正規化、第 409)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.exists_rat_normalization") 3,
    .citation "[ABC3]" "separation_step_not_mem(f₀(β) ∉ S′、§9-984)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.separation_step_not_mem") 4,
    .citation "[ABC3]" "exists_descend_data(最大層の点と f₀、第 416)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.exists_descend_data") 3,
    .citation "[mathlib]" "Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly(埋め込みの値域)"
      (.inMathlib "Algebra.IsAlgebraic.range_eval_eq_rootSet_minpoly") 3,
    .implicitStep
      ("★★★★★★★測定(2026-08-29): 原文 p.5 の帰納段(正規化 → f₀ を取る → " ++
       "f₀(β) ∉ S′ → S を S′ に、β を f₀(β) に置き換える)を**1 本の定理**にした。" ++
       "★組み立てに要ったのは 7 つ: 正規化(第 409)、測度の Möbius 不変性(§9-983 と本ファイルの " ++
       "dSum_mobius)、Gal-安定性の Möbius 不変性(§9-985)、最大層の点(第 416)、" ++
       "分離(§9-984)、**多項式の像での Gal-安定性**(本ファイルの conjSet_aeval_subset)、" ++
       "測度の降下(第 408)") 8,
    .implicitStep
      ("★★★これで Lemma 2.4 に残るのは**配管 1 つだけ**である" ++
       "——合成を ℙ¹ の有理写像として組み立てること" ++
       "(原文の f(x) ∈ ℚ(x) は Möbius と多項式の交互合成。Separation.lean の P1C が受け皿)。" ++
       "★β を運ぶ帰納の骨格は §9-986 に、降下段は本ファイルにある") 7 ]

end ABC3.Found.NCBelyi
