/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.MobiusRedDeg
import ABC3.Found.NCBelyi.ConjStable
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★`ℚ`-Möbius は `Gal`-安定性を保つ（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> applying an automorphism (with rational coeﬃcients!) as in Lemma 2.3 and then

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— 配管の (1)

`Found/NCBelyi/SeparationStep.lean` の測定で、`Lemma 2.4` の (b) に残るのは
**配管 2 つ**だけになった:

1. ★**`ℚ`-Möbius が `Gal`-安定性を保つこと** ← ★★**本ファイル**
2. 合成を `ℙ¹` の有理写像として組み立てること

★測度の側（`redDeg`／`maxRedDeg`）は `MobiusRedDeg.lean` で済んでいる。
★★本ファイルは**共役の側**を取る:

    **`conjSet (σ x) ⊆ σ (conjSet x)`**   （`σ x ≔ ν/(x−λ) + μ`、`λ, ν, μ ∈ ℚ`、`ν ≠ 0`）

★★★これで `S` が `Gal`-安定なら `σ(S)` も `Gal`-安定である
——`separation_step` の仮定がそろう。

## ★機構 —— 変換多項式

`p ≔ minpoly ℚ x`、`d ≔ deg p` として

    `q(Y) ≔ ∑_{k≤d} p_k · (ν + λ(Y−μ))^k · (Y−μ)^{d−k}`   （`mobiusTransform`）

と置くと、`y ≠ μ` で

    **`q(y) = (y−μ)^d · p(ν/(y−μ) + λ)`**

である（`aeval_mobiusTransform`）。★したがって

* `q(σx) = (σx−μ)^d · p(x) = 0` ——`ν/(σx−μ) + λ = x` だから
* `q(μ) = p_d · ν^d = ν^d ≠ 0` ——`p` はモニックだから（`eval_mobiusTransform_at_mu`）

★★前者から `minpoly ℚ (σx) ∣ q`、後者から `q ≠ 0` かつ **`μ` は `q` の根でない**。
★★★ゆえに `z ∈ conjSet (σx)` なら `z ≠ μ` かつ `p(ν/(z−μ)+λ) = 0`
——すなわち `σ⁻¹ z ∈ conjSet x` であり、`z = σ(σ⁻¹ z)` である。

## ★★これで `Lemma 2.4` に残るのは配管 1 つ

★**合成を `ℙ¹` の有理写像として組み立てること**（原文の `f(x) ∈ ℚ(x)` は
Möbius と多項式の交互合成。`Separation.lean` の `P1C` が受け皿）。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★★★変換多項式 -/

/-- ★★★★★★★★★★**Möbius で引き戻した多項式**。

    `q(Y) ≔ ∑_{k≤d} p_k · (ν + λ(Y−μ))^k · (Y−μ)^{d−k}`

★`y ≠ μ` で `q(y) = (y−μ)^d · p(ν/(y−μ) + λ)` になる（分母を払った形）。 -/
noncomputable def mobiusTransform (p : ℚ[X]) (lam nu mu : ℚ) : ℚ[X] :=
  ∑ k ∈ Finset.range (p.natDegree + 1),
    C (p.coeff k) * (C nu + C lam * (X - C mu)) ^ k * (X - C mu) ^ (p.natDegree - k)

/-- ★★★★★★**`q(y) = (y−μ)^d · p(ν/(y−μ) + λ)`**（`y ≠ μ`）。 -/
theorem aeval_mobiusTransform (p : ℚ[X]) (lam nu mu : ℚ) (y : ℂ) (hy : y ≠ (mu : ℂ)) :
    aeval y (mobiusTransform p lam nu mu)
      = (y - (mu : ℂ)) ^ p.natDegree * aeval ((nu : ℂ) / (y - (mu : ℂ)) + (lam : ℂ)) p := by
  have ht0 : y - (mu : ℂ) ≠ 0 := sub_ne_zero.mpr hy
  rw [mobiusTransform, map_sum]
  rw [Polynomial.aeval_eq_sum_range (p := p) (x := (nu : ℂ) / (y - (mu : ℂ)) + (lam : ℂ))]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun k hk => ?_)
  have hkd : k ≤ p.natDegree := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
  have hsplit : (y - (mu : ℂ)) ^ p.natDegree
      = (y - (mu : ℂ)) ^ k * (y - (mu : ℂ)) ^ (p.natDegree - k) := by
    rw [← pow_add]; congr 1; omega
  have hz : ((nu : ℂ) / (y - (mu : ℂ)) + (lam : ℂ))
      = ((nu : ℂ) + (lam : ℂ) * (y - (mu : ℂ))) / (y - (mu : ℂ)) := by field_simp
  simp only [map_mul, map_pow, map_add, map_sub, aeval_C, aeval_X, Algebra.smul_def,
    eq_ratCast]
  rw [hz, div_pow, hsplit]
  field_simp

/-- ★★★★★**`q(μ) = p_d · ν^d`** —— `p` がモニックなら `ν^d ≠ 0`。 -/
theorem eval_mobiusTransform_at_mu (p : ℚ[X]) (lam nu mu : ℚ) :
    (mobiusTransform p lam nu mu).eval mu = p.coeff p.natDegree * nu ^ p.natDegree := by
  rw [mobiusTransform, Polynomial.eval_finsetSum]
  rw [Finset.sum_eq_single p.natDegree]
  · simp
  · intro k hk hne
    have hkd : k < p.natDegree := by
      have := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
      omega
    have hne0 : p.natDegree - k ≠ 0 := by omega
    simp [hne0]
  · intro h
    exact absurd (Finset.self_mem_range_succ p.natDegree) h

/-! ## ★★★★★★★★★★★★★★★★★★★★★★共役は Möbius で写る -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**`conjSet (σ x) ⊆ σ (conjSet x)`**。

原文 (NCBelyi p.5):
> applying an automorphism (with rational coeﬃcients!) as in Lemma 2.3 and then

★★これが原文の『**with rational coefficients!**』のもう一つの帰結である
——`Gal`-安定性が正規化で壊れない。
★（測度が壊れないことは `MobiusRedDeg.lean` で取った。） -/
theorem conjSet_mobius_subset (lam nu mu : ℚ) (hnu : nu ≠ 0) (x : ℂ)
    (hint : IsIntegral ℚ x) (hx : x ≠ (lam : ℂ)) :
    conjSet ((nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ))
      ⊆ (conjSet x).image (fun w => (nu : ℂ) * ((w - (lam : ℂ))⁻¹) + (mu : ℂ)) := by
  classical
  have hmonic : (minpoly ℚ x).Monic := minpoly.monic hint
  have hnuC : ((nu : ℂ)) ≠ 0 := by exact_mod_cast hnu
  have hxl : x - (lam : ℂ) ≠ 0 := sub_ne_zero.mpr hx
  set q : ℚ[X] := mobiusTransform (minpoly ℚ x) lam nu mu with hq
  set y : ℂ := (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ) with hy
  have hymu : y - (mu : ℂ) = (nu : ℂ) * ((x - (lam : ℂ))⁻¹) := by rw [hy]; ring
  have hyne : y ≠ (mu : ℂ) := by
    intro h
    rw [sub_eq_zero.mpr h] at hymu
    exact (mul_ne_zero hnuC (inv_ne_zero hxl)) hymu.symm
  have hqmu : q.eval mu = nu ^ (minpoly ℚ x).natDegree := by
    rw [hq, eval_mobiusTransform_at_mu, hmonic.coeff_natDegree, one_mul]
  have hq0 : q ≠ 0 := by
    intro hc
    rw [hc] at hqmu
    simp only [Polynomial.eval_zero] at hqmu
    exact (pow_ne_zero _ hnu) hqmu.symm
  have hqy : aeval y q = 0 := by
    rw [hq, aeval_mobiusTransform (minpoly ℚ x) lam nu mu y hyne, hymu]
    have harg : (nu : ℂ) / ((nu : ℂ) * ((x - (lam : ℂ))⁻¹)) + (lam : ℂ) = x := by
      field_simp
      ring
    rw [harg, show (aeval x) (minpoly ℚ x) = 0 from minpoly.aeval ℚ x, mul_zero]
  have hyint : IsIntegral ℚ y := isIntegral_of_aeval_eq_zero hq0 hqy
  have hdvd : minpoly ℚ y ∣ q := minpoly.dvd ℚ y hqy
  intro z hz
  have hzq : aeval z q = 0 := by
    obtain ⟨r, hr⟩ := hdvd
    rw [hr, map_mul, (mem_conjSet_iff hyint).1 hz, zero_mul]
  have hzmu : z ≠ (mu : ℂ) := by
    intro hc
    rw [hc] at hzq
    have h1 : ((mu : ℂ)) = algebraMap ℚ ℂ mu := by simp
    rw [h1, Polynomial.aeval_algebraMap_apply] at hzq
    have h2 : q.eval mu = 0 := by simpa using hzq
    rw [hqmu] at h2
    exact (pow_ne_zero _ hnu) h2
  have hzmu0 : z - (mu : ℂ) ≠ 0 := sub_ne_zero.mpr hzmu
  have hkey := aeval_mobiusTransform (minpoly ℚ x) lam nu mu z hzmu
  rw [← hq, hzq] at hkey
  have hw : aeval ((nu : ℂ) / (z - (mu : ℂ)) + (lam : ℂ)) (minpoly ℚ x) = 0 := by
    rcases mul_eq_zero.1 hkey.symm with h | h
    · exact absurd h (pow_ne_zero _ hzmu0)
    · exact h
  refine Finset.mem_image.2 ⟨(nu : ℂ) / (z - (mu : ℂ)) + (lam : ℂ),
    (mem_conjSet_iff hint).2 hw, ?_⟩
  rw [add_sub_cancel_right]
  field_simp
  ring

/-- ★★★★★★★★★★★★★★★★★★**`ℚ`-Möbius は `Gal`-安定性を保つ**。

★これで `exists_rat_normalization`（`|α| ≤ 1`、`|β| ≥ 4`）で正規化したあとも
`separation_step` の仮定（`IsConjStable`）がそろう。 -/
theorem isConjStable_mobius (lam nu mu : ℚ) (hnu : nu ≠ 0) (S : Finset ℂ)
    (hint : ∀ x ∈ S, IsIntegral ℚ x) (hlam : ∀ x ∈ S, x ≠ (lam : ℂ))
    (hstab : IsConjStable S) :
    IsConjStable (S.image (fun x => (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ))) := by
  classical
  intro y hy z hz
  obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
  have hsub := conjSet_mobius_subset lam nu mu hnu x (hint x hxS) (hlam x hxS) hz
  obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hsub
  exact Finset.mem_image.2 ⟨w, hstab x hxS hw, rfl⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def mobiusTransform.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(Möbius で引き戻した多項式)",
    sectionId := "ncbelyi-lemma-2-4" }

def conjSet_mobius_subset.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(conjSet (σ x) ⊆ σ (conjSet x))",
    sectionId := "ncbelyi-lemma-2-4" }

def isConjStable_mobius.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ-Möbius は Gal-安定性を保つ)",
    sectionId := "ncbelyi-lemma-2-4" }

def isConjStable_mobius.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "mem_conjSet_iff(共役 ⟺ 最小多項式の根、第 410)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.mem_conjSet_iff") 2,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 原文の『with rational coefficients!』の帰結は 2 つある。" ++
       "(1) 測度(redDeg／maxRedDeg)が壊れない——MobiusRedDeg.lean。" ++
       "(2) ★**Gal-安定性が壊れない**——本ファイル。" ++
       "★機構は変換多項式 q(Y) ≔ ∑ p_k (ν+λ(Y−μ))^k (Y−μ)^{d−k} であり、" ++
       "q(σx) = 0(⟹ minpoly(σx) ∣ q)と q(μ) = ν^d ≠ 0(⟹ μ は根でない)から出る") 6,
    .implicitStep
      ("★★★これで Lemma 2.4 に残るのは**配管 1 つ**である" ++
       "——合成を ℙ¹ の有理写像として組み立てること" ++
       "(原文の f(x) ∈ ℚ(x) は Möbius と多項式の交互合成。Separation.lean の P1C が受け皿)。" ++
       "★数値の核(SeparationStep.lean)・測度(MobiusRedDeg.lean)・" ++
       "Gal-安定性(本ファイル)はすべて取れた") 6 ]

end ABC3.Found.NCBelyi
