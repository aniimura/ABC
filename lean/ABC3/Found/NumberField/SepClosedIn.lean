/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.FieldTheory.SeparableClosure
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# 分離閉な部分体の上では分離既約多項式は既約のまま(鎖 `sec6items` の `thm62-i-ext`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★言明

`K₁ ⊆ K₂` で `K₁` が `K₂` の中で**分離閉**(`separableClosure K₁ K₂ = ⊥`)なら、
`K₁` 上の**分離的で既約**な多項式は `K₂` 上でも既約である。
★系として `[L₁·K₂ : K₂] = [L₁ : K₁]`(原文 `Theorem 6.2, (i)` の仮定 (b)(c) の帰結)。

## ★★★中身は 4 段

| 段 | 中身 |
|---|---|
| 1 | 分離的な多項式の根は分離的(`isSeparable_of_aeval_eq_zero`) |
| 2 | `K₂[X]` のモニック因子の係数は代数閉包で分離的 —— **Vieta** |
| 3 | 分離閉性から係数は `K₁` に落ちる(`exists_map_of_monic_dvd`) |
| 4 | `minpoly K₂ α` を持ち上げて `minpoly K₁ α = f` と次数を比べる |
-/

namespace ABC3.Found.NF

open Polynomial

variable {K₁ K₂ : Type*} [Field K₁] [Field K₂] [Algebra K₁ K₂]

/-! ## ★1. 分離的な多項式の根は分離的 -/

/-- ★**分離的な多項式の根は分離的**。 -/
theorem isSeparable_of_aeval_eq_zero {Ω : Type*} [Field Ω] [Algebra K₁ Ω] {f : K₁[X]}
    (hsep : f.Separable) {r : Ω} (hr : aeval r f = 0) : IsSeparable K₁ r := by
  haveI : IsIntegral K₁ r := IsAlgebraic.isIntegral ⟨f, hsep.ne_zero, hr⟩
  exact hsep.of_dvd (minpoly.dvd K₁ r hr)

/-! ## ★2. モニック因子は `K₁` から来る -/

/-- ★★★**分離閉なら、`f` のモニック因子は `K₁[X]` から来る**。

★中身は Vieta —— 因子の根は `f` の根なので `K₁` 上分離的、
したがって係数は `separableClosure K₁ Ω` に入る。
`K₁` が `K₂` の中で分離閉なので、`K₂` の元でもある係数は `K₁` に落ちる。 -/
theorem exists_map_of_monic_dvd (hsc : separableClosure K₁ K₂ = ⊥)
    {f : K₁[X]} (hsep : f.Separable) {g : K₂[X]} (hg : g.Monic)
    (hdvd : g ∣ f.map (algebraMap K₁ K₂)) :
    ∃ h : K₁[X], h.Monic ∧ h.map (algebraMap K₁ K₂) = g := by
  classical
  set Ω := AlgebraicClosure K₂
  set E : IntermediateField K₁ Ω := separableClosure K₁ Ω with hE
  -- `g` の係数は `E` に入る
  have hcoeff : ∀ n, algebraMap K₂ Ω (g.coeff n) ∈ E := by
    have hgΩ : (g.map (algebraMap K₂ Ω)).Monic := hg.map _
    have hsplit : (g.map (algebraMap K₂ Ω)).Splits := IsAlgClosed.splits _
    have hprod := hsplit.eq_prod_roots_of_monic hgΩ
    have hroot : ∀ r ∈ (g.map (algebraMap K₂ Ω)).roots, r ∈ E := by
      intro r hr
      have hr0 : (g.map (algebraMap K₂ Ω)).eval r = 0 := (mem_roots'.mp hr).2
      obtain ⟨c, hc⟩ := hdvd
      have hfr : aeval r f = 0 := by
        have : (f.map (algebraMap K₁ Ω)).eval r = 0 := by
          rw [show (algebraMap K₁ Ω) = (algebraMap K₂ Ω).comp (algebraMap K₁ K₂) from
              IsScalarTower.algebraMap_eq K₁ K₂ Ω,
            ← map_map, hc, Polynomial.map_mul, eval_mul, hr0, zero_mul]
        simpa [aeval_def, eval_map] using this
      exact isSeparable_of_aeval_eq_zero hsep hfr
    -- 根がすべて `E` にあるので、積の係数も `E` に入る
    have hlift : g.map (algebraMap K₂ Ω) ∈ lifts (algebraMap E Ω) := by
      rw [hprod]
      refine Subsemiring.multiset_prod_mem _ _ (fun p hp => ?_)
      obtain ⟨r, hr, rfl⟩ := Multiset.mem_map.mp hp
      exact ⟨X - C ⟨r, hroot r hr⟩, by simp⟩
    intro n
    obtain ⟨y, hy⟩ := (lifts_iff_coeff_lifts _).mp hlift n
    have : algebraMap K₂ Ω (g.coeff n) = algebraMap E Ω y := by
      rw [hy, coeff_map]
    rw [this]
    exact y.2
  -- 分離閉性から `K₁` に落ちる
  have hK₁ : ∀ n, g.coeff n ∈ Set.range (algebraMap K₁ K₂) := by
    intro n
    have hsep' : IsSeparable K₁ (g.coeff n) := by
      have h1 : IsSeparable K₁ (algebraMap K₂ Ω (g.coeff n)) := hcoeff n
      rwa [IsSeparable, minpoly.algebraMap_eq
        (FaithfulSMul.algebraMap_injective K₂ Ω) (g.coeff n)] at h1
    have : g.coeff n ∈ separableClosure K₁ K₂ := hsep'
    rw [hsc] at this
    exact IntermediateField.mem_bot.mp this
  have hgl : g ∈ lifts (algebraMap K₁ K₂) := (lifts_iff_coeff_lifts _).mpr hK₁
  obtain ⟨h, hmap, -, hmon⟩ := lifts_and_degree_eq_and_monic hgl hg
  exact ⟨h, hmon, hmap⟩

/-! ## ★3. 既約性は保たれる -/

/-- ★★★★★**分離閉な部分体の上では分離既約多項式は既約のまま**。

★`α` を代数閉包の根として `minpoly K₂ α` を `K₁[X]` へ持ち上げると
`minpoly K₁ α = f` が割るので、次数が一致する。 -/
theorem irreducible_map_of_separableClosure_eq_bot (hsc : separableClosure K₁ K₂ = ⊥)
    {f : K₁[X]} (hmon : f.Monic) (hf : Irreducible f) (hsep : f.Separable) :
    Irreducible (f.map (algebraMap K₁ K₂)) := by
  classical
  set Ω := AlgebraicClosure K₂
  have hdeg0 : 0 < f.degree := hf.degree_pos
  obtain ⟨α, hα⟩ : ∃ α : Ω, aeval α f = 0 :=
    IsAlgClosed.exists_aeval_eq_zero Ω f (by exact_mod_cast hdeg0.ne')
  have hfm : f = minpoly K₁ α := minpoly.eq_of_irreducible_of_monic hf hα hmon
  haveI : IsIntegral K₁ α := ⟨f, hmon, by simpa [aeval_def] using hα⟩
  haveI : IsIntegral K₂ α := IsIntegral.tower_top this
  set m : K₂[X] := minpoly K₂ α with hm
  have hmmon : m.Monic := minpoly.monic ‹IsIntegral K₂ α›
  have hαm : aeval α (f.map (algebraMap K₁ K₂)) = 0 := by
    rw [aeval_map_algebraMap]
    exact hα
  have hdvd : m ∣ f.map (algebraMap K₁ K₂) := minpoly.dvd K₂ α hαm
  obtain ⟨h, hmonh, hmap⟩ := exists_map_of_monic_dvd hsc hsep hmmon hdvd
  have hαh : aeval α h = 0 := by
    have : aeval α (h.map (algebraMap K₁ K₂)) = 0 := by rw [hmap]; exact minpoly.aeval K₂ α
    rwa [aeval_map_algebraMap] at this
  have hfh : f ∣ h := hfm ▸ minpoly.dvd K₁ α hαh
  have hdh : f.degree ≤ h.degree := Polynomial.degree_le_of_dvd hfh hmonh.ne_zero
  have hdm : h.degree = m.degree := by
    rw [← hmap, Polynomial.degree_map_eq_of_injective
      (FaithfulSMul.algebraMap_injective K₁ K₂)]
  have hdfm : (f.map (algebraMap K₁ K₂)).degree = f.degree :=
    Polynomial.degree_map_eq_of_injective (FaithfulSMul.algebraMap_injective K₁ K₂) f
  have hle : m.degree ≤ f.degree := by
    rw [← hdfm]
    exact Polynomial.degree_le_of_dvd hdvd (hmon.map _).ne_zero
  have heq : m.degree = (f.map (algebraMap K₁ K₂)).degree := by
    rw [hdfm]
    exact le_antisymm hle (hdm ▸ hdh)
  obtain ⟨c, hc⟩ := hdvd
  have hfd : (f.map (algebraMap K₁ K₂)).degree ≠ ⊥ := by
    simpa using (hmon.map (algebraMap K₁ K₂)).ne_zero
  have hcdeg : c.degree = 0 := by
    have hmul := Polynomial.degree_mul (p := m) (q := c)
    rw [← hc, heq] at hmul
    cases ha : (f.map (algebraMap K₁ K₂)).degree with
    | bot => exact absurd ha hfd
    | coe n =>
      cases hb : c.degree with
      | bot =>
        exfalso
        have hc0 : c = 0 := Polynomial.degree_eq_bot.mp hb
        rw [hc0, mul_zero] at hc
        exact (hmon.map (algebraMap K₁ K₂)).ne_zero hc
      | coe k =>
        rw [ha, hb] at hmul
        have hnk : n = n + k := by exact_mod_cast hmul
        have hk0 : k = 0 := by omega
        simp [hk0]
  have hcu : c = 1 := by
    have hcm : c.Monic := by
      have h1 : (f.map (algebraMap K₁ K₂)).leadingCoeff = 1 := (hmon.map _).leadingCoeff
      rw [hc, Polynomial.leadingCoeff_mul, hmmon.leadingCoeff, one_mul] at h1
      exact h1
    exact hcm.eq_one_of_isUnit (Polynomial.isUnit_iff_degree_eq_zero.mpr hcdeg)
  rw [hc, hcu, mul_one]
  exact minpoly.irreducible ‹IsIntegral K₂ α›

/-! ## ★4. 最小多項式は変わらない —— 次数の保存 -/

/-- ★★★★★**分離閉なら最小多項式は変わらない**。 -/
theorem minpoly_eq_map_of_separableClosure_eq_bot (hsc : separableClosure K₁ K₂ = ⊥)
    {Ω : Type*} [Field Ω] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
    {α : Ω} (hint : IsIntegral K₁ α) (hs : IsSeparable K₁ α) :
    minpoly K₂ α = (minpoly K₁ α).map (algebraMap K₁ K₂) := by
  refine (minpoly.eq_of_irreducible_of_monic ?_ ?_ ?_).symm
  · exact irreducible_map_of_separableClosure_eq_bot hsc (minpoly.monic hint)
      (minpoly.irreducible hint) hs
  · rw [aeval_map_algebraMap]
    exact minpoly.aeval K₁ α
  · exact (minpoly.monic hint).map _

open IntermediateField in
/-- ★★★★★★**`[K₂(α) : K₂] = [K₁(α) : K₁]`** —— 原文 `Theorem 6.2, (i)` の
仮定 (b)(c) から `[L₂:K₂] = [L₁:K₁]` を出す段。 -/
theorem finrank_adjoin_eq_of_separableClosure_eq_bot (hsc : separableClosure K₁ K₂ = ⊥)
    {Ω : Type*} [Field Ω] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
    {α : Ω} (hint : IsIntegral K₁ α) (hs : IsSeparable K₁ α) :
    Module.finrank K₂ K₂⟮α⟯ = Module.finrank K₁ K₁⟮α⟯ := by
  rw [IntermediateField.adjoin.finrank (hint.tower_top), IntermediateField.adjoin.finrank hint,
    minpoly_eq_map_of_separableClosure_eq_bot hsc hint hs,
    Polynomial.natDegree_map_eq_of_injective (FaithfulSMul.algebraMap_injective K₁ K₂)]

open IntermediateField in
/-- ★★★★★★**原文 `Theorem 6.2, (i)` の言明** —— `L₂ := L₁·K₂` は `[L₂:K₂] = [L₁:K₁]`。

★原始元定理で `L₁ = K₁⟮α⟯` と書き、`finrank_adjoin_eq_of_separableClosure_eq_bot` を当てる。 -/
theorem finrank_adjoin_coe_eq_of_separableClosure_eq_bot
    (hsc : separableClosure K₁ K₂ = ⊥)
    {Ω : Type*} [Field Ω] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
    (L : IntermediateField K₁ Ω) [FiniteDimensional K₁ L] [Algebra.IsSeparable K₁ L] :
    Module.finrank K₂ (IntermediateField.adjoin K₂ (L : Set Ω)) = Module.finrank K₁ L := by
  obtain ⟨β, hβ⟩ := Field.exists_primitive_element K₁ L
  have hLeq : L = K₁⟮(β : Ω)⟯ := by
    have h := congrArg IntermediateField.lift hβ
    rwa [IntermediateField.lift_adjoin_simple, IntermediateField.lift_top, eq_comm] at h
  have hmem : (β : Ω) ∈ IntermediateField.adjoin K₂ (L : Set Ω) :=
    IntermediateField.subset_adjoin K₂ _ β.2
  have hAeq : IntermediateField.adjoin K₂ (L : Set Ω) = K₂⟮(β : Ω)⟯ := by
    refine le_antisymm (IntermediateField.adjoin_le_iff.mpr ?_)
      (IntermediateField.adjoin_simple_le_iff.mpr hmem)
    intro x hx
    have hx' : x ∈ K₁⟮(β : Ω)⟯ := hLeq ▸ hx
    have hle : K₁⟮(β : Ω)⟯ ≤ (K₂⟮(β : Ω)⟯).restrictScalars K₁ :=
      IntermediateField.adjoin_simple_le_iff.mpr
        (IntermediateField.mem_adjoin_simple_self K₂ _)
    exact hle hx'
  have hint : IsIntegral K₁ (β : Ω) :=
    (IsIntegral.of_finite K₁ β).map (IntermediateField.val L)
  have hsep : IsSeparable K₁ (β : Ω) := by
    have h2 : minpoly K₁ ((IntermediateField.val L) β) = minpoly K₁ β :=
      minpoly.algebraMap_eq (IntermediateField.val L).injective β
    show (minpoly K₁ ((IntermediateField.val L) β)).Separable
    rw [h2]
    exact Algebra.IsSeparable.isSeparable K₁ β
  rw [hAeq, finrank_adjoin_eq_of_separableClosure_eq_bot hsc hint hsep]
  exact congrArg (fun M : IntermediateField K₁ Ω => Module.finrank K₁ M) hLeq.symm

/-! ## ★5. 射の側の材料(`thm62-i-Dfun` が使う) -/

/-- ★★★**`K₁` 上の共役は `K₂` 上でも共役** —— `minpoly K₂ β` の根は
`minpoly K₁ β` の根と同じ。

★★これが `thm62-i-Dfun`(関手 `𝒟₁ → 𝒟₂` の**射の側**)の要である ——
`K₁`-代数射 `σ : L → M` に対し `σβ` は `minpoly K₁ β` の根なので、
`minpoly K₂ β` の根でもあり、したがって `K₂⟮β⟯ → K₂⟮σβ⟯` が定まる。 -/
theorem aeval_minpoly_extend (hsc : separableClosure K₁ K₂ = ⊥)
    {Ω : Type*} [Field Ω] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
    {β : Ω} (hint : IsIntegral K₁ β) (hs : IsSeparable K₁ β)
    {γ : Ω} (hγ : aeval γ (minpoly K₁ β) = 0) :
    aeval γ (minpoly K₂ β) = 0 := by
  rw [minpoly_eq_map_of_separableClosure_eq_bot hsc hint hs, aeval_map_algebraMap]
  exact hγ

/-- ★逆向き —— `minpoly K₂ β` の根は `minpoly K₁ β` の根。 -/
theorem aeval_minpoly_descend (hsc : separableClosure K₁ K₂ = ⊥)
    {Ω : Type*} [Field Ω] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
    {β : Ω} (hint : IsIntegral K₁ β) (hs : IsSeparable K₁ β)
    {γ : Ω} (hγ : aeval γ (minpoly K₂ β) = 0) :
    aeval γ (minpoly K₁ β) = 0 := by
  rw [minpoly_eq_map_of_separableClosure_eq_bot hsc hint hs, aeval_map_algebraMap] at hγ
  exact hγ

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.2, (i)` の「`K₁` が `K₂` の中で分離閉 ⟹ 次数が保たれる」。 -/
def irreducible_map_of_separableClosure_eq_bot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 分離閉な部分体の上では分離既約多項式は既約のまま",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.NF
