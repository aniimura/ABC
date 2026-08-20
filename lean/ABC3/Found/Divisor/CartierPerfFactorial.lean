/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierMonoprime
import ABC3.Found.FrdI.Def24

/-!
# `Φ(L)` は perf-factorial

★★★★★[FrdI] `Example 6.1` の

> one verifies immediately that Φ(L) is perf-factorial

を閉じる。`Φ := Γ ∩ ℤ≥0[S]`(`Γ` は Cartier 因子の群、`Q`-Cartier)について
`IsPerfFactorialWith (effSub Γ) (iotaAt hQ)` を組む。

## ★鍵は 1 本の定理

族 `ι p` を「`Φ^pf → ℚ[S]`(`CartierPf.lean` の `pfCoeff`)の `p` 成分」とすると、

**`factorMap_iotaAt` —— `factorMap ι a p = ι p a`**

すなわち★**原文の `sup(Bound^p_{0}(a))` は、`a` の `p` 成分そのものになる**。
これが言えると (c1)(c2)(c3) がすべて `pfCoeff` の加法性・単射性に落ちる。

証明は 2 方向:
* `≤` —— `Bound` の元は `a` 以下だから係数も以下(`ratCoeffAt_mono`)。
* `≥` —— ★**`a` の `p` 成分だけを取り出した元が実際に `Φ^pf` にある**
  (`exists_pPart`。`pfCoeffHom_range` が `ℚ≥0[S]` 全体を像に持つことによる)、
  しかもそれは `pCarrierPf p` に入る(`mem_pCarrierPf_of_pfCoeff_single`。
  ここで `Q`-Cartier の生成元 `d` を使い、`n = c.den · d` 倍すると
  `Φ` の元 `d·c.num·s` になる)。

## ★`realScale` は空虚

`M_p ≃+ ℕ`(`CartierMonoprime.lean`)であり、**`ℕ` と `ℝ≥0` は同型でない**
(`ℝ≥0` は可除、`ℕ` は違う)ので、`ℝ`-monoprime の仮定は偽である。
-/

namespace ABC3.Found.FrdI

open Finsupp

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}

/-! ## ★1. 素点ごとの係数と `ι` -/

/-- ★`Φ^pf` の元の、素点 `p` での係数(有理数)。 -/
noncomputable def ratCoeffAt (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (x : Pf (effSub Γ)) : ℚ :=
  pfCoeff x (effSubPrimeEquiv hQ p)

theorem ratCoeffAt_nonneg (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (x : Pf (effSub Γ)) : 0 ≤ ratCoeffAt hQ p x :=
  Finsupp.le_def.mp (pfCoeff_nonneg x) _

theorem ratCoeffAt_add (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (x y : Pf (effSub Γ)) :
    ratCoeffAt hQ p (x + y) = ratCoeffAt hQ p x + ratCoeffAt hQ p y := by
  show (pfCoeffHom (x + y)) _ = _
  rw [map_add]
  rfl

/-- ★★**`ι p`** —— 原文の `M^rlf_p ≃ ℝ≥0` の同一視。 -/
noncomputable def iotaAt (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (x : Pf (effSub Γ)) : NNReal :=
  Real.toNNReal ((ratCoeffAt hQ p x : ℚ) : ℝ)

theorem iotaAt_add (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) (x y : Pf (effSub Γ)) :
    iotaAt hQ p (x + y) = iotaAt hQ p x + iotaAt hQ p y := by
  have hx : (0:ℝ) ≤ ((ratCoeffAt hQ p x : ℚ) : ℝ) := by
    exact_mod_cast ratCoeffAt_nonneg hQ p x
  have hy : (0:ℝ) ≤ ((ratCoeffAt hQ p y : ℚ) : ℝ) := by
    exact_mod_cast ratCoeffAt_nonneg hQ p y
  rw [iotaAt, iotaAt, iotaAt, ratCoeffAt_add]
  push_cast
  exact Real.toNNReal_add hx hy

theorem ratCoeffAt_of_iotaAt_eq (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    {x y : Pf (effSub Γ)} (h : iotaAt hQ p x = iotaAt hQ p y) :
    ratCoeffAt hQ p x = ratCoeffAt hQ p y := by
  have hx : (0:ℝ) ≤ ((ratCoeffAt hQ p x : ℚ) : ℝ) := by
    exact_mod_cast ratCoeffAt_nonneg hQ p x
  have hy : (0:ℝ) ≤ ((ratCoeffAt hQ p y : ℚ) : ℝ) := by
    exact_mod_cast ratCoeffAt_nonneg hQ p y
  have h2 := congrArg (fun t : NNReal => (t : ℝ)) h
  simp only [iotaAt] at h2
  rw [Real.coe_toNNReal _ hx, Real.coe_toNNReal _ hy] at h2
  exact_mod_cast h2

theorem iotaAt_le_of_ratCoeffAt_le (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    {x y : Pf (effSub Γ)} (h : ratCoeffAt hQ p x ≤ ratCoeffAt hQ p y) :
    iotaAt hQ p x ≤ iotaAt hQ p y := by
  rw [iotaAt, iotaAt]
  exact Real.toNNReal_mono (by exact_mod_cast h)

/-- ★`MLe` は係数を増やす。 -/
theorem ratCoeffAt_mono (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    {x y : Pf (effSub Γ)} (h : MLe x y) : ratCoeffAt hQ p x ≤ ratCoeffAt hQ p y := by
  obtain ⟨c, rfl⟩ := h
  rw [ratCoeffAt_add]
  exact le_add_of_nonneg_right (ratCoeffAt_nonneg hQ p c)

theorem iotaAt_eq_zero_iff (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (x : Pf (effSub Γ)) : iotaAt hQ p x = 0 ↔ ratCoeffAt hQ p x = 0 := by
  rw [iotaAt, Real.toNNReal_eq_zero]
  have hnn := ratCoeffAt_nonneg hQ p x
  constructor
  · intro h
    have hz : ((ratCoeffAt hQ p x : ℚ) : ℝ) = 0 := le_antisymm h (by exact_mod_cast hnn)
    exact_mod_cast hz
  · intro h
    rw [h]
    simp

/-! ## ★2. `pCarrierPf p` の同定 -/

theorem pfCoeff_of (b : effSub Γ) :
    pfCoeffHom (Pf.of b) = intToRatFs ((b : effSub Γ) : S →₀ ℤ) := by
  show pfCoeff (Pf.mk b 1) = _
  rw [pfCoeff_mk]
  simp

/-- ★★`pCarrierPf p` の元は、係数が `s` にしか乗らない。 -/
theorem support_pfCoeff_of_mem_pCarrierPf (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    {x : Pf (effSub Γ)} (hx : x ∈ pCarrierPf (effSub Γ) p) :
    (pfCoeffHom x : S →₀ ℚ).support ⊆ {effSubPrimeEquiv hQ p} := by
  classical
  rcases hx with ⟨n, b, hb, hnx⟩ | hx0
  · have hbsupp : ((b : effSub Γ) : S →₀ ℤ).support = {effSubPrimeEquiv hQ p} :=
      (mem_primeCarrier_iff hQ p b).mp hb
    have h1 : pfCoeffHom (((n : ℕ+) : ℕ) • x) = pfCoeffHom (Pf.of b) := by rw [hnx]
    rw [map_nsmul, pfCoeff_of] at h1
    intro t ht
    have hne : (pfCoeffHom x : S →₀ ℚ) t ≠ 0 := Finsupp.mem_support_iff.mp ht
    have hnz : ((n : ℕ) : ℚ) ≠ 0 := by
      have := n.2
      positivity
    have h2 := congrArg (fun v : S →₀ ℚ => v t) h1
    simp only [Finsupp.smul_apply, nsmul_eq_mul, intToRatFs_apply] at h2
    have h3 : ((((b : effSub Γ) : S →₀ ℤ) t : ℤ) : ℚ) ≠ 0 := by
      rw [← h2]
      exact mul_ne_zero hnz hne
    have h4 : ((b : effSub Γ) : S →₀ ℤ) t ≠ 0 := by exact_mod_cast h3
    rw [← hbsupp]
    exact Finsupp.mem_support_iff.mpr h4
  · rw [Set.mem_singleton_iff] at hx0
    rw [hx0, map_zero]
    simp

/-- ★★係数が `s` にしか乗らない `Φ^pf` の元は `pCarrierPf p` に入る。

★ここで `Q`-Cartier の生成元 `d` を使う —— `c = num/den` に対し
`n = den·d` 倍すると `Φ` の元 `d·num·s` になる。 -/
theorem mem_pCarrierPf_of_pfCoeff_single (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    {x : Pf (effSub Γ)} {c : ℚ} (hc : 0 ≤ c)
    (hx : (pfCoeffHom x : S →₀ ℚ) = single (effSubPrimeEquiv hQ p) c) :
    x ∈ pCarrierPf (effSub Γ) p := by
  classical
  set s := effSubPrimeEquiv hQ p with hs
  rcases eq_or_lt_of_le hc with hc0 | hcpos
  · refine Or.inr ?_
    refine Set.mem_singleton_iff.mpr (pfCoeffHom_injective ?_)
    rw [hx, map_zero, ← hc0]
    simp
  refine Or.inl ?_
  obtain ⟨d, hd, hdiv⟩ := exists_pos_gen_sGrpAt hQ s
  have hnum : 0 < c.num := Rat.num_pos.mpr hcpos
  have hdcnum : 0 < d * c.num := mul_pos hd hnum
  have hbmem : (single s (d * c.num) : S →₀ ℤ) ∈ effSub Γ := by
    refine mem_effSub.mpr ⟨(hdiv _).mpr ⟨c.num, rfl⟩, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simp only [Finsupp.single_eq_same]
      omega
    · simp [hst]
  have hbsupp : ((⟨_, hbmem⟩ : effSub Γ) : S →₀ ℤ).support = {s} := by
    show (single s (d * c.num) : S →₀ ℤ).support = {s}
    exact Finsupp.support_single_ne_zero s (by omega)
  have hnpos : 0 < c.den * d.toNat := by
    have h1 : 0 < c.den := c.pos
    have h2 : 0 < d.toNat := by omega
    exact Nat.mul_pos h1 h2
  refine ⟨⟨c.den * d.toNat, hnpos⟩, ⟨_, hbmem⟩,
    (mem_primeCarrier_iff hQ p _).mpr hbsupp, ?_⟩
  refine pfCoeffHom_injective ?_
  rw [map_nsmul, hx, pfCoeff_of]
  show ((c.den * d.toNat : ℕ)) • (single s c : S →₀ ℚ) = intToRatFs (single s (d * c.num))
  have hden : ((c.den : ℚ)) ≠ 0 := by
    have := c.pos
    positivity
  have hnd : (c.num : ℚ) = c * (c.den : ℚ) := (div_eq_iff hden).mp (Rat.num_div_den c)
  have hdt : ((d.toNat : ℤ)) = d := Int.toNat_of_nonneg (le_of_lt hd)
  have hdtQ : ((d.toNat : ℕ) : ℚ) = ((d : ℤ) : ℚ) := by
    exact_mod_cast congrArg (fun z : ℤ => (z : ℚ)) hdt
  ext t
  rcases eq_or_ne s t with rfl | hst
  · simp only [Finsupp.smul_apply, Finsupp.single_eq_same, nsmul_eq_mul, intToRatFs_apply]
    push_cast
    rw [hdtQ, hnd]
    ring
  · simp [hst]

/-- ★★**`a` の `p` 成分だけを取り出した元が `Φ^pf` にあり、`a` 以下**。 -/
theorem exists_pPart (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) (a : Pf (effSub Γ)) :
    ∃ x : Pf (effSub Γ), (pfCoeffHom x : S →₀ ℚ)
        = single (effSubPrimeEquiv hQ p) (ratCoeffAt hQ p a)
      ∧ MLe x a := by
  classical
  set s := effSubPrimeEquiv hQ p with hs
  set c := ratCoeffAt hQ p a with hc
  obtain ⟨x, hx⟩ := single_mem_range_pfCoeff hQ s (ratCoeffAt_nonneg hQ p a)
  refine ⟨x, hx, ?_⟩
  have hdiff : (0 : S →₀ ℚ) ≤ (pfCoeffHom a : S →₀ ℚ) - single s c := by
    refine Finsupp.le_def.mpr (fun t => ?_)
    simp only [Finsupp.coe_zero, Pi.zero_apply, Finsupp.sub_apply]
    rcases eq_or_ne s t with rfl | hst
    · simp only [Finsupp.single_eq_same, hc, ratCoeffAt]
      simp only [sub_nonneg]
      exact le_refl _
    · rw [Finsupp.single_apply, if_neg hst, sub_zero]
      exact Finsupp.le_def.mp (pfCoeff_nonneg a) t
  have hw : (pfCoeffHom a : S →₀ ℚ) - single s c ∈ Set.range (pfCoeffHom (Γ := Γ)) := by
    rw [pfCoeffHom_range hQ]
    exact hdiff
  obtain ⟨w, hwv⟩ := hw
  refine ⟨w, pfCoeffHom_injective ?_⟩
  rw [map_add, hx, hwv]
  abel

/-! ## ★3. `factorMap ι a p = ι p a` -/

theorem iotaAt_le_of_mem_bound (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ))
    (a : Pf (effSub Γ)) {x : Pf (effSub Γ)}
    (hx : x ∈ Bound (Pf (effSub Γ)) (pCarrierPf (effSub Γ) p) a) :
    iotaAt hQ p x ≤ iotaAt hQ p a :=
  iotaAt_le_of_ratCoeffAt_le hQ p (ratCoeffAt_mono hQ p hx.2)

theorem bddAbove_iotaAt (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) (a : Pf (effSub Γ)) :
    BddAbove (iotaAt hQ p '' Bound (Pf (effSub Γ)) (pCarrierPf (effSub Γ) p) a) := by
  refine ⟨iotaAt hQ p a, ?_⟩
  rintro _ ⟨x, hx, rfl⟩
  exact iotaAt_le_of_mem_bound hQ p a hx

/-- ★★★★**因子分解写像は「その素点の係数を取る」だけ**。

★原文の `sup(Bound^p_{0}(a))` が、`a` の `p` 成分そのものになる。 -/
theorem factorMap_iotaAt (hQ : IsQCartierSubgroup Γ) (a : Pf (effSub Γ))
    (p : Prime (effSub Γ)) : factorMap (iotaAt hQ) a p = iotaAt hQ p a := by
  refine le_antisymm ?_ ?_
  · exact boundSup_le _ (zero_mem_pCarrierPf (M := effSub Γ) p) a
      (fun x hx => iotaAt_le_of_mem_bound hQ p a hx)
  · obtain ⟨x, hxv, hxle⟩ := exists_pPart hQ p a
    have hxc : x ∈ pCarrierPf (effSub Γ) p :=
      mem_pCarrierPf_of_pfCoeff_single hQ p (ratCoeffAt_nonneg hQ p a) hxv
    have heq : ratCoeffAt hQ p x = ratCoeffAt hQ p a := by
      show (pfCoeffHom x : S →₀ ℚ) (effSubPrimeEquiv hQ p) = _
      rw [hxv, Finsupp.single_eq_same]
    have hiota : iotaAt hQ p x = iotaAt hQ p a := by rw [iotaAt, iotaAt, heq]
    rw [← hiota]
    exact le_boundSup _ _ a (bddAbove_iotaAt hQ p a) ⟨hxc, hxle⟩

/-! ## ★4. 各条 -/

/-- ★★`ι p` は `pCarrierPf p` の上で単射。 -/
theorem iotaAt_injOn (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) :
    Set.InjOn (iotaAt hQ p) (pCarrierPf (effSub Γ) p) := by
  classical
  intro x hx y hy h
  have hxs := support_pfCoeff_of_mem_pCarrierPf hQ p hx
  have hys := support_pfCoeff_of_mem_pCarrierPf hQ p hy
  have hc := ratCoeffAt_of_iotaAt_eq hQ p h
  refine pfCoeffHom_injective ?_
  ext t
  rcases eq_or_ne (effSubPrimeEquiv hQ p) t with rfl | hst
  · exact hc
  · have hx0 : (pfCoeffHom x : S →₀ ℚ) t = 0 := by
      by_contra hne
      exact hst (Finset.mem_singleton.mp (hxs (Finsupp.mem_support_iff.mpr hne))).symm
    have hy0 : (pfCoeffHom y : S →₀ ℚ) t = 0 := by
      by_contra hne
      exact hst (Finset.mem_singleton.mp (hys (Finsupp.mem_support_iff.mpr hne))).symm
    rw [hx0, hy0]

theorem factorMap_iotaAt_add (hQ : IsQCartierSubgroup Γ) (a b : Pf (effSub Γ)) :
    factorMap (iotaAt hQ) (a + b) = factorMap (iotaAt hQ) a + factorMap (iotaAt hQ) b := by
  funext p
  rw [factorMap_iotaAt hQ (a + b) p]
  show iotaAt hQ p (a + b) = factorMap (iotaAt hQ) a p + factorMap (iotaAt hQ) b p
  rw [factorMap_iotaAt hQ a p, factorMap_iotaAt hQ b p, iotaAt_add]

theorem factorMap_iotaAt_injective (hQ : IsQCartierSubgroup Γ) :
    Function.Injective (factorMap (iotaAt (Γ := Γ) hQ)) := by
  intro a b h
  refine pfCoeffHom_injective ?_
  ext t
  set p := (effSubPrimeEquiv hQ).symm t with hp
  have hpt : effSubPrimeEquiv hQ p = t := (effSubPrimeEquiv hQ).apply_symm_apply t
  have h1 : iotaAt hQ p a = iotaAt hQ p b := by
    have h2 := congrFun h p
    rwa [factorMap_iotaAt hQ a p, factorMap_iotaAt hQ b p] at h2
  have h3 := ratCoeffAt_of_iotaAt_eq hQ p h1
  rw [ratCoeffAt, ratCoeffAt, hpt] at h3
  exact h3

theorem factorMap_iotaAt_mem (hQ : IsQCartierSubgroup Γ) (a : Pf (effSub Γ))
    (p : Prime (effSub Γ)) :
    factorMap (iotaAt hQ) a p ∈ iotaAt hQ p '' pCarrierPf (effSub Γ) p := by
  obtain ⟨x, hxv, _⟩ := exists_pPart hQ p a
  have hxc : x ∈ pCarrierPf (effSub Γ) p :=
    mem_pCarrierPf_of_pfCoeff_single hQ p (ratCoeffAt_nonneg hQ p a) hxv
  refine ⟨x, hxc, ?_⟩
  have heq : ratCoeffAt hQ p x = ratCoeffAt hQ p a := by
    show (pfCoeffHom x : S →₀ ℚ) (effSubPrimeEquiv hQ p) = _
    rw [hxv, Finsupp.single_eq_same]
  rw [factorMap_iotaAt hQ a p, iotaAt, iotaAt, heq]

/-- ★★★[FrdI] `Definition 2.4, (i)` の条件 **(d)**。 -/
theorem supp_field (hQ : IsQCartierSubgroup Γ) (a : RlfFactor (effSub Γ))
    (hmem : ∀ p, a p ∈ iotaAt hQ p '' pCarrierPf (effSub Γ) p)
    (b : Pf (effSub Γ)) (hsub : Supp a ⊆ Supp (factorMap (iotaAt hQ) b)) :
    a ∈ Set.range (factorMap (iotaAt (Γ := Γ) hQ)) := by
  classical
  set q : Prime (effSub Γ) → ℚ := fun p => ratCoeffAt hQ p (hmem p).choose with hqdef
  have hqnn : ∀ p, 0 ≤ q p := fun p => ratCoeffAt_nonneg hQ p _
  have haq : ∀ p, a p = Real.toNNReal ((q p : ℚ) : ℝ) := by
    intro p
    have h := (hmem p).choose_spec.2
    rw [← h, iotaAt]
  have hqne : ∀ p, q p ≠ 0 → a p ≠ 0 := by
    intro p hqp hap
    rw [haq p] at hap
    have hz : q p = 0 := by
      have h1 : iotaAt hQ p (hmem p).choose = 0 := by rw [iotaAt]; exact hap
      exact (iotaAt_eq_zero_iff hQ p _).mp h1
    exact hqp hz
  set f : S → ℚ := fun t => q ((effSubPrimeEquiv hQ).symm t) with hfdef
  have hf : ∀ t : S, f t ≠ 0 → t ∈ (pfCoeffHom b : S →₀ ℚ).support := by
    intro t ht
    set p := (effSubPrimeEquiv hQ).symm t with hp
    have hap : a p ≠ 0 := hqne p ht
    have hb : factorMap (iotaAt hQ) b p ≠ 0 := hsub hap
    rw [factorMap_iotaAt hQ b p] at hb
    have hr : ratCoeffAt hQ p b ≠ 0 := fun hz => hb ((iotaAt_eq_zero_iff hQ p b).mpr hz)
    refine Finsupp.mem_support_iff.mpr ?_
    have hpt : effSubPrimeEquiv hQ p = t := (effSubPrimeEquiv hQ).apply_symm_apply t
    rw [← hpt]
    exact hr
  set v : S →₀ ℚ := Finsupp.onFinset (pfCoeffHom b : S →₀ ℚ).support f hf with hv
  have hvapp : ∀ t, v t = f t := fun t => Finsupp.onFinset_apply
  have hvnn : (0 : S →₀ ℚ) ≤ v := by
    refine Finsupp.le_def.mpr (fun t => ?_)
    rw [hvapp t]
    simpa using hqnn _
  have hxr : v ∈ Set.range (pfCoeffHom (Γ := Γ)) := by
    rw [pfCoeffHom_range hQ]
    exact hvnn
  obtain ⟨x, hx⟩ := hxr
  refine ⟨x, ?_⟩
  funext p
  rw [factorMap_iotaAt hQ x p, iotaAt, haq p]
  congr 2
  show ratCoeffAt hQ p x = q p
  have h1 : ratCoeffAt hQ p x = v (effSubPrimeEquiv hQ p) := by
    show (pfCoeffHom x : S →₀ ℚ) _ = _
    rw [hx]
  rw [h1, hvapp, hfdef]
  simp

/-! ## ★5. `realScale` は空虚 -/

/-- ★`ℕ` と `ℝ≥0` は加法モノイドとして同型でない —— `ℝ≥0` は可除、`ℕ` は違う。 -/
theorem not_nonempty_natEquivNNReal : ¬ Nonempty (ℕ ≃+ NNReal) := by
  rintro ⟨e⟩
  set y : ℕ := e.symm ((e 1) / 2) with hy
  have h2 : (2 : NNReal) ≠ 0 := two_ne_zero
  have hhalf : (2 : ℕ) • ((e 1) / 2) = e 1 := by
    rw [two_nsmul, ← two_mul]
    field_simp
  have hy2 : (2 : ℕ) • y = 1 := by
    rw [hy, ← map_nsmul, hhalf, AddEquiv.symm_apply_apply]
  simp only [smul_eq_mul] at hy2
  omega

theorem nonempty_Mp_equiv_nat (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) :
    Nonempty (Mp (effSub Γ) p ≃+ ℕ) := by
  rw [Mp_eq_suppAtSub hQ p]
  exact nonempty_suppAtSubEquivNat hQ _

/-- ★★`Φ` の `M_p` は **ℝ-monoprime ではない**(`≃+ ℕ` なので)——
`realScale` の仮定は空虚である。 -/
theorem not_isLambdaMonoprime_real (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) :
    ¬ IsLambdaMonoprime (Mp (effSub Γ) p) MonoidType.real := by
  rintro ⟨er⟩
  obtain ⟨en⟩ := nonempty_Mp_equiv_nat hQ p
  exact not_nonempty_natEquivNNReal ⟨en.symm.trans er⟩

/-! ## ★6. まとめ -/

/-- ★★★★★**[FrdI] `Example 6.1`** —— `Φ(L)` は **perf-factorial** である
(族 `ι` を明示した形)。

原文の「one verifies immediately that `Φ(L)` is perf-factorial」。 -/
theorem isPerfFactorialWith_effSub (hQ : IsQCartierSubgroup Γ) :
    IsPerfFactorialWith (effSub Γ) (iotaAt hQ) where
  divisorial := isDivisorial_effSub Γ
  monoprimeAt := isMonoprime_Mp_effSub hQ
  embedAdd p x y _ _ := iotaAt_add hQ p x y
  embedInj := iotaAt_injOn hQ
  embedMono p _ _ _ _ h := iotaAt_le_of_ratCoeffAt_le hQ p (ratCoeffAt_mono hQ p h)
  bounded a p := bddAbove_iotaAt hQ p a
  factorAdd := factorMap_iotaAt_add hQ
  factorInj := factorMap_iotaAt_injective hQ
  factorMem := factorMap_iotaAt_mem hQ
  supp := supp_field hQ
  realScale p h := absurd h (not_isLambdaMonoprime_real hQ p)

/-- ★★★★★**`Φ(L)` は perf-factorial**。 -/
theorem isPerfFactorial_effSub (hQ : IsQCartierSubgroup Γ) :
    IsPerfFactorial (effSub Γ) :=
  ⟨iotaAt hQ, isPerfFactorialWith_effSub hQ⟩

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の「`Φ(L)` is perf-factorial」。 -/
def isPerfFactorial_effSub.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L) は perf-factorial",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
