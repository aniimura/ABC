/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPf
import ABC3.Found.Divisor.CartierPerfFactorial

/-!
# 実係数の `Φ` は perf-factorial —— `realScale` が**空虚でない**場合

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★幾何版との違いはただ 1 つ

`CartierPerfFactorial.lean`(`Example 6.1`)では `M_p ≃+ ℕ` がつねに成り立つので、
`Definition 2.4, (i)` の `realScale`(**`ℝ`-monoprime な素点では正の実数倍で閉じる**)は
**空虚**だった(`not_isLambdaMonoprime_real`)。

★算術(`Example 6.3`)ではアルキメデス素点で `M_p ≃+ ℝ≥0` なので、
`realScale` は**本当に確かめなければならない**。
★確かめる中身は「局所群が `ℝ` 全体なら `ι p` の像も `ℝ≥0` 全体」であり、
`IsLocallyMonoprimeR` の 2 分岐がそのまま効く。

## ★因子分解写像の同定

幾何版と同じく **`factorMap ι a p = ι p a`**(原文の `sup(Bound^p_0(a))` が
`a` の `p` 成分そのもの)が鍵である。★`≥` 向きに要るのは
「`a` の `p` 成分だけを取り出した元が `Φ^pf` にある」ことで、
それには `Γ` が**座標ごとに閉じている**こと(`IsCoordwiseR`)を使う。
`Example 6.3` の `arithDivGroup L` は条件が座標ごとなので満たす。
-/

namespace ABC3.Found.Divisor

open Finsupp ABC3.Found.FrdI

open scoped NNReal

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℝ)}

/-! ## ★0. 座標ごとに閉じていること -/

/-- ★★**`Γ` は座標ごとに閉じている** —— `y ∈ Γ` なら各成分だけ取り出しても `Γ` に入る。

★`Example 6.3` の `arithDivGroup L` は「各非アルキメデス成分が `log(N v)` の整数倍」
という**座標ごとの条件**で定まっているので、これを満たす。 -/
def IsCoordwiseR (Γ : AddSubgroup (S →₀ ℝ)) : Prop :=
  ∀ y ∈ Γ, ∀ s : S, single s (y s) ∈ Γ

theorem single_mem_satQ (hC : IsCoordwiseR Γ) {x : S →₀ ℝ} (hx : x ∈ satQ Γ) (s : S) :
    single s (x s) ∈ satQ Γ := by
  obtain ⟨a, ha⟩ := hx
  refine ⟨a, ?_⟩
  have h1 : ((a : ℕ)) • (single s (x s) : S →₀ ℝ)
      = single s ((((a : ℕ)) • x) s) := by
    ext t
    rcases eq_or_ne s t with rfl | hst
    · simp
    · simp [hst]
  rw [h1]
  exact hC _ ha s

/-- ★`Φ^pf` の像はちょうど `Γ^ℚ ∩ ℝ≥0[S]`。 -/
theorem pfCoeffRHom_range :
    Set.range (pfCoeffRHom (Γ := Γ)) = {x : S →₀ ℝ | x ∈ effR (satQ Γ)} := by
  refine Set.Subset.antisymm ?_ ?_
  · rintro _ ⟨x, rfl⟩
    exact pfCoeffR_mem x
  · intro x hx
    exact pfCoeffRHom_surjective hx

/-! ## ★1. 素点ごとの係数と `ι` -/

/-- ★`Φ^pf` の元の、素点 `p` での係数(実数)。 -/
noncomputable def realCoeffAt (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    (x : Pf (effR Γ)) : ℝ :=
  (pfCoeffRHom x : S →₀ ℝ) (effRPrimeEquiv hG p)

theorem realCoeffAt_nonneg (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (x : Pf (effR Γ)) :
    0 ≤ realCoeffAt hG p x :=
  Finsupp.le_def.mp (pfCoeffR_nonneg x) _

theorem realCoeffAt_add (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (x y : Pf (effR Γ)) :
    realCoeffAt hG p (x + y) = realCoeffAt hG p x + realCoeffAt hG p y := by
  show (pfCoeffRHom (x + y)) _ = _
  rw [map_add]
  rfl

/-- ★★**`ι p`** —— 原文の `M^rlf_p ≃ ℝ≥0` の同一視。 -/
noncomputable def iotaR (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    (x : Pf (effR Γ)) : ℝ≥0 :=
  Real.toNNReal (realCoeffAt hG p x)

theorem iotaR_add (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (x y : Pf (effR Γ)) :
    iotaR hG p (x + y) = iotaR hG p x + iotaR hG p y := by
  rw [iotaR, iotaR, iotaR, realCoeffAt_add]
  exact Real.toNNReal_add (realCoeffAt_nonneg hG p x) (realCoeffAt_nonneg hG p y)

theorem realCoeffAt_of_iotaR_eq (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    {x y : Pf (effR Γ)} (h : iotaR hG p x = iotaR hG p y) :
    realCoeffAt hG p x = realCoeffAt hG p y := by
  have h2 := congrArg (fun t : ℝ≥0 => (t : ℝ)) h
  simp only [iotaR] at h2
  rwa [Real.coe_toNNReal _ (realCoeffAt_nonneg hG p x),
    Real.coe_toNNReal _ (realCoeffAt_nonneg hG p y)] at h2

theorem iotaR_le_of_realCoeffAt_le (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    {x y : Pf (effR Γ)} (h : realCoeffAt hG p x ≤ realCoeffAt hG p y) :
    iotaR hG p x ≤ iotaR hG p y :=
  Real.toNNReal_mono h

/-- ★`MLe` は係数を増やす。 -/
theorem realCoeffAt_mono (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    {x y : Pf (effR Γ)} (h : MLe x y) : realCoeffAt hG p x ≤ realCoeffAt hG p y := by
  obtain ⟨c, rfl⟩ := h
  rw [realCoeffAt_add]
  exact le_add_of_nonneg_right (realCoeffAt_nonneg hG p c)

theorem iotaR_eq_zero_iff (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (x : Pf (effR Γ)) :
    iotaR hG p x = 0 ↔ realCoeffAt hG p x = 0 := by
  rw [iotaR, Real.toNNReal_eq_zero]
  have hnn := realCoeffAt_nonneg hG p x
  exact ⟨fun h => le_antisymm h hnn, fun h => h.le⟩

/-! ## ★2. `pCarrierPf p` の同定 -/

theorem pfCoeffR_of (b : effR Γ) :
    (pfCoeffRHom (Pf.of b) : S →₀ ℝ) = ((b : effR Γ) : S →₀ ℝ) := by
  show pfCoeffR (Pf.mk b 1) = _
  rw [pfCoeffR_mk]
  simp

/-- ★★`pCarrierPf p` の元は、係数が `s` にしか乗らない。 -/
theorem support_pfCoeffR_of_mem_pCarrierPf (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    {x : Pf (effR Γ)} (hx : x ∈ pCarrierPf (effR Γ) p) :
    (pfCoeffRHom x : S →₀ ℝ).support ⊆ {effRPrimeEquiv hG p} := by
  classical
  rcases hx with ⟨n, b, hb, hnx⟩ | hx0
  · have hbsupp : ((b : effR Γ) : S →₀ ℝ).support = {effRPrimeEquiv hG p} :=
      (mem_primeCarrier_iff_R hG p b).mp hb
    have h1 : pfCoeffRHom (((n : ℕ+) : ℕ) • x) = pfCoeffRHom (Pf.of b) := by rw [hnx]
    rw [map_nsmul, pfCoeffR_of] at h1
    intro t ht
    have hne : (pfCoeffRHom x : S →₀ ℝ) t ≠ 0 := Finsupp.mem_support_iff.mp ht
    have hnz : ((n : ℕ) : ℝ) ≠ 0 := by
      have := n.2
      positivity
    have h2 := congrArg (fun v : S →₀ ℝ => v t) h1
    simp only [Finsupp.smul_apply, nsmul_eq_mul] at h2
    have h3 : ((b : effR Γ) : S →₀ ℝ) t ≠ 0 := by
      rw [← h2]
      exact mul_ne_zero hnz hne
    rw [← hbsupp]
    exact Finsupp.mem_support_iff.mpr h3
  · rw [Set.mem_singleton_iff] at hx0
    rw [hx0, map_zero]
    simp

/-- ★★係数が `s` にしか乗らない `Φ^pf` の元は `pCarrierPf p` に入る。

★幾何版と違って `Q`-Cartier の生成元は要らない —— `single s c` が
`Γ^ℚ` に入っていること(`pfCoeffR_mem_satQ`)から直接 `n` が取れる。 -/
theorem mem_pCarrierPf_of_pfCoeffR_single (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    {x : Pf (effR Γ)} {c : ℝ} (hc : 0 ≤ c)
    (hx : (pfCoeffRHom x : S →₀ ℝ) = single (effRPrimeEquiv hG p) c) :
    x ∈ pCarrierPf (effR Γ) p := by
  classical
  set s := effRPrimeEquiv hG p with hs
  rcases eq_or_lt_of_le hc with hc0 | hcpos
  · refine Or.inr ?_
    refine Set.mem_singleton_iff.mpr (pfCoeffRHom_injective ?_)
    rw [hx, map_zero, ← hc0]
    simp
  refine Or.inl ?_
  -- ★`single s c ∈ Γ^ℚ` から分母 `n` を取る
  have hmemQ : (single s c : S →₀ ℝ) ∈ satQ Γ := by
    rw [← hx]
    exact pfCoeffR_mem_satQ x
  obtain ⟨n, hn⟩ := hmemQ
  have hsmul : ((n : ℕ)) • (single s c : S →₀ ℝ) = single s (((n : ℕ) : ℝ) * c) := by
    ext t
    rcases eq_or_ne s t with rfl | hst
    · simp
    · simp [hst]
  rw [hsmul] at hn
  have hncpos : 0 < ((n : ℕ) : ℝ) * c := by
    have : (0 : ℝ) < ((n : ℕ) : ℝ) := by
      have := n.2
      positivity
    exact mul_pos this hcpos
  have hbmem : (single s (((n : ℕ) : ℝ) * c) : S →₀ ℝ) ∈ effR Γ := by
    refine mem_effR.mpr ⟨hn, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simpa using hncpos.le
    · simp [hst]
  have hbsupp : ((⟨_, hbmem⟩ : effR Γ) : S →₀ ℝ).support = {s} := by
    show (single s (((n : ℕ) : ℝ) * c) : S →₀ ℝ).support = {s}
    exact Finsupp.support_single_ne_zero s hncpos.ne'
  refine ⟨n, ⟨_, hbmem⟩, (mem_primeCarrier_iff_R hG p _).mpr hbsupp, ?_⟩
  refine pfCoeffRHom_injective ?_
  rw [map_nsmul, hx, pfCoeffR_of]
  show ((n : ℕ)) • (single s c : S →₀ ℝ) = single s (((n : ℕ) : ℝ) * c)
  exact hsmul

/-- ★★**`a` の `p` 成分だけを取り出した元が `Φ^pf` にあり、`a` 以下**。 -/
theorem exists_pPartR (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    (a : Pf (effR Γ)) :
    ∃ x : Pf (effR Γ), (pfCoeffRHom x : S →₀ ℝ)
        = single (effRPrimeEquiv hG p) (realCoeffAt hG p a)
      ∧ MLe x a := by
  classical
  set s := effRPrimeEquiv hG p with hs
  set c := realCoeffAt hG p a with hc
  -- ★`single s c` は `Γ^ℚ` の有効元
  have hsingle : (single s c : S →₀ ℝ) ∈ effR (satQ Γ) := by
    refine mem_effR.mpr ⟨?_, fun t => ?_⟩
    · have h1 : ((pfCoeffRHom a : S →₀ ℝ) s) = c := rfl
      have := single_mem_satQ hC (pfCoeffR_mem_satQ a) s
      rwa [h1] at this
    · rcases eq_or_ne s t with rfl | hst
      · simpa using realCoeffAt_nonneg hG p a
      · simp [hst]
  obtain ⟨x, hx⟩ := pfCoeffRHom_surjective hsingle
  refine ⟨x, hx, ?_⟩
  have hdiff : (pfCoeffRHom a : S →₀ ℝ) - single s c ∈ effR (satQ Γ) := by
    refine mem_effR.mpr ⟨(satQ Γ).sub_mem (pfCoeffR_mem_satQ a) (mem_effR.mp hsingle).1,
      fun t => ?_⟩
    simp only [Finsupp.sub_apply]
    rcases eq_or_ne s t with rfl | hst
    · simp only [Finsupp.single_eq_same, hc]
      show (0:ℝ) ≤ (pfCoeffRHom a : S →₀ ℝ) s - (pfCoeffRHom a : S →₀ ℝ) s
      simp
    · rw [Finsupp.single_apply, if_neg hst, sub_zero]
      exact Finsupp.le_def.mp (pfCoeffR_nonneg a) t
  obtain ⟨w, hwv⟩ := pfCoeffRHom_surjective hdiff
  refine ⟨w, pfCoeffRHom_injective ?_⟩
  rw [map_add, hx, hwv]
  abel

/-! ## ★3. `factorMap ι a p = ι p a` -/

theorem iotaR_le_of_mem_bound (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    (a : Pf (effR Γ)) {x : Pf (effR Γ)}
    (hx : x ∈ Bound (Pf (effR Γ)) (pCarrierPf (effR Γ) p) a) :
    iotaR hG p x ≤ iotaR hG p a :=
  iotaR_le_of_realCoeffAt_le hG p (realCoeffAt_mono hG p hx.2)

theorem bddAbove_iotaR (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (a : Pf (effR Γ)) :
    BddAbove (iotaR hG p '' Bound (Pf (effR Γ)) (pCarrierPf (effR Γ) p) a) := by
  refine ⟨iotaR hG p a, ?_⟩
  rintro _ ⟨x, hx, rfl⟩
  exact iotaR_le_of_mem_bound hG p a hx

/-- ★★★★**因子分解写像は「その素点の係数を取る」だけ**。 -/
theorem factorMap_iotaR (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ) (a : Pf (effR Γ))
    (p : Prime (effR Γ)) : factorMap (iotaR hG) a p = iotaR hG p a := by
  refine le_antisymm ?_ ?_
  · exact boundSup_le _ (zero_mem_pCarrierPf (M := effR Γ) p) a
      (fun x hx => iotaR_le_of_mem_bound hG p a hx)
  · obtain ⟨x, hxv, hxle⟩ := exists_pPartR hC hG p a
    have hxc : x ∈ pCarrierPf (effR Γ) p :=
      mem_pCarrierPf_of_pfCoeffR_single hG p (realCoeffAt_nonneg hG p a) hxv
    have heq : realCoeffAt hG p x = realCoeffAt hG p a := by
      show (pfCoeffRHom x : S →₀ ℝ) (effRPrimeEquiv hG p) = _
      rw [hxv, Finsupp.single_eq_same]
    have hiota : iotaR hG p x = iotaR hG p a := by rw [iotaR, iotaR, heq]
    rw [← hiota]
    exact le_boundSup _ _ a (bddAbove_iotaR hG p a) ⟨hxc, hxle⟩

/-! ## ★4. 各条 -/

theorem iotaR_injOn (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) :
    Set.InjOn (iotaR hG p) (pCarrierPf (effR Γ) p) := by
  classical
  intro x hx y hy h
  have hxs := support_pfCoeffR_of_mem_pCarrierPf hG p hx
  have hys := support_pfCoeffR_of_mem_pCarrierPf hG p hy
  have hc := realCoeffAt_of_iotaR_eq hG p h
  refine pfCoeffRHom_injective ?_
  ext t
  rcases eq_or_ne (effRPrimeEquiv hG p) t with rfl | hst
  · exact hc
  · have hx0 : (pfCoeffRHom x : S →₀ ℝ) t = 0 := by
      by_contra hne
      exact hst (Finset.mem_singleton.mp (hxs (Finsupp.mem_support_iff.mpr hne))).symm
    have hy0 : (pfCoeffRHom y : S →₀ ℝ) t = 0 := by
      by_contra hne
      exact hst (Finset.mem_singleton.mp (hys (Finsupp.mem_support_iff.mpr hne))).symm
    rw [hx0, hy0]

theorem factorMap_iotaR_add (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ)
    (a b : Pf (effR Γ)) :
    factorMap (iotaR hG) (a + b) = factorMap (iotaR hG) a + factorMap (iotaR hG) b := by
  funext p
  rw [factorMap_iotaR hC hG (a + b) p]
  show iotaR hG p (a + b) = factorMap (iotaR hG) a p + factorMap (iotaR hG) b p
  rw [factorMap_iotaR hC hG a p, factorMap_iotaR hC hG b p, iotaR_add]

theorem factorMap_iotaR_injective (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ) :
    Function.Injective (factorMap (iotaR (Γ := Γ) hG)) := by
  intro a b h
  refine pfCoeffRHom_injective ?_
  ext t
  set p := (effRPrimeEquiv hG).symm t with hp
  have hpt : effRPrimeEquiv hG p = t := (effRPrimeEquiv hG).apply_symm_apply t
  have h1 : iotaR hG p a = iotaR hG p b := by
    have h2 := congrFun h p
    rwa [factorMap_iotaR hC hG a p, factorMap_iotaR hC hG b p] at h2
  have h3 := realCoeffAt_of_iotaR_eq hG p h1
  rw [realCoeffAt, realCoeffAt, hpt] at h3
  exact h3

theorem factorMap_iotaR_mem (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ) (a : Pf (effR Γ))
    (p : Prime (effR Γ)) :
    factorMap (iotaR hG) a p ∈ iotaR hG p '' pCarrierPf (effR Γ) p := by
  obtain ⟨x, hxv, _⟩ := exists_pPartR hC hG p a
  have hxc : x ∈ pCarrierPf (effR Γ) p :=
    mem_pCarrierPf_of_pfCoeffR_single hG p (realCoeffAt_nonneg hG p a) hxv
  refine ⟨x, hxc, ?_⟩
  have heq : realCoeffAt hG p x = realCoeffAt hG p a := by
    show (pfCoeffRHom x : S →₀ ℝ) (effRPrimeEquiv hG p) = _
    rw [hxv, Finsupp.single_eq_same]
  rw [factorMap_iotaR hC hG a p, iotaR, iotaR, heq]

/-- ★★★[FrdI] `Definition 2.4, (i)` の条件 **(d)**。 -/
theorem supp_fieldR (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ) (a : RlfFactor (effR Γ))
    (hmem : ∀ p, a p ∈ iotaR hG p '' pCarrierPf (effR Γ) p)
    (b : Pf (effR Γ)) (hsub : Supp a ⊆ Supp (factorMap (iotaR hG) b)) :
    a ∈ Set.range (factorMap (iotaR (Γ := Γ) hG)) := by
  classical
  set q : Prime (effR Γ) → ℝ := fun p => realCoeffAt hG p (hmem p).choose with hqdef
  have hqnn : ∀ p, 0 ≤ q p := fun p => realCoeffAt_nonneg hG p _
  have haq : ∀ p, a p = Real.toNNReal (q p) := by
    intro p
    have h := (hmem p).choose_spec.2
    rw [← h, iotaR]
  have hqne : ∀ p, q p ≠ 0 → a p ≠ 0 := by
    intro p hqp hap
    rw [haq p] at hap
    have h1 : iotaR hG p (hmem p).choose = 0 := by rw [iotaR]; exact hap
    exact hqp ((iotaR_eq_zero_iff hG p _).mp h1)
  set f : S → ℝ := fun t => q ((effRPrimeEquiv hG).symm t) with hfdef
  have hf : ∀ t : S, f t ≠ 0 → t ∈ (pfCoeffRHom b : S →₀ ℝ).support := by
    intro t ht
    set p := (effRPrimeEquiv hG).symm t with hp
    have hap : a p ≠ 0 := hqne p ht
    have hb : factorMap (iotaR hG) b p ≠ 0 := hsub hap
    rw [factorMap_iotaR hC hG b p] at hb
    have hr : realCoeffAt hG p b ≠ 0 := fun hz => hb ((iotaR_eq_zero_iff hG p b).mpr hz)
    refine Finsupp.mem_support_iff.mpr ?_
    have hpt : effRPrimeEquiv hG p = t := (effRPrimeEquiv hG).apply_symm_apply t
    rw [← hpt]
    exact hr
  set v : S →₀ ℝ := Finsupp.onFinset (pfCoeffRHom b : S →₀ ℝ).support f hf with hv
  have hvapp : ∀ t, v t = f t := fun t => Finsupp.onFinset_apply
  have hvnn : ∀ t, (0 : ℝ) ≤ v t := by
    intro t
    rw [hvapp t]
    exact hqnn _
  -- ★`v` は座標ごとに `Γ^ℚ` の元の成分なので、有限和として `Γ^ℚ` に入る
  have hvsingle : ∀ t : S, (single t (v t) : S →₀ ℝ) ∈ satQ Γ := by
    intro t
    set p := (effRPrimeEquiv hG).symm t with hp
    have hpt : effRPrimeEquiv hG p = t := (effRPrimeEquiv hG).apply_symm_apply t
    have h1 : v t = (pfCoeffRHom (hmem p).choose : S →₀ ℝ) t := by
      rw [hvapp t, hfdef]
      show q p = _
      rw [hqdef]
      show realCoeffAt hG p (hmem p).choose = _
      rw [realCoeffAt, hpt]
    rw [h1]
    exact single_mem_satQ hC (pfCoeffR_mem_satQ _) t
  have hvmem : v ∈ effR (satQ Γ) := by
    refine mem_effR.mpr ⟨?_, hvnn⟩
    have hsum : ∑ t ∈ v.support, single t (v t) = v := Finsupp.sum_single v
    have hall := AddSubgroup.sum_mem (satQ Γ) (fun t (_ : t ∈ v.support) => hvsingle t)
    rwa [hsum] at hall
  obtain ⟨x, hx⟩ := pfCoeffRHom_surjective hvmem
  refine ⟨x, ?_⟩
  funext p
  rw [factorMap_iotaR hC hG x p, iotaR, haq p]
  congr 1
  show realCoeffAt hG p x = q p
  have h1 : realCoeffAt hG p x = v (effRPrimeEquiv hG p) := by
    show (pfCoeffRHom x : S →₀ ℝ) _ = _
    rw [hx]
  rw [h1, hvapp, hfdef]
  simp

/-! ## ★5. `realScale` —— アルキメデス素点では**空虚でない** -/

/-- ★★★★★**局所群が `ℝ` 全体なら `ι p` の像は `ℝ≥0` 全体**。

★★これが `Example 6.1`(幾何)には現れない条で、`Example 6.3` のアルキメデス素点で
`realScale` を**実際に確かめる**中身である。 -/
theorem image_iotaR_univ_of_full (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ))
    (hfull : ∀ c : ℝ, single (effRPrimeEquiv hG p) c ∈ Γ) (r : ℝ≥0) :
    r ∈ iotaR hG p '' pCarrierPf (effR Γ) p := by
  classical
  set s := effRPrimeEquiv hG p with hs
  have hmem : (single s (r : ℝ) : S →₀ ℝ) ∈ effR Γ := by
    refine mem_effR.mpr ⟨hfull _, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simpa using r.2
    · simp [hst]
  refine ⟨Pf.of ⟨_, hmem⟩, ?_, ?_⟩
  · refine mem_pCarrierPf_of_pfCoeffR_single hG p r.2 ?_
    rw [pfCoeffR_of]
    show (single s (r : ℝ) : S →₀ ℝ) = single (effRPrimeEquiv hG p) (r : ℝ)
    rw [hs]
  · rw [iotaR, realCoeffAt]
    have h1 : (pfCoeffRHom (Pf.of (⟨_, hmem⟩ : effR Γ)) : S →₀ ℝ) = single s (r : ℝ) := by
      rw [pfCoeffR_of]
    rw [h1, Finsupp.single_eq_same]
    exact Real.toNNReal_coe

/-- ★★`M_p` が `ℝ`-monoprime なら、その素点の局所群は `ℝ` 全体である。

★`IsLocallyMonoprimeR` の 2 分岐のうち離散の方は `M_p ≃+ ℕ` を与えるが、
`ℕ` と `ℝ≥0` は同型でない(`not_nonempty_natEquivNNReal`)。 -/
theorem full_of_isLambdaMonoprime_real (hL : IsLocallyMonoprimeR Γ) (p : Prime (effR Γ))
    (h : IsLambdaMonoprime (Mp (effR Γ) p) MonoidType.real) :
    ∀ c : ℝ, single (effRPrimeEquiv (isGenSubgroupR_of_isLocallyMonoprimeR hL) p) c ∈ Γ := by
  set hG := isGenSubgroupR_of_isLocallyMonoprimeR hL with hGdef
  rcases hL (effRPrimeEquiv hG p) with ⟨d, hd, hspec⟩ | hfull
  · exfalso
    obtain ⟨er⟩ := h
    have hn : Nonempty (Mp (effR Γ) p ≃+ ℕ) := by
      rw [Mp_eq_suppAtSubR hG p]
      exact nonempty_suppAtSubR_equiv_nat hd hspec
    obtain ⟨en⟩ := hn
    exact not_nonempty_natEquivNNReal ⟨en.symm.trans er⟩
  · exact hfull

/-! ## ★6. まとめ -/

/-- ★★★★★★**実係数の `Φ` は perf-factorial**(族 `ι` を明示した形)。

★★`realScale` が**空虚でない**のが幾何版との違いである ——
アルキメデス素点では `M_p ≃+ ℝ≥0` なので、正の実数倍で閉じることを
`image_iotaR_univ_of_full` で実際に示す。 -/
theorem isPerfFactorialWith_effR (hC : IsCoordwiseR Γ) (hL : IsLocallyMonoprimeR Γ) :
    IsPerfFactorialWith (effR Γ) (iotaR (isGenSubgroupR_of_isLocallyMonoprimeR hL)) where
  divisorial := isDivisorial_effR Γ
  monoprimeAt := isMonoprime_Mp_effR hL
  embedAdd p x y _ _ := iotaR_add _ p x y
  embedInj := iotaR_injOn _
  embedMono p _ _ _ _ h := iotaR_le_of_realCoeffAt_le _ p (realCoeffAt_mono _ p h)
  bounded a p := bddAbove_iotaR _ p a
  factorAdd := factorMap_iotaR_add hC _
  factorInj := factorMap_iotaR_injective hC _
  factorMem := factorMap_iotaR_mem hC _
  supp := supp_fieldR hC _
  realScale p h r _ x _ := by
    have hfull := full_of_isLambdaMonoprime_real hL p h
    exact image_iotaR_univ_of_full _ p hfull (r * x)

/-- ★★★★★★**実係数の `Φ` は perf-factorial**。 -/
theorem isPerfFactorial_effR (hC : IsCoordwiseR Γ) (hL : IsLocallyMonoprimeR Γ) :
    IsPerfFactorial (effR Γ) :=
  ⟨_, isPerfFactorialWith_effR hC hL⟩

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の「`Φ(L)` は perf-factorial」(実係数の骨格)。 -/
def isPerfFactorial_effR.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) は perf-factorial(実係数の骨格)",
    sectionId := "frdi-example-6-3" }

/-! ## ★★`Φ(L)` は primary 元で生成される —— non-dilating の入力(2026-08-21)

`Theorem 6.4, (i)` / `Theorem 6.2, (iii)` の non-dilating は
`isNonDilating_of_primary_sharp`(`MonoidTransport.lean`)に 2 つ入力すれば出る:

| 入力 | ここで閉じるか |
|---|---|
| `hgen` : primary 元で生成される | ★**本節で閉じる**(`IsCoordwiseR` から) |
| `hfix` : 自己同型は素点を素点へ(係数 1 で)移す | ★まだ —— `Datum` に条項が無い |

★★**`IsCoordwiseR` がちょうど `hgen` を与える** —— 座標ごとに閉じているので、
`x = Σ_{s ∈ supp x} single s (x s)` の各項が `Φ(L)` に入る。
★各項が primary であることは「台が 1 点なら primary」で、
中身は**実数の Archimedes 性**(`exists_nat_gt`)だけである。 -/

/-- ★★**台が 1 点の元は primary**。

★`b ≼ a` なら `b` の台も `{s}` に入る(非負だから)。逆向きは
`m·b ≥ a` となる `m` を Archimedes 性で取るだけ。 -/
theorem isPrimaryElt_effR_single {s : S} {r : ℝ} (hr : 0 < r)
    (hmem : Finsupp.single s r ∈ effR Γ) :
    IsPrimaryElt (⟨Finsupp.single s r, hmem⟩ : effR Γ) := by
  refine ⟨?_, ?_⟩
  · intro h0
    have h : Finsupp.single s r = (0 : S →₀ ℝ) := congrArg Subtype.val h0
    rw [Finsupp.single_eq_zero] at h
    exact hr.ne' h
  · rintro b hb0 ⟨n, hn, c, hbc⟩
    have hcoe : (b : S →₀ ℝ) + (c : S →₀ ℝ) = n • (Finsupp.single s r : S →₀ ℝ) :=
      congrArg Subtype.val hbc
    have hbnn : ∀ t, 0 ≤ (b : S →₀ ℝ) t := (mem_effR.mp b.2).2
    have hcnn : ∀ t, 0 ≤ (c : S →₀ ℝ) t := (mem_effR.mp c.2).2
    have hbt : ∀ t, t ≠ s → (b : S →₀ ℝ) t = 0 := by
      intro t ht
      have h1 := congrArg (fun f : S →₀ ℝ => f t) hcoe
      have hz : (Finsupp.single s r : S →₀ ℝ) t = 0 := Finsupp.single_eq_of_ne ht
      simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
      rw [hz, smul_zero] at h1
      linarith [hbnn t, hcnn t]
    have hbsupp : (b : S →₀ ℝ) = Finsupp.single s ((b : S →₀ ℝ) s) := by
      refine Finsupp.ext fun t => ?_
      rcases eq_or_ne t s with rfl | ht
      · rw [Finsupp.single_eq_same]
      · have hz : (Finsupp.single s ((b : S →₀ ℝ) s) : S →₀ ℝ) t = 0 :=
          Finsupp.single_eq_of_ne ht
        rw [hbt t ht, hz]
    have hbs : 0 < (b : S →₀ ℝ) s := by
      rcases lt_or_eq_of_le (hbnn s) with h | h
      · exact h
      · refine absurd (Subtype.ext ?_) hb0
        rw [hbsupp, ← h, Finsupp.single_zero]
        rfl
    obtain ⟨m, hm⟩ := exists_nat_gt (r / (b : S →₀ ℝ) s)
    have hrm : r < (m : ℝ) * (b : S →₀ ℝ) s := (div_lt_iff₀ hbs).mp hm
    have hm0 : 0 < m := by
      rcases Nat.eq_zero_or_pos m with h | h
      · rw [h] at hrm; norm_num at hrm; linarith
      · exact h
    have hmb : (m : ℕ) • (b : S →₀ ℝ) = Finsupp.single s ((m : ℝ) * (b : S →₀ ℝ) s) := by
      rw [hbsupp, Finsupp.smul_single, Finsupp.single_eq_same, nsmul_eq_mul]
    have hdmem : Finsupp.single s ((m : ℝ) * (b : S →₀ ℝ) s - r) ∈ effR Γ := by
      refine mem_effR.mpr ⟨?_, ?_⟩
      · have h1 : Finsupp.single s ((m : ℝ) * (b : S →₀ ℝ) s - r)
            = (m : ℕ) • (b : S →₀ ℝ) - Finsupp.single s r := by
          rw [hmb, ← Finsupp.single_sub]
        rw [h1]
        exact Γ.sub_mem (Γ.nsmul_mem (mem_effR.mp b.2).1 m) (mem_effR.mp hmem).1
      · intro t
        rcases eq_or_ne t s with rfl | ht
        · rw [Finsupp.single_eq_same]; linarith
        · have hz : (Finsupp.single s ((m : ℝ) * (b : S →₀ ℝ) s - r) : S →₀ ℝ) t = 0 :=
            Finsupp.single_eq_of_ne ht
          rw [hz]
    refine ⟨m, hm0, ⟨⟨_, hdmem⟩, Subtype.ext ?_⟩⟩
    show Finsupp.single s r + Finsupp.single s ((m : ℝ) * (b : S →₀ ℝ) s - r)
      = (m : ℕ) • (b : S →₀ ℝ)
    rw [← Finsupp.single_add, hmb]
    congr 1
    ring

/-- ★座標ごとに切り出した元。 -/
noncomputable def effRSingle (hc : IsCoordwiseR Γ) (x : effR Γ) (s : S) : effR Γ :=
  ⟨Finsupp.single s ((x : S →₀ ℝ) s),
    mem_effR.mpr ⟨hc (x : S →₀ ℝ) (mem_effR.mp x.2).1 s, by
      intro t
      rcases eq_or_ne t s with rfl | h
      · rw [Finsupp.single_eq_same]; exact (mem_effR.mp x.2).2 t
      · rw [Finsupp.single_eq_of_ne h]⟩⟩

/-- ★元は座標ごとの和。 -/
theorem effR_eq_sum_single (hc : IsCoordwiseR Γ) (x : effR Γ) :
    x = ∑ s ∈ ((x : S →₀ ℝ)).support, effRSingle hc x s := by
  refine Subtype.ext ?_
  rw [AddSubmonoid.coe_finsetSum]
  exact (Finsupp.sum_single (x : S →₀ ℝ)).symm

/-- ★★★★**`Φ(L)` は primary 元で生成される** —— non-dilating の入力 `hgen`。 -/
theorem closure_primary_effR_eq_top (hc : IsCoordwiseR Γ) :
    AddSubmonoid.closure {a : effR Γ | IsPrimaryElt a} = ⊤ := by
  refine eq_top_iff.mpr fun x _ => ?_
  rw [effR_eq_sum_single hc x]
  refine sum_mem fun s hs => AddSubmonoid.subset_closure ?_
  have hpos : 0 < (x : S →₀ ℝ) s :=
    lt_of_le_of_ne ((mem_effR.mp x.2).2 s) (Ne.symm (Finsupp.mem_support_iff.mp hs))
  exact isPrimaryElt_effR_single hpos _

/-- ★★★locator —— `Example 6.3` の `Φ(L)` が primary 元で生成されること
(`Theorem 6.4, (i)` の non-dilating の入力)。 -/
def closure_primary_effR_eq_top.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) は primary 元(台が 1 点の元)で生成される",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
