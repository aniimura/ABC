/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithDivisor
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem
import Mathlib.NumberTheory.NumberField.ClassNumber

/-!
# `δ_A : Pic_Φ(A) ≅ ℝ`(鎖 `sec6items` の `thm64-i-dirichlet`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> isomorphism of groups

## ★★主張

`A ∈ Ob(𝒞^rlf)` なので `Φ^rlf(L)` は**全素点上の実係数**である。
`Pic_Φ(A) = Φ^gp(A)/Φ^birat(A)` で `Φ^birat` は主因子の張る部分空間だから、
示すべきは

    Submodule.span ℝ (Set.range (arithDiv : Lˣ → (ArithPlace L →₀ ℝ)))
      = LinearMap.ker arithDegreeLin

であり、これから `Pic_Φ(A) ≅ ℝ`(`arithPicIso`)が出る。

## ★★★中身は「Dirichlet 単数定理 + 類数有限」

| 段 | 中身 |
|---|---|
| ⊆ | 積公式(`arithDegree_arithDiv`、在庫) |
| ⊇ の第 1 段 | ★**有限素点を消す** —— `𝔭^h`(`h` = 類数)は単項だから、`ord_𝔭 > 0` で他の有限素点では `0` の元が取れる |
| ⊇ の第 2 段 | 残るのは無限素点にしか台がない次数 0 の因子 |
| ⊇ の第 3 段 | ★★**Dirichlet** —— mathlib の `unitLattice_span_eq_top` |

★★橋は 1 本だけ ——
`arithPlaceLog x (inr w) = −mult(w)·log(w x)` は mathlib の
`logEmbedding_component` の**符号違い**である(`toLogSpace_arithDiv_unit`)。

★★★`w₀`(mathlib が落とす素点)の成分は「次数 0」で決まるので、
`Θ : D ↦ (w ↦ −D(inr w))` は「無限素点にしか台がない次数 0 の因子」の上で**単射**。
これで `span ℝ (unitLattice) = ⊤` がそのまま使える。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_gen_pow_classNumber` | `𝔭^h` は単項 |
| `ordFin_eq_zero_of_notMem` / `ordFin_pos_of_mem` | `ord` と素イデアルの所属 |
| `exists_ordFin_pos_zero_elsewhere` | ★1 つの有限素点だけを台にもつ主因子 |
| `toLogSpace_arithDiv_unit` | ★`Θ(div u) = logEmbedding u` |
| `eq_zero_of_inf_deg_theta` | ★`Θ` は単射(次数 0・無限素点上) |
| `map_toLogSpace_unitSpan` | ★★Dirichlet(`Θ` の像が全体) |
| `single_inl_mem_sup` | ★有限素点の δ は `span + 無限素点` に入る |
| `principalSpan_eq_ker` | ★★★★**主定理** |
| `arithPicIso` | ★★★★★**`Pic_Φ(A) ≅ ℝ`** |
-/

namespace ABC3.Found.Divisor

open _root_.NumberField _root_.NumberField.Units
open _root_.NumberField.Units.dirichletUnitTheorem IsDedekindDomain Ideal
open scoped _root_.NumberField Classical

variable {L : Type*} [Field L] [NumberField L]

/-! ## ★1. 類数 —— `𝔭^h` は単項 -/

/-- ★★`𝔭^h`(`h` = 類数)は単項。 -/
theorem exists_gen_pow_classNumber (𝔭 : HeightOneSpectrum (𝓞 L)) :
    ∃ y : 𝓞 L, y ≠ 0 ∧ Ideal.span {y} = 𝔭.asIdeal ^ (classNumber L) := by
  have hne0 : 𝔭.asIdeal ≠ 0 := by rw [Ideal.zero_eq_bot]; exact 𝔭.ne_bot
  have hmem : 𝔭.asIdeal ∈ nonZeroDivisors (Ideal (𝓞 L)) := mem_nonZeroDivisors_of_ne_zero hne0
  have hpowmem : 𝔭.asIdeal ^ (classNumber L) ∈ nonZeroDivisors (Ideal (𝓞 L)) := pow_mem hmem _
  have h1 : (ClassGroup.mk0 ⟨𝔭.asIdeal, hmem⟩) ^ (classNumber L) = 1 := by
    rw [NumberField.classNumber]
    exact pow_card_eq_one
  have h2 : ClassGroup.mk0 (⟨𝔭.asIdeal ^ (classNumber L), hpowmem⟩) = 1 := by
    rw [show (⟨𝔭.asIdeal ^ (classNumber L), hpowmem⟩ : ↥(nonZeroDivisors (Ideal (𝓞 L))))
        = ⟨𝔭.asIdeal, hmem⟩ ^ (classNumber L) from Subtype.ext rfl, map_pow, h1]
  obtain ⟨y, hy⟩ := (ClassGroup.mk0_eq_one_iff hpowmem).mp h2
  refine ⟨y, ?_, hy.symm⟩
  intro hy0
  rw [hy0] at hy
  have hbot : 𝔭.asIdeal ^ (classNumber L) = ⊥ := by
    rw [hy]
    exact Submodule.span_zero_singleton _
  exact (pow_ne_zero (classNumber L) hne0) (by rw [Ideal.zero_eq_bot]; exact hbot)

/-! ## ★2. `ord` と素イデアルの所属 -/

theorem one_lt_absNorm_real (v : HeightOneSpectrum (𝓞 L)) :
    (1 : ℝ) < (Ideal.absNorm v.asIdeal : ℝ) := by
  exact_mod_cast NumberField.HeightOneSpectrum.one_lt_absNorm v

/-- ★`y ∉ 𝔭` なら `ord_𝔭(y) = 0`。 -/
theorem ordFin_eq_zero_of_notMem (v : HeightOneSpectrum (𝓞 L)) (y : 𝓞 L) (hy : y ≠ 0)
    (hnot : y ∉ v.asIdeal) : ordFin v (y : L) = 0 := by
  have hyL : (y : L) ≠ 0 := RingOfIntegers.coe_ne_zero_iff.mpr hy
  have hnorm : (NumberField.FinitePlace.mk v) ((y : L)) = 1 := by
    rw [NumberField.FinitePlace.mk_apply]
    exact le_antisymm (NumberField.FinitePlace.norm_le_one L v y)
      (not_lt.mp ((NumberField.FinitePlace.norm_lt_one_iff_mem L v y).not.mpr hnot))
  have h := finitePlace_eq_absNorm_zpow v (y : L) hyL
  rw [hnorm] at h
  have hb := one_lt_absNorm_real (L := L) v
  have hinj := zpow_right_injective₀ (by linarith : (0:ℝ) < (Ideal.absNorm v.asIdeal : ℝ))
    (by linarith : (Ideal.absNorm v.asIdeal : ℝ) ≠ 1)
    (a₁ := -(ordFin v (y:L))) (a₂ := 0) (by simpa using h.symm)
  omega

/-- ★`y ∈ 𝔭` なら `ord_𝔭(y) > 0`。 -/
theorem ordFin_pos_of_mem (v : HeightOneSpectrum (𝓞 L)) (y : 𝓞 L) (hy : y ≠ 0)
    (hmem : y ∈ v.asIdeal) : 0 < ordFin v (y : L) := by
  have hyL : (y : L) ≠ 0 := RingOfIntegers.coe_ne_zero_iff.mpr hy
  have hnorm : (NumberField.FinitePlace.mk v) ((y : L)) < 1 := by
    rw [NumberField.FinitePlace.mk_apply]
    exact (NumberField.FinitePlace.norm_lt_one_iff_mem L v y).mpr hmem
  have h := finitePlace_eq_absNorm_zpow v (y : L) hyL
  rw [h] at hnorm
  have hb := one_lt_absNorm_real (L := L) v
  have := (zpow_lt_one_iff_right₀ hb).mp hnorm
  omega

/-- ★★★**1 つの有限素点だけを台にもつ主因子** ——
`𝔭^h`(`h` = 類数)が単項であることから。 -/
theorem exists_ordFin_pos_zero_elsewhere (𝔭 : HeightOneSpectrum (𝓞 L)) :
    ∃ x : Lˣ, 0 < ordFin 𝔭 (x : L) ∧
      ∀ 𝔮 : HeightOneSpectrum (𝓞 L), 𝔮 ≠ 𝔭 → ordFin 𝔮 (x : L) = 0 := by
  obtain ⟨y, hy0, hspan⟩ := exists_gen_pow_classNumber (L := L) 𝔭
  have hh : classNumber L ≠ 0 := (NumberField.classNumber_pos L).ne'
  have hyL : (y : L) ≠ 0 := RingOfIntegers.coe_ne_zero_iff.mpr hy0
  have hmem : y ∈ 𝔭.asIdeal := by
    have h1 : Ideal.span {y} ≤ 𝔭.asIdeal := by
      rw [hspan]
      exact Ideal.pow_le_self hh
    exact h1 (Ideal.mem_span_singleton_self y)
  refine ⟨Units.mk0 (y : L) hyL, ordFin_pos_of_mem 𝔭 y hy0 hmem, fun 𝔮 hne => ?_⟩
  refine ordFin_eq_zero_of_notMem 𝔮 y hy0 (fun hmem2 => hne ?_)
  have h2 : 𝔭.asIdeal ^ (classNumber L) ≤ 𝔮.asIdeal := by
    rw [← hspan]
    exact (Ideal.span_singleton_le_iff_mem _).mpr hmem2
  have h3 : 𝔮.asIdeal ∣ 𝔭.asIdeal ^ (classNumber L) := Ideal.dvd_iff_le.mpr h2
  haveI := 𝔮.isPrime
  have h4 : Prime 𝔮.asIdeal := (Ideal.prime_iff_isPrime (by simpa using 𝔮.ne_bot)).mpr 𝔮.isPrime
  have h5 : 𝔮.asIdeal ∣ 𝔭.asIdeal := h4.dvd_of_dvd_pow h3
  have h6 : 𝔭.asIdeal ≤ 𝔮.asIdeal := Ideal.dvd_iff_le.mp h5
  haveI := 𝔭.isMaximal
  exact HeightOneSpectrum.ext (((𝔭.isMaximal).eq_of_le (𝔮.isPrime.ne_top) h6).symm)

/-! ## ★3. 単数の因子 -/

/-- ★`(𝓞 L)ˣ → Lˣ`。 -/
noncomputable def unitToLUnits (u : (𝓞 L)ˣ) : Lˣ :=
  Units.map (algebraMap (𝓞 L) L).toMonoidHom u

@[simp] theorem unitToLUnits_val (u : (𝓞 L)ˣ) : ((unitToLUnits u : Lˣ) : L) = (u : 𝓞 L) := rfl

theorem unit_notMem (v : HeightOneSpectrum (𝓞 L)) (u : (𝓞 L)ˣ) : (u : 𝓞 L) ∉ v.asIdeal := by
  intro hmem
  exact v.isPrime.ne_top (Ideal.eq_top_of_isUnit_mem _ hmem u.isUnit)

/-- ★★単数の因子は無限素点にしか台がない。 -/
theorem arithDiv_unit_finite_zero (v : HeightOneSpectrum (𝓞 L)) (u : (𝓞 L)ˣ) :
    arithDiv (unitToLUnits u) (Sum.inl (NumberField.FinitePlace.mk v)) = 0 := by
  have hne : ((u : 𝓞 L) : L) ≠ 0 := RingOfIntegers.coe_ne_zero_iff.mpr u.ne_zero
  show arithPlaceLog ((unitToLUnits u : Lˣ) : L) (Sum.inl (NumberField.FinitePlace.mk v)) = 0
  rw [unitToLUnits_val, arithPlaceLog_finite v _ hne,
    ordFin_eq_zero_of_notMem v (u : 𝓞 L) u.ne_zero (unit_notMem v u)]
  simp

/-! ## ★4. `Θ` —— 無限素点への射影 -/

/-- ★無限素点にしか台がない因子。 -/
def infSupported (L : Type*) [Field L] [NumberField L] :
    Submodule ℝ (ArithPlace L →₀ ℝ) :=
  Finsupp.supported ℝ ℝ (Set.range (Sum.inr : InfinitePlace L → ArithPlace L))

/-- ★`Θ D w = −D(inr w)`(`w ≠ w₀`)。 -/
noncomputable def toLogSpace :
    (ArithPlace L →₀ ℝ) →ₗ[ℝ] NumberField.Units.dirichletUnitTheorem.logSpace L where
  toFun D := fun w => -(D (Sum.inr (w : InfinitePlace L)))
  map_add' D E := by funext w; simp; ring
  map_smul' c D := by funext w; simp

@[simp] theorem toLogSpace_apply (D : ArithPlace L →₀ ℝ)
    (w : {w : InfinitePlace L // w ≠ w₀}) :
    toLogSpace D w = -(D (Sum.inr (w : InfinitePlace L))) := rfl

/-- ★★**`Θ(div u) = logEmbedding u`**(単数について)—— **符号違いだけ**。 -/
theorem toLogSpace_arithDiv_unit (u : (𝓞 L)ˣ) :
    toLogSpace (arithDiv (unitToLUnits u))
      = NumberField.Units.logEmbedding L (Additive.ofMul u) := by
  funext w
  rw [toLogSpace_apply, logEmbedding_component]
  show -(arithPlaceLog ((unitToLUnits u : Lˣ) : L) (Sum.inr (w : InfinitePlace L))) = _
  show -(-((w : InfinitePlace L).mult : ℝ) *
    Real.log ((w : InfinitePlace L) ((unitToLUnits u : Lˣ) : L))) = _
  rw [unitToLUnits_val]
  ring

/-- ★★**`Θ` は「無限素点にしか台がない次数 0 の因子」の上で単射**。

★`w₀` の成分は次数 0 で決まる。 -/
theorem eq_zero_of_inf_deg_theta (E : ArithPlace L →₀ ℝ) (hW : E ∈ infSupported L)
    (hdeg : arithDegree E = 0) (hΘ : toLogSpace E = 0) : E = 0 := by
  have hWs : (↑E.support : Set (ArithPlace L))
      ⊆ Set.range (Sum.inr : InfinitePlace L → ArithPlace L) :=
    (Finsupp.mem_supported ℝ E).mp hW
  have h1 : ∀ w : InfinitePlace L, w ≠ w₀ → E (Sum.inr w) = 0 := by
    intro w hw
    have h := congrFun hΘ ⟨w, hw⟩
    have h2 : -(E (Sum.inr w)) = 0 := h
    linarith
  have h2 : ∀ w : FinitePlace L, E (Sum.inl w) = 0 := by
    intro w
    by_contra hne
    obtain ⟨w', hw'⟩ := hWs (Finsupp.mem_support_iff.mpr hne)
    cases hw'
  have hsub : E.support ⊆ {Sum.inr (w₀ : InfinitePlace L)} := by
    intro v hv
    rw [Finset.mem_singleton]
    rcases v with w | w
    · exact absurd (Finsupp.mem_support_iff.mp hv) (by rw [h2 w]; simp)
    · by_cases hw : w = w₀
      · rw [hw]
      · exact absurd (Finsupp.mem_support_iff.mp hv) (by rw [h1 w hw]; simp)
  have hval : E (Sum.inr (w₀ : InfinitePlace L)) = 0 := by
    have hd : arithDegree E
        = ∑ v ∈ ({Sum.inr (w₀ : InfinitePlace L)} : Finset (ArithPlace L)), E v :=
      Finsupp.sum_of_support_subset _ hsub _ (fun _ _ => rfl)
    rw [Finset.sum_singleton] at hd
    rw [← hd, hdeg]
  refine Finsupp.ext fun v => ?_
  rcases v with w | w
  · exact h2 w
  · by_cases hw : w = w₀
    · rw [hw]; exact hval
    · exact h1 w hw

/-! ## ★5. 主因子の張る部分空間 -/

/-- ★★**算術次数の線形版**。 -/
noncomputable def arithDegreeLin : (ArithPlace L →₀ ℝ) →ₗ[ℝ] ℝ :=
  Finsupp.lsum ℝ (fun _ : ArithPlace L => LinearMap.id)

@[simp] theorem arithDegreeLin_apply (D : ArithPlace L →₀ ℝ) :
    arithDegreeLin D = arithDegree D := rfl

/-- ★主因子の張る部分空間(`Φ^birat`)。 -/
noncomputable def principalSpan (L : Type*) [Field L] [NumberField L] :
    Submodule ℝ (ArithPlace L →₀ ℝ) :=
  Submodule.span ℝ (Set.range (fun x : Lˣ => arithDiv x))

/-- ★単数の主因子の張る部分空間。 -/
noncomputable def unitSpan (L : Type*) [Field L] [NumberField L] :
    Submodule ℝ (ArithPlace L →₀ ℝ) :=
  Submodule.span ℝ (Set.range (fun u : (𝓞 L)ˣ => arithDiv (unitToLUnits u)))

theorem arithDiv_unit_mem_infSupported (u : (𝓞 L)ˣ) :
    arithDiv (unitToLUnits u) ∈ infSupported L := by
  refine (Finsupp.mem_supported ℝ _).mpr ?_
  intro v hv
  rcases v with w | w
  · exfalso
    rw [Finset.mem_coe, Finsupp.mem_support_iff] at hv
    refine hv ?_
    have := arithDiv_unit_finite_zero (L := L) w.maximalIdeal u
    rwa [NumberField.FinitePlace.mk_maximalIdeal] at this
  · exact ⟨w, rfl⟩

theorem unitSpan_le_infSupported : unitSpan L ≤ infSupported L := by
  rw [unitSpan]
  refine Submodule.span_le.mpr ?_
  rintro _ ⟨u, rfl⟩
  exact arithDiv_unit_mem_infSupported u

theorem unitSpan_le_principalSpan : unitSpan L ≤ principalSpan L := by
  rw [unitSpan, principalSpan]
  refine Submodule.span_le.mpr ?_
  rintro _ ⟨u, rfl⟩
  exact Submodule.subset_span ⟨unitToLUnits u, rfl⟩

/-- ★★★**Dirichlet 単数定理** —— `Θ(単数の主因子の張る空間) = ⊤`。 -/
theorem map_toLogSpace_unitSpan : Submodule.map toLogSpace (unitSpan L) = ⊤ := by
  rw [unitSpan, Submodule.map_span]
  have himg : (toLogSpace : (ArithPlace L →₀ ℝ) →ₗ[ℝ] _) ''
      (Set.range (fun u : (𝓞 L)ˣ => arithDiv (unitToLUnits u)))
      = Set.range (fun u : (𝓞 L)ˣ => NumberField.Units.logEmbedding L (Additive.ofMul u)) := by
    ext g
    constructor
    · rintro ⟨_, ⟨u, rfl⟩, rfl⟩
      exact ⟨u, (toLogSpace_arithDiv_unit u).symm⟩
    · rintro ⟨u, rfl⟩
      exact ⟨_, ⟨u, rfl⟩, toLogSpace_arithDiv_unit u⟩
  rw [himg]
  have hlat : (NumberField.Units.unitLattice L : Set (logSpace L))
      = Set.range (fun u : (𝓞 L)ˣ => NumberField.Units.logEmbedding L (Additive.ofMul u)) := by
    ext g
    constructor
    · rintro ⟨x, -, rfl⟩
      exact ⟨Additive.toMul x, rfl⟩
    · rintro ⟨u, rfl⟩
      exact ⟨Additive.ofMul u, trivial, rfl⟩
  rw [← hlat]
  exact unitLattice_span_eq_top L

/-- ★★**主因子は次数 0**(積公式)。 -/
theorem principalSpan_le_ker : principalSpan L ≤ LinearMap.ker arithDegreeLin := by
  rw [principalSpan]
  refine Submodule.span_le.mpr ?_
  rintro _ ⟨x, rfl⟩
  simp only [SetLike.mem_coe, LinearMap.mem_ker, arithDegreeLin_apply]
  exact arithDegree_arithDiv x

/-- ★★★★**無限素点にしか台がない次数 0 の因子は単数の主因子で張られる**(Dirichlet)。 -/
theorem mem_unitSpan_of_inf_of_deg_zero (E : ArithPlace L →₀ ℝ) (hW : E ∈ infSupported L)
    (hdeg : arithDegree E = 0) : E ∈ unitSpan L := by
  have h1 : toLogSpace E ∈ Submodule.map toLogSpace (unitSpan L) := by
    rw [map_toLogSpace_unitSpan]; trivial
  obtain ⟨s, hs, hΘs⟩ := h1
  have hsW : s ∈ infSupported L := unitSpan_le_infSupported hs
  have hsdeg : arithDegree s = 0 :=
    principalSpan_le_ker (unitSpan_le_principalSpan hs)
  have hzero : s - E = 0 := by
    refine eq_zero_of_inf_deg_theta (s - E) (Submodule.sub_mem _ hsW hW) ?_ ?_
    · have hh : arithDegreeLin (s - E) = arithDegreeLin s - arithDegreeLin E := map_sub _ _ _
      simpa [arithDegreeLin_apply, hsdeg, hdeg] using hh
    · rw [map_sub, hΘs, sub_self]
  rw [← sub_eq_zero.mp hzero]
  exact hs

/-! ## ★6. 有限素点を消す -/

theorem arithDiv_inl (x : Lˣ) (w : FinitePlace L) :
    arithDiv x (Sum.inl w)
      = (ordFin w.maximalIdeal (x : L) : ℝ)
        * Real.log (Ideal.absNorm w.maximalIdeal.asIdeal : ℝ) := by
  have h := arithPlaceLog_finite w.maximalIdeal (x : L) x.ne_zero
  rw [NumberField.FinitePlace.mk_maximalIdeal] at h
  exact h

theorem maximalIdeal_injective {w w' : FinitePlace L} (h : w.maximalIdeal = w'.maximalIdeal) :
    w = w' := by
  rw [← NumberField.FinitePlace.mk_maximalIdeal w, ← NumberField.FinitePlace.mk_maximalIdeal w', h]

/-- ★★★**有限素点の `δ` は「主因子の空間 + 無限素点の空間」に入る**。 -/
theorem single_inl_mem_sup (w : FinitePlace L) (c : ℝ) :
    (Finsupp.single (Sum.inl w) c : ArithPlace L →₀ ℝ)
      ∈ principalSpan L ⊔ infSupported L := by
  obtain ⟨x, hpos, hzero⟩ := exists_ordFin_pos_zero_elsewhere (L := L) w.maximalIdeal
  set c₀ : ℝ := (ordFin w.maximalIdeal (x : L) : ℝ)
    * Real.log (Ideal.absNorm w.maximalIdeal.asIdeal : ℝ) with hc₀
  have hc₀pos : 0 < c₀ := by
    refine mul_pos ?_ (log_absNorm_pos (L := L) w.maximalIdeal)
    exact_mod_cast hpos
  have hc₀ne : c₀ ≠ 0 := hc₀pos.ne'
  have hzero' : ∀ w' : FinitePlace L, w' ≠ w → arithDiv x (Sum.inl w') = 0 := by
    intro w' hw'
    rw [arithDiv_inl, hzero w'.maximalIdeal (fun h => hw' (maximalIdeal_injective h))]
    simp
  have hrest : arithDiv x - (Finsupp.single (Sum.inl w) c₀ : ArithPlace L →₀ ℝ)
      ∈ infSupported L := by
    refine (Finsupp.mem_supported ℝ _).mpr ?_
    intro v hv
    rcases v with w' | w'
    · exfalso
      rw [Finset.mem_coe, Finsupp.mem_support_iff] at hv
      refine hv ?_
      by_cases hw' : w' = w
      · subst hw'
        rw [Finsupp.sub_apply, Finsupp.single_eq_same, arithDiv_inl, ← hc₀, sub_self]
      · rw [Finsupp.sub_apply, Finsupp.single_eq_of_ne (by simp [hw']), hzero' w' hw', sub_zero]
    · exact ⟨w', rfl⟩
  have hmem0 : (Finsupp.single (Sum.inl w) c₀ : ArithPlace L →₀ ℝ)
      ∈ principalSpan L ⊔ infSupported L := by
    have h1 : arithDiv x ∈ principalSpan L ⊔ infSupported L :=
      Submodule.mem_sup_left (Submodule.subset_span ⟨x, rfl⟩)
    have h2 : arithDiv x - (Finsupp.single (Sum.inl w) c₀ : ArithPlace L →₀ ℝ)
        ∈ principalSpan L ⊔ infSupported L := Submodule.mem_sup_right hrest
    have h3 := Submodule.sub_mem _ h1 h2
    simpa using h3
  have hsm : (Finsupp.single (Sum.inl w) c : ArithPlace L →₀ ℝ)
      = (c / c₀) • (Finsupp.single (Sum.inl w) c₀ : ArithPlace L →₀ ℝ) := by
    rw [Finsupp.smul_single]
    congr 1
    rw [smul_eq_mul]
    field_simp
  rw [hsm]
  exact Submodule.smul_mem _ _ hmem0

/-! ## ★7. 主定理 -/

theorem sup_principalSpan_infSupported : principalSpan L ⊔ infSupported L = ⊤ := by
  rw [eq_top_iff, ← Finsupp.span_single_eq_top (R := ℝ) (M := ℝ) (α := ArithPlace L)]
  refine Submodule.span_le.mpr ?_
  rintro D ⟨v, c, rfl⟩
  rcases v with w | w
  · exact single_inl_mem_sup w c
  · refine Submodule.mem_sup_right ((Finsupp.mem_supported ℝ _).mpr ?_)
    intro u hu
    have h := Finsupp.support_single_subset hu
    rw [Finset.mem_singleton] at h
    exact ⟨w, h ▸ rfl⟩

/-- ★★★★★★**[FrdI] Theorem 6.4, (i)** —— 主因子の張る空間は次数 0 の因子の全体。

原文 (FrdI p.114):
> isomorphism of groups

★★これが `δ_A : Pic_Φ(A) ≅ ℝ` の中身である。 -/
theorem principalSpan_eq_ker : principalSpan L = LinearMap.ker arithDegreeLin := by
  refine le_antisymm principalSpan_le_ker ?_
  intro D hD
  have htop : D ∈ principalSpan L ⊔ infSupported L := by
    rw [sup_principalSpan_infSupported]; trivial
  obtain ⟨s, hs, e, he, rfl⟩ := Submodule.mem_sup.mp htop
  have hsdeg : arithDegree s = 0 := principalSpan_le_ker hs
  have hedeg : arithDegree e = 0 := by
    have h1 : arithDegreeLin (s + e) = 0 := LinearMap.mem_ker.mp hD
    rw [map_add] at h1
    simp only [arithDegreeLin_apply] at h1 hsdeg ⊢
    linarith
  exact Submodule.add_mem _ hs
    (unitSpan_le_principalSpan (mem_unitSpan_of_inf_of_deg_zero e he hedeg))

theorem arithDegreeLin_surjective :
    Function.Surjective (arithDegreeLin : (ArithPlace L →₀ ℝ) →ₗ[ℝ] ℝ) := by
  intro c
  obtain ⟨w⟩ : Nonempty (InfinitePlace L) := inferInstance
  refine ⟨Finsupp.single (Sum.inr w) c, ?_⟩
  rw [arithDegreeLin_apply, arithDegree_apply, Finsupp.sum_single_index]
  rfl

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)** —— `δ_A : Pic_Φ(A) ≅ ℝ`。 -/
noncomputable def arithPicIso :
    ((ArithPlace L →₀ ℝ) ⧸ principalSpan L) ≃ₗ[ℝ] ℝ :=
  (Submodule.quotEquivOfEq _ _ principalSpan_eq_ker).trans
    (arithDegreeLin.quotKerEquivOfSurjective arithDegreeLin_surjective)

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Theorem 6.4, (i)` の `δ_A : Pic_Φ(A) ≅ ℝ`。 -/
def arithPicIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — δ_A : Pic_Φ(A) ≅ ℝ(Dirichlet 単数定理)",
    sectionId := "frdi-thm-6-4" }

def arithPicIso.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "NumberField.Units.dirichletUnitTheorem.unitLattice_span_eq_top(Dirichlet 単数定理)"
      (.inMathlib "NumberField.Units.dirichletUnitTheorem.unitLattice_span_eq_top") 114,
    .citation "[mathlib]" "NumberField.classNumber(類数の有限性)"
      (.inMathlib "NumberField.classNumber") 114,
    .citation "[ABC3]" "arithDegree_arithDiv(積公式)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDegree_arithDiv") 114,
    .implicitStep
      "★原文は「an immediate consequence of the well-known Dirichlet unit theorem」と 1 行で書く" 114 ]

end ABC3.Found.Divisor
