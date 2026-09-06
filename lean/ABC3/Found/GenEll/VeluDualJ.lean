/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluLatticePoint
import ABC3.Found.GenEll.VeluEllipticField
import ABC3.Found.GenEll.VeluEllipticNF
import ABC3.Found.GenEll.PointDescentFinite
import ABC3.Found.GenEll.VeluImage
import ABC3.Found.GaloisRep.TorsionCard
import ABC3.Found.GenEll.LatticeScale
import ABC3.Found.GenEll.LatticeDiscNonzero
import ABC3.Meta.Claim

/-!
# 第 1442 ブロック —— **双対同種の `j`（解析側）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——**対偶の道（道 C）の②**

`Skeleton/GenEll/VeluSemistable.lean` に残る `sorry` は

  `p ∣ l` かつ `0 ≤ jExp p E` のとき `0 ≤ jExp p (E/⟨Q⟩)`

である。★順方向（`v(c₄′) ≥ 4S`）は形式群が本体なので、**対偶**で行く:

    jExp p E′ < 0
      → E′ の `l`-巡回商は jExp < 0          …… ①
      → その商が E と同じ j を持つ（双対同種） …… ②
      → jExp p E < 0 で矛盾

★★★本ファイルは**②の解析側**である。すなわち `ℂ` 上で

    `Λ ⊂ Λ′` が指数 `l` なら、`Λ′ ⊂ (1/l)Λ` も指数 `l` であり、
    `(1/l)Λ` は `Λ` と相似だから `j` が等しい

を、**`veluQuotientFull` の語彙で**述べる。

## ★★得られるもの

| 定理 | 内容 |
|---|---|
| `latticeCurve_congr_lattice` | ☆束が同じなら曲線も同じ |
| `latticeCurve_j_scalePair` | ★相似な束は同じ `j` |
| `latticeDivPair` | ☆`(1/l)Λ` の周期対 |
| `lattice_dual_step` | ★★★★核心の格子計算（`Λ′ + ℤz₂ = (1/l)Λ`） |
| `exists_dual_elt` | ★★★★★双対の元がとれる |
| `latticeCurve_eq_veluQuotientFull_nsmul` | ★★`Λ′` を**指定して**の商の等式 |
| `exists_dual_veluQuotientFull_lattice` | ★★★★★★**②の格子側（`j` の一致）** |
| `exists_dual_veluQuot_complex` | ★★★★★★**②の `ℂ` 側** |
| `exists_rhPoint_of_nsmul_eq_zero` | ★★★代数閉体の間で `n`-捩れは全射 |
| `exists_dual_veluQuot_j_numberField` | ★★★★★★★★**②（数体の上、有限次拡大を経て）** |

## ★☆残っているもの——①（`jExp < 0` は `l`-巡回商で保たれる）

本ファイルは②だけである。対偶の道を閉じるには**①**が要る:

    jExp P X < 0 → jExp P (X/⟨Q′⟩) < 0

☆在庫の `semistableAt_veluQuot_badPrime_all`（第 1436）は
**`SemistableAt p X` を要求する**ので、そのままでは使えない
（`jExp < 0` だけでは加法還元でありうる——`c₄` が単元とは限らない）。
★①には「潜在的乗法還元は 2 次拡大で乗法還元になる」（`j` が同じ曲線は
2 次拡大で同型）が要る。**これは新しい節点である。**
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair Finset ABC3.Meta

open scoped Classical

/-! ## ☆束が同じなら曲線も同じ -/

/-- ☆`G` は束だけで決まる。 -/
theorem G_congr_lattice {P P' : PeriodPair} (h : P.lattice = P'.lattice) (n : ℕ) :
    P.G n = P'.G n := by
  have key : ∀ S T : Submodule ℤ ℂ, S = T →
      (∑' l : ↥S, ((l : ℂ) ^ n)⁻¹) = ∑' l : ↥T, ((l : ℂ) ^ n)⁻¹ := by
    rintro S T rfl; rfl
  exact key _ _ h

/-- ☆**束が同じなら `latticeCurve` も同じ**——`g₂`・`g₃` は束の関数だから。 -/
theorem latticeCurve_congr_lattice {P P' : PeriodPair} (h : P.lattice = P'.lattice) :
    latticeCurve P = latticeCurve P' := by
  simp only [latticeCurve, PeriodPair.g₂, PeriodPair.g₃, G_congr_lattice h]

/-! ## ★相似な束は同じ `j` -/

/-- ★**相似な束の曲線は同じ `j` を持つ**——`g₂³` も `Δ` も `c⁻¹²` 倍だから。 -/
theorem latticeCurve_j_scalePair (P : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    haveI := isElliptic_latticeCurve' (scalePair P c hc)
    haveI := isElliptic_latticeCurve' P
    (latticeCurve (scalePair P c hc)).j = (latticeCurve P).j := by
  haveI := isElliptic_latticeCurve' (scalePair P c hc)
  haveI := isElliptic_latticeCurve' P
  rw [latticeCurve_j, latticeCurve_j, g₂_scalePair, latticeDisc_scalePair]
  have h12 : (c ^ 12 : ℂ) ≠ 0 := pow_ne_zero 12 hc
  have h4 : (c ^ 4 : ℂ) ≠ 0 := pow_ne_zero 4 hc
  have hD : latticeDisc P ≠ 0 := latticeDisc_ne_zero P
  field_simp

/-! ## ☆`(1/l)Λ` の周期対 -/

/-- ☆**`(1/l)Λ` の周期対**——`scalePair` の `c = 1/l` の場合。 -/
noncomputable def latticeDivPair (P : PeriodPair) {l : ℕ} (hl : 0 < l) : PeriodPair :=
  scalePair P ((l : ℂ)⁻¹) (inv_ne_zero (Nat.cast_ne_zero.2 hl.ne'))

/-- ☆`(1/l)Λ` の元は `l` 倍すると `Λ` に入る（そしてその逆）。 -/
theorem mem_latticeDivPair_iff (P : PeriodPair) {l : ℕ} (hl : 0 < l) (x : ℂ) :
    x ∈ (latticeDivPair P hl).lattice ↔ (l : ℂ) * x ∈ P.lattice := by
  have hl0 : ((l : ℂ)) ≠ 0 := Nat.cast_ne_zero.2 hl.ne'
  rw [latticeDivPair, scalePair_lattice]
  constructor
  · rintro ⟨y, hy, rfl⟩
    have : (l : ℂ) * ((l : ℂ)⁻¹ * y) = y := by field_simp
    rw [LinearMap.mulLeft_apply, this]
    exact hy
  · intro h
    exact ⟨(l : ℂ) * x, h, by rw [LinearMap.mulLeft_apply]; field_simp⟩

/-- ☆`Λ ⊆ (1/l)Λ`。 -/
theorem le_latticeDivPair (P : PeriodPair) {l : ℕ} (hl : 0 < l) :
    P.lattice ≤ (latticeDivPair P hl).lattice := by
  intro x hx
  rw [mem_latticeDivPair_iff]
  have : (l : ℂ) * x = (l : ℤ) • x := by push_cast [zsmul_eq_mul]; ring
  rw [this]
  exact P.lattice.smul_mem _ hx

/-! ## ★★核心の格子計算 -/

/-- ☆`ℝ`-独立なら `ℤ` 係数の関係式も自明。 -/
theorem intCast_eq_zero_of_indep {u v : ℂ} (h : LinearIndependent ℝ ![u, v]) {p q : ℤ}
    (hpq : (p : ℂ) * u + (q : ℂ) * v = 0) : p = 0 ∧ q = 0 := by
  have hr : ((p : ℝ)) • u + ((q : ℝ)) • v = 0 := by
    simp only [Complex.real_smul]
    push_cast
    exact hpq
  obtain ⟨h1, h2⟩ := LinearIndependent.pair_iff.1 h _ _ hr
  exact ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩

/-- ★★★★**双対の一段**——`u, v` が `ℝ`-独立、`z₀ = a·u + b·v` で `l ∤ a` のとき

* `v ∉ ⟨l·u, l·v⟩ + ℤz₀`（したがって `v` の像は位数 `l` である）
* `(⟨l·u, l·v⟩ + ℤz₀) + ℤv = ⟨u, v⟩`

☆前者は `Λ′` の外にあること、後者は「もう一段上げると `(1/l)Λ` になる」ことである。 -/
theorem lattice_dual_step {u v : ℂ} (hindep : LinearIndependent ℝ ![u, v])
    {l : ℕ} (hl : l.Prime) {a b : ℤ} (ha : ¬ ((l : ℤ) ∣ a)) {z₀ : ℂ}
    (hz₀ : z₀ = (a : ℂ) * u + (b : ℂ) * v) :
    v ∉ Submodule.span ℤ ({(l : ℂ) * u, (l : ℂ) * v} : Set ℂ)
        ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (Submodule.span ℤ ({(l : ℂ) * u, (l : ℂ) * v} : Set ℂ)
          ⊔ Submodule.span ℤ ({z₀} : Set ℂ)) ⊔ Submodule.span ℤ ({v} : Set ℂ)
        = Submodule.span ℤ ({u, v} : Set ℂ) := by
  have hlprime : Prime ((l : ℕ) : ℤ) := Nat.prime_iff_prime_int.1 hl
  have hcop : IsCoprime ((l : ℕ) : ℤ) a := hlprime.coprime_iff_not_dvd.2 ha
  constructor
  · intro hmem
    rw [Submodule.mem_sup] at hmem
    obtain ⟨y, hy, w, hw, hsum⟩ := hmem
    obtain ⟨m, n, hmn⟩ := Submodule.mem_span_pair.1 hy
    obtain ⟨k, hk⟩ := Submodule.mem_span_singleton.1 hw
    have h1 : y = (m : ℂ) * ((l : ℂ) * u) + (n : ℂ) * ((l : ℂ) * v) := by
      rw [← hmn]; simp [zsmul_eq_mul]
    have h2 : w = (k : ℂ) * ((a : ℂ) * u + (b : ℂ) * v) := by
      rw [← hk, hz₀, zsmul_eq_mul]
    have hz : ((m * l + k * a : ℤ) : ℂ) * u + ((n * l + k * b - 1 : ℤ) : ℂ) * v = 0 := by
      push_cast
      linear_combination hsum - h1 - h2
    obtain ⟨e1, e2⟩ := intCast_eq_zero_of_indep hindep hz
    have hdvd : ((l : ℕ) : ℤ) ∣ k * a := ⟨-m, by linarith⟩
    have hlk : ((l : ℕ) : ℤ) ∣ k := (hlprime.dvd_mul.1 hdvd).resolve_right ha
    obtain ⟨c, rfl⟩ := hlk
    have hdvd1 : ((l : ℕ) : ℤ) ∣ 1 := ⟨n + c * b, by linarith⟩
    have hle : ((l : ℕ) : ℤ) ≤ 1 := Int.le_of_dvd one_pos hdvd1
    have := hl.two_le
    omega
  · apply le_antisymm
    · refine sup_le (sup_le ?_ ?_) ?_
      · rw [Submodule.span_le, Set.insert_subset_iff, Set.singleton_subset_iff]
        exact ⟨Submodule.mem_span_pair.2 ⟨(l : ℤ), 0, by simp [zsmul_eq_mul]⟩,
          Submodule.mem_span_pair.2 ⟨0, (l : ℤ), by simp [zsmul_eq_mul]⟩⟩
      · rw [Submodule.span_le, Set.singleton_subset_iff, hz₀]
        exact Submodule.mem_span_pair.2 ⟨a, b, by simp [zsmul_eq_mul]⟩
      · rw [Submodule.span_le, Set.singleton_subset_iff]
        exact Submodule.subset_span (by simp)
    · set M := (Submodule.span ℤ ({(l : ℂ) * u, (l : ℂ) * v} : Set ℂ)
        ⊔ Submodule.span ℤ ({z₀} : Set ℂ)) ⊔ Submodule.span ℤ ({v} : Set ℂ) with hM
      have hvM : v ∈ M := by
        rw [hM]
        exact Submodule.mem_sup_right (Submodule.subset_span (by simp))
      have hz₀M : z₀ ∈ M := by
        rw [hM]
        exact Submodule.mem_sup_left (Submodule.mem_sup_right
          (Submodule.subset_span (by simp)))
      have hluM : (l : ℂ) * u ∈ M := by
        rw [hM]
        exact Submodule.mem_sup_left (Submodule.mem_sup_left
          (Submodule.subset_span (by simp)))
      have haM : (a : ℂ) * u ∈ M := by
        have hbv : (b : ℤ) • v ∈ M := M.smul_mem _ hvM
        have : (a : ℂ) * u = z₀ - (b : ℤ) • v := by
          rw [hz₀, zsmul_eq_mul]; ring
        rw [this]
        exact M.sub_mem hz₀M hbv
      have huM : u ∈ M := by
        obtain ⟨s, t, hst⟩ := hcop
        have hstC : (s : ℂ) * ((l : ℕ) : ℂ) + (t : ℂ) * (a : ℂ) = 1 := by
          have := congrArg (fun z : ℤ => (z : ℂ)) hst
          push_cast at this
          exact this
        have hu : (s : ℤ) • ((l : ℂ) * u) + (t : ℤ) • ((a : ℂ) * u) = u := by
          simp only [zsmul_eq_mul]
          linear_combination u * hstC
        rw [← hu]
        exact M.add_mem (M.smul_mem _ hluM) (M.smul_mem _ haM)
      rw [Submodule.span_le, Set.insert_subset_iff, Set.singleton_subset_iff]
      exact ⟨huM, hvM⟩

/-- ☆`ℝ`-独立は入れ替えても保たれる。 -/
theorem linearIndependent_pair_swap {u v : ℂ} (h : LinearIndependent ℝ ![u, v]) :
    LinearIndependent ℝ ![v, u] := by
  refine LinearIndependent.pair_iff.2 (fun s t hst => ?_)
  obtain ⟨h1, h2⟩ := LinearIndependent.pair_iff.1 h t s (by rw [add_comm]; exact hst)
  exact ⟨h2, h1⟩

/-- ★★★★★★**双対の元がとれる**——`Λ` の外の `l`-捩れ元 `z₀` に対し、

* `z₂ ∉ Λ + ℤz₀`
* `l·z₂ ∈ Λ`
* `(Λ + ℤz₀) + ℤz₂ = (1/l)Λ`

なる `z₂` がある。☆これが「もう一段上げると `(1/l)Λ` になる」ことである。 -/
theorem exists_dual_elt (P : PeriodPair) {l : ℕ} (hl : l.Prime) {z₀ : ℂ}
    (hlz : (l : ℂ) * z₀ ∈ P.lattice) (hz₀ : z₀ ∉ P.lattice) :
    ∃ z₂ : ℂ, z₂ ∉ P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (l : ℂ) * z₂ ∈ P.lattice ∧
      (P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ)) ⊔ Submodule.span ℤ ({z₂} : Set ℂ)
        = (latticeDivPair P hl.pos).lattice := by
  have hl0 : ((l : ℂ)) ≠ 0 := Nat.cast_ne_zero.2 hl.pos.ne'
  have hlw₁ : (l : ℂ) * (((l : ℂ))⁻¹ * P.ω₁) = P.ω₁ := by field_simp
  have hlw₂ : (l : ℂ) * (((l : ℂ))⁻¹ * P.ω₂) = P.ω₂ := by field_simp
  have hlatP : P.lattice
      = Submodule.span ℤ ({(l : ℂ) * (((l : ℂ))⁻¹ * P.ω₁),
          (l : ℂ) * (((l : ℂ))⁻¹ * P.ω₂)} : Set ℂ) := by
    rw [hlw₁, hlw₂]; rfl
  have hlatPs : P.lattice
      = Submodule.span ℤ ({(l : ℂ) * (((l : ℂ))⁻¹ * P.ω₂),
          (l : ℂ) * (((l : ℂ))⁻¹ * P.ω₁)} : Set ℂ) := by
    rw [hlw₁, hlw₂, Set.pair_comm]; rfl
  have hlatD : (latticeDivPair P hl.pos).lattice
      = Submodule.span ℤ ({((l : ℂ))⁻¹ * P.ω₁, ((l : ℂ))⁻¹ * P.ω₂} : Set ℂ) := rfl
  have hlatDs : (latticeDivPair P hl.pos).lattice
      = Submodule.span ℤ ({((l : ℂ))⁻¹ * P.ω₂, ((l : ℂ))⁻¹ * P.ω₁} : Set ℂ) := by
    rw [hlatD, Set.pair_comm]
  have hindep : LinearIndependent ℝ ![((l : ℂ))⁻¹ * P.ω₁, ((l : ℂ))⁻¹ * P.ω₂] :=
    (latticeDivPair P hl.pos).indep
  have hω₁ : P.ω₁ ∈ P.lattice := PeriodPair.mem_lattice.2 ⟨1, 0, by push_cast; ring⟩
  have hω₂ : P.ω₂ ∈ P.lattice := PeriodPair.mem_lattice.2 ⟨0, 1, by push_cast; ring⟩
  obtain ⟨a, b, hab⟩ := PeriodPair.mem_lattice.1 hlz
  have hz₀eq : z₀ = (a : ℂ) * (((l : ℂ))⁻¹ * P.ω₁) + (b : ℂ) * (((l : ℂ))⁻¹ * P.ω₂) := by
    refine mul_left_cancel₀ hl0 ?_
    field_simp
    linear_combination -hab
  by_cases ha : ((l : ℤ) ∣ a)
  · -- ☆`l ∣ a` なら `l ∤ b`（さもないと `z₀ ∈ Λ`）
    have hb : ¬ ((l : ℤ) ∣ b) := by
      intro hb'
      apply hz₀
      obtain ⟨a', rfl⟩ := ha
      obtain ⟨b', rfl⟩ := hb'
      refine PeriodPair.mem_lattice.2 ⟨a', b', ?_⟩
      have h := hab
      push_cast at h
      refine mul_left_cancel₀ hl0 ?_
      linear_combination h
    obtain ⟨hnot, hsup⟩ := lattice_dual_step (u := ((l : ℂ))⁻¹ * P.ω₂)
      (v := ((l : ℂ))⁻¹ * P.ω₁) (linearIndependent_pair_swap hindep) hl hb (z₀ := z₀)
      (by rw [hz₀eq]; ring)
    refine ⟨((l : ℂ))⁻¹ * P.ω₁, ?_, ?_, ?_⟩
    · rwa [hlatPs]
    · rw [hlw₁]; exact hω₁
    · rw [hlatPs, hlatDs]; exact hsup
  · obtain ⟨hnot, hsup⟩ := lattice_dual_step (u := ((l : ℂ))⁻¹ * P.ω₁)
      (v := ((l : ℂ))⁻¹ * P.ω₂) hindep hl ha (z₀ := z₀) hz₀eq
    refine ⟨((l : ℂ))⁻¹ * P.ω₂, ?_, ?_, ?_⟩
    · rwa [hlatP]
    · rw [hlw₂]; exact hω₂
    · rw [hlatP, hlatD]; exact hsup

/-! ## ★★`Λ′` を指定しての商の等式 -/

/-- ☆`Λ` の外の `l`-捩れ元の像は位数ちょうど `l`（`l` は素数）。 -/
theorem addOrderOf_uniformMap_eq (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {l : ℕ}
    (hl : l.Prime) {z : ℂ} (hlz : (l : ℂ) * z ∈ P.lattice) (hz : z ∉ P.lattice) :
    addOrderOf (uniformMap P hΔ z) = l := by
  have hzero : l • uniformMap P hΔ z = 0 := by
    have hmap : l • uniformMap P hΔ z = uniformMap P hΔ ((l : ℕ) • z) := by
      rw [← uniformHom_apply, ← uniformHom_apply, map_nsmul]
    rw [hmap, uniformMap_eq_zero_iff]
    have hc : (l : ℕ) • z = (l : ℂ) * z := by rw [nsmul_eq_mul]
    rw [hc]
    exact hlz
  have hne : uniformMap P hΔ z ≠ 0 := by
    rw [Ne, uniformMap_eq_zero_iff]
    exact hz
  have hdvd := addOrderOf_dvd_of_nsmul_eq_zero hzero
  rcases hl.eq_one_or_self_of_dvd _ hdvd with h1 | h2
  · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h1) hne
  · exact h2

/-- ★★**`Λ′` を指定しての Vélu の商の等式**——
`isElliptic_veluQuotientFull_nsmul_lattice`（第 1333）の中身を、
楕円性ではなく**曲線の等式**として取り出した形である。 -/
theorem latticeCurve_eq_veluQuotientFull_nsmul (P P' : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q)
    (hP' : P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ)) :
    latticeCurve P' = veluQuotientFull (latticeCurve P)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) := by
  classical
  have hord := intCast_mul_mem_lattice_iff P hΔ hQ hz₀
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ := exists_velu_rep P P' z₀ l hl hP' hord
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  have heq := latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)
  have hnot : ∀ k ∈ (range l).erase 0, ((k : ℂ)) * z₀ ∉ P.lattice := by
    intro k hk hmem
    have hk0 : k ≠ 0 := Finset.ne_of_mem_erase hk
    have hkl : k < l := Finset.mem_range.1 (Finset.mem_of_mem_erase hk)
    have hmem' : ((k : ℤ) : ℂ) * z₀ ∈ P.lattice := by push_cast; exact hmem
    have hd := (hord (k : ℤ)).1 hmem'
    have hle' : (l : ℤ) ≤ (k : ℤ) :=
      Int.le_of_dvd (by exact_mod_cast Nat.pos_of_ne_zero hk0) hd
    omega
  have hTe : T.erase 0 = ((range l).erase 0).image (fun k : ℕ => (k : ℂ) * z₀) := by
    ext w
    simp only [hTdef, Finset.mem_erase, Finset.mem_image, Finset.mem_range]
    constructor
    · rintro ⟨hw0, k, hkl, rfl⟩
      refine ⟨k, ⟨?_, hkl⟩, rfl⟩
      rintro rfl
      exact hw0 (by simp)
    · rintro ⟨k, ⟨hk0, hkl⟩, rfl⟩
      refine ⟨?_, k, hkl, rfl⟩
      intro hz
      exact hnot k (Finset.mem_erase.2 ⟨hk0, Finset.mem_range.2 hkl⟩)
        (hz ▸ P.lattice.zero_mem)
  have hset : ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))
      = (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)) := by
    rw [hTe, Finset.image_image, ← hz₀]
    exact image_pointCoords_nsmul_eq P hΔ z₀ l hnot
  rw [hset]
  exact heq

/-! ## ★★★★★★双対同種の `j`（解析側） -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**双対同種の `j`（解析側）**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`ℂ` 上の格子曲線 `E = latticeCurve P` と位数 `l`（素数）の点 `Q` に対し、

* `E′ ≔ latticeCurve P′` は `E/⟨Q⟩`（Vélu の商）であり、
* `E′` の上に位数 `l` の点 `Q′` があって
* `E′/⟨Q′⟩ = latticeCurve P″` の `j` は `E` の `j` に等しい

☆これが**双対同種**である——`Λ ⊂ Λ′ ⊂ (1/l)Λ` で、`(1/l)Λ` は `Λ` と相似だから。 -/
theorem exists_dual_veluQuotientFull_lattice (P : PeriodPair)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : l.Prime) (hQ : addOrderOf Q = l) :
    ∃ (P' P'' : PeriodPair) (Q' : (latticeCurve P').toAffine.Point),
      latticeCurve P' = veluQuotientFull (latticeCurve P)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) ∧
        addOrderOf Q' = l ∧
        latticeCurve P'' = veluQuotientFull (latticeCurve P')
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))) ∧
        (latticeCurve P'').j = (latticeCurve P).j := by
  classical
  have hΔ : latticeDisc P ≠ 0 := latticeDisc_ne_zero P
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl.pos hQ
  have hΔ' : latticeDisc P' ≠ 0 := latticeDisc_ne_zero P'
  have hord := intCast_mul_mem_lattice_iff P hΔ hQ hz₀
  have hlz : (l : ℂ) * z₀ ∈ P.lattice := by
    have hd := (hord (l : ℤ)).2 dvd_rfl
    push_cast at hd
    exact hd
  have hz₀not : z₀ ∉ P.lattice := by
    intro hc
    have hd := (hord 1).1 (by push_cast; simpa using hc)
    have hle := Int.le_of_dvd one_pos hd
    have := hl.two_le
    omega
  obtain ⟨z₂, hz₂not, hlz₂, hsup⟩ := exists_dual_elt P hl hlz hz₀not
  have hlz₂' : (l : ℂ) * z₂ ∈ P'.lattice := by
    rw [hP']; exact Submodule.mem_sup_left hlz₂
  have hz₂not' : z₂ ∉ P'.lattice := by rw [hP']; exact hz₂not
  refine ⟨P', latticeDivPair P hl.pos, uniformMap P' hΔ' z₂,
    latticeCurve_eq_veluQuotientFull_nsmul P P' hΔ hl.pos hQ hz₀ hP',
    addOrderOf_uniformMap_eq P' hΔ' hl hlz₂' hz₂not', ?_, ?_⟩
  · refine latticeCurve_eq_veluQuotientFull_nsmul P' (latticeDivPair P hl.pos) hΔ' hl.pos
      (addOrderOf_uniformMap_eq P' hΔ' hl hlz₂' hz₂not') rfl ?_
    rw [hP']
    exact hsup.symm
  · exact latticeCurve_j_scalePair P ((l : ℂ))⁻¹ _

/-! ## ★★★★★★`ℂ` 上の任意の楕円曲線へ -/

/-- ☆曲線が等しければ `⟨Q⟩` による商も等しい（`isElliptic_velu_congr_curve` の等式版）。 -/
theorem veluQuot_congr_curve {F : Type} [Field F] {W₁ W₂ : WeierstrassCurve F}
    (h : W₁ = W₂) (Q : W₁.toAffine.Point) (l : ℕ) :
    veluQuotientFull W₁ (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))
      = veluQuotientFull W₂
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • (h ▸ Q : W₂.toAffine.Point)))) := by
  subst h; rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**双対同種の `j`（`ℂ` 上の任意の楕円曲線で）**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`ℂ` 上の楕円曲線 `V` と位数 `l`（素数）の点 `Q` に対し、
`V′ = V/⟨Q⟩` の上に位数 `l` の点 `Q′` があって `V′/⟨Q′⟩` の `j` は `V` の `j` に等しい。

☆結論を「`Z` がその商に等しいなら `Z.j = V.j`」の形で述べるのは、
商の `IsElliptic` インスタンスを存在束縛の中で作らないためである。 -/
theorem exists_dual_veluQuot_complex (V : WeierstrassCurve ℂ) [hV : V.IsElliptic]
    {Q : V.toAffine.Point} {l : ℕ} (hl : l.Prime) (hQ : addOrderOf Q = l)
    [hV' : (veluQuotientFull V
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic] :
    ∃ Q' : (veluQuotientFull V
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).toAffine.Point,
      addOrderOf Q' = l ∧
      ∀ (Z : WeierstrassCurve ℂ) [Z.IsElliptic],
        Z = veluQuotientFull (veluQuotientFull V
            (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))) →
        Z.j = V.j := by
  classical
  obtain ⟨P, C, hCP⟩ := exists_periodPair_of_isElliptic V hV
  haveI hCVell : (C • V).IsElliptic := by rw [hCP]; infer_instance
  have h2 : (2 : ℂ) ≠ 0 := two_ne_zero
  have hQC : addOrderOf (vcPoint C V Q) = l := by rw [addOrderOf_vcPoint, hQ]
  have hQP : addOrderOf ((hCP ▸ vcPoint C V Q : (latticeCurve P).toAffine.Point)) = l := by
    rw [addOrderOf_congr_curve, hQC]
  obtain ⟨P', P'', Q'ₗ, heq1, hQ'ₗ, heq2, hj⟩ :=
    exists_dual_veluQuotientFull_lattice P hl hQP
  have hCV' : C • (veluQuotientFull V
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
      = veluQuotientFull (C • V)
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C V Q))) :=
    veluQuotientFull_vcPoint_eq C V _ hQ h2 rfl
  have hP'eq : latticeCurve P'
      = C • (veluQuotientFull V
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
    rw [heq1, hCV', veluQuot_congr_curve hCP (vcPoint C V Q) l]
  obtain ⟨R, hR⟩ := vcPoint_surjective C (veluQuotientFull V
    (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hP'eq ▸ Q'ₗ : (C • (veluQuotientFull V
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))).toAffine.Point)
  haveI hCV'ell : (C • (veluQuotientFull V
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))).IsElliptic := by
    rw [← hP'eq]; infer_instance
  have hRord : addOrderOf R = l := by
    have hv := addOrderOf_vcPoint C (veluQuotientFull V
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) R
    rw [hR] at hv
    rw [← hv, addOrderOf_congr_curve hP'eq Q'ₗ, hQ'ₗ]
  refine ⟨R, hRord, ?_⟩
  intro Z hZ hZeq
  have hstep : C • Z = latticeCurve P'' := by
    rw [hZeq, veluQuotientFull_vcPoint_eq C _ _ hRord h2 rfl, hR,
      ← veluQuot_congr_curve hP'eq Q'ₗ l, ← heq2]
  haveI hZC : (C • Z).IsElliptic := by rw [hstep]; infer_instance
  calc Z.j = (C • Z).j := (WeierstrassCurve.variableChange_j Z C).symm
    _ = (latticeCurve P'').j := j_congr_curve hstep
    _ = (latticeCurve P).j := hj
    _ = (C • V).j := j_congr_curve hCP.symm
    _ = V.j := WeierstrassCurve.variableChange_j V C

/-! ## ★★数体への降下の道具 -/

/-- ☆`j` は環準同型で写る。 -/
theorem j_map {F K : Type} [Field F] [Field K] (φ : F →+* K) (W : WeierstrassCurve F)
    [W.IsElliptic] [(W.map φ).IsElliptic] : (W.map φ).j = φ W.j := by
  rw [j_eq_inv_Delta_mul, j_eq_inv_Delta_mul, WeierstrassCurve.map_Δ, WeierstrassCurve.map_c₄,
    map_mul, map_inv₀, map_pow]

/-- ★★★**代数閉体の間では `n`-捩れは `rhPoint` で全射**——★**無条件**。

☆`torsionMap_bijective`（第 1271）と同じ論法を `rhPoint` の語彙で書いた形である
（単射＋両側の個数が `n²`）。 -/
theorem exists_rhPoint_of_nsmul_eq_zero {F K : Type} [Field F] [Field K]
    [IsAlgClosed F] [IsAlgClosed K] (φ : F →+* K) (W : WeierstrassCurve F) [W.IsElliptic]
    [(W.map φ).IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (hcF : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hcK : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : K) ≠ 0)
    (P : (W.map φ).toAffine.Point) (hP : n • P = 0) :
    ∃ R : W.toAffine.Point, n • R = 0 ∧ rhPoint φ W R = P := by
  classical
  haveI hfF : Finite {R : W.toAffine.Point // n • R = 0} :=
    (ABC3.Found.GaloisRep.finite_torsion W n hn (hcF n hn le_rfl)).to_subtype
  haveI hfK : Finite {R : (W.map φ).toAffine.Point // n • R = 0} :=
    (ABC3.Found.GaloisRep.finite_torsion (W.map φ) n hn (hcK n hn le_rfl)).to_subtype
  have hbij : Function.Bijective
      (fun R : {R : W.toAffine.Point // n • R = 0} =>
        (⟨rhPoint φ W (R : W.toAffine.Point), by
          rw [← rhPoint_nsmul, R.2, rhPoint_zero]⟩ :
            {R : (W.map φ).toAffine.Point // n • R = 0})) := by
    refine (Nat.bijective_iff_injective_and_card _).2 ⟨?_, ?_⟩
    · intro a b hab
      exact Subtype.ext (rhPoint_injective φ W (Subtype.ext_iff.1 hab))
    · rw [ABC3.Found.GaloisRep.torsion_card W W.isUnit_Δ.ne_zero n hn hcF,
        ABC3.Found.GaloisRep.torsion_card (W.map φ) (W.map φ).isUnit_Δ.ne_zero n hn hcK]
  obtain ⟨R, hR⟩ := hbij.2 ⟨P, hP⟩
  exact ⟨(R : W.toAffine.Point), R.2, congrArg Subtype.val hR⟩

/-! ## ★★★★★★★★双対同種の `j`（数体の上、有限次拡大を経て） -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**双対同種の `j`（数体の上）**——★**無条件**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

数体 `L` 上の楕円曲線 `E` と位数 `l`（素数）の点 `Q`、`E′ = E/⟨Q⟩` に対し、
**`L` の有限次拡大 `M`** と `E′ ⊗ M` の位数 `l` の点 `Q′` があって

    j( (E′ ⊗ M)/⟨Q′⟩ ) = j( E ⊗ M )

☆すなわち**双対同種**が `M` の上で Vélu の形で書ける。

★道は 3 段:

1. `ℂ` へ埋め込んで一意化し、`(1/l)Λ` を取る（`exists_dual_veluQuot_complex`）
2. `ℂ` の `l`-捩れ点を `L̄` へ引き戻す（`exists_rhPoint_of_nsmul_eq_zero`）
3. `L̄` の点を有限次拡大へ降ろす（`exists_finite_point_descent`、第 1207）

☆結論を「`Z` がその商に等しいなら」の形で述べるのは、
商の `IsElliptic` インスタンスを存在束縛の中で作らないためである。 -/
theorem exists_dual_veluQuot_j_numberField (L : Type) [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic] {l : ℕ} (hl : l.Prime)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    ∃ M : IntermediateField L (AlgebraicClosure L), FiniteDimensional L M ∧
      letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
      ∃ Q' : (E'.baseChange (M : Type)).toAffine.Point, addOrderOf Q' = l ∧
        ∀ (Z : WeierstrassCurve (M : Type)) [Z.IsElliptic],
          Z = veluQuotientFull (E'.baseChange (M : Type))
            (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))) →
          Z.j = (E.map (algebraMap L (M : Type))).j := by
  classical
  letI : Algebra L ℂ := (embedComplex L).toAlgebra
  obtain ⟨f⟩ : Nonempty (AlgebraicClosure L →ₐ[L] ℂ) := ⟨IsAlgClosed.lift⟩
  -- ☆底変換の楕円性
  haveI hEbar : (E.baseChange (AlgebraicClosure L)).IsElliptic := by
    show (E.map (algebraMap L (AlgebraicClosure L))).IsElliptic; infer_instance
  haveI hE'bar : (E'.baseChange (AlgebraicClosure L)).IsElliptic := by
    show (E'.map (algebraMap L (AlgebraicClosure L))).IsElliptic; infer_instance
  haveI hEC : (E.map (algebraMap L ℂ)).IsElliptic := inferInstance
  haveI hE'C : (E'.map (algebraMap L ℂ)).IsElliptic := inferInstance
  -- ☆`L̄ → ℂ` は `L` の上の写像
  have hmapbar : ∀ W : WeierstrassCurve L,
      (W.baseChange (AlgebraicClosure L)).map (f : AlgebraicClosure L →+* ℂ)
        = W.map (algebraMap L ℂ) := by
    intro W
    have hcomp : (f : AlgebraicClosure L →+* ℂ).comp (algebraMap L (AlgebraicClosure L))
        = algebraMap L ℂ := RingHom.ext (fun x => f.commutes x)
    show (W.map (algebraMap L (AlgebraicClosure L))).map (f : AlgebraicClosure L →+* ℂ)
      = W.map (algebraMap L ℂ)
    rw [← hcomp]; rfl
  haveI hE'barC : ((E'.baseChange (AlgebraicClosure L)).map
      (f : AlgebraicClosure L →+* ℂ)).IsElliptic := by rw [hmapbar]; infer_instance
  -- ★段 1: `ℂ` の側で双対の点を作る
  have hQC : addOrderOf (rhPoint (algebraMap L ℂ) E Q) = l := by rw [addOrderOf_rhPoint, hQ]
  have hquot := veluQuotientFull_baseChange (algebraMap L ℂ) E E' hQ hE'
  haveI hquotEll : (veluQuotientFull (E.map (algebraMap L ℂ))
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • rhPoint (algebraMap L ℂ) E Q)))).IsElliptic := by
    rw [← hquot]; infer_instance
  obtain ⟨Q'C, hQ'Cord, hjC⟩ :=
    exists_dual_veluQuot_complex (E.map (algebraMap L ℂ)) hl hQC
  -- ☆`E′ ⊗ ℂ` の点に移す
  have hQ2ord : addOrderOf (castPoint hquot.symm Q'C) = l := by
    rw [addOrderOf_castPoint, hQ'Cord]
  have hQ3ord : addOrderOf (castPoint (hmapbar E').symm (castPoint hquot.symm Q'C)) = l := by
    rw [addOrderOf_castPoint, hQ2ord]
  have hQ3zero : l • (castPoint (hmapbar E').symm (castPoint hquot.symm Q'C)) = 0 := by
    have h := addOrderOf_nsmul_eq_zero (castPoint (hmapbar E').symm (castPoint hquot.symm Q'C))
    rwa [hQ3ord] at h
  -- ★段 2: `L̄` へ引き戻す
  have hchar : ∀ (F : Type) [Field F] [CharZero F], ∀ k : ℕ, 1 ≤ k → k ≤ l → (k : F) ≠ 0 := by
    intro F _ _ k hk _
    exact Nat.cast_ne_zero.2 (by omega)
  obtain ⟨R, hR0, hRmap⟩ := exists_rhPoint_of_nsmul_eq_zero (f : AlgebraicClosure L →+* ℂ)
    (E'.baseChange (AlgebraicClosure L)) l hl.pos (hchar (AlgebraicClosure L)) (hchar ℂ)
    (castPoint (hmapbar E').symm (castPoint hquot.symm Q'C)) hQ3zero
  have hRord : addOrderOf R = l := by
    rw [← addOrderOf_rhPoint (f : AlgebraicClosure L →+* ℂ) (E'.baseChange (AlgebraicClosure L)),
      hRmap, hQ3ord]
  -- ★段 3: 有限次拡大へ降ろす
  obtain ⟨M, hfin, Q_M, hQM, himg⟩ := exists_finite_point_descent E' hl R hRord
  refine ⟨M, hfin, Q_M, hQM, ?_⟩
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  intro Z hZ hZeq
  -- ☆`M → ℂ`
  set ψ : (M : Type) →+* ℂ :=
    (f : AlgebraicClosure L →+* ℂ).comp (algebraMap (M : Type) (AlgebraicClosure L)) with hψdef
  have hmapM : ∀ W : WeierstrassCurve L, (W.map (algebraMap L (M : Type))).map ψ
      = W.map (algebraMap L ℂ) := by
    intro W
    have hcomp : ψ.comp (algebraMap L (M : Type)) = algebraMap L ℂ := by
      refine RingHom.ext (fun x => ?_)
      show f ((algebraMap (M : Type) (AlgebraicClosure L)) ((algebraMap L (M : Type)) x))
        = algebraMap L ℂ x
      rw [← IsScalarTower.algebraMap_apply L (M : Type) (AlgebraicClosure L) x]
      exact f.commutes x
    show (W.map (algebraMap L (M : Type))).map ψ = W.map (algebraMap L ℂ)
    rw [← hcomp]; rfl
  haveI hEMell : (E.map (algebraMap L (M : Type))).IsElliptic := inferInstance
  haveI hE'Mell : (E'.baseChange (M : Type)).IsElliptic := by
    show (E'.map (algebraMap L (M : Type))).IsElliptic; infer_instance
  haveI hZmapell : (Z.map ψ).IsElliptic := inferInstance
  haveI hEMCell : ((E.map (algebraMap L (M : Type))).map ψ).IsElliptic := inferInstance
  -- ☆点集合が `ℂ` の側で一致する
  have hsetC : ((((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q_M))).image
        (fun q : (↥M) × (↥M) => (ψ q.1, ψ q.2)))
      = ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'C)) := by
    have h1 : ((((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q_M))).image
          (fun q : (↥M) × (↥M) => (ψ q.1, ψ q.2)))
        = (((range l).erase 0).image (fun k : ℕ => pointCoords (k • R))).image
          (fun q : AlgebraicClosure L × AlgebraicClosure L =>
            ((f : AlgebraicClosure L →+* ℂ) q.1, (f : AlgebraicClosure L →+* ℂ) q.2)) := by
      have hψapp : ∀ x : (M : Type), ψ x
          = (f : AlgebraicClosure L →+* ℂ) (algebraMap (M : Type) (AlgebraicClosure L) x) :=
        fun _ => rfl
      rw [← himg]
      simp only [Finset.image_image, Function.comp_def, hψapp]
    rw [h1, ← image_pointCoords_rhPoint_nsmul (f : AlgebraicClosure L →+* ℂ)
      (E'.baseChange (AlgebraicClosure L)) hRord, hRmap]
    refine Finset.image_congr (fun k _ => ?_)
    rw [← castPoint_nsmul, pointCoords_castPoint, ← castPoint_nsmul, pointCoords_castPoint]
  -- ☆`Z ⊗ ℂ` を `ℂ` 側の商に一致させる
  have hZC : Z.map ψ = veluQuotientFull (veluQuotientFull (E.map (algebraMap L ℂ))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • rhPoint (algebraMap L ℂ) E Q))))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'C))) := by
    have hmapM' : (E'.baseChange (M : Type)).map ψ = E'.map (algebraMap L ℂ) := hmapM E'
    rw [hZeq, veluQuotientFull_map, hsetC, hmapM', hquot]
  have hjZ := hjC (Z.map ψ) hZC
  -- ☆`ψ` は単射
  have hψinj : Function.Injective ψ := ψ.injective
  refine hψinj ?_
  have e1 : ψ Z.j = (Z.map ψ).j := (j_map ψ Z).symm
  have e2 : ψ (E.map (algebraMap L (M : Type))).j
      = ((E.map (algebraMap L (M : Type))).map ψ).j := (j_map ψ _).symm
  rw [e1, e2, hjZ]
  exact (j_congr_curve (hmapM E)).symm

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_dual_veluQuotientFull_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(双対同種の j——格子側。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_dual_veluQuotientFull_lattice.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_isogeny_periodPair(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_isogeny_periodPair") 1,
    .citation "[ABC3]" "latticeDisc_scalePair(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeDisc_scalePair") 1,
    .implicitStep
      ("★★★★★**2026-09-06（第 1442）**——`Λ ⊂ Λ′` が指数 `l` なら " ++
       "`Λ′ ⊂ (1/l)Λ` も指数 `l` であり、`(1/l)Λ` は `Λ` と相似だから `j` が等しい。" ++
       "☆これが原文が括弧で畳んだ「同種なので自動」の**双対の側**である。") 17 ]

def exists_dual_veluQuot_complex.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(双対同種の j——ℂ 上の任意の楕円曲線で。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_dual_veluQuot_j_numberField.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(双対同種の j——数体の上、有限次拡大を経て。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_dual_veluQuot_j_numberField.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_dual_veluQuot_complex(本ファイル、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_dual_veluQuot_complex") 1,
    .citation "[ABC3]" "exists_finite_point_descent(第 1207、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finite_point_descent") 1,
    .citation "[ABC3]" "torsion_card(証明済み——代数閉体では E[n] は n² 個)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.torsion_card") 1,
    .implicitStep
      ("★★★★★**2026-09-06（第 1442）**——`ℂ` で作った双対の点を " ++
       "`L̄` へ引き戻し（個数の勘定だけ）、有限次拡大へ降ろす。" ++
       "☆分体多項式も Néron モデルも要らない。") 17 ]

def exists_rhPoint_of_nsmul_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(代数閉体の間では n-捩れは rhPoint で全射。★無条件)",
    sectionId := "genell-thm-3-8" }

def latticeCurve_j_scalePair.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(相似な束の曲線は同じ j を持つ。★無条件)",
    sectionId := "genell-prop-3-4" }

def lattice_dual_step.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(双対の一段——(Λ + ℤz₀) + ℤv = (1/l)Λ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_dual_elt.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(双対の元がとれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeCurve_eq_veluQuotientFull_nsmul.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ′ を指定しての Vélu の商の等式。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
