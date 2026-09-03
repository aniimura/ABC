/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Uniformization.GroupIso

/-!
# 一様化 —— `Λ′` の基底と行列式・巡回の場合の代表系・代表系の対称性・Laurent

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Λ′` の基底と行列式 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**指数 `l` の格子 `Λ′ = Λ + ℤz₀` の基底と行列式**

`l·z₀ = a·ω₁ + b·ω₂`（`gcd(a, b, l) = 1`）のとき、`Λ′` の基底 `ω₁′, ω₂′` と
整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★構成は初等的である。`h = gcd(a,b)`、`a = h a₁`、`b = h b₁` として
`a₁p + b₁q = 1` を取り

    η₁ ≔ a₁ω₁ + b₁ω₂,   η₂ ≔ −qω₁ + pω₂

とすると `(η₁, η₂)` は `Λ` の基底（行列式 `a₁p + b₁q = 1`）で `l z₀ = h η₁`。
`gcd(h, l) = 1`（`gcd(a,b,l) = 1` から）なので `xh + yl = 1` が取れ、

    ω₁′ ≔ η₁/l,   ω₂′ ≔ η₂

とすれば `z₀ = h·ω₁′`・`ω₁′ = x·z₀ + y·η₁` となって
`Λ′ = ℤω₁′ + ℤω₂′`。行列式は

    (pl)·a₁ − (−b₁)·(ql) = l·(pa₁ + b₁q) = l

★★★★☆**これで `Lemma 3.5` の格子側——「位数 `l` の巡回部分群 ↔
指数 `l` の格子」——が完全に閉じた。**
☆`htFalt_isogeny_le_of_analytic_minimal`（第 617）が要求する
`h₁`・`h₂`・`hdet` はこの `A, B, C, D` そのものである。 -/
theorem exists_lattice_basis_of_cyclic (P : PeriodPair) (z₀ : ℂ) (l : ℕ) (hl : 0 < l)
    (a b : ℤ) (hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂)
    (hgcd : Nat.gcd (Int.gcd a b) l = 1) :
    ∃ (ω₁' ω₂' : ℂ) (A B C D : ℤ),
      P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂' ∧
      P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂' ∧
      (A * D - B * C).natAbs = l ∧
      Submodule.span ℤ ({ω₁', ω₂'} : Set ℂ)
        = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  have hlC : (l : ℂ) ≠ 0 := Nat.cast_ne_zero.2 hl.ne'
  by_cases hab : Int.gcd a b = 0
  · -- `a = b = 0` の退化した場合。`l = 1`・`z₀ = 0`。
    have ha : a = 0 := (Int.gcd_eq_zero_iff.1 hab).1
    have hb : b = 0 := (Int.gcd_eq_zero_iff.1 hab).2
    have hl1 : l = 1 := by simpa [hab] using hgcd
    have hz0 : z₀ = 0 := by
      have := hz
      rw [ha, hb] at this
      simp only [Int.cast_zero, zero_mul, add_zero] at this
      exact (mul_eq_zero.1 this).resolve_left hlC
    refine ⟨P.ω₁, P.ω₂, 1, 0, 0, 1, by push_cast; ring, by push_cast; ring, by
      simp [hl1], ?_⟩
    rw [hz0, Submodule.span_zero_singleton, sup_bot_eq]
    rfl
  · -- 本体。`h = gcd(a,b) ≠ 0`。
    set h : ℤ := (Int.gcd a b : ℤ) with hh
    have hhpos : 0 < h := by
      simpa [hh] using Nat.pos_of_ne_zero hab
    have hhne : h ≠ 0 := hhpos.ne'
    set a₁ : ℤ := a / h with ha₁
    set b₁ : ℤ := b / h with hb₁
    have hae : a = h * a₁ := by
      rw [ha₁, Int.mul_ediv_cancel' (Int.gcd_dvd_left a b)]
    have hbe : b = h * b₁ := by
      rw [hb₁, Int.mul_ediv_cancel' (Int.gcd_dvd_right a b)]
    have hg1 : Int.gcd a₁ b₁ = 1 := by
      rw [ha₁, hb₁, hh]
      exact Int.gcd_div_gcd_div_gcd (Nat.pos_of_ne_zero hab)
    -- `a₁ p + b₁ q = 1`
    set p : ℤ := Int.gcdA a₁ b₁ with hp
    set q : ℤ := Int.gcdB a₁ b₁ with hq
    have hbez1 : a₁ * p + b₁ * q = 1 := by
      have := Int.gcd_eq_gcd_ab a₁ b₁
      rw [hg1] at this
      simpa [hp, hq] using this.symm
    -- `gcd(h, l) = 1` から `x h + y l = 1`
    have hghl : Int.gcd h (l : ℤ) = 1 := by
      simpa [hh] using hgcd
    set x : ℤ := Int.gcdA h (l : ℤ) with hx
    set y : ℤ := Int.gcdB h (l : ℤ) with hy
    have hbez2 : h * x + (l : ℤ) * y = 1 := by
      have := Int.gcd_eq_gcd_ab h (l : ℤ)
      rw [hghl] at this
      simpa [hx, hy] using this.symm
    -- 新しい基底
    set η₁ : ℂ := (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ with hη₁
    set η₂ : ℂ := (-q : ℤ) * P.ω₁ + (p : ℂ) * P.ω₂ with hη₂
    set ω₁' : ℂ := η₁ / (l : ℂ) with hω₁'
    have hlω : (l : ℂ) * ω₁' = η₁ := by
      rw [hω₁']; field_simp
    clear_value ω₁'
    have hlz : (l : ℂ) * z₀ = (h : ℂ) * η₁ := by
      rw [hz, hη₁, hae, hbe]
      push_cast
      ring
    have hz0 : z₀ = (h : ℂ) * ω₁' := by
      have : (l : ℂ) * z₀ = (l : ℂ) * ((h : ℂ) * ω₁') := by
        rw [hlz, ← hlω]; ring
      exact mul_left_cancel₀ hlC this
    have hbez1C : (a₁ : ℂ) * (p : ℂ) + (b₁ : ℂ) * (q : ℂ) = 1 := by
      exact_mod_cast congrArg (fun n : ℤ => (n : ℂ)) hbez1
    have hbez2C : (h : ℂ) * (x : ℂ) + (l : ℂ) * (y : ℂ) = 1 := by
      exact_mod_cast congrArg (fun n : ℤ => (n : ℂ)) hbez2
    have hω1 : P.ω₁ = ((p * l : ℤ) : ℂ) * ω₁' + ((-b₁ : ℤ) : ℂ) * η₂ := by
      rw [hη₂]
      push_cast
      have hlω' : (l : ℂ) * ω₁' = (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ := by
        rw [hlω, hη₁]
      linear_combination (-(p : ℂ)) * hlω' - P.ω₁ * hbez1C
    have hω2 : P.ω₂ = ((q * l : ℤ) : ℂ) * ω₁' + ((a₁ : ℤ) : ℂ) * η₂ := by
      rw [hη₂]
      push_cast
      have hlω' : (l : ℂ) * ω₁' = (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ := by
        rw [hlω, hη₁]
      linear_combination (-(q : ℂ)) * hlω' - P.ω₂ * hbez1C
    refine ⟨ω₁', η₂, p * l, -b₁, q * l, a₁, hω1, hω2, ?_, ?_⟩
    · have : (p * l) * a₁ - (-b₁) * (q * l) = (l : ℤ) * (a₁ * p + b₁ * q) := by ring
      rw [this, hbez1, mul_one, Int.natAbs_natCast]
    · -- `span {ω₁′, η₂} = Λ ⊔ span {z₀}`
      have hη₁mem : η₁ ∈ P.lattice := by
        rw [PeriodPair.mem_lattice]
        exact ⟨a₁, b₁, hη₁.symm⟩
      have hη₂mem : η₂ ∈ P.lattice := by
        rw [PeriodPair.mem_lattice]
        exact ⟨-q, p, hη₂.symm⟩
      refine le_antisymm ?_ ?_
      · -- `ω₁′ = x·z₀ + y·η₁`
        have hval : ω₁' = (x : ℂ) * z₀ + (y : ℂ) * η₁ := by
          rw [hz0, ← hlω]
          linear_combination (-ω₁') * hbez2C
        have hm1 : ω₁' ∈ P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
          rw [hval]
          refine Submodule.add_mem _ ?_ ?_
          · exact Submodule.mem_sup_right (by
              rw [show (x : ℂ) * z₀ = x • z₀ by simp [zsmul_eq_mul]]
              exact Submodule.smul_mem _ _ (Submodule.mem_span_singleton_self z₀))
          · exact Submodule.mem_sup_left (by
              rw [show (y : ℂ) * η₁ = y • η₁ by simp [zsmul_eq_mul]]
              exact Submodule.smul_mem _ _ hη₁mem)
        have hm2 : η₂ ∈ P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) :=
          Submodule.mem_sup_left hη₂mem
        refine Submodule.span_le.2 ?_
        intro w hw
        simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hw
        rcases hw with h1 | h1
        · rw [SetLike.mem_coe, h1]; exact hm1
        · rw [SetLike.mem_coe, h1]; exact hm2
      · refine sup_le ?_ ?_
        · intro w hw
          obtain ⟨m, n, rfl⟩ := PeriodPair.mem_lattice.1 hw
          have h1 : P.ω₁ ∈ Submodule.span ℤ ({ω₁', η₂} : Set ℂ) := by
            rw [Submodule.mem_span_pair]
            exact ⟨p * l, -b₁, by rw [hω1]; push_cast; simp [zsmul_eq_mul]⟩
          have h2 : P.ω₂ ∈ Submodule.span ℤ ({ω₁', η₂} : Set ℂ) := by
            rw [Submodule.mem_span_pair]
            exact ⟨q * l, a₁, by rw [hω2]; push_cast; simp [zsmul_eq_mul]⟩
          refine Submodule.add_mem _ ?_ ?_
          · rw [show (m : ℂ) * P.ω₁ = m • P.ω₁ by simp [zsmul_eq_mul]]
            exact Submodule.smul_mem _ _ h1
          · rw [show (n : ℂ) * P.ω₂ = n • P.ω₂ by simp [zsmul_eq_mul]]
            exact Submodule.smul_mem _ _ h2
        · rw [Submodule.span_le, Set.singleton_subset_iff, SetLike.mem_coe, hz0]
          rw [show (h : ℂ) * ω₁' = h • ω₁' by simp [zsmul_eq_mul]]
          exact Submodule.smul_mem _ _ (Submodule.subset_span (by simp))

def exists_lattice_basis_of_cyclic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(指数 l の格子 Λ′ = Λ + ℤz₀ の基底と行列式 = l。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★**位数がちょうど `l` なら `gcd(a, b, l) = 1`**。

`l z₀ = a ω₁ + b ω₂` のとき、`g ≔ gcd(a, b, l)` とすると
`(l/g)·z₀ = (a/g)ω₁ + (b/g)ω₂ ∈ Λ`、すなわち `(l/g)·Q = 0`。
`Q` の位数が `l` なら `l ∣ l/g` なので `g = 1`。 -/
theorem gcd_eq_one_of_addOrderOf (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) {a b : ℤ}
    (hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂) :
    Nat.gcd (Int.gcd a b) l = 1 := by
  set g : ℕ := Nat.gcd (Int.gcd a b) l with hg
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_right _ hl
  have hgC : (g : ℂ) ≠ 0 := Nat.cast_ne_zero.2 hgpos.ne'
  obtain ⟨n, hn⟩ := (Nat.gcd_dvd_right (Int.gcd a b) l : g ∣ l)
  have hnpos : 0 < n := by
    rcases Nat.eq_zero_or_pos n with hc | hc
    · rw [hc, Nat.mul_zero] at hn; omega
    · exact hc
  have hga : (g : ℤ) ∣ a :=
    dvd_trans (Int.natCast_dvd_natCast.2 (Nat.gcd_dvd_left _ _)) (Int.gcd_dvd_left a b)
  have hgb : (g : ℤ) ∣ b :=
    dvd_trans (Int.natCast_dvd_natCast.2 (Nat.gcd_dvd_left _ _)) (Int.gcd_dvd_right a b)
  have hA : ((a / (g : ℤ) : ℤ) : ℂ) * (g : ℂ) = (a : ℂ) := by
    exact_mod_cast congrArg (fun t : ℤ => (t : ℂ)) (Int.ediv_mul_cancel hga)
  have hB : ((b / (g : ℤ) : ℤ) : ℂ) * (g : ℂ) = (b : ℂ) := by
    exact_mod_cast congrArg (fun t : ℤ => (t : ℂ)) (Int.ediv_mul_cancel hgb)
  have hnz : (n : ℂ) * z₀
      = ((a / (g : ℤ) : ℤ) : ℂ) * P.ω₁ + ((b / (g : ℤ) : ℤ) : ℂ) * P.ω₂ := by
    refine mul_left_cancel₀ hgC ?_
    calc (g : ℂ) * ((n : ℂ) * z₀) = (l : ℂ) * z₀ := by
          rw [hn]; push_cast; ring
      _ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂ := hz
      _ = (g : ℂ) * (((a / (g : ℤ) : ℤ) : ℂ) * P.ω₁
            + ((b / (g : ℤ) : ℤ) : ℂ) * P.ω₂) := by
          rw [← hA, ← hB]; ring
  have hnzmem : (n : ℂ) * z₀ ∈ P.lattice :=
    PeriodPair.mem_lattice.2 ⟨a / (g : ℤ), b / (g : ℤ), hnz.symm⟩
  have hnQ : n • Q = 0 := by
    have hzero : uniformMap P hΔ ((n : ℂ) * z₀) = 0 := uniformMap_of_mem P hΔ hnzmem
    rw [show ((n : ℂ) * z₀) = n • z₀ by simp [nsmul_eq_mul], ← uniformHom_apply,
      map_nsmul, uniformHom_apply, hz₀] at hzero
    exact hzero
  have hdvd : l ∣ n := hQ ▸ addOrderOf_dvd_of_nsmul_eq_zero hnQ
  have hle : l ≤ n := Nat.le_of_dvd hnpos hdvd
  have hne : n ≤ l := by
    rw [hn]
    exact Nat.le_mul_of_pos_left n hgpos
  have hnl : n = l := le_antisymm hne hle
  have hn' : l = g * l := by
    conv_lhs => rw [hn]
    rw [hnl]
  refine Nat.eq_of_mul_eq_mul_right hl ?_
  rw [← hn', Nat.one_mul]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（格子側）——位数 `l` の点から指数 `l` の格子を作る**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`（`Φ(z₀) = Q`）と
`Λ′ = Λ + ℤz₀` の基底 `ω₁′, ω₂′` と整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

★★★★★☆**これが `htFalt_isogeny_le_of_analytic_minimal`（第 617）の
`h₁`・`h₂`・`hdet` そのものである。**

★部品:

* 全射（第 663）で `z₀` を取る
* `l·Q = 0` と核 `= Λ`（第 663）で `l z₀ ∈ Λ`、すなわち `l z₀ = aω₁ + bω₂`
* 位数がちょうど `l` だから `gcd(a, b, l) = 1`（本ブロック）
* Hermite 標準形（第 666）で基底と行列式 -/
theorem exists_isogeny_lattice_basis (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ ω₁' ω₂' : ℂ) (A B C D : ℤ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂' ∧
      P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂' ∧
      (A * D - B * C).natAbs = l ∧
      Submodule.span ℤ ({ω₁', ω₂'} : Set ℂ)
        = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  obtain ⟨z₀, hz₀⟩ := uniformMap_surjective P hΔ Q
  have hlQ : l • Q = 0 := hQ ▸ addOrderOf_nsmul_eq_zero Q
  have hlz : ((l : ℂ) * z₀) ∈ P.lattice := by
    refine (uniformMap_eq_zero_iff P hΔ _).1 ?_
    rw [show ((l : ℂ) * z₀) = l • z₀ by simp [nsmul_eq_mul], ← uniformHom_apply,
      map_nsmul, uniformHom_apply, hz₀, hlQ]
  obtain ⟨a, b, hab⟩ := PeriodPair.mem_lattice.1 hlz
  have hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂ := hab.symm
  obtain ⟨ω₁', ω₂', A, B, C, D, h1, h2, hdet, hspan⟩ :=
    exists_lattice_basis_of_cyclic P z₀ l hl a b hz
      (gcd_eq_one_of_addOrderOf P hΔ hl hQ hz₀ hz)
  exact ⟨z₀, ω₁', ω₂', A, B, C, D, hz₀, h1, h2, hdet, hspan⟩

def gcd_eq_one_of_addOrderOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数がちょうど l なら gcd(a, b, l) = 1。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isogeny_lattice_basis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子側——位数 l の点から指数 l の格子を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★**基底変換は ℝ-線型独立性を保つ**。

`ω₁ = Aω₁′ + Bω₂′`・`ω₂ = Cω₁′ + Dω₂′`・`AD − BC ≠ 0` なら
`ω₁′, ω₂′` も ℝ 上独立。

★逆行列は `Δω₁′ = Dω₁ − Bω₂`・`Δω₂′ = −Cω₁ + Aω₂`（`Δ = AD − BC`）。
`rω₁′ + tω₂′ = 0` に `Δ` を掛けて `ω₁, ω₂` の独立性に帰着させる。 -/
theorem linearIndependent_of_basis_change (P : PeriodPair) {ω₁' ω₂' : ℂ} {A B C D : ℤ}
    (h1 : P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂')
    (h2 : P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂')
    (hdet : A * D - B * C ≠ 0) :
    LinearIndependent ℝ ![ω₁', ω₂'] := by
  have hP := LinearIndependent.pair_iff.1 P.indep
  have hΔR : ((A : ℝ) * D - B * C) ≠ 0 := by
    have : ((A * D - B * C : ℤ) : ℝ) ≠ 0 := Int.cast_ne_zero.2 hdet
    push_cast at this
    exact this
  refine LinearIndependent.pair_iff.2 ?_
  intro r t hrt
  have hrt' : (r : ℂ) * ω₁' + (t : ℂ) * ω₂' = 0 := by
    simpa [Complex.real_smul] using hrt
  have hkey : (r * D - t * C : ℝ) • P.ω₁ + (-(r * B) + t * A : ℝ) • P.ω₂ = 0 := by
    simp only [Complex.real_smul]
    push_cast
    rw [h1, h2]
    linear_combination ((A : ℂ) * (D : ℂ) - (B : ℂ) * (C : ℂ)) * hrt'
  obtain ⟨e1, e2⟩ := hP _ _ hkey
  constructor
  · have hr : r * ((A : ℝ) * D - B * C) = 0 := by
      linear_combination (A : ℝ) * e1 + (C : ℝ) * e2
    exact (mul_eq_zero.1 hr).resolve_right hΔR
  · have ht : t * ((A : ℝ) * D - B * C) = 0 := by
      linear_combination (B : ℝ) * e1 + (D : ℝ) * e2
    exact (mul_eq_zero.1 ht).resolve_right hΔR

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（格子側・完成形）——位数 `l` の点から `PeriodPair` `P′` を作る**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`（`Φ(z₀) = Q`）と
**周期対** `P′`（格子は `Λ′ = Λ + ℤz₀`）と整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

★★★★★☆**これで `htFalt_isogeny_le_of_analytic_minimal`（第 617）が要求する
`P′`・`h₁`・`h₂`・`hdet` がすべて揃った。**
☆残るのは `α`（`u′ = α·u`）——代数的な同種写像 `E → E/H` と
この解析的な `Λ ⊆ Λ′` を突き合わせることである。 -/
theorem exists_isogeny_periodPair (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  obtain ⟨z₀, ω₁', ω₂', A, B, C, D, hz₀, h1, h2, hdet, hspan⟩ :=
    exists_isogeny_lattice_basis P hΔ hl hQ
  have hdet0 : A * D - B * C ≠ 0 := by
    intro hc
    rw [hc] at hdet
    simp at hdet
    omega
  exact ⟨z₀, ⟨ω₁', ω₂', linearIndependent_of_basis_change P h1 h2 hdet0⟩,
    A, B, C, D, hz₀, h1, h2, hdet, hspan⟩

def linearIndependent_of_basis_change.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(基底変換は ℝ-線型独立性を保つ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isogeny_periodPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子側・完成形——位数 l の点から PeriodPair P′ を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★巡回の場合の代表系 -/

/-- ★`|k| < l` で `l ∣ k` なら `k = 0`。 -/
theorem eq_zero_of_dvd_of_abs_lt {l k : ℤ} (hl : 0 < l) (h : l ∣ k)
    (h1 : -l < k) (h2 : k < l) : k = 0 := by
  obtain ⟨c, rfl⟩ := h
  rcases lt_trichotomy c 0 with hc | hc | hc
  · nlinarith
  · simp [hc]
  · nlinarith

open scoped Classical in
/-- ★★★★★★★★★★★★★★`k·z₀ ∈ Λ ⟺ l ∣ k`——`Q` の位数がちょうど `l` のとき。 -/
theorem intCast_mul_mem_lattice_iff (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) (k : ℤ) :
    (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k := by
  rw [← uniformMap_eq_zero_iff P hΔ,
    show ((k : ℂ) * z₀) = k • z₀ by simp [zsmul_eq_mul],
    ← uniformHom_apply, map_zsmul, uniformHom_apply, hz₀, ← hQ,
    addOrderOf_dvd_iff_zsmul_eq_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**巡回の場合の代表系 `T = {0, z₀, 2z₀, …, (l−1)z₀}`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`Λ′ = Λ + ℤz₀` で `z₀` の `Λ` を法とした位数がちょうど `l` のとき、
`T` は `Λ′/Λ` の代表系である（第 601 `weierstrassP_eq_velu_of_rep` の仮説）。

★★★`p = y + n z₀ ∈ Λ′` に対し `w₀ ≔ ((−n) mod l)·z₀` を取れば
`p + w₀ ∈ Λ`。一意性は `|k − m| < l` かつ `l ∣ k − m` から。 -/
theorem exists_velu_rep (P P' : PeriodPair) (z₀ : ℂ) (l : ℕ) (hl : 0 < l)
    (hP' : P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ))
    (hord : ∀ k : ℤ, (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k) :
    ∃ T : Finset ℂ, (0 : ℂ) ∈ T ∧ T.card = l ∧ (∀ w ∈ T, w ∈ P'.lattice) ∧
      (∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
        ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
      ∧ T = (Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀) := by
  have hlZ : (0 : ℤ) < (l : ℤ) := by exact_mod_cast hl
  have hinj : ∀ k ∈ Finset.range l, ∀ k' ∈ Finset.range l,
      (k : ℂ) * z₀ = (k' : ℂ) * z₀ → k = k' := by
    intro k hk k' hk' he
    simp only [Finset.mem_range] at hk hk'
    have h0 : (((k : ℤ) - (k' : ℤ) : ℤ) : ℂ) * z₀ ∈ P.lattice := by
      have hz : (((k : ℤ) - (k' : ℤ) : ℤ) : ℂ) * z₀ = 0 := by
        push_cast
        linear_combination he
      rw [hz]
      exact P.lattice.zero_mem
    have hd := (hord _).1 h0
    have hzero := eq_zero_of_dvd_of_abs_lt hlZ hd (by omega) (by omega)
    omega
  refine ⟨(Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀), ?_, ?_, ?_, ?_, rfl⟩
  · exact Finset.mem_image.2 ⟨0, Finset.mem_range.2 hl, by simp⟩
  · have hio : Set.InjOn (fun k : ℕ => (k : ℂ) * z₀) ↑(Finset.range l) := by
      intro k hk k' hk' he
      exact hinj k (by simpa using hk) k' (by simpa using hk') he
    rw [Finset.card_image_of_injOn hio, Finset.card_range]
  · intro w hw
    obtain ⟨k, -, rfl⟩ := Finset.mem_image.1 hw
    rw [hP']
    refine Submodule.mem_sup_right ?_
    rw [show ((k : ℂ) * z₀) = (k : ℤ) • z₀ by simp [zsmul_eq_mul]]
    exact Submodule.smul_mem _ _ (Submodule.mem_span_singleton_self z₀)
  · intro p hp
    rw [hP', Submodule.mem_sup] at hp
    obtain ⟨y, hy, w, hw, rfl⟩ := hp
    obtain ⟨n, rfl⟩ := Submodule.mem_span_singleton.1 hw
    set r : ℤ := (-n) % (l : ℤ) with hrdef
    have hr0 : 0 ≤ r := Int.emod_nonneg _ hlZ.ne'
    have hrl : r < (l : ℤ) := Int.emod_lt_of_pos _ hlZ
    set m : ℕ := r.toNat with hmdef
    have hmr : (m : ℤ) = r := Int.toNat_of_nonneg hr0
    have hml : m < l := by omega
    have hdvd : (l : ℤ) ∣ (n + (m : ℤ)) := by
      rw [hmr]
      refine ⟨-((-n) / (l : ℤ)), ?_⟩
      have hq := Int.emod_add_ediv (-n) (l : ℤ)
      linear_combination hq
    refine ⟨(m : ℂ) * z₀, Finset.mem_image.2 ⟨m, Finset.mem_range.2 hml, rfl⟩, ?_, ?_⟩
    · have he : y + (n • z₀) + (m : ℂ) * z₀ = y + ((n + (m : ℤ) : ℤ) : ℂ) * z₀ := by
        push_cast [zsmul_eq_mul]
        ring
      rw [he]
      exact P.lattice.add_mem hy ((hord _).2 hdvd)
    · intro v hv hvne hmem
      obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hv
      rw [Finset.mem_range] at hk
      have he : y + (n • z₀) + (k : ℂ) * z₀ = y + ((n + (k : ℤ) : ℤ) : ℂ) * z₀ := by
        push_cast [zsmul_eq_mul]
        ring
      rw [he] at hmem
      have hk2 : ((n + (k : ℤ) : ℤ) : ℂ) * z₀ ∈ P.lattice := by
        have hd2 := P.lattice.sub_mem hmem hy
        simpa using hd2
      have hsub : (l : ℤ) ∣ ((k : ℤ) - (m : ℤ)) := by
        have h2 := (hord (n + (k : ℤ))).1 hk2
        have h3 : (l : ℤ) ∣ ((n + (k : ℤ)) - (n + (m : ℤ))) := dvd_sub h2 hdvd
        have h4 : (n + (k : ℤ)) - (n + (m : ℤ)) = (k : ℤ) - (m : ℤ) := by ring
        rwa [h4] at h3
      have hzero : ((k : ℤ) - (m : ℤ)) = 0 :=
        eq_zero_of_dvd_of_abs_lt hlZ hsub (by omega) (by omega)
      exact hvne (by rw [show k = m by omega])

def exists_velu_rep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(巡回の場合の代表系 T = {0, z₀, …, (l−1)z₀}。★無条件)",
    sectionId := "genell-lemma-3-5" }

def intCast_mul_mem_lattice_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(k·z₀ ∈ Λ ⟺ l ∣ k——位数がちょうど l のとき。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（解析側・完成形）——位数 `l` の点から Vélu の公式まで**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`・周期対 `P′`・整数 `A, B, C, D`・
代表系 `T`（`|T| = l`）が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

    ℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)

★★★★★★☆**これで「位数 `l` の点 → 指数 `l` の格子 → その `℘` は
`Λ` の `℘` の Vélu 和で書ける」が一本につながった。**

☆残るのは、この解析的な等式を代数的な Vélu の商
（`Found/GenEll/Velu.lean` の `veluQuotient`）と突き合わせて
`α`（`u′ = α·u`）を作ることである。 -/
theorem exists_velu_formula_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ) (T : Finset ℂ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (0 : ℂ) ∈ T ∧ T.card = l ∧
      (∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  exact ⟨z₀, P', A, B, C, D, T, hz₀, h1, h2, hdet, hP', h0T, hcard,
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z⟩

def exists_velu_formula_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側・完成形——位数 l の点から Vélu の公式まで。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★代表系の対称性 -/

/-- ★★★★★★★★★★★★**代表系の元は `Λ` を法として相異なる**。 -/
theorem rep_sub_mem_lattice_imp_eq (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {w v : ℂ} (hw : w ∈ T) (hv : v ∈ T) (hd : w - v ∈ P.lattice) : w = v := by
  obtain ⟨w₀, hw₀T, hw₀Λ, hw₀u⟩ := hrep (-v) (neg_mem (hT v hv))
  have hv0 : v = w₀ := by
    by_contra hc
    exact hw₀u v hv hc (by simpa using P.lattice.zero_mem)
  have hw0 : w = w₀ := by
    by_contra hc
    exact hw₀u w hw hc (by rw [show -v + w = w - v by ring]; exact hd)
  rw [hw0, hv0]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系の上での `℘′` の和は消える**

    `Σ_{w ∈ T∖{0}} ℘′_Λ(w) = 0`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`w ↦ ν(w)`（`w + ν w ∈ Λ` を満たす唯一の代表元、つまり `ν w ≡ −w`）は
`T∖{0}` の対合であり、`℘′(ν w) = ℘′(−w) = −℘′(w)`。
したがって `S = Σ ℘′(ν w) = −S`、すなわち `S = 0`。

★★☆これが Vélu の `ω`-正規化（第 586-593 の代数側 `velu_omega_gen`）の解析版である。 -/
theorem sum_derivWeierstrassP_rep_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    ∑ w ∈ T.erase 0, P.derivWeierstrassP w = 0 := by
  classical
  have huniq : ∀ {w v : ℂ}, w ∈ T → v ∈ T → w - v ∈ P.lattice → w = v :=
    fun hw hv hd => rep_sub_mem_lattice_imp_eq P P' T hT hrep hw hv hd
  have hex : ∀ w ∈ T, ∃ v, v ∈ T ∧ w + v ∈ P.lattice := by
    intro w hw
    obtain ⟨v, hv, hv2, -⟩ := hrep w (hT w hw)
    exact ⟨v, hv, hv2⟩
  choose! ν hνT hνΛ using hex
  have hνe : ∀ w ∈ T.erase 0, ν w ∈ T.erase 0 := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have hw0 : w ≠ 0 := Finset.ne_of_mem_erase hw
    refine Finset.mem_erase.2 ⟨?_, hνT w hw'⟩
    intro hc
    refine hw0 (huniq hw' h0T ?_)
    have hz := hνΛ w hw'
    rw [hc, add_zero] at hz
    simpa using hz
  have hinvol : ∀ w ∈ T.erase 0, ν (ν w) = w := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have h1 := hνΛ w hw'
    have h2 := hνΛ (ν w) (hνT w hw')
    refine huniq (hνT (ν w) (hνT w hw')) hw' ?_
    have hd := P.lattice.sub_mem h2 h1
    rw [show ν w + ν (ν w) - (w + ν w) = ν (ν w) - w by ring] at hd
    exact hd
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0, ν w = ν v → w = v := by
    intro w hw v hv he
    rw [← hinvol w hw, ← hinvol v hv, he]
  have hodd : ∀ w ∈ T, P.derivWeierstrassP (ν w) = -P.derivWeierstrassP w := by
    intro w hw
    have hl : w + ν w ∈ P.lattice := hνΛ w hw
    have he : ν w = -w + (w + ν w) := by ring
    rw [he, P.derivWeierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.derivWeierstrassP_neg]
  have hinjOn : Set.InjOn ν ↑(T.erase 0) := fun w hw v hv he =>
    hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he
  have himg : (T.erase 0).image ν = T.erase 0 :=
    Finset.eq_of_subset_of_card_le
      (fun v hv => by
        obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hv
        exact hνe w hw)
      (le_of_eq (Finset.card_image_of_injOn hinjOn).symm)
  have h1 : ∑ v ∈ T.erase 0, P.derivWeierstrassP v
      = ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w) := by
    conv_lhs => rw [← himg]
    exact Finset.sum_image (fun w hw v hv he => hinj w hw v hv he)
  have h2 : ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w)
      = -∑ w ∈ T.erase 0, P.derivWeierstrassP w := by
    have hc : ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w)
        = ∑ w ∈ T.erase 0, (-P.derivWeierstrassP w) :=
      Finset.sum_congr rfl (fun w hw => hodd w (Finset.mem_of_mem_erase hw))
    rw [hc, Finset.sum_neg_distrib]
  have hSS : ∑ w ∈ T.erase 0, P.derivWeierstrassP w
      = -∑ w ∈ T.erase 0, P.derivWeierstrassP w := h1.trans h2
  have h3 : (2 : ℂ) * ∑ w ∈ T.erase 0, P.derivWeierstrassP w = 0 := by
    linear_combination hSS
  simpa using h3

def rep_sub_mem_lattice_imp_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の元は Λ を法として相異なる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_derivWeierstrassP_rep_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の上での ℘′ の和は消える——Vélu の ω-正規化の解析版。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★Laurent 比較の入口 -/

/-- ★★★★★★★★★★**代表系の `0` 以外の元は格子の外**——`w ≡ 0` なら `w = 0` だから。 -/
theorem rep_notMem_lattice (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {w : ℂ} (hw : w ∈ T.erase 0) : w ∉ P.lattice := by
  intro hc
  refine Finset.ne_of_mem_erase hw ?_
  refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw) h0T ?_
  simpa using hc

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`℘_{Λ′} − ℘_Λ` は `Λ` の `℘` の平行移動の和**

    `℘_{Λ′}(z) − ℘_Λ(z) = Σ_{w ∈ T∖{0}} (℘_Λ(z + w) − ℘_Λ(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 669 の Vélu の公式から `w = 0` の項を分離しただけ。
☆左辺は原点の `z⁻²` の極が打ち消し合っており、右辺は原点で解析的
（`T∖{0}` の元は格子の外）。**これが Laurent 係数の比較の入口である**——
`z²` の係数から `g₂′ − g₂`、`z⁴` の係数から `g₃′ − g₃` が出る。 -/
theorem weierstrassP_sub_eq_sum (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) (z : ℂ) :
    P'.weierstrassP z - P.weierstrassP z
      = ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) := by
  rw [hvelu z]
  simp only [veluAnalyticX, veluAnalyticC]
  rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (z + w)) h0T]
  simp only [add_zero]
  rw [Finset.sum_sub_distrib]
  ring

/-- ★★★★★★★★★★★★**平行移動の和は原点で解析的**。 -/
theorem analyticAt_veluShiftSum (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => ∑ w ∈ T.erase 0,
      (P.weierstrassP (z + w) - P.weierstrassP w)) 0 := by
  have h : AnalyticAt ℂ (∑ w ∈ T.erase 0,
      fun z : ℂ => P.weierstrassP (z + w) - P.weierstrassP w) 0 := by
    refine Finset.analyticAt_sum _ ?_
    intro w hw
    refine AnalyticAt.sub ?_ analyticAt_const
    refine shifted_analyticAt P 0 w ?_
    rw [zero_add]
    exact rep_notMem_lattice P P' T h0T hT hrep hw
  refine h.congr ?_
  filter_upwards with z
  simp [Finset.sum_apply]

/-- ★★★★★★★★★★★★★★★★
**`℘[Λ′ − 0] − ℘[Λ − 0]` は平行移動の和に等しい**——`z⁻²` が消える形。

☆左辺は mathlib の `iteratedDeriv_weierstrassPExcept_self` で Taylor 係数が
`sumInvPow`（＝`G n`、したがって `g₂`・`g₃`）で書けている。 -/
theorem weierstrassPExcept_sub_eq_sum (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) (z : ℂ) :
    P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z
      = ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) := by
  rw [← weierstrassP_sub_invSq P' z, ← weierstrassP_sub_invSq P z,
    ← weierstrassP_sub_eq_sum P P' T h0T hvelu z]
  ring

def weierstrassP_sub_eq_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘_{Λ′} − ℘_Λ は ℘_Λ の平行移動の和——Laurent 比較の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassPExcept_sub_eq_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘[Λ′−0] − ℘[Λ−0] は平行移動の和——z⁻² が消える形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def rep_notMem_lattice.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の 0 以外の元は格子の外。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
