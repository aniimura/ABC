import ABC3.Found.GenEll.OrdvIntegral
import ABC3.Found.GenEll.HeightADiv
import ABC3.Found.GenEll.NorthcottNF
import Mathlib.RingTheory.Ideal.Norm.AbsNorm

/-!
# [GenEll] Proposition 1.4, (iv) の基底 —— 分母の評価(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★何を取るか

`Found/GenEll/NorthcottNF.lean` は**整元**についての Northcott 性を取った。
残っていたのは**分母の評価**である:

> `x ∈ K` の高さが `B` 以下なら、**`B` 以下の自然数 `d`** があって `d·x` は整元。

★これで `{x : K | H(x) ≤ B}` の有限性が出る。

## ★機構 —— 分母を潰すイデアル

`𝔡_x ≝ ∏_v v^{max(−ord_v(x), 0)}` と置く(積は `ord_v(x) ≠ 0` の有限個の `v` 上)。

- `a ∈ 𝔡_x` なら `ord_v(a) ≥ max(−ord_v(x),0) ≥ −ord_v(x)` なので `a·x` は整元
- `Ideal.absNorm I ∈ I`(mathlib)より `d ≝ N(𝔡_x)` が使える
- `N(𝔡_x) = ∏_v q_v^{max(−ord_v(x),0)} = ∏ᶠ_v max(|x|_v, 1) ≤ H(x)`

★★**最後の等式が `HeightADiv.lean` の内容そのものである**
——`max(|x|_v,1) = q_v^{max(−ord_v(x),0)}`。
高さの有限素点側が、まさに分母イデアルのノルムである。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★分母を潰すイデアル -/

/-- ★`v` での分母の指数 `max(−ord_v(x), 0)`。 -/
noncomputable def denExp (x : Kˣ) (v : FinitePlace K) : ℕ := (max (-(ordv v x)) 0).toNat

theorem denExp_cast (x : Kˣ) (v : FinitePlace K) :
    ((denExp x v : ℕ) : ℤ) = max (-(ordv v x)) 0 :=
  Int.toNat_of_nonneg (le_max_right _ _)

theorem neg_ordv_le_denExp (x : Kˣ) (v : FinitePlace K) :
    -(ordv v x) ≤ ((denExp x v : ℕ) : ℤ) := by
  rw [denExp_cast]
  exact le_max_left _ _

theorem denExp_eq_zero_of_ordv_eq_zero {x : Kˣ} {v : FinitePlace K} (h : ordv v x = 0) :
    denExp x v = 0 := by
  simp [denExp, h]

open scoped Classical in
/-- ★★**分母を潰すイデアル** `𝔡_x ≝ ∏_v v^{max(−ord_v(x),0)}`。 -/
noncomputable def denIdeal (x : Kˣ) : Ideal (𝓞 K) :=
  ∏ v ∈ (ordv_finite_support x).toFinset, v.asIdeal ^ denExp x v

open scoped Classical in
theorem denIdeal_ne_bot (x : Kˣ) : denIdeal x ≠ ⊥ := by
  rw [denIdeal]
  refine Finset.prod_ne_zero_iff.2 fun v _ => ?_
  exact pow_ne_zero _ v.ne_bot

open scoped Classical in
/-- ★`𝔡_x` の元は、各素点で指数以上の `ord` を持つ。 -/
theorem denExp_le_ordv_of_mem {x : Kˣ} {a : 𝓞 K} (ha : a ≠ 0)
    (hmem : a ∈ denIdeal x) (u : Kˣ) (hu : (u : K) = algebraMap (𝓞 K) K a)
    (v : FinitePlace K) :
    ((denExp x v : ℕ) : ℤ) ≤ ordv v u := by
  classical
  by_cases hv : v ∈ (ordv_finite_support x).toFinset
  · -- `v^{e_v}` は積の因子なので `𝔡_x ≤ v^{e_v}`
    have hdvd : v.asIdeal ^ denExp x v ∣ denIdeal x := by
      rw [denIdeal]
      exact Finset.dvd_prod_of_mem _ hv
    have hle : denIdeal x ≤ v.asIdeal ^ denExp x v := Ideal.le_of_dvd hdvd
    have hmem' : a ∈ v.asIdeal ^ denExp x v := hle hmem
    -- `v^{e_v} ∣ (a)`
    have hspan : (Ideal.span {a} : Ideal (𝓞 K)) ≤ v.asIdeal ^ denExp x v := by
      rw [Ideal.span_le, Set.singleton_subset_iff]
      exact hmem'
    have hdvd2 : v.asIdeal ^ denExp x v ∣ (Ideal.span {a} : Ideal (𝓞 K)) :=
      Ideal.dvd_iff_le.2 hspan
    -- `count` へ
    have hspan_ne : (Ideal.span {a} : Ideal (𝓞 K)) ≠ 0 := by
      simpa [Ideal.zero_eq_bot, Ideal.span_singleton_eq_bot] using ha
    have hirr : Irreducible (Associates.mk v.asIdeal) := by
      rw [Associates.irreducible_mk]
      exact (Ideal.prime_of_isPrime v.ne_bot v.isPrime).irreducible
    have hmk_ne : Associates.mk (Ideal.span {a} : Ideal (𝓞 K)) ≠ 0 := by
      simpa [Associates.mk_eq_zero] using hspan_ne
    have hcount : denExp x v
        ≤ (Associates.mk v.asIdeal).count (Associates.mk (Ideal.span {a})).factors := by
      rw [← Associates.prime_pow_dvd_iff_le hmk_ne hirr, ← Associates.mk_pow]
      exact Associates.mk_dvd_mk.2 hdvd2
    rw [ordv_algebraMap_eq_count v a ha u hu]
    exact_mod_cast hcount
  · -- 台の外では指数 0、`a` は整元なので `ord ≥ 0`
    have h0 : ordv v x = 0 := by
      simpa [Set.Finite.mem_toFinset] using hv
    rw [denExp_eq_zero_of_ordv_eq_zero h0]
    simpa using ordv_algebraMap_nonneg v a ha u hu

open scoped Classical in
/-- ★★**`𝔡_x` の元を掛ければ整元になる**。 -/
theorem mul_isIntegral_of_mem_denIdeal {x : Kˣ} {a : 𝓞 K} (ha : a ≠ 0)
    (hmem : a ∈ denIdeal x) :
    (algebraMap (𝓞 K) K a * (x : K)) ∈ (algebraMap (𝓞 K) K).range := by
  classical
  have hane : algebraMap (𝓞 K) K a ≠ 0 := by
    simpa using (FaithfulSMul.algebraMap_injective (𝓞 K) K).ne ha
  set u : Kˣ := Units.mk0 (algebraMap (𝓞 K) K a) hane with hudef
  have hu : (u : K) = algebraMap (𝓞 K) K a := rfl
  have key : ((u * x : Kˣ) : K) = algebraMap (𝓞 K) K a * (x : K) := rfl
  rw [← key]
  refine mem_range_algebraMap_of_ordv_nonneg (u * x) (fun v => ?_)
  rw [ordv_mul]
  have h1 : ((denExp x v : ℕ) : ℤ) ≤ ordv v u := denExp_le_ordv_of_mem ha hmem u hu v
  have h2 : -(ordv v x) ≤ ((denExp x v : ℕ) : ℤ) := neg_ordv_le_denExp x v
  linarith

/-! ## ★★分母イデアルのノルムは高さの有限素点側である -/

/-- ★**`max(|x|_w, 1) = q_v^{max(−ord_v(x),0)}`**。

★`q_v > 1` なので、指数の `max` が値の `max` に翻訳される。 -/
theorem max_finitePlace_eq_pow (w : NumberField.FinitePlace K) (x : Kˣ) :
    max (w ((x : K))) 1
      = (residueCard (NumberField.FinitePlace.maximalIdeal w) : ℝ)
          ^ denExp x (NumberField.FinitePlace.maximalIdeal w) := by
  set v := NumberField.FinitePlace.maximalIdeal w with hv
  have hq1 : (1 : ℝ) < (residueCard v : ℝ) := by
    have : 1 < residueCard v := by
      simpa [residueCard] using NumberField.HeightOneSpectrum.one_lt_absNorm v
    exact_mod_cast this
  have hq0 : (0 : ℝ) < (residueCard v : ℝ) := lt_trans zero_lt_one hq1
  rw [finitePlace_apply_units w x, denExp, ← zpow_natCast,
    Int.toNat_of_nonneg (le_max_right (-(ordv v x)) 0)]
  rcases le_or_gt 0 (-(ordv v x)) with h | h
  · rw [max_eq_left h, max_eq_left]
    exact one_le_zpow₀ hq1.le h
  · rw [max_eq_right h.le, max_eq_right, zpow_zero]
    rw [← zpow_zero ((residueCard v : ℝ))]
    exact (zpow_le_zpow_iff_right₀ hq1).2 h.le

open scoped Classical in
/-- ★★**`N(𝔡_x) = ∏ᶠ_v max(|x|_v, 1)`** —— 分母イデアルのノルムは高さの有限素点側。 -/
theorem absNorm_denIdeal (x : Kˣ) :
    ((Ideal.absNorm (denIdeal x) : ℕ) : ℝ)
      = ∏ᶠ w : NumberField.FinitePlace K, max (w ((x : K))) 1 := by
  classical
  -- 右辺の添字を `HeightOneSpectrum` に付け替える
  have hR : ∏ᶠ w : NumberField.FinitePlace K, max (w ((x : K))) 1
      = ∏ᶠ v : FinitePlace K, (residueCard v : ℝ) ^ denExp x v := by
    rw [finprod_congr (fun w => max_finitePlace_eq_pow w x)]
    exact finprod_comp_equiv (M := ℝ)
      (e := NumberField.FinitePlace.equivHeightOneSpectrum (K := K))
      (f := fun v : FinitePlace K => (residueCard v : ℝ) ^ denExp x v)
  rw [hR]
  -- `finprod` を有限積に落とす
  have hsupp : (Function.mulSupport
      (fun v : FinitePlace K => (residueCard v : ℝ) ^ denExp x v))
      ⊆ (ordv_finite_support x).toFinset := by
    intro v hv
    rw [Function.mem_mulSupport] at hv
    have he : denExp x v ≠ 0 := by
      intro hcon
      exact hv (by rw [hcon, pow_zero])
    have h0 : ordv v x ≠ 0 := by
      intro hcon
      exact he (denExp_eq_zero_of_ordv_eq_zero hcon)
    simpa [Set.Finite.mem_toFinset] using h0
  rw [finprod_eq_prod_of_mulSupport_subset _ hsupp]
  -- 左辺を計算する
  rw [denIdeal, map_prod]
  push_cast
  refine Finset.prod_congr rfl fun v _ => ?_
  rw [map_pow]
  norm_cast

open scoped Classical in
/-- ★★**`N(𝔡_x) ≤ H(x)`** —— 分母は高さで抑えられる。

★高さの無限素点側の因子はすべて 1 以上なので、有限素点側だけ残せる。 -/
theorem absNorm_denIdeal_le_mulHeight₁ (x : Kˣ) :
    ((Ideal.absNorm (denIdeal x) : ℕ) : ℝ) ≤ Height.mulHeight₁ ((x : K)) := by
  classical
  rw [absNorm_denIdeal x, NumberField.mulHeight₁_eq]
  have hP : (0 : ℝ) ≤ ∏ᶠ w : NumberField.FinitePlace K, max (w ((x : K))) 1 :=
    le_trans zero_le_one (one_le_finprod (fun w => le_max_right _ _))
  have hA : (1 : ℝ) ≤ ∏ v : InfinitePlace K, max (v ((x : K))) 1 ^ v.mult :=
    Finset.one_le_prod (fun v _ => one_le_pow₀ (le_max_right _ _))
  nth_rewrite 1 [← one_mul (∏ᶠ w : NumberField.FinitePlace K, max (w ((x : K))) 1)]
  exact mul_le_mul_of_nonneg_right hA hP

/-! ## ★★★数体上の Northcott 性(全体) -/

open scoped Classical in
/-- ★★★**`{y : K | H(y) ≤ B}` は有限である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★mathlib は `Northcott (mulHeight₁ (K := K))` の**基底 instance を持たない**
(2026-08-17 実測)。本定理がそれである。

★機構: `d ≝ N(𝔡_y) ≤ H(y) ≤ B` を取ると `d·y` は**整元**で、
その共役はすべて `B²` 以下。整元の側は
`NumberField.Embeddings.finite_of_norm_le` で有限、`d` は `⌊B⌋₊` 以下で有限。
`y = (d·y)/d` なので有限集合の像に入る。 -/
theorem finite_mulHeight₁_le (B : ℝ) :
    {y : K | Height.mulHeight₁ y ≤ B}.Finite := by
  classical
  have hfin1 : {z : K | IsIntegral ℤ z ∧ ∀ φ : K →+* ℂ, ‖φ z‖ ≤ B * B}.Finite :=
    NumberField.Embeddings.finite_of_norm_le K ℂ (B * B)
  have hfin2 : (Set.Iic (⌊B⌋₊)).Finite := Set.finite_Iic _
  refine Set.Finite.subset
    (Set.Finite.insert 0 (Set.Finite.image (fun p : K × ℕ => p.1 / (p.2 : K))
      (Set.Finite.prod hfin1 hfin2))) ?_
  intro y hy
  simp only [Set.mem_setOf_eq] at hy
  rcases eq_or_ne y 0 with rfl | hy0
  · exact Set.mem_insert _ _
  refine Set.mem_insert_of_mem _ ?_
  set u : Kˣ := Units.mk0 y hy0 with hu
  set d : ℕ := Ideal.absNorm (denIdeal u) with hd
  have hdle : (d : ℝ) ≤ B := le_trans (absNorm_denIdeal_le_mulHeight₁ u) hy
  have hdne : d ≠ 0 := by
    rw [hd]
    simpa [Ideal.absNorm_eq_zero_iff] using denIdeal_ne_bot u
  have hdpos : 0 < d := Nat.pos_of_ne_zero hdne
  -- `d` は `𝔡_u` の元
  have hmem : ((d : ℕ) : 𝓞 K) ∈ denIdeal u := Ideal.absNorm_mem (denIdeal u)
  have hdne' : ((d : ℕ) : 𝓞 K) ≠ 0 := Nat.cast_ne_zero.2 hdne
  have hint := mul_isIntegral_of_mem_denIdeal (x := u) hdne' hmem
  rw [map_natCast] at hint
  obtain ⟨z, hz⟩ := hint
  -- `d·y` は整元で共役が `B²` 以下
  have hzy : algebraMap (𝓞 K) K z = (d : K) * y := hz
  have hzint : IsIntegral ℤ ((d : K) * y) := by
    rw [← hzy]
    exact NumberField.RingOfIntegers.isIntegral_coe z
  have hnorm : ∀ φ : K →+* ℂ, ‖φ ((d : K) * y)‖ ≤ B * B := by
    intro φ
    have hy1 : ‖φ y‖ ≤ B := by
      calc ‖φ y‖ = (InfinitePlace.mk φ) y := (InfinitePlace.apply φ y).symm
        _ ≤ Height.mulHeight₁ y := infinitePlace_le_mulHeight₁ _ y
        _ ≤ B := hy
    have hBnn : (0 : ℝ) ≤ B := le_trans (Nat.cast_nonneg d) hdle
    calc ‖φ ((d : K) * y)‖ = (d : ℝ) * ‖φ y‖ := by
          rw [map_mul, norm_mul, map_natCast, Complex.norm_natCast]
      _ ≤ B * B := mul_le_mul hdle hy1 (norm_nonneg _) hBnn
  refine ⟨((d : K) * y, d), ⟨⟨hzint, hnorm⟩, ?_⟩, ?_⟩
  · exact Nat.le_floor hdle
  · have hdK : (d : K) ≠ 0 := Nat.cast_ne_zero.2 hdne
    field_simp

/-- ★★★**`Northcott (mulHeight₁ (K := K))` の基底 instance**。

★★mathlib は `Height/Northcott.lean` に「`mulHeight₁` で Northcott なら
`logHeight₁` でも」という**条件つき instance を 1 つ持つだけ**で、
基底 instance を 1 つも持たない(2026-08-17 実測)。**これがそれである。**

★これが入ると mathlib 側の条件つき instance が発火して、
`Northcott (logHeight₁ (K := K))` も自動で付く。

★`Found/GenEll/NorthcottRat.lean` の `northcott_mulHeight₁_rat`(`ℚ` の場合)は
本 instance の特別な場合になった。 -/
instance northcott_mulHeight₁ : Northcott (Height.mulHeight₁ (K := K)) where
  finite_le B := finite_mulHeight₁_le B

/-! ## ★出典の紐付け(`.src`) -/

def finite_mulHeight₁_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(数体上の Northcott 性——高さ理論の側のみ)",
    sectionId := "genell-prop-1-4" }

def denIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(分母の評価)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
