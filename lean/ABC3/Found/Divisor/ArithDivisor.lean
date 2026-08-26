/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import Mathlib.NumberTheory.NumberField.Units.Basic

/-!
# 算術因子 —— `deg^arith` が `B(L)` の像で消えること(`Example 6.3`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113–114。

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group

原文 (FrdI p.114):
> vanishes on the

## ★★何を閉じたか

依存グラフ(`ResearchPaper/frdi-decomposition.json`)の鎖 `arith` のうち
**`deg-vanish`** ——「`deg^arith` が `B(L)` の像で消える」—— を閉じた。
★原文は `Example 6.3` の末尾で「one verifies immediately that `deg^arith_L`
vanishes on the image of `B(L)`」と一言で畳むが、その中身は**積公式**である。

★★mathlib に `NumberField.prod_abs_eq_one` が在る:

  `(∏ w : InfinitePlace, w x ^ w.mult) * ∏ᶠ (w : FinitePlace), w x = 1`

両辺の**対数**を取れば `Σ_v log|x|_v = 0` になり、
`deg^arith(div(x)) = −Σ_v log|x|_v = 0` が出る。

## ★設計 —— 対数を先に取る

原文は `ord(F_v)` を `F_v^×/O_v^×` として定め、非アルキメデスなら `≅ ℤ`、
アルキメデスなら `≅ ℝ` と述べ、`deg^arith` を素点ごとに
「`−[F_v:ℝ]·log λ`」「`log #(O_v/(λ))`」で定める。

★★**本ファイルは `deg^arith` の値の側を直接扱う** —— すなわち各素点の成分を
`−log|x|_v`(無限素点は重複度 `mult` つき)と置く。
非アルキメデスでは `log|x|_v = −ord_v(x)·log(Nv)` なので、
これは原文の `ord_v(x)·log(Nv)` と**同じ数**である。

★★★**逸脱(記録)**: したがって本ファイルの `ArithPlace L →₀ ℝ` は
原文の `Φ(L)^gp = ⊕_v ord(F_v)`(有限素点成分が `ℤ`)を
**`deg^arith` で見た像**に相当する。`ord(F_v) ≅ ℤ` の同一視そのものは
鎖 `arith` の節点 `ord-mon` であり、まだ実装していない
(`Skeleton/Divisor/ArithDivisor.lean` に statement がある)。
-/

namespace ABC3.Found.Divisor

open NumberField
open scoped Classical

variable {L : Type*} [Field L] [NumberField L]

/-! ## ★1. 素点と `−log|x|_v` -/

/-- ★**`V(L)`** —— 有限素点と無限素点の合併。 -/
abbrev ArithPlace (L : Type*) [Field L] [NumberField L] : Type _ :=
  FinitePlace L ⊕ InfinitePlace L

/-- ★各素点での `−log|x|_v`(無限素点は重複度 `mult` つき)。

★非アルキメデスでは `= ord_v(x)·log(Nv)`、
アルキメデスでは `= −[F_v:ℝ]·log|x|_v` である。 -/
noncomputable def arithPlaceLog (x : L) : ArithPlace L → ℝ :=
  Sum.elim (fun w : FinitePlace L => -Real.log (w x))
    (fun w : InfinitePlace L => -(w.mult : ℝ) * Real.log (w x))

/-- ★`|x|_v ≠ 1` となる有限素点は有限個(mathlib)。 -/
noncomputable def arithFinSupport {x : L} (hx : x ≠ 0) : Finset (FinitePlace L) :=
  (FinitePlace.hasFiniteMulSupport hx).toFinset

theorem prod_arithFinSupport {x : L} (hx : x ≠ 0) :
    ∏ᶠ (w : FinitePlace L), w x = ∏ w ∈ arithFinSupport hx, w x :=
  finprod_eq_prod_of_mulSupport_toFinset_subset _ (FinitePlace.hasFiniteMulSupport hx)
    (Finset.Subset.refl _)

/-- ★★★**積公式の対数版** —— `Σ_v log|x|_v = 0`。 -/
theorem sum_log_place_eq_zero {x : L} (hx : x ≠ 0) :
    (∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w x))
      + ∑ w ∈ arithFinSupport hx, Real.log (w x) = 0 := by
  have hpf := NumberField.prod_abs_eq_one hx
  rw [prod_arithFinSupport hx] at hpf
  have h1 : ∀ w : InfinitePlace L, w x ≠ 0 := fun w => (InfinitePlace.pos_iff.mpr hx).ne'
  have h2 : ∀ w : FinitePlace L, w x ≠ 0 := fun w => (FinitePlace.pos_iff.mpr hx).ne'
  have hA : (∏ w : InfinitePlace L, (w x) ^ w.mult) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun w _ => pow_ne_zero _ (h1 w))
  have hB : (∏ w ∈ arithFinSupport hx, w x) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun w _ => h2 w)
  have hlog := congrArg Real.log hpf
  rw [Real.log_mul hA hB, Real.log_one] at hlog
  rw [Real.log_prod (f := fun w : InfinitePlace L => (w x) ^ w.mult)
      (fun w _ => pow_ne_zero _ (h1 w)),
    Real.log_prod (f := fun w : FinitePlace L => w x) (fun w _ => h2 w)] at hlog
  simpa [Real.log_pow] using hlog

/-! ## ★2. 算術因子 `div(x)` -/

theorem finite_support_arithPlaceLog {x : L} (hx : x ≠ 0) :
    (Function.support (arithPlaceLog x)).Finite := by
  refine Set.Finite.subset
    ((Set.Finite.image Sum.inl (FinitePlace.hasFiniteMulSupport hx)).union
      (Set.finite_range (Sum.inr : InfinitePlace L → ArithPlace L))) ?_
  rintro (w | w) hw
  · refine Or.inl ⟨w, ?_, rfl⟩
    simp only [Function.mem_mulSupport]
    intro h
    exact hw (by simp [arithPlaceLog, h])
  · exact Or.inr ⟨w, rfl⟩

/-- ★★**`x ∈ L^×` の算術因子** —— 各素点の成分が `−log|x|_v`。 -/
noncomputable def arithDiv (x : Lˣ) : ArithPlace L →₀ ℝ :=
  Finsupp.ofSupportFinite (arithPlaceLog (x : L)) (finite_support_arithPlaceLog x.ne_zero)

@[simp] theorem arithDiv_apply (x : Lˣ) (v : ArithPlace L) :
    arithDiv x v = arithPlaceLog (x : L) v := rfl

theorem arithPlaceLog_mul {x y : L} (hx : x ≠ 0) (hy : y ≠ 0) (v : ArithPlace L) :
    arithPlaceLog (x * y) v = arithPlaceLog x v + arithPlaceLog y v := by
  rcases v with w | w
  · show -Real.log (w (x * y)) = -Real.log (w x) + -Real.log (w y)
    rw [map_mul, Real.log_mul (FinitePlace.pos_iff.mpr hx).ne' (FinitePlace.pos_iff.mpr hy).ne']
    ring
  · show -(w.mult : ℝ) * Real.log (w (x * y))
      = -(w.mult : ℝ) * Real.log (w x) + -(w.mult : ℝ) * Real.log (w y)
    rw [map_mul, Real.log_mul (InfinitePlace.pos_iff.mpr hx).ne'
      (InfinitePlace.pos_iff.mpr hy).ne']
    ring

/-- ★★**`B(L) → Φ(L)^gp` は準同型**。 -/
theorem arithDiv_mul (x y : Lˣ) : arithDiv (x * y) = arithDiv x + arithDiv y := by
  refine Finsupp.ext (fun v => ?_)
  show arithPlaceLog ((x : L) * (y : L)) v = arithPlaceLog (x : L) v + arithPlaceLog (y : L) v
  exact arithPlaceLog_mul x.ne_zero y.ne_zero v

/-- ★★`B(L) = L^× → Φ(L)^gp` を群準同型として。 -/
noncomputable def arithDivHom : Lˣ →* Multiplicative (ArithPlace L →₀ ℝ) where
  toFun x := Multiplicative.ofAdd (arithDiv x)
  map_one' := by
    refine congrArg Multiplicative.ofAdd (Finsupp.ext (fun v => ?_))
    show arithPlaceLog ((1 : Lˣ) : L) v = 0
    rcases v with w | w <;> simp [arithPlaceLog]
  map_mul' x y := congrArg Multiplicative.ofAdd (arithDiv_mul x y)

/-! ## ★3. 算術次数と、その `B(L)` 上での消滅 -/

/-- ★**算術次数** `deg^arith : Φ(L)^gp → ℝ` —— 成分の総和。 -/
noncomputable def arithDegree : (ArithPlace L →₀ ℝ) →+ ℝ :=
  Finsupp.liftAddHom (fun _ : ArithPlace L => (AddMonoidHom.id ℝ))

@[simp] theorem arithDegree_apply (d : ArithPlace L →₀ ℝ) :
    arithDegree d = d.sum (fun _ r => r) := rfl

/-- ★台を含む有限集合。 -/
noncomputable def arithBigSupport {x : L} (hx : x ≠ 0) : Finset (ArithPlace L) :=
  (arithFinSupport hx).image Sum.inl ∪ (Finset.univ : Finset (InfinitePlace L)).image Sum.inr

theorem support_subset_arithBigSupport (x : Lˣ) :
    (arithDiv x).support ⊆ arithBigSupport x.ne_zero := by
  intro v hv
  rw [Finsupp.mem_support_iff, arithDiv_apply] at hv
  rcases v with w | w
  · refine Finset.mem_union_left _ (Finset.mem_image.mpr ⟨w, ?_, rfl⟩)
    simp only [arithFinSupport, Set.Finite.mem_toFinset, Function.mem_mulSupport]
    intro h
    exact hv (by simp [arithPlaceLog, h])
  · exact Finset.mem_union_right _ (Finset.mem_image.mpr ⟨w, Finset.mem_univ w, rfl⟩)

/-- ★★★★★**[FrdI] Example 6.3** —— **`deg^arith` は `B(L)` の像で消える**。

原文 (FrdI p.114):
> image of B(L).

★★これは**積公式**そのものである(`NumberField.prod_abs_eq_one` の対数版)。 -/
theorem arithDegree_arithDiv (x : Lˣ) : arithDegree (arithDiv x) = 0 := by
  have hx : (x : L) ≠ 0 := x.ne_zero
  have hsum : arithDegree (arithDiv x)
      = ∑ v ∈ arithBigSupport hx, arithPlaceLog (x : L) v := by
    rw [arithDegree_apply]
    exact Finsupp.sum_of_support_subset _ (support_subset_arithBigSupport x) _ (fun _ _ => rfl)
  rw [hsum, arithBigSupport]
  have hdisj : Disjoint ((arithFinSupport hx).image (Sum.inl : FinitePlace L → ArithPlace L))
      ((Finset.univ : Finset (InfinitePlace L)).image Sum.inr) := by
    rw [Finset.disjoint_left]
    rintro v hv1 hv2
    obtain ⟨a, _, rfl⟩ := Finset.mem_image.mp hv1
    obtain ⟨b, _, hb⟩ := Finset.mem_image.mp hv2
    cases hb
  rw [Finset.sum_union hdisj,
    Finset.sum_image (fun a _ b _ h => Sum.inl_injective h),
    Finset.sum_image (fun a _ b _ h => Sum.inr_injective h)]
  have h := sum_log_place_eq_zero hx
  show (∑ w ∈ arithFinSupport hx, -Real.log (w (x : L)))
      + ∑ w : InfinitePlace L, -(w.mult : ℝ) * Real.log (w (x : L)) = 0
  have h1 : (∑ w ∈ arithFinSupport hx, -Real.log (w (x : L)))
      = -∑ w ∈ arithFinSupport hx, Real.log (w (x : L)) := by
    rw [← Finset.sum_neg_distrib]
  have h2 : ∑ w : InfinitePlace L, -(w.mult : ℝ) * Real.log (w (x : L))
      = -∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w (x : L)) := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl (fun w _ => by ring)
  rw [h1, h2]
  linarith [h]

/-! ## ★4. 有限素点の `ord` との一致

原文 (FrdI p.114):
> the natural logarithm of the cardinality of the finite set Ov/(λ) [where Ov is the

★★原文は非アルキメデスの `deg^arith` を「`log #(O_v/(λ))`」と定める。
`#(O_v/(λ)) = (N v)^{ord_v(λ)}` なので、これは `ord_v(λ)·log(N v)` であり、
本ファイルの `arithPlaceLog x (inl w) = −log|x|_w` と**一致する**。
★これでファイル冗頭に記した逸脱（対数を先に取る形を採ったこと）の
**非アルキメデス側が埋まる**（鎖 `arith` の `ord-mon` の半分）。 -/

/-- ★**有限素点での `ord`** —— `v.valuation` の指数の符号を反転したもの。 -/
noncomputable def ordFin (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L)) (x : L) : ℤ :=
  if h : v.valuation L x = 0 then 0 else -(WithZero.unzero h).toAdd

/-- ★★**`|x|_v = (N v)^{-ord_v(x)}`**。 -/
theorem finitePlace_eq_absNorm_zpow (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L))
    (x : L) (hx : x ≠ 0) :
    ((FinitePlace.mk v) x : ℝ) = (Ideal.absNorm v.asIdeal : ℝ) ^ (-(ordFin v x)) := by
  have hval : v.valuation L x ≠ 0 := by
    simp only [ne_eq, Valuation.zero_iff]
    exact hx
  rw [FinitePlace.mk_apply, FinitePlace.norm_embedding,
    NumberField.RingOfIntegers.HeightOneSpectrum.adicAbv_def]
  rw [WithZeroMulInt.toNNReal_neg_apply
    (NumberField.RingOfIntegers.HeightOneSpectrum.absNorm_ne_zero v) hval]
  rw [ordFin, dif_neg hval, neg_neg]
  push_cast
  rfl

/-- ★★★**原文の `log #(O_v/(λ))` との一致** ——
`arithPlaceLog x (inl w) = ord_v(x)·log(N v)`。 -/
theorem arithPlaceLog_finite (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L))
    (x : L) (hx : x ≠ 0) :
    arithPlaceLog x (Sum.inl (FinitePlace.mk v))
      = (ordFin v x : ℝ) * Real.log (Ideal.absNorm v.asIdeal : ℝ) := by
  show -Real.log ((FinitePlace.mk v) x) = _
  rw [finitePlace_eq_absNorm_zpow v x hx, Real.log_zpow]
  push_cast
  ring

/-- ★`log(N v) > 0`（`N v > 1` は mathlib）。 -/
theorem log_absNorm_pos (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L)) :
    (0:ℝ) < Real.log (Ideal.absNorm v.asIdeal : ℝ) := by
  refine Real.log_pos ?_
  exact_mod_cast NumberField.RingOfIntegers.HeightOneSpectrum.one_lt_absNorm v

/-- ★★**`ord_v` は加法的** —— `arithPlaceLog` の加法性を `log(N v)` で割る。 -/
theorem ordFin_mul (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L)) {x y : L}
    (hx : x ≠ 0) (hy : y ≠ 0) : ordFin v (x * y) = ordFin v x + ordFin v y := by
  have hlog := log_absNorm_pos (L := L) v
  have h := arithPlaceLog_mul hx hy (Sum.inl (FinitePlace.mk v))
  rw [arithPlaceLog_finite v _ (mul_ne_zero hx hy), arithPlaceLog_finite v x hx,
    arithPlaceLog_finite v y hy] at h
  have h2 : ((ordFin v (x*y) : ℝ)) = (ordFin v x : ℝ) + (ordFin v y : ℝ) := by
    field_simp at h
    linarith [h]
  exact_mod_cast h2

/-- ★★★**一次元子の存在** —— `ord_v(π) = 1` なる `π` がある。

★mathlib の `HeightOneSpectrum.valuation_exists_uniformizer`。
★★これと `ordFin_mul` で **`ord(F_v) ≅ ℤ`**（鎖 `arith` の `ord-mon`）が出る。 -/
theorem exists_ordFin_eq_one (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L)) :
    ∃ π : L, π ≠ 0 ∧ ordFin v π = 1 := by
  obtain ⟨π, hπ⟩ := v.valuation_exists_uniformizer L
  have hne : v.valuation L π ≠ 0 := by rw [hπ]; simp
  refine ⟨π, fun h => hne (by rw [h, map_zero]), ?_⟩
  rw [ordFin, dif_neg hne]
  have h2 : (WithZero.unzero hne) = Multiplicative.ofAdd (-1 : ℤ) := by
    rw [← WithZero.coe_inj, WithZero.coe_unzero, hπ]
    rfl
  rw [h2]
  rfl

/-! ## ★5. `Φ(F)` と `Φ(F)^gp`（鎖 `arith` の `arith-phi`）

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group

★★対数の模型では、`ord(O_v^▷)` は
* 非アルキメデス: `ℕ·log(N v)`（`≅ ℕ`）
* アルキメデス: `ℝ≥0`
に対応する。そこで `Φ(F)` を「各成分が非負で、非アルキメデス成分は
`log(N v)` の非負整数倍」として定める。
★`Φ(F)^gp` は「整数倍（符号は問わない）」である。 -/

/-- ★★**`Φ(F)^gp`** —— 非アルキメデス成分が `log(N v)` の整数倍であるもの。 -/
def arithDivGroup (L : Type*) [Field L] [NumberField L] :
    AddSubgroup (ArithPlace L →₀ ℝ) where
  carrier := {d | ∀ w : FinitePlace L, ∃ n : ℤ,
    d (Sum.inl w) = (n : ℝ) * Real.log (Ideal.absNorm (FinitePlace.maximalIdeal w).asIdeal : ℝ)}
  zero_mem' := fun w => ⟨0, by simp⟩
  add_mem' {a b} ha hb := by
    intro w
    obtain ⟨m, hm⟩ := ha w
    obtain ⟨n, hn⟩ := hb w
    refine ⟨m + n, ?_⟩
    show a (Sum.inl w) + b (Sum.inl w) = _
    rw [hm, hn]
    push_cast
    ring
  neg_mem' {a} ha := by
    intro w
    obtain ⟨m, hm⟩ := ha w
    refine ⟨-m, ?_⟩
    show -a (Sum.inl w) = _
    rw [hm]
    push_cast
    ring

/-- ★★**`Φ(F)`** —— `Φ(F)^gp` のうち各成分が非負なもの（有効算術因子）。 -/
def arithEff (L : Type*) [Field L] [NumberField L] :
    AddSubmonoid (ArithPlace L →₀ ℝ) where
  carrier := {d | d ∈ arithDivGroup L ∧ ∀ v, 0 ≤ d v}
  zero_mem' := ⟨(arithDivGroup L).zero_mem, by simp⟩
  add_mem' {a b} ha hb :=
    ⟨(arithDivGroup L).add_mem ha.1 hb.1, fun v => by
      rw [Finsupp.add_apply]
      exact add_nonneg (ha.2 v) (hb.2 v)⟩

/-- ★★★**`B(L) → Φ(L)^gp`** —— 有理函数の因子は実際に `Φ(L)^gp` に入る。

★非アルキメデス成分が `ord_v(x)·log(N v)` であること（`arithPlaceLog_finite`）から直ちに出る。 -/
theorem arithDiv_mem_arithDivGroup (x : Lˣ) : arithDiv x ∈ arithDivGroup L := by
  intro w
  refine ⟨ordFin (FinitePlace.maximalIdeal w) (x : L), ?_⟩
  have h := arithPlaceLog_finite (FinitePlace.maximalIdeal w) (x : L) x.ne_zero
  rw [FinitePlace.mk_maximalIdeal] at h
  exact h

/-! ## ★6. `𝓞^×(A) = 𝓞^▷(A) = μ(L)`（鎖 `arith` の `mu-units`）

原文 (FrdI p.113):
> to Spec(L) ∈Ob(B(G)0), we have

★★`div(x) = 0` ならすべての素点で `|x|_v = 1` であり、
非アルキメデス側から `x ∈ (𝓞 L)ˣ`、アルキメデス側と mathlib の
`NumberField.Units.mem_torsion` から `x` は 1 の冪根になる。 -/

/-- ★`ord_v(x) = 0` なら `v.valuation x = 1`。 -/
theorem valuation_eq_one_of_ordFin_eq_zero (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L))
    {x : L} (hx : x ≠ 0) (h : ordFin v x = 0) : v.valuation L x = 1 := by
  have hne : v.valuation L x ≠ 0 := by
    simp only [ne_eq, Valuation.zero_iff]; exact hx
  rw [ordFin, dif_neg hne, neg_eq_zero] at h
  rw [← WithZero.coe_unzero hne]
  rw [show WithZero.unzero hne = 1 from
    Multiplicative.toAdd.injective (by simpa using h)]
  rfl

/-- ★`div(x) = 0` ならどの有限素点でも `ord_v(x) = 0`。 -/
theorem ordFin_eq_zero_of_arithDiv_eq_zero {x : Lˣ} (h : arithDiv x = 0)
    (v : IsDedekindDomain.HeightOneSpectrum (𝓞 L)) : ordFin v (x : L) = 0 := by
  have h1 : arithPlaceLog (x : L) (Sum.inl (FinitePlace.mk v)) = 0 := by
    have := congrFun (congrArg (fun d : ArithPlace L →₀ ℝ => (d : ArithPlace L → ℝ)) h)
      (Sum.inl (FinitePlace.mk v))
    simpa using this
  rw [arithPlaceLog_finite v _ x.ne_zero] at h1
  have hlog := log_absNorm_pos (L := L) v
  have h2 : (ordFin v (x:L) : ℝ) = 0 := by
    rcases mul_eq_zero.mp h1 with h3 | h3
    · exact h3
    · exact absurd h3 hlog.ne'
  exact_mod_cast h2

/-- ★`div(x) = 0` ならどの無限素点でも `|x|_w = 1`。 -/
theorem infinitePlace_eq_one_of_arithDiv_eq_zero {x : Lˣ} (h : arithDiv x = 0)
    (w : InfinitePlace L) : w (x : L) = 1 := by
  have h1 : arithPlaceLog (x : L) (Sum.inr w) = 0 := by
    have := congrFun (congrArg (fun d : ArithPlace L →₀ ℝ => (d : ArithPlace L → ℝ)) h)
      (Sum.inr w)
    simpa using this
  have h2 : -(w.mult : ℝ) * Real.log (w (x : L)) = 0 := h1
  have hm : (w.mult : ℝ) ≠ 0 := by
    exact_mod_cast (InfinitePlace.mult_pos (w := w)).ne'
  have h3 : Real.log (w (x : L)) = 0 := by
    rcases mul_eq_zero.mp h2 with h4 | h4
    · exact absurd (neg_eq_zero.mp h4) hm
    · exact h4
  have hpos : 0 < w (x : L) := InfinitePlace.pos_iff.mpr x.ne_zero
  have := Real.exp_log hpos
  rw [h3, Real.exp_zero] at this
  exact this.symm

/-- ★★★★★**[FrdI] Example 6.3** —— `𝓞^×(A) = 𝓞^▷(A) = μ(L)`。

★`div(x) = 0` なら `x` は 1 の冪根である。
★非アルキメデス側で `x ∈ (𝓞 L)ˣ` を出し、アルキメデス側と
mathlib の `NumberField.Units.mem_torsion` で有限位数にする。 -/
theorem exists_pow_eq_one_of_arithDiv_eq_zero (x : Lˣ) (h : arithDiv x = 0) :
    ∃ n : ℕ, 0 < n ∧ (x : L) ^ n = 1 := by
  have hv : ∀ v : IsDedekindDomain.HeightOneSpectrum (𝓞 L), v.valuation L (x : L) = 1 :=
    fun v => valuation_eq_one_of_ordFin_eq_zero v x.ne_zero
      (ordFin_eq_zero_of_arithDiv_eq_zero h v)
  have hvi : ∀ v : IsDedekindDomain.HeightOneSpectrum (𝓞 L),
      v.valuation L ((x : L)⁻¹) = 1 := by
    intro v
    rw [map_inv₀, hv v, inv_one]
  obtain ⟨a, ha⟩ := IsDedekindDomain.HeightOneSpectrum.mem_integers_of_valuation_le_one
    (R := 𝓞 L) L (x : L) (fun v => le_of_eq (hv v))
  obtain ⟨b, hb⟩ := IsDedekindDomain.HeightOneSpectrum.mem_integers_of_valuation_le_one
    (R := 𝓞 L) L ((x : L)⁻¹) (fun v => le_of_eq (hvi v))
  have hinj : Function.Injective (algebraMap (𝓞 L) L) :=
    IsFractionRing.injective (𝓞 L) L
  have hab : a * b = 1 := by
    refine hinj ?_
    rw [map_mul, ha, hb, map_one]
    field_simp
  have hba : b * a = 1 := by rw [mul_comm]; exact hab
  set u : (𝓞 L)ˣ := ⟨a, b, hab, hba⟩ with hu
  have hua : algebraMap (𝓞 L) L (u : 𝓞 L) = (x : L) := ha
  have htor : u ∈ _root_.NumberField.Units.torsion L := by
    refine (_root_.NumberField.Units.mem_torsion (K := L) (x := u)).mpr (fun w => ?_)
    have := infinitePlace_eq_one_of_arithDiv_eq_zero h w
    rwa [← hua] at this
  obtain ⟨n, hn, hpow⟩ := isOfFinOrder_iff_pow_eq_one.mp htor
  refine ⟨n, hn, ?_⟩
  have := congrArg (fun t : (𝓞 L)ˣ => algebraMap (𝓞 L) L (t : 𝓞 L)) hpow
  simpa [hua, map_pow] using this

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.3` の素点の集合 `V(F)`。 -/
def ArithPlace.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 112, item := "Example 6.3 — V(F)(素点の集合)",
    sectionId := "frdi-example-6-3" }

/-- ★locator —— `Example 6.3` の `B(F) = F^× → Φ(F)^gp`。 -/
def arithDivHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — B(F) = F^× → Φ(F)^gp",
    sectionId := "frdi-example-6-3" }

/-- ★locator —— `Example 6.3` の `deg^arith_L`。 -/
def arithDegree.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — deg^arith_L : Φ(L)^gp → ℝ",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の「`deg^arith` は `B(L)` の像で消える」(=積公式)。 -/
def arithDegree_arithDiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114, item := "Example 6.3 — deg^arith は B(L) の像で消える",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の非アルキメデス成分が `log #(O_v/(λ))` と一致すること。 -/
def arithPlaceLog_finite.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 非アルキメデスで deg^arith = ord_v·log(Nv)",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の `ord(F_v) ≅ ℤ`（一次元子の存在と加法性）。 -/
def exists_ordFin_eq_one.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — ord(F_v) ≅ ℤ（非アルキメデス）",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の `Φ(F)` と `Φ(F)^gp`。 -/
def arithDivGroup.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(F) と Φ(F)^gp",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の `𝓞^×(A) = 𝓞^▷(A) = μ(L)`。 -/
def exists_pow_eq_one_of_arithDiv_eq_zero.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — O^×(A) = O^▷(A) = μ(L)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
