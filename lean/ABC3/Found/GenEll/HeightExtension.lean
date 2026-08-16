import ABC3.Found.GenEll.BaseChange
import ABC3.Found.GenEll.HeightADiv
import Mathlib.NumberTheory.RamificationInertia.Valuation

/-!
# [GenEll] §1 —— **高さの拡大公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★何を取るか

> **`x ∈ K ⊆ L` について `logHeight₁_L(x) = [L:K] · logHeight₁_K(x)`**

★★mathlib は**高さの体拡大に関する公式を持たない**(2026-08-17 実測——
`Mathlib/NumberTheory/Height/` を `algebraMap` / `extension` / `finrank` で検索して 0 件)。

★これが無いと「**絶対**高さ」——体に依らない量——が定義できない。
原文 `Proposition 1.4, (iv)` が `X(ℚ̄)^{≤d}` を走らせるのは絶対高さの話であり、
その基礎が本ファイルである。

## ★★機構 —— 今朝の 2 つの相対和がそのまま要る

- **アルキメデス側**: `Σ_{w|v} mult_w = mult_v·[L:K]`
  (`InfinitePlaceRel.lean` の `sum_mult_comap_eq`。★mathlib は絶対版しか持たない)
- **非アルキメデス側**: `Σ_{w|v} e(w|v)·f(w|v) = [L:K]`
  (`FinitePlaceRel.lean` の `sum_ramification_inertia_hos`。
  ★mathlib の基本等式を `HeightOneSpectrum` に翻訳したもの)

★★**`degNormalized` の底変換不変性のために作った道具が、そのまま高さに効いた。**
どちらも「素点の族を `[L:K]` 倍に数える」という同じ形だからである。

## ★局所の関係

- `w(x) = v(x)^{e(w|v)·f(w|v)}`(`x ∈ K`)
  ——`ord_w(x) = e·ord_v(x)`(mathlib の
  `IsDedekindDomain.HeightOneSpectrum.valuation_liesOver`、2026 年に入ったもの)と
  `q_w = q_v^{f}`(`Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver`)から。
- `w(x) = (w|_K)(x)`(アルキメデス、`InfinitePlace.comap` の定義)
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

section HeightExtension

variable (K L : Type*) [Field K] [NumberField K] [Field L] [NumberField L] [Algebra K L]

/-- `x ∈ Kˣ` の `Lˣ` への像。 -/
noncomputable def unitsExt (x : Kˣ) : Lˣ := Units.map (algebraMap K L).toMonoidHom x

@[simp] theorem unitsExt_coe (x : Kˣ) :
    ((unitsExt K L x : Lˣ) : L) = algebraMap K L (x : K) := rfl

/-- ★**`ord_w(x) = e(w|v)·ord_v(x)`**(`x ∈ K`、`v = w` の引き戻し)。

★mathlib の `valuation_liesOver` の `toAdd` を取っただけ。
`LiesOver` インスタンスは `hosComap_liesOver`(`FinitePlaceRel.lean`)が供給する。 -/
theorem ordv_unitsExt (W : FinitePlace L) (x : Kˣ) :
    ordv W (unitsExt K L x)
      = (ramIdxOver K L W : ℤ) * ordv (hosComap K L W) x := by
  set V := hosComap K L W with hV
  have hval : (V.valuation K ((x : K))) ^ (V.asIdeal.ramificationIdx W.asIdeal)
      = W.valuation L (algebraMap K L ((x : K))) :=
    IsDedekindDomain.HeightOneSpectrum.valuation_liesOver L V W ((x : K))
  have hneV : V.valuation K ((x : K)) ≠ 0 :=
    (Valuation.ne_zero_iff (V.valuation K)).2 (Units.ne_zero x)
  have hneW : W.valuation L (((unitsExt K L x : Lˣ) : L)) ≠ 0 :=
    (Valuation.ne_zero_iff (W.valuation L)).2 (Units.ne_zero _)
  have hunz : WithZero.unzero hneW
      = (WithZero.unzero hneV) ^ (V.asIdeal.ramificationIdx W.asIdeal) := by
    apply WithZero.coe_injective
    rw [WithZero.coe_unzero, WithZero.coe_pow, WithZero.coe_unzero, hval]
    rfl
  simp only [ordv, hunz, ramIdxOver, hV]
  rw [toAdd_pow]
  ring

/-- ★**`w(x) = v(x)^{e·f}`**(`x ∈ K`)。

★`ord` の関係(分岐指数)と剰余体の大きさ(剰余次数)が**両方**効く。 -/
theorem posLog_finitePlace_unitsExt (w : NumberField.FinitePlace L) (x : Kˣ) :
    Real.posLog (w (((unitsExt K L x : Lˣ) : L)))
      = ((ramIdxOver K L (NumberField.FinitePlace.maximalIdeal w) *
            inertiaOver K L (NumberField.FinitePlace.maximalIdeal w) : ℕ) : ℝ)
          * (((max (-(ordv (hosComap K L (NumberField.FinitePlace.maximalIdeal w)) x)) 0 : ℤ) : ℝ)
              * Real.log (residueCard
                  (hosComap K L (NumberField.FinitePlace.maximalIdeal w)))) := by
  set W := NumberField.FinitePlace.maximalIdeal w with hW
  have hq : (0 : ℝ) ≤ Real.log (residueCard (hosComap K L W)) :=
    (log_residueCard_pos (hosComap K L W)).le
  have hef : (0 : ℝ) ≤ (ramIdxOver K L W : ℝ) := Nat.cast_nonneg _
  have hin : (0 : ℝ) ≤ (inertiaOver K L W : ℝ) := Nat.cast_nonneg _
  have hef2 : (0 : ℝ) ≤ (ramIdxOver K L W : ℝ) * (inertiaOver K L W : ℝ) :=
    mul_nonneg hef hin
  rw [Real.posLog, log_finitePlace_apply w (unitsExt K L x), ordv_unitsExt K L W x,
    log_residueCard_over K L W]
  push_cast
  rcases le_or_gt (0 : ℝ) (-(ordv (hosComap K L W) x : ℝ)) with h | h
  · rw [max_eq_left h, max_eq_right ?_]
    · ring
    · have hac : (0 : ℝ) ≤ -(ordv (hosComap K L W) x : ℝ)
          * Real.log (residueCard (hosComap K L W)) := mul_nonneg h hq
      nlinarith [mul_nonneg hef2 hac]
  · rw [max_eq_right h.le, max_eq_left ?_]
    · ring
    · have hac : -(ordv (hosComap K L W) x : ℝ)
          * Real.log (residueCard (hosComap K L W)) ≤ 0 :=
        mul_nonpos_iff.mpr (Or.inr ⟨h.le, hq⟩)
      nlinarith [mul_nonneg hef2 (neg_nonneg.2 hac)]

/-! ## ★有限素点側 -/

open scoped Classical in
/-- ★★**有限素点側は `[L:K]` 倍になる**。

★`Σ_{w|v} e(w|v)·f(w|v) = [L:K]`(`FinitePlaceRel.lean`)そのものである。 -/
theorem finsum_posLog_extension (x : Kˣ) :
    ∑ᶠ w : NumberField.FinitePlace L, Real.posLog (w (((unitsExt K L x : Lˣ) : L)))
      = (Module.finrank K L : ℝ)
        * ∑ᶠ v : NumberField.FinitePlace K, Real.posLog (v ((x : K))) := by
  classical
  -- 右辺を `HeightOneSpectrum` の言葉へ
  have hRHS : ∑ᶠ v : NumberField.FinitePlace K, Real.posLog (v ((x : K)))
      = ∑ᶠ V : FinitePlace K, ((max (-(ordv V x)) 0 : ℤ) : ℝ) * Real.log (residueCard V) := by
    rw [finsum_congr (fun v => posLog_finitePlace_apply v x)]
    exact finsum_comp_equiv (M := ℝ)
      (e := NumberField.FinitePlace.equivHeightOneSpectrum (K := K))
      (f := fun V : FinitePlace K => ((max (-(ordv V x)) 0 : ℤ) : ℝ) * Real.log (residueCard V))
  -- 左辺を `HeightOneSpectrum` の言葉へ
  have hLHS : ∑ᶠ w : NumberField.FinitePlace L, Real.posLog (w (((unitsExt K L x : Lˣ) : L)))
      = ∑ᶠ W : FinitePlace L,
          ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
            * (((max (-(ordv (hosComap K L W) x)) 0 : ℤ) : ℝ)
                * Real.log (residueCard (hosComap K L W))) := by
    rw [finsum_congr (fun w => posLog_finitePlace_unitsExt K L w x)]
    exact finsum_comp_equiv (M := ℝ)
      (e := NumberField.FinitePlace.equivHeightOneSpectrum (K := L))
      (f := fun W : FinitePlace L => ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
        * (((max (-(ordv (hosComap K L W) x)) 0 : ℤ) : ℝ)
            * Real.log (residueCard (hosComap K L W))))
  rw [hLHS, hRHS]
  -- 有限和へ落とす
  set S : Finset (FinitePlace K) := (ordv_finite_support x).toFinset with hS
  set T : Finset (FinitePlace L) := S.biUnion (fun V => hosFiber K L V) with hT
  have hGsupp : (Function.support
      (fun V : FinitePlace K =>
        ((max (-(ordv V x)) 0 : ℤ) : ℝ) * Real.log (residueCard V))) ⊆ ↑S := by
    intro V hV
    rw [Function.mem_support] at hV
    have h0 : ordv V x ≠ 0 := by
      intro hcon
      exact hV (by simp [hcon])
    simpa [hS, Set.Finite.mem_toFinset] using h0
  have hFsupp : (Function.support
      (fun W : FinitePlace L => ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
        * (((max (-(ordv (hosComap K L W) x)) 0 : ℤ) : ℝ)
            * Real.log (residueCard (hosComap K L W))))) ⊆ ↑T := by
    intro W hW
    rw [Function.mem_support] at hW
    have h0 : ordv (hosComap K L W) x ≠ 0 := by
      intro hcon
      exact hW (by simp [hcon])
    rw [Finset.mem_coe, hT, Finset.mem_biUnion]
    exact ⟨hosComap K L W, by simpa [hS, Set.Finite.mem_toFinset] using h0,
      (mem_hosFiber K L).2 rfl⟩
  rw [finsum_eq_finsetSum_of_support_subset _ hFsupp,
    finsum_eq_finsetSum_of_support_subset _ hGsupp, hT,
    Finset.sum_biUnion (hosFiber_pairwiseDisjoint K L S), Finset.mul_sum]
  refine Finset.sum_congr rfl fun V _ => ?_
  have hterm : ∀ W ∈ hosFiber K L V,
      ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
          * (((max (-(ordv (hosComap K L W) x)) 0 : ℤ) : ℝ)
              * Real.log (residueCard (hosComap K L W)))
        = ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
            * (((max (-(ordv V x)) 0 : ℤ) : ℝ) * Real.log (residueCard V)) := by
    intro W hW
    rw [mem_hosFiber] at hW
    rw [hW]
  rw [Finset.sum_congr rfl hterm, ← Finset.sum_mul]
  have hsum : ∑ W ∈ hosFiber K L V, ((ramIdxOver K L W * inertiaOver K L W : ℕ) : ℝ)
      = ((Module.finrank K L : ℕ) : ℝ) := by
    rw [← Nat.cast_sum]
    congr 1
    rw [← sum_ramification_inertia_hos K L V]
    refine Finset.sum_congr rfl fun W hW => ?_
    rw [mem_hosFiber] at hW
    rw [ramIdxOver, inertiaOver, hW]
  rw [hsum]

/-! ## ★無限素点側 -/

open scoped Classical in
/-- ★★**無限素点側も `[L:K]` 倍になる**。

★`Σ_{w|v} mult_w = mult_v·[L:K]`(`InfinitePlaceRel.lean`)そのものである。 -/
theorem sum_arc_extension (x : Kˣ) :
    ∑ w : InfinitePlace L, (w.mult : ℝ) * Real.posLog (w (((unitsExt K L x : Lˣ) : L)))
      = (Module.finrank K L : ℝ)
        * ∑ v : InfinitePlace K, (v.mult : ℝ) * Real.posLog (v ((x : K))) := by
  classical
  have hcomap : ∀ w : InfinitePlace L,
      w (((unitsExt K L x : Lˣ) : L)) = (w.comap (algebraMap K L)) ((x : K)) := by
    intro w; rfl
  rw [Finset.sum_congr rfl (fun w _ => by rw [hcomap w])]
  rw [← Finset.sum_fiberwise_of_maps_to
    (g := fun w : InfinitePlace L => w.comap (algebraMap K L))
    (t := (Finset.univ : Finset (InfinitePlace K)))
    (fun w _ => Finset.mem_univ _)]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun v _ => ?_
  have hinner : ∑ w ∈ Finset.univ.filter
      (fun w : InfinitePlace L => w.comap (algebraMap K L) = v),
        (w.mult : ℝ) * Real.posLog ((w.comap (algebraMap K L)) ((x : K)))
      = (∑ w ∈ Finset.univ.filter
          (fun w : InfinitePlace L => w.comap (algebraMap K L) = v), (w.mult : ℝ))
        * Real.posLog (v ((x : K))) := by
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun w hw => ?_
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hw
    rw [hw]
  rw [hinner]
  have hsum : (∑ w ∈ Finset.univ.filter
      (fun w : InfinitePlace L => w.comap (algebraMap K L) = v), (w.mult : ℝ))
      = ((InfinitePlace.mult v : ℕ) : ℝ) * (Module.finrank K L : ℝ) := by
    have h := sum_mult_comap_eq K L v
    have hcast : ((∑ w ∈ Finset.univ.filter
        (fun w : InfinitePlace L => w.comap (algebraMap K L) = v),
          InfinitePlace.mult w : ℕ) : ℝ)
        = ∑ w ∈ Finset.univ.filter
            (fun w : InfinitePlace L => w.comap (algebraMap K L) = v), (w.mult : ℝ) := by
      push_cast; ring
    rw [← hcast, h]
    push_cast
    ring
  rw [hsum]
  ring

/-! ## ★★★高さの拡大公式 -/

open scoped Classical in
/-- ★★★**`logHeight₁_L(x) = [L:K]·logHeight₁_K(x)`**(`x ∈ K ⊆ L`)。

★★mathlib は高さの体拡大の公式を持たない(2026-08-17 実測)。
★これが無いと「**絶対**高さ」——体に依らない量——が定義できない。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.
-/
theorem logHeight₁_extension (x : Kˣ) :
    Height.logHeight₁ (((unitsExt K L x : Lˣ) : L))
      = (Module.finrank K L : ℝ) * Height.logHeight₁ ((x : K)) := by
  classical
  rw [NumberField.logHeight₁_eq, NumberField.logHeight₁_eq, mul_add,
    sum_arc_extension K L x, finsum_posLog_extension K L x]

/-- ★**乗法版** —— `H_L(x) = H_K(x)^{[L:K]}`。 -/
theorem mulHeight₁_extension (x : Kˣ) :
    Height.mulHeight₁ (((unitsExt K L x : Lˣ) : L))
      = Height.mulHeight₁ ((x : K)) ^ (Module.finrank K L) := by
  have h := logHeight₁_extension K L x
  rw [Height.logHeight₁_eq_log_mulHeight₁, Height.logHeight₁_eq_log_mulHeight₁] at h
  have h1 : (0 : ℝ) < Height.mulHeight₁ (((unitsExt K L x : Lˣ) : L)) :=
    Height.mulHeight₁_pos _
  have h2 : (0 : ℝ) < Height.mulHeight₁ ((x : K)) := Height.mulHeight₁_pos _
  have := congrArg Real.exp h
  rwa [Real.exp_log h1, ← Real.log_pow, Real.exp_log (by positivity)] at this

/-! ## ★★★絶対高さ —— 体に依らない量 -/

end HeightExtension

/-- ★★**絶対対数高さ** `h(x) ≝ logHeight₁_K(x) / [K:ℚ]`。

★原文が `X(ℚ̄)` の上で高さを語れるのは、この量が**体の取り方に依らない**からである。
`degNormalized`(`ArithDiv.lean`)と**まったく同じ正規化**である。 -/
noncomputable def logAbsHeight (K : Type*) [Field K] [NumberField K] (x : K) : ℝ :=
  Height.logHeight₁ x / (Module.finrank ℚ K : ℝ)

section AbsHeight

variable (K L : Type*) [Field K] [NumberField K] [Field L] [NumberField L] [Algebra K L]

/-- ★★★**絶対高さは体拡大で不変である**。

> `h_L(x) = h_K(x)`(`x ∈ K ⊆ L`)

★★これが `Definition 1.2, (i)` の `ht : X(ℚ̄) → ℝ` が well-defined である理由の、
**高さ理論の側の中身**である。

★`Found/GenEll/BaseChange.lean` の `degNormalized_baseChange`(算術因子の次数)と
**同じ形**である——どちらも `[L:K]` が約分される。 -/
theorem logAbsHeight_extension (x : Kˣ) :
    logAbsHeight L (((unitsExt K L x : Lˣ) : L)) = logAbsHeight K ((x : K)) := by
  haveI : FiniteDimensional K L := Module.Finite.of_restrictScalars_finite ℚ K L
  have htower : Module.finrank ℚ K * Module.finrank K L = Module.finrank ℚ L :=
    Module.finrank_mul_finrank ℚ K L
  have hKL : (0 : ℝ) < (Module.finrank K L : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := K) (M := L)
  have hK : (0 : ℝ) < (Module.finrank ℚ K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := K)
  rw [logAbsHeight, logAbsHeight, logHeight₁_extension K L x, ← htower]
  push_cast
  field_simp

end AbsHeight

/-! ## ★出典の紐付け(`.src`) -/

def ordv_unitsExt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さの拡大公式——局所の関係)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
