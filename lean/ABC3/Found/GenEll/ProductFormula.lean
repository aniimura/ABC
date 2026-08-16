import ABC3.Found.GenEll.ArithDivHom
import Mathlib.NumberTheory.NumberField.ProductFormula

/-!
# [GenEll] §1 地の文 —— 主算術因子の次数は `0`(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

## ★★これが `APic(Spec(O_F)) → ℝ` の中身である

原文は「`deg_F` は準同型 `APic(Spec(O_F)) → ℝ` を定める」と書く。
★**そう言えるためには `deg(APrc(F)) = 0` が要る**——すなわち

> **`deg(ADiv(f)) = 0`(すべての `f ∈ Fˣ`)**

でなければ商の上に降りない。★**それが積公式である。**

## ★積公式との対応

mathlib: `NumberField.prod_abs_eq_one`

```
(∏ w : InfinitePlace K, w x ^ w.mult) * ∏ᶠ w : FinitePlace K, w x = 1
```

対数を取ると `Σ_w mult_w·log(w x) + Σᶠ_v log(v x) = 0`。
ここで有限素点では **`v x = q_v^{−ord_v(x)}`** なので `log(v x) = −ord_v(x)·log q_v`。
したがって

```
Σ_v ord_v(x)·log q_v = Σ_w mult_w·log(w x)
```

であり、`deg(ADiv(f)) = Σ_v ord_v log q_v − Σ_w mult_w log|f|_w = 0`。★これが証明の全部である。

## ★★符号がここで効く(2026-08-17 の訂正)

★`ordv` の符号が逆だと `deg = −2Σ_w mult_w log|f|_w` になって **0 にならない**。
`ArithDiv.lean` の `ordv` の符号バグは、**まさにこの定理を取りにいって発覚した**。
`ordv_algebraMap_eq_count` がその規約を機械的に固定している。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {F : Type*} [Field F] [NumberField F]

/-! ## ★段 1 —— 有限素点の値を `ord_v` で書く -/

/-- ★**`w(f) = q_v^{−ord_v(f)}`**(`v = w.maximalIdeal`)。

mathlib の `adicAbv` は `toNNReal` 経由で `q_v ^ (unzero).toAdd` を返し、
我々の `ordv` はその **`toAdd` の符号を反転したもの**である。 -/
theorem finitePlace_apply_units (w : NumberField.FinitePlace F) (f : Fˣ) :
    w ((f : F))
      = (residueCard (NumberField.FinitePlace.maximalIdeal w) : ℝ)
          ^ (-(ordv (NumberField.FinitePlace.maximalIdeal w) f)) := by
  set v := NumberField.FinitePlace.maximalIdeal w with hv
  have hne : (v.valuation F) ((f : F)) ≠ 0 :=
    (Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)
  rw [← NumberField.FinitePlace.norm_embedding_eq w ((f : F)),
    NumberField.FinitePlace.norm_embedding',
    WithZeroMulInt.toNNReal_neg_apply _ hne]
  simp only [ordv, neg_neg, residueCard]
  push_cast
  ring

/-! ## ★段 2 —— 対数を取る -/

theorem log_finitePlace_apply (w : NumberField.FinitePlace F) (f : Fˣ) :
    Real.log (w ((f : F)))
      = -((ordv (NumberField.FinitePlace.maximalIdeal w) f : ℤ) : ℝ)
          * Real.log (residueCard (NumberField.FinitePlace.maximalIdeal w)) := by
  rw [finitePlace_apply_units w f, Real.log_zpow]
  push_cast
  ring

/-! ## ★段 3 —— 有限素点の和を `Finsupp.sum` に繋ぐ -/

/-- ★**有限側の次数は、有限素点の対数和の符号違いである**。 -/
theorem finsum_log_finitePlace (f : Fˣ) :
    ∑ᶠ w : NumberField.FinitePlace F, Real.log (w ((f : F)))
      = -((principalADiv f).fin.sum
            (fun v n => (n : ℝ) * Real.log (residueCard v))) := by
  classical
  have h1 : ∀ w : NumberField.FinitePlace F, Real.log (w ((f : F)))
      = -(((ordv (NumberField.FinitePlace.maximalIdeal w) f : ℤ) : ℝ)
            * Real.log (residueCard (NumberField.FinitePlace.maximalIdeal w))) := by
    intro w
    rw [log_finitePlace_apply w f]
    ring
  rw [finsum_congr h1]
  rw [finsum_neg_distrib]
  congr 1
  -- 添字を `HeightOneSpectrum` に付け替える
  have h2 : ∑ᶠ w : NumberField.FinitePlace F,
      (((ordv (NumberField.FinitePlace.maximalIdeal w) f : ℤ) : ℝ)
        * Real.log (residueCard (NumberField.FinitePlace.maximalIdeal w)))
      = ∑ᶠ v : FinitePlace F, ((ordv v f : ℤ) : ℝ) * Real.log (residueCard v) :=
    finsum_comp_equiv (M := ℝ)
      (e := NumberField.FinitePlace.equivHeightOneSpectrum (K := F))
      (f := fun v : FinitePlace F => ((ordv v f : ℤ) : ℝ) * Real.log (residueCard v))
  rw [h2]
  -- `finsum` を `Finsupp.sum` に落とす
  simp only [principalADiv, ADiv.fin]
  rw [Finsupp.onFinset_sum _ (fun a => by simp)]
  refine finsum_eq_finsetSum_of_support_subset _ ?_
  intro v hv
  rw [Function.mem_support] at hv
  have h0 : ordv v f ≠ 0 := by
    intro hcon
    apply hv
    rw [hcon]
    push_cast
    ring
  simpa using h0

/-! ## ★段 4 —— 積公式から `deg = 0` -/

/-- ★★**[GenEll] §1 地の文 —— `deg(ADiv(f)) = 0`**。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

★**これが「`deg_F` が `APic(Spec(O_F)) → ℝ` を定める」ことの中身である**
——主算術因子の次数が `0` でなければ、商の上に降りない。 -/
theorem deg_principalADiv_eq_zero (f : Fˣ) : deg (principalADiv f) = (0 : ℝ) := by
  classical
  have hx : ((f : F)) ≠ 0 := Units.ne_zero f
  -- 積公式
  have hpf := NumberField.prod_abs_eq_one (K := F) hx
  have harchpos : ∀ w : InfinitePlace F, (0 : ℝ) < w ((f : F)) ^ w.mult := by
    intro w
    exact pow_pos ((InfinitePlace.pos_iff).2 hx) _
  have hfinpos : ∀ w : NumberField.FinitePlace F, (0 : ℝ) < w ((f : F)) :=
    fun w => (NumberField.FinitePlace.pos_iff).2 hx
  have harch : (0 : ℝ) < ∏ w : InfinitePlace F, w ((f : F)) ^ w.mult :=
    Finset.prod_pos fun w _ => harchpos w
  have hfin : (0 : ℝ) < ∏ᶠ w : NumberField.FinitePlace F, w ((f : F)) := by
    rcases lt_trichotomy (∏ᶠ w : NumberField.FinitePlace F, w ((f : F))) 0 with h | h | h
    · have := mul_neg_of_pos_of_neg harch h
      rw [hpf] at this
      norm_num at this
    · rw [h, mul_zero] at hpf
      norm_num at hpf
    · exact h
  -- 対数を取る
  have hlog : Real.log (∏ w : InfinitePlace F, w ((f : F)) ^ w.mult)
      + Real.log (∏ᶠ w : NumberField.FinitePlace F, w ((f : F))) = 0 := by
    rw [← Real.log_mul harch.ne' hfin.ne', hpf, Real.log_one]
  -- アルキメデス側
  have harchlog : Real.log (∏ w : InfinitePlace F, w ((f : F)) ^ w.mult)
      = ∑ w : InfinitePlace F, (w.mult : ℝ) * Real.log (w ((f : F))) := by
    rw [Real.log_prod (fun w _ => (harchpos w).ne')]
    exact Finset.sum_congr rfl fun w _ => Real.log_pow _ _
  -- 有限側
  have hfinlog : Real.log (∏ᶠ w : NumberField.FinitePlace F, w ((f : F)))
      = ∑ᶠ w : NumberField.FinitePlace F, Real.log (w ((f : F))) :=
    Real.log_finprod hfinpos
  rw [harchlog, hfinlog, finsum_log_finitePlace f] at hlog
  -- `deg` を展開して突き合わせる
  have harcsum : (principalADiv f).arc.sum (fun _ r => r)
      = -∑ w : InfinitePlace F, (w.mult : ℝ) * Real.log (w ((f : F))) := by
    simp only [principalADiv, ADiv.arc]
    rw [Finsupp.onFinset_sum _ (fun a => rfl), ← Finset.sum_neg_distrib]
  rw [deg, harcsum]
  linarith [hlog]

/-! ## ★★段 5 —— `deg` は `APic(Spec(O_F))` に降りる

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

★**ここが原文の一文の中身である。** -/

/-- **算術 Picard 群** `APic(Spec(O_F)) ≝ ADiv(F)/APrc(F)`。

★これは `Definition 1.1` の `APic(X)` の **`X = Spec(O_F)` の場合**である
(一般の `X` には `X^arc` とエルミート計量が要り、mathlib に無い)。 -/
abbrev APicOF (F : Type*) [Field F] [NumberField F] := ADiv F ⧸ APrc F

/-- ★`deg` は主算術因子の上で消える(積公式)。 -/
theorem deg_eq_zero_of_mem_APrc (a : ADiv F) (ha : a ∈ APrc F) : deg a = 0 := by
  obtain ⟨f, rfl⟩ := (mem_APrc_iff a).1 ha
  exact deg_principalADiv_eq_zero f

/-- ★★**原文「determines a homomorphism `APic(Spec(O_F)) → ℝ`」**。

`deg` が `APrc(F)` の上で消えるので、商の上に降りる。
★**降りることが言えて初めて `deg` は「算術直線束の次数」になる。** -/
noncomputable def degAPic : APicOF F →+ ℝ :=
  QuotientAddGroup.lift (APrc F) degHom deg_eq_zero_of_mem_APrc

@[simp] theorem degAPic_mk (a : ADiv F) :
    degAPic (QuotientAddGroup.mk' (APrc F) a) = deg a := rfl

/-- ★正規化した版も同様に降りる(`deg` を `[F:ℚ]` で割るだけ)。 -/
noncomputable def degNormalizedAPic : APicOF F →+ ℝ :=
  (AddMonoidHom.mulRight ((Module.finrank ℚ F : ℝ))⁻¹).comp (degAPic (F := F))

@[simp] theorem degNormalizedAPic_mk (a : ADiv F) :
    degNormalizedAPic (QuotientAddGroup.mk' (APrc F) a) = degNormalized a := by
  simp only [degNormalizedAPic, AddMonoidHom.coe_comp, Function.comp_apply, degAPic_mk,
    AddMonoidHom.mulRight_apply, degNormalized]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def deg_principalADiv_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

def degAPic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

end ABC3.Found.GenEll
