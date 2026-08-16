import ABC3.Found.GenEll.ArithDiv

/-!
# [GenEll] §1 地の文 —— `APrc(F)` が**主算術因子の像そのもの**であること(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> An arithmetic divisor on F is defined to be a finite formal sum

## ★★自分で残した宿題を埋める

`Found/GenEll/ArithDiv.lean` の `APrc` の docstring はこう書いていた:

> ★`principalADiv` が準同型であることは**まだ示していない**。示せば
> `AddMonoidHom.range` と一致することが言える。

★本ファイルがそれである。`principalADiv : Fˣ → ADiv(F)` が
**乗法から加法への準同型**であることを示し、
`APrc(F)`(生成する部分群)が**像そのもの**であることを結論する。

★★**なぜこれが要るか。** 原文は「the principal arithmetic divisors determine a subgroup
`APrc(F) ⊆ ADiv(F)`」と書く。★**「生成する部分群」と「像」が一致することは自明ではない**
——一致するのは `principalADiv` が準同型だからである。
`ArithDiv.lean` では**安全側に倒して「生成する部分群」と定義した**ので、
ここで一致を示して初めて原文の読みに追いつく。

## ★これが Arakelov 側で持つ意味

`APic(Spec(O_F)) = ADiv(F)/APrc(F)` が**算術 Picard 群**である。
★`APrc` が像だと分かったので、商は `ADiv(F)/im(principalADiv)` として扱える
——これは `Definition 1.1` の `APic(X)` の **`X = Spec(O_F)` の場合**にあたる。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {F : Type*} [Field F] [NumberField F]

/-! ## ★`ord_v` は加法的 -/

/-- ★**`ord_v(fg) = ord_v(f) + ord_v(g)`**。

付値が乗法的であることを `WithZero` の外へ運ぶ。 -/
theorem ordv_mul (v : FinitePlace F) (f g : Fˣ) :
    ordv v (f * g) = ordv v f + ordv v g := by
  have hf : v.valuation F ((f : F)) ≠ 0 :=
    (Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero f)
  have hg : v.valuation F ((g : F)) ≠ 0 :=
    (Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero g)
  have hfg : v.valuation F (((f * g : Fˣ) : F)) ≠ 0 :=
    (Valuation.ne_zero_iff (v.valuation F)).2 (Units.ne_zero (f * g))
  have hval : v.valuation F (((f * g : Fˣ) : F))
      = v.valuation F ((f : F)) * v.valuation F ((g : F)) := by
    push_cast
    exact map_mul _ _ _
  have key : WithZero.unzero hfg = WithZero.unzero hf * WithZero.unzero hg := by
    apply WithZero.coe_injective
    rw [WithZero.coe_unzero, WithZero.coe_mul, WithZero.coe_unzero, WithZero.coe_unzero]
    exact hval
  simp only [ordv]
  rw [key, toAdd_mul]
  ring

/-! ## ★`principalADiv` は準同型 -/

/-- ★**`ADiv(fg) = ADiv(f) + ADiv(g)`** —— 原文の `APrc(F)` が部分群である根拠。

有限側は `ord_v` の加法性、アルキメデス側は `log` の乗法→加法。 -/
theorem principalADiv_mul (f g : Fˣ) :
    principalADiv (f * g) = principalADiv f + principalADiv g := by
  refine Prod.ext ?_ ?_
  · ext v
    simp only [principalADiv, Finsupp.onFinset_apply, Prod.fst_add, Finsupp.add_apply]
    exact ordv_mul v f g
  · ext w
    have hf : w ((f : F)) ≠ 0 := by
      simpa using (map_ne_zero_iff _ w.injective).2 (Units.ne_zero f)
    have hg : w ((g : F)) ≠ 0 := by
      simpa using (map_ne_zero_iff _ w.injective).2 (Units.ne_zero g)
    have hw : w (((f * g : Fˣ) : F)) = w ((f : F)) * w ((g : F)) := by
      push_cast
      exact map_mul _ _ _
    rw [show ((principalADiv f + principalADiv g : ADiv F)).2
        = (principalADiv f).2 + (principalADiv g).2 from rfl]
    simp only [principalADiv, Finsupp.onFinset_apply, Finsupp.add_apply]
    rw [hw, Real.log_mul hf hg]
    ring

@[simp] theorem principalADiv_one : principalADiv (1 : Fˣ) = (0 : ADiv F) := by
  have h := principalADiv_mul (1 : Fˣ) (1 : Fˣ)
  rw [one_mul] at h
  have h2 : (0 : ADiv F) + principalADiv (1 : Fˣ)
      = principalADiv (1 : Fˣ) + principalADiv (1 : Fˣ) := by rw [zero_add]; exact h
  exact (add_right_cancel h2).symm

theorem principalADiv_inv (f : Fˣ) : principalADiv f⁻¹ = -principalADiv f := by
  have h := principalADiv_mul f f⁻¹
  rw [mul_inv_cancel, principalADiv_one] at h
  exact eq_neg_of_add_eq_zero_right h.symm

/-! ## ★`APrc(F)` は像そのもの -/

/-- 主算術因子の**像**を部分群として直に取ったもの。 -/
noncomputable def principalRange (F : Type*) [Field F] [NumberField F] :
    AddSubgroup (ADiv F) where
  carrier := Set.range (principalADiv : Fˣ → ADiv F)
  zero_mem' := ⟨1, principalADiv_one⟩
  add_mem' := by
    rintro _ _ ⟨f, rfl⟩ ⟨g, rfl⟩
    exact ⟨f * g, principalADiv_mul f g⟩
  neg_mem' := by
    rintro _ ⟨f, rfl⟩
    exact ⟨f⁻¹, principalADiv_inv f⟩

/-- ★★**原文の読みに追いついた** —— `APrc(F)`(生成する部分群)は**像そのもの**である。

★`ArithDiv.lean` は安全側に倒して「生成する部分群」と定義していた。
一致するのは `principalADiv` が準同型だからであり、それが `principalADiv_mul` である。 -/
theorem APrc_eq_principalRange : APrc F = principalRange F := by
  rw [APrc]
  exact AddSubgroup.closure_eq (principalRange F)

/-- ★言い換え: `APrc(F)` の元は**主算術因子そのもの**である(和を取る必要がない)。 -/
theorem mem_APrc_iff (a : ADiv F) :
    a ∈ APrc F ↔ ∃ f : Fˣ, principalADiv f = a := by
  rw [APrc_eq_principalRange]
  exact Iff.rfl

/-! ## ★出典の紐付け(`.src`) -/

def principalADiv_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

def APrc_eq_principalRange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(算術因子 ADiv(F))",
    sectionId := "genell-adiv" }

end ABC3.Found.GenEll
