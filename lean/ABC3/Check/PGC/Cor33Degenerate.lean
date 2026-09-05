import ABC3.Check.PGC.Def32Degenerate
import ABC3.Check.PGC.Theorem42Degenerate
import ABC3.Found.PGC.AbsGalUnitsSurjective

/-!
# [pGC] Corollary 3.3 の**旧**形は偽だった——`ρ` と `ρ'` が無関係だから

原文 (pGC p.6):

> Given a continuous E[Γ_K]-module V of E-dimension 1, the issue of whether or not V is
> uniformizing can be determined entirely group-theoretically from the filtered group Γ_K.

原文は**ひとつの** `V` を `α` で移して見比べる。ところが旧形は

```
∀ {K K'} (_α) (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ),
  IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ'
```

と、`ρ` と `ρ'` を**互いに無関係な**パラメータで受け取っていた。すると
`K' = K`・`_α = 恒等` を取って「どんな2つの表現も uniformizing 性が一致する」
と言っていることになる。

## 反例(`cor_3_3_statement_false`、`sorry` 無し)

一致しない2つの `ρ` を実際に作る:

* **真になる側**: `E := K`・`ι := id`・`I := U_K` として、
  `Γ_K ↠ 𝒪_K^×`(`Found/PGC/AbsGalUnitsSurjective.lean::
  exists_surjective_absGal_units`、Lubin-Tate 相互律を無条件化したもの)から
  `ρ := Γ_K ↠ 𝒪_K^× ↪ K^×` と `toGal := その全射の切断(選択)`を取れば
  `ρ (toGal x) = x` がすべての `x ∈ U_K` で成り立つ
  (`isUniformizing_unitsHomField`)。
* **偽になる側**: 自明な表現 `ρ' = 1`
  (`Check/PGC/Def32Degenerate.lean::not_isUniformizing_one`)。

★**本日構築した相互律の半分がそのまま反例の材料になった**——「真になる `ρ`」
を作るには `U_K` の開部分群の上で `ρ∘toGal` が体準同型と一致する必要があり、
それは相互律相当の構成を要する。それが手に入ったので初めて反証できた。

## 修理

`ρ'` を `α` で `ρ` と結ぶ条件 `∀ g, ρ' (α.equiv g) = ρ g` を課した
(原文どおり「ひとつの `V` を移して見る」)。`K' = K`・`α = 恒等` では
`ρ' = ρ` が強制されるので、上の反例は塞がる。

★これで「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例は
4 つ目(`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・本件)。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- `‖x‖ = 1` なら `x` は `𝒪_K` の単数。 -/
noncomputable def unitOfNormOne (K : PAdicLocalField p) (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) :
    (𝒪[K.carrier])ˣ := by
  have hmem : (x : K.carrier) ∈ 𝒪[K.carrier] := by
    rw [Valuation.mem_integer_iff]
    have hv : Valued.v (x : K.carrier) = (‖(x : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    have h1 : ‖(x : K.carrier)‖ ≤ 1 := le_of_eq x.2
    exact_mod_cast h1
  refine (?_ : IsUnit (⟨(x : K.carrier), hmem⟩ : 𝒪[K.carrier])).unit
  rw [Valued.integer.isUnit_iff_norm_eq_one]
  exact x.2

/-- `Γ_K ↠ 𝒪_K^×` の全射を1つ選ぶ(第 962 で無条件になったもの)。 -/
noncomputable def unitsHom (K : PAdicLocalField p) : K.absGal →* (𝒪[K.carrier])ˣ :=
  Classical.choose (exists_surjective_absGal_units K)

theorem unitsHom_surjective (K : PAdicLocalField p) : Function.Surjective (unitsHom K) :=
  Classical.choose_spec (exists_surjective_absGal_units K)

/-- `Γ_K →* K^×`(像は `U_K`)。 -/
noncomputable def unitsHomField (K : PAdicLocalField p) : K.absGal →* (K.carrier)ˣ :=
  (Units.map (Subring.subtype 𝒪[K.carrier] : 𝒪[K.carrier] →* K.carrier)).comp (unitsHom K)

/-- 全射の切断(選択公理)——`toGal` の役をこれで務める。 -/
noncomputable def toGalChoice (K : PAdicLocalField p)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) : K.absGal :=
  Classical.choose (unitsHom_surjective K (unitOfNormOne K x))

theorem unitsHomField_toGalChoice (K : PAdicLocalField p)
    (x : {y : K.carrier // ‖y‖ = (1 : ℝ)}) :
    ((unitsHomField K (toGalChoice K x) : (K.carrier)ˣ) : K.carrier) = (x : K.carrier) := by
  have h := Classical.choose_spec (unitsHom_surjective K (unitOfNormOne K x))
  show (((unitsHom K (toGalChoice K x) : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) = _
  rw [show unitsHom K (toGalChoice K x) = unitOfNormOne K x from h]
  rfl

/-- `U_K = {‖x‖ = 1}` は開(超距離)。 -/
theorem isOpen_normOne (K : PAdicLocalField p) : IsOpen {x : K.carrier | ‖x‖ = 1} := by
  rw [Metric.isOpen_iff]
  intro x hx
  refine ⟨1, one_pos, ?_⟩
  intro y hy
  simp only [Metric.mem_ball, dist_eq_norm] at hy
  simp only [Set.mem_setOf_eq] at hx ⊢
  have hne : ‖y - x‖ ≠ ‖x‖ := by rw [hx]; exact ne_of_lt hy
  have hyx : y = (y - x) + x := by ring
  rw [hyx, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier) hne, hx]
  exact max_eq_right (le_of_lt (by rw [← hx] at hy ⊢; exact hy))

/-- **★★★★★真になる `ρ` は存在する**——`Γ_K ↠ 𝒪_K^×` を使う。 -/
theorem isUniformizing_unitsHomField (K : PAdicLocalField p) :
    IsUniformizing K K.carrier (toGalChoice K) (unitsHomField K) := by
  refine ⟨{x : K.carrier | ‖x‖ = 1}, subset_rfl, isOpen_normOne K, by simp, ?_, ?_,
    RingHom.id K.carrier, ?_⟩
  · intro a ha b hb
    simp only [Set.mem_setOf_eq] at *
    rw [norm_mul, ha, hb, one_mul]
  · intro a ha
    simp only [Set.mem_setOf_eq] at *
    rw [norm_inv, ha, inv_one]
  · intro x hx
    exact unitsHomField_toGalChoice K ⟨x, hx⟩

/-- **★★★★★★[pGC] Corollary 3.3 の旧形は偽**——`ρ` と `ρ'` が無関係なので、
`K' = K`・`α = 恒等` で「どんな2つの表現も uniformizing 性が一致する」に
なってしまう。真になる `ρ`(上)と偽になる `ρ' = 1` を並べれば落ちる。 -/
theorem cor_3_3_statement_false (p : ℕ) [Fact p.Prime] (K₀ : PAdicLocalField p) :
    ¬ (∀ (RF : RamificationFiltration p) (E : Type) (_ : Field E) (_ : Algebra ℚ_[p] E)
        (toGal : ∀ K : PAdicLocalField p, {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal)
        {K K' : PAdicLocalField p}
        (_α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
        (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ),
        IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ') := by
  intro h
  have key := (h (degenerateRF p) K₀.carrier K₀.isField K₀.isAlgebra toGalChoice
    (idFilteredIso (filteredGroupOf (degenerateRF p) K₀)) (unitsHomField K₀) 1).mp
    (isUniformizing_unitsHomField K₀)
  exact not_isUniformizing_one K₀ K₀.carrier (toGalChoice K₀) key

end ABC3.Check.PGC
