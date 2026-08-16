import ABC3.Found.FrdI.Def31

/-!
# [FrdI] Proposition 3.3 —— Base-identity Pre-steps and Units

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.59–p.61。

原文 (FrdI p.59):
> (Base-identity Pre-steps and Units)

## ★この命題の規模(測定)

**5 条、主張は 8**:

| 条 | # | 内容 | 状態 |
|---|---|---|---|
| (i) | 1 | `𝒟` が Frobenius-slim なら `𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` の `1` の像は base-identity pre-step | 未 |
| (i) | 2 | 逆に Frobenius-normalized ＋ `A` が Frobenius-trivial なら、どの base-identity pre-step もそう現れる | 未 |
| (ii) | 3 | `α₁ ≈ α₂`(unit-equivalent)⟺ `degFr`・`Div`・`Base` が一致 —— **必要性** | ★**ここで実装** |
| (ii) | 4 | 同 —— **十分性** | ★**ここで実装** |
| (iii) | 5 | `𝒞^istr → 𝒞^un-tr` は full かつ本質的全射、圏同値 ⟺ `𝒞^istr` が unit-trivial 型 | 未 |
| (iv) | 6 | `𝒞^un-tr → 𝔽_Φ` は忠実かつ本質的全射で、`𝒞^un-tr` に Frobenioid の構造を与える | 未 |
| (iv) | 7 | `𝒞^un-tr` の射の類型は `𝒞^istr` の射の類型から来る | 未 |
| (v) | 8 | `𝒞 → 𝔽_Φ` が圏同値 ⟺ `𝒞` が Aut-ample・unit-trivial・base-trivial 型 | 未 |

★**十分性の証明は原文どおり 4 手**(`Definition 1.3, (iv), (a)` の 3 分解、
`Definition 1.3, (ii)` の Frobenius 型射の本質的一意性、
`Definition 1.3, (i), (c)` の pull-back の圏同値、`Definition 1.3, (vi)` の
faithfulness up to units)であり、`Remark 2.7.2` の分解と同じ形をしている。

★**(i)(iii)(iv)(v) が入るまで `.src` は付けない。**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(ii) の必要性

原文 (FrdI p.59):
> three conditions are satisfied: (a) degFr(

★**理由は「`𝒪^×` の元は LB-invertible な base-identity linear 自己射」**
(原文どおり)——`Div = 0`、`degFr = 1`、`Base = 𝟙` なので、
3 つの量のどれにも寄与しない。
-/

/-- ★★★**[FrdI] Proposition 3.3, (ii) の必要性** —— unit-equivalent な 2 射は
`Frobenius 次数`・`Div`・`Base` がすべて一致する。

★すなわち **`𝔽_Φ` へ同じ射に写る**。 -/
theorem prop_3_3_ii_necessity {A B : C} {α₁ α₂ : A ⟶ B}
    (h : IsUnitEquivalent P α₁ α₂) :
    P.degFr α₁ = P.degFr α₂ ∧ P.Div α₁ = P.Div α₂ ∧ P.Base α₁ = P.Base α₂ := by
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ := h
  haveI : IsIso ((δ : Cc ⟶ Cc)) := isIso_of_mem_otimes P hδ
  -- `δ` の 3 つの量
  have hdb : P.Base ((δ : Cc ⟶ Cc)) = 𝟙 _ := by
    have h : P.Base ((δ : Cc ⟶ Cc)) = P.Base (𝟙 Cc) := hδ.1.1
    rwa [P.Base_id] at h
  have hdd : P.degFr ((δ : Cc ⟶ Cc)) = 1 := hδ.1.2
  have hdv : P.Div ((δ : Cc ⟶ Cc)) = 0 := isIsometric_of_isIso P _
  -- `δ ≫ β` は `β` と 3 つの量が一致する
  have hdegβ : P.degFr (((δ : Cc ⟶ Cc)) ≫ β) = P.degFr β := by
    rw [P.degFr_comp, hdd, mul_one]
  have hbaseβ : P.Base (((δ : Cc ⟶ Cc)) ≫ β) = P.Base β := by
    rw [P.Base_comp, hdb, Category.id_comp]
  have hdivβ : P.Div (((δ : Cc ⟶ Cc)) ≫ β) = P.Div β := by
    rw [P.Div_comp, hdb, MonoidOn.map_id, hdv]
    simp
  refine ⟨?_, ?_, ?_⟩
  · rw [h₁, h₂, P.degFr_comp, P.degFr_comp, hdegβ]
  · rw [h₁, h₂, P.Div_comp, P.Div_comp, hdegβ, hdivβ]
  · rw [h₁, h₂, P.Base_comp, P.Base_comp, hbaseβ]

/-! ## ★(ii) の十分性

原文 (FrdI p.61):
> uniqueness of morphisms of Frobenius type of a given Frobenius degree [cf. Defi-

★**原文どおり 4 手**:
1. `Definition 1.3, (iv), (a)`(`arbFactor`)で `α₁ = γ₁ ≫ β₁ ≫ δ₁`、
   `α₂ = γ₂ ≫ β₂ ≫ δ₂`(Frobenius 型 ≫ pre-step ≫ pull-back)と分ける
2. **条件 (a)** と `Definition 1.3, (ii)`(`frobDegUniq`)で `γ` を**共通に揃える**
3. **条件 (c)** と `Definition 1.3, (i), (c)`(pull-back の普遍性)で `δ` を**共通に揃える**
   ——★**ここで pull-back の全単射性を 2 度使い、`ε` が同型であることを出す**
4. **条件 (b)** で `Div β₁ = Div β₂`(★`γ`・`δ` が LB-invertible ＝ `Div = 0` なので
   `Div α = Φ.map (Base γ) (Div β)` に潰れ、`Φ.map` の単射性で消える)、
   `Definition 1.3, (vi)`(`faithfulUpToUnits`)で `β₂ = β₁ ≫ ζ`
-/

/-- ★★★**[FrdI] Proposition 3.3, (ii) の十分性** —— `degFr`・`Div`・`Base` が
一致する co-objective な 2 射は **unit-equivalent** である。

★`𝒞` が isotropic 型であることを使う(原文が「`𝒞` を `𝒞^istr` に置き換えてよい」
と述べているとおり)。 -/
theorem prop_3_3_ii_sufficiency (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} {α₁ α₂ : A ⟶ B}
    (hdeg : P.degFr α₁ = P.degFr α₂) (hdiv : P.Div α₁ = P.Div α₂)
    (hbase : P.Base α₁ = P.Base α₂) : IsUnitEquivalent P α₁ α₂ := by
  -- 手 1: 3 分解
  obtain ⟨X₁, Y₁, γ₁, β₁, δ₁, hfac₁, hγ₁, hβ₁, hδ₁⟩ := Fc.arbFactor α₁
  obtain ⟨X₂, Y₂, γ₂, β₂, δ₂, hfac₂, hγ₂, hβ₂, hδ₂⟩ := Fc.arbFactor α₂
  -- 補助: `γ` は次数を担い、`β`・`δ` は linear
  have hdβ₁ : P.degFr β₁ = 1 := hβ₁.1
  have hdβ₂ : P.degFr β₂ = 1 := hβ₂.1
  have hdδ₁ : P.degFr δ₁ = 1 := (Fc.pullBackLB δ₁ hδ₁).2
  have hdδ₂ : P.degFr δ₂ = 1 := (Fc.pullBackLB δ₂ hδ₂).2
  have hdegγ : P.degFr γ₁ = P.degFr γ₂ := by
    have e₁ : P.degFr α₁ = P.degFr γ₁ := by
      rw [hfac₁, P.degFr_comp, P.degFr_comp, hdβ₁, hdδ₁]; simp
    have e₂ : P.degFr α₂ = P.degFr γ₂ := by
      rw [hfac₂, P.degFr_comp, P.degFr_comp, hdβ₂, hdδ₂]; simp
    rw [← e₁, ← e₂, hdeg]
  -- 手 2: `γ` を揃える
  obtain ⟨e, hei, hee⟩ := Fc.frobDegUniq A X₁ X₂ γ₁ γ₂ hγ₁ hγ₂ hdegγ
  haveI := hei
  have hβ₂' : IsPreStep P (e ≫ β₂) := IsPreStep.comp P (isPreStep_of_isIso P e) hβ₂
  have hfac₂' : α₂ = γ₁ ≫ (e ≫ β₂) ≫ δ₂ := by
    rw [hfac₂, ← hee]; simp
  -- 底の同型たち
  haveI hbγ : IsIso (P.Base γ₁) := hγ₁.2
  haveI hbβ₁ : IsIso (P.Base β₁) := hβ₁.2
  haveI hbβ₂' : IsIso (P.Base (e ≫ β₂)) := hβ₂'.2
  -- 手 3: `δ` を揃える
  set u : (P.toElem.obj Y₁).base ⟶ (P.toElem.obj Y₂).base :=
    inv (P.Base β₁) ≫ P.Base (e ≫ β₂) with hu
  haveI : IsIso u := by rw [hu]; infer_instance
  have hbaseδ : P.Base δ₁ = u ≫ P.Base δ₂ := by
    have h := hbase
    rw [hfac₁, hfac₂', P.Base_comp, P.Base_comp, P.Base_comp, P.Base_comp] at h
    have h2 : P.Base β₁ ≫ P.Base δ₁ = P.Base (e ≫ β₂) ≫ P.Base δ₂ :=
      (cancel_epi (P.Base γ₁)).mp h
    rw [hu, Category.assoc, ← h2, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨-, hsurj₂⟩ := hδ₂ Y₁
  obtain ⟨ε, hε⟩ := hsurj₂ ⟨(δ₁, u), hbaseδ⟩
  have hε' := congrArg Subtype.val hε
  have hεδ : ε ≫ δ₂ = δ₁ := congrArg Prod.fst hε'
  have hεb : P.Base ε = u := congrArg Prod.snd hε'
  obtain ⟨-, hsurj₁⟩ := hδ₁ Y₂
  obtain ⟨hinj₁, -⟩ := hδ₁ Y₁
  have hbaseδ' : P.Base δ₂ = inv u ≫ P.Base δ₁ := by
    rw [hbaseδ, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨ζ₀, hζ₀⟩ := hsurj₁ ⟨(δ₂, inv u), hbaseδ'⟩
  have hζ₀' := congrArg Subtype.val hζ₀
  have hζδ : ζ₀ ≫ δ₁ = δ₂ := congrArg Prod.fst hζ₀'
  have hζb : P.Base ζ₀ = inv u := congrArg Prod.snd hζ₀'
  haveI : IsIso ε := by
    refine ⟨ζ₀, ?_, ?_⟩
    · refine hinj₁ (Subtype.ext (Prod.ext ?_ ?_))
      · show (ε ≫ ζ₀) ≫ δ₁ = 𝟙 Y₁ ≫ δ₁
        rw [Category.assoc, hζδ, hεδ, Category.id_comp]
      · show P.Base (ε ≫ ζ₀) = P.Base (𝟙 Y₁)
        rw [P.Base_comp, hεb, hζb, IsIso.hom_inv_id, P.Base_id]
    · obtain ⟨hinj₂, -⟩ := hδ₂ Y₂
      refine hinj₂ (Subtype.ext (Prod.ext ?_ ?_))
      · show (ζ₀ ≫ ε) ≫ δ₂ = 𝟙 Y₂ ≫ δ₂
        rw [Category.assoc, hεδ, hζδ, Category.id_comp]
      · show P.Base (ζ₀ ≫ ε) = P.Base (𝟙 Y₂)
        rw [P.Base_comp, hζb, hεb, IsIso.inv_hom_id, P.Base_id]
  -- `β₁'' := β₁ ≫ ε`
  have hβ₁'' : IsPreStep P (β₁ ≫ ε) := IsPreStep.comp P hβ₁ (isPreStep_of_isIso P ε)
  have hfac₁' : α₁ = γ₁ ≫ (β₁ ≫ ε) ≫ δ₂ := by
    rw [hfac₁, Category.assoc, hεδ]
  -- 底が一致する
  have hbb : P.Base (β₁ ≫ ε) = P.Base (e ≫ β₂) := by
    rw [P.Base_comp, hεb, hu, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  -- 手 4: `Div` が一致する
  have hDivγ : P.Div γ₁ = 0 := hγ₁.1.2
  have hDivδ₂ : P.Div δ₂ = 0 := (Fc.pullBackLB δ₂ hδ₂).1.2
  have hDivFac : ∀ β : X₁ ⟶ Y₂,
      P.Div (γ₁ ≫ β ≫ δ₂) = Φ.map (P.Base γ₁) (P.Div β) := by
    intro β
    rw [P.Div_comp, P.Div_comp, hDivδ₂, hdδ₂, hDivγ]
    simp
  have hdivβ : P.Div (β₁ ≫ ε) = P.Div (e ≫ β₂) := by
    refine Φ.map_injective (P.Base γ₁) ?_
    rw [← hDivFac _, ← hDivFac _, ← hfac₁', ← hfac₂', hdiv]
  -- `faithfulUpToUnits`
  obtain ⟨ζ, hζm, hζe⟩ := Fc.faithfulUpToUnits (e ≫ β₂) (β₁ ≫ ε)
    hbb.symm hdivβ.symm
    (isCoAngular_of_isotropic_dom P Fc (hiso X₁) _) hβ₂'
    (isCoAngular_of_isotropic_dom P Fc (hiso X₁) _) hβ₁''
  -- 組み立て
  refine ⟨Y₂, γ₁ ≫ (β₁ ≫ ε), δ₂, ζ, hζm, ?_, ?_⟩
  · rw [hfac₁']
    simp
  · rw [hfac₂', hζe]
    simp

/-- ★★★**[FrdI] Proposition 3.3, (ii)** —— 両向きを合わせた形。 -/
theorem prop_3_3_ii (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (α₁ α₂ : A ⟶ B) :
    IsUnitEquivalent P α₁ α₂ ↔
      (P.degFr α₁ = P.degFr α₂ ∧ P.Div α₁ = P.Div α₂ ∧ P.Base α₁ = P.Base α₂) :=
  ⟨prop_3_3_ii_necessity P,
   fun ⟨hd, hv, hb⟩ => prop_3_3_ii_sufficiency P Fc hiso hd hv hb⟩

/-- ★**`𝔽_Φ` へ同じ射に写る**という原文の言い方そのもの。 -/
theorem prop_3_3_ii_toElem {A B : C} {α₁ α₂ : A ⟶ B} (h : IsUnitEquivalent P α₁ α₂) :
    P.toElem.map α₁ = P.toElem.map α₂ := by
  obtain ⟨hd, hv, hb⟩ := prop_3_3_ii_necessity P h
  exact ElemFrobCat.Hom.ext hb hv hd

end ABC3.Found.FrdI
