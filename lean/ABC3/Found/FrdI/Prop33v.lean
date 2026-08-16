import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.Prop33

/-!
# [FrdI] Proposition 3.3, (v) —— `𝒞 → 𝔽_Φ` が圏同値になる条件

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.60–p.61。

原文 (FrdI p.60):
> is an equivalence of categories if and only if the Frobenioid C is of Aut-ample,

★**`𝒞 → 𝔽_Φ` が圏同値 ⟺ `𝒞` が Aut-ample・unit-trivial・base-trivial 型**。

## ★原文の証明

**必要性**は `Proposition 1.5, (i), (ii)`(`𝔽_Φ` がその 3 型である)から。

**十分性**は 4 手:
1. base-trivial なら **isotropic 型**(isotropic hull は base-isomorphism だから)
2. **忠実性** —— `Proposition 3.3, (ii)` で unit-equivalent に落ち、
   `𝒪^× = ⊥` なので 2 射は等しい
3. **本質的全射性** —— `Definition 1.3, (i), (a)`(`baseSurj`)
4. ★**充満性** —— 原文は「`Definition 1.3, (iv), (a)` の分解、
   `Definition 1.3, (ii)` の次数、`Definition 1.3, (i), (c)` の pull-back の圏同値により、
   `𝒪^▷(−)` への全射性に帰着する」と言う。
   ★**その全射性が `Proposition 2.5, (i)`**(`prop_2_5_i_surjective`)であり、
   base-trivial ⟹ metrically trivial(`isMetricallyTrivial_of_isBaseTrivial`)で
   仮定が合う。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★手 1 —— base-trivial なら isotropic 型 -/

include P in
/-- ★**base-trivial 型 ⟹ isotropic 型**。

★isotropic hull は pre-step ゆえ base-isomorphism なので、base-trivial 性で
`A ≅ A^istr`。同型で isotropic 性が移る。 -/
theorem isOfIsotropicType_of_baseTrivial (Fc : FrobenioidCore P)
    (hbt : IsOfBaseTrivialType P) : IsOfIsotropicType P := by
  intro A
  obtain ⟨Ai, ι, hι⟩ := Fc.isotropicHullExists A
  haveI : IsIso (P.Base ι) := hι.2.1.2
  obtain ⟨e⟩ := hbt A Ai ⟨asIso (P.Base ι)⟩
  exact isIsotropic_of_iso P e.symm hι.2.2.1

/-! ## ★手 2 —— 忠実性 -/

include P in
/-- ★★**忠実性** —— base-trivial ＋ unit-trivial なら `𝒞 → 𝔽_Φ` は忠実。 -/
theorem toElem_faithful_of (Fc : FrobenioidCore P)
    (hbt : IsOfBaseTrivialType P) (hut : IsOfUnitTrivialType P) : P.toElem.Faithful := by
  haveI hiso := isOfIsotropicType_of_baseTrivial P Fc hbt
  refine ⟨fun {A B} {f g} h => ?_⟩
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ :=
    prop_3_3_ii_sufficiency P Fc hiso
      (congrArg ElemFrobCat.Hom.deg h) (congrArg ElemFrobCat.Hom.div h)
      (congrArg ElemFrobCat.Hom.base h)
  have hd1 : (δ : Cc ⟶ Cc) = 𝟙 Cc := by
    have hb : δ ∈ (⊥ : Submonoid (End Cc)) := by rw [← hut Cc]; exact hδ
    simpa using hb
  rw [h₁, h₂, hd1]
  simp

/-! ## ★手 3 —— 本質的全射性 -/

include P in
/-- ★**本質的全射性** —— `Definition 1.3, (i), (a)`(`baseSurj`)そのもの。

★`𝔽_Φ` の対象は底だけで決まり、底の同型 `e` から `⟨e, 0, 1⟩` が同型を与える。 -/
theorem toElem_essSurj_of (Fc : FrobenioidCore P) : P.toElem.EssSurj := by
  refine ⟨fun X => ?_⟩
  obtain ⟨A, -, ⟨e⟩⟩ := Fc.baseSurj X.base
  haveI : IsIso ((⟨e.hom, 0, 1⟩ : P.toElem.obj A ⟶ X)) :=
    (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, isAddUnit_zero, rfl⟩
  exact ⟨A, ⟨asIso ((⟨e.hom, 0, 1⟩ : P.toElem.obj A ⟶ X))⟩⟩

/-! ## ★手 4 —— 充満性

★**3 つの部品**を組む: base-identity な Frobenius 型自己射(次数を担う)、
`𝒪^▷` の元(`Div` を担う)、底だけを動かす射(`Base` を担う)。
-/

include P in
/-- ★**部品 1** —— base-trivial ＋ Aut-ample なら、各次数に
**base-identity な Frobenius 型自己射**がある。

★`Definition 1.3, (ii)` で次数 `n` の Frobenius 型射を取り、base-trivial で
終域を `A` に戻し、Aut-ample でずれた底の自己同型を打ち消す。 -/
theorem frobEndo_of (Fc : FrobenioidCore P) (hbt : IsOfBaseTrivialType P)
    (haa : IsOfAutAmpleType P) (A : C) (n : ℕ+) :
    ∃ γ : A ⟶ A, IsBaseIdentity P γ ∧ IsFrobeniusType P γ ∧ P.degFr γ = n := by
  obtain ⟨A', γ₀, hγ₀, hd₀⟩ := Fc.frobDegSurj A n
  haveI hb₀ : IsIso (P.Base γ₀) := hγ₀.2
  obtain ⟨e⟩ := hbt A A' ⟨asIso (P.Base γ₀)⟩
  have hγ₁ : IsFrobeniusType P (γ₀ ≫ e.hom) :=
    IsFrobeniusType.comp P Fc hγ₀ (isFrobeniusType_of_isIso P e.hom)
  haveI hb₁ : IsIso (P.Base (γ₀ ≫ e.hom)) := hγ₁.2
  obtain ⟨u, hui, hub⟩ := haa A (P.Base (γ₀ ≫ e.hom)) hb₁
  haveI := hui
  obtain ⟨w, hw1, hw2⟩ : ∃ w : A ⟶ A, (u : A ⟶ A) ≫ w = 𝟙 A ∧ w ≫ (u : A ⟶ A) = 𝟙 A :=
    ⟨inv (u : A ⟶ A), IsIso.hom_inv_id _, IsIso.inv_hom_id _⟩
  haveI hwi : IsIso w := ⟨(u : A ⟶ A), hw2, hw1⟩
  refine ⟨(γ₀ ≫ e.hom) ≫ w, ?_, ?_, ?_⟩
  · show P.Base ((γ₀ ≫ e.hom) ≫ w) = P.Base (𝟙 A)
    rw [P.Base_comp, P.Base_id, ← hub, ← P.Base_comp, hw1, P.Base_id]
  · exact IsFrobeniusType.comp P Fc hγ₁ (isFrobeniusType_of_isIso P w)
  · rw [P.degFr_comp, degFr_of_isIso P w, one_mul, P.degFr_comp,
      degFr_of_isIso P e.hom, one_mul, hd₀]

include P in
/-- ★**部品 3** —— base-trivial ＋ Aut-ample なら、`𝒟` の任意の射 `b` を
**底に持つ等長・linear な射**がある。

★`Definition 1.3, (i), (c)` で `b` を pull-back 射の底として実現し、
base-trivial で始域を `A` に戻し、Aut-ample でずれを打ち消す。 -/
theorem baseHom_of (Fc : FrobenioidCore P) (hbt : IsOfBaseTrivialType P)
    (haa : IsOfAutAmpleType P) (A B : C)
    (b : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) :
    ∃ α : A ⟶ B, P.Base α = b ∧ P.Div α = 0 ∧ P.degFr α = 1 := by
  obtain ⟨X, π, hπ, θ, hπb⟩ := plBk_realize P Fc B b
  obtain ⟨eX⟩ := hbt A X ⟨θ.symm⟩
  have hbc : IsIso (P.Base eX.inv ≫ θ.hom) := by
    haveI : IsIso (P.Base eX.inv) := isBaseIsomorphism_of_isIso P eX.inv
    infer_instance
  obtain ⟨u, hui, hub⟩ := haa A (P.Base eX.inv ≫ θ.hom) hbc
  haveI := hui
  obtain ⟨w, hw1, hw2⟩ : ∃ w : A ⟶ A, (u : A ⟶ A) ≫ w = 𝟙 A ∧ w ≫ (u : A ⟶ A) = 𝟙 A :=
    ⟨inv (u : A ⟶ A), IsIso.hom_inv_id _, IsIso.inv_hom_id _⟩
  haveI hwi : IsIso w := ⟨(u : A ⟶ A), hw2, hw1⟩
  have key : P.Base w ≫ (P.Base eX.inv ≫ θ.hom) = 𝟙 _ := by
    rw [← hub, ← P.Base_comp, hw2, P.Base_id]
  refine ⟨w ≫ eX.inv ≫ π, ?_, ?_, ?_⟩
  · rw [P.Base_comp, P.Base_comp, hπb, ← Category.assoc (P.Base eX.inv),
      ← Category.assoc, key, Category.id_comp]
  · have h1 : P.Div π = 0 := (Fc.pullBackLB π hπ).1.2
    have h2 : P.Div eX.inv = 0 := isIsometric_of_isIso P eX.inv
    have h3 : P.Div w = 0 := isIsometric_of_isIso P w
    rw [P.Div_comp, P.Div_comp, h1, h2, h3]
    simp
  · rw [P.degFr_comp, P.degFr_comp, degFr_of_isIso P w, degFr_of_isIso P eX.inv,
      (Fc.pullBackLB π hπ).2]
    simp

include P in
/-- ★★**充満性** —— base-trivial ＋ Aut-ample なら `𝒞 → 𝔽_Φ` は充満。

★`ψ = (b, d, n)` に対し `φ := γ ≫ β ≫ α`(次数・`Div`・底をそれぞれ担う 3 部品)。 -/
theorem toElem_full_of (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hbt : IsOfBaseTrivialType P) (haa : IsOfAutAmpleType P) : P.toElem.Full := by
  refine ⟨fun {A B} ψ => ?_⟩
  obtain ⟨γ, hγb, hγf, hγd⟩ := frobEndo_of P Fc hbt haa A ψ.deg
  obtain ⟨α, hαb, hαd, hαdeg⟩ := baseHom_of P Fc hbt haa A B ψ.base
  obtain ⟨β, hβ⟩ := prop_2_5_i_surjective P G
    (isMetricallyTrivial_of_isBaseTrivial P (hbt A)) (haa A) ψ.div
  have hγb' : P.Base γ = 𝟙 _ := by
    have h : P.Base γ = P.Base (𝟙 A) := hγb
    rwa [P.Base_id] at h
  have hβb : P.Base ((β : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((β : End A) : A ⟶ A) = P.Base (𝟙 A) := β.2.1
    rwa [P.Base_id] at h
  refine ⟨γ ≫ ((β : End A) : A ⟶ A) ≫ α, ?_⟩
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · show P.Base (γ ≫ ((β : End A) : A ⟶ A) ≫ α) = ψ.base
    rw [P.Base_comp, P.Base_comp, hγb', hβb, hαb, Category.id_comp, Category.id_comp]
  · show P.Div (γ ≫ ((β : End A) : A ⟶ A) ≫ α) = ψ.div
    rw [P.Div_comp, P.Div_comp, hαd, hβb, MonoidOn.map_id, hγb',
      show P.Div γ = 0 from hγf.1.2, hαdeg]
    simpa using hβ
  · show P.degFr (γ ≫ ((β : End A) : A ⟶ A) ≫ α) = ψ.deg
    rw [P.degFr_comp, P.degFr_comp, hαdeg,
      show P.degFr ((β : End A) : A ⟶ A) = 1 from β.2.2, hγd]
    simp

/-! ## ★本体 -/

include P in
/-- ★★★**[FrdI] Proposition 3.3, (v) の `⟸`** —— Aut-ample・unit-trivial・
base-trivial 型なら `𝒞 → 𝔽_Φ` は圏同値。 -/
theorem prop_3_3_v_mpr (Fc : FrobenioidCore P) (G : Frobenioid P)
    (haa : IsOfAutAmpleType P) (hut : IsOfUnitTrivialType P)
    (hbt : IsOfBaseTrivialType P) : P.toElem.IsEquivalence := by
  haveI := toElem_faithful_of P Fc hbt hut
  haveI := toElem_full_of P Fc G hbt haa
  haveI := toElem_essSurj_of P Fc
  exact ⟨inferInstance, inferInstance, ‹_›⟩

include P in
/-- ★★★**[FrdI] Proposition 3.3, (v) の `⟹`** —— `𝒞 → 𝔽_Φ` が圏同値なら
`𝒞` は Aut-ample・unit-trivial・base-trivial 型。

★**`𝔽_Φ` 側でそれぞれが成り立つ**(`Proposition 1.5, (i), (ii)`)ことを、
圏同値が `Base`・`Div`・`degFr` を保つことを通して移す。 -/
theorem prop_3_3_v_mp (heq : P.toElem.IsEquivalence) :
    IsOfAutAmpleType P ∧ IsOfUnitTrivialType P ∧ IsOfBaseTrivialType P := by
  haveI := heq
  refine ⟨?_, ?_, ?_⟩
  · intro A g hg
    haveI := hg
    haveI : IsIso ((⟨g, 0, 1⟩ : P.toElem.obj A ⟶ P.toElem.obj A)) :=
      (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, isAddUnit_zero, rfl⟩
    obtain ⟨φ, hφ⟩ := P.toElem.map_surjective (⟨g, 0, 1⟩ : P.toElem.obj A ⟶ P.toElem.obj A)
    haveI : IsIso (P.toElem.map φ) := by rw [hφ]; infer_instance
    exact ⟨φ, isIso_of_reflects_iso φ P.toElem, congrArg ElemFrobCat.Hom.base hφ⟩
  · intro A
    refine le_antisymm (fun δ hδ => ?_) bot_le
    haveI : IsIso ((δ : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso _).mp hδ.2
    have hdb : P.Base ((δ : A ⟶ A)) = 𝟙 _ := by
      have h : P.Base ((δ : A ⟶ A)) = P.Base (𝟙 A) := hδ.1.1
      rwa [P.Base_id] at h
    have hdd : P.degFr ((δ : A ⟶ A)) = 1 := hδ.1.2
    have hdv : P.Div ((δ : A ⟶ A)) = 0 := isIsometric_of_isIso P _
    have hmap : P.toElem.map ((δ : A ⟶ A)) = P.toElem.map (𝟙 A) := by
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · show P.Base _ = P.Base (𝟙 A)
        rw [hdb, P.Base_id]
      · show P.Div _ = P.Div (𝟙 A)
        rw [hdv, P.Div_id]
      · show P.degFr _ = P.degFr (𝟙 A)
        rw [hdd, P.degFr_id]
    exact Submonoid.mem_bot.mpr (P.toElem.map_injective hmap)
  · intro A Dd hbi
    obtain ⟨e⟩ := hbi
    haveI : IsIso ((⟨e.inv, 0, 1⟩ : P.toElem.obj Dd ⟶ P.toElem.obj A)) :=
      (ElemFrobCat.isIso_iff _).mpr ⟨inferInstance, isAddUnit_zero, rfl⟩
    obtain ⟨φ, hφ⟩ := P.toElem.map_surjective
      (⟨e.inv, 0, 1⟩ : P.toElem.obj Dd ⟶ P.toElem.obj A)
    haveI : IsIso (P.toElem.map φ) := by rw [hφ]; infer_instance
    haveI : IsIso φ := isIso_of_reflects_iso φ P.toElem
    exact ⟨asIso φ⟩

include P in
/-- ★★★**[FrdI] Proposition 3.3, (v)** —— 両向き。 -/
theorem prop_3_3_v (Fc : FrobenioidCore P) (G : Frobenioid P) :
    P.toElem.IsEquivalence ↔
      (IsOfAutAmpleType P ∧ IsOfUnitTrivialType P ∧ IsOfBaseTrivialType P) :=
  ⟨prop_3_3_v_mp P, fun ⟨haa, hut, hbt⟩ => prop_3_3_v_mpr P Fc G haa hut hbt⟩

end ABC3.Found.FrdI
