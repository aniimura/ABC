import ABC3.Found.FrdI.Prop114
import Mathlib.CategoryTheory.InducedCategory

/-!
# [FrdI] Proposition 2.2 —— The Functor `𝒪^▷(−)`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.45–p.46。

原文 (FrdI p.45):
> Proposition 2.2.

原文 (FrdI p.45):
> objects are the objects of Cistr and whose morphisms are given as follows:

## ★`𝒟*` の正体

原文の
> HomD∗(A, B)
> = HomD(AD, BD)

は ★**mathlib の `InducedCategory`** そのものである ——
対象は `Ob(𝒞^istr)`、射は `Base` で送った先の `𝒟` の射。
★**したがって `𝒟* → 𝒟` が充満忠実であることは構成から出る**(`inducedFunctor`)。
原文が「fully faithfulness follows formally from the definition of the category D∗」
と書くのは、この意味である。

## ★この命題の規模(測定)

**4 条 (i)–(iv)**、主張の数は **9**:

| 条 | 主張の数 | 内容 |
|---|---|---|
| (i) | 1 | `𝒟* → 𝒟` が圏同値 |
| (ii) | 4 | 反変関手の存在 / 一意性 / (a) linear の場合 / (b) pre-step の場合 |
| (iii) | 3 | `𝒪^×` が部分関手 / `𝒪^▷(−)^±` と一致 / `Div` が関手的準同型で `𝒪^▷(A)^char ↪ Φ(A)` |
| (iv) | 1 | isotropic hull が `𝒪^▷(A) ↪ 𝒪^▷(A^istr)` を誘導 |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★`𝒟*` —— 対象は `𝒞^istr`、射は `𝒟` から

原文 (FrdI p.45):
> [where A, B ∈Ob(Cistr); AD
-/

/-- ★`𝒞^istr` の対象を `𝒟` へ送る写像。 -/
def istrBase (A : Istr P) : D := (P.toElem.obj A.obj).base

/-- ★★**`𝒟*`** —— 対象は `Ob(𝒞^istr)`、`Hom_{𝒟*}(A,B) = Hom_𝒟(A_𝒟, B_𝒟)`。

★★**mathlib の `InducedCategory` そのもの**なので、`𝒟* → 𝒟`(`inducedFunctor`)の
充満忠実性は**構成から出る**。★原文が
「fully faithfulness follows formally from the definition of the category D∗」
と書くのは、この意味である。

★**型シノニムを別名で包むとインスタンス解決が `Istr P` 側に落ちる**ので、
ここでは `InducedCategory` を直接使う。 -/
abbrev dStarToD : InducedCategory D (istrBase P) ⥤ D := inducedFunctor (istrBase P)

/-! ## ★(i) —— `𝒟* → 𝒟` は圏同値

原文 (FrdI p.46):
> tion 1.3, (i), (a) [i.e., applied to the Frobenioid Cistr — cf. Proposition 1.9, (v)],

★**充満忠実性は `InducedCategory` の構成から**、
★**本質的全射性は `Definition 1.3, (i), (a)` を `𝒞^istr` に当てる**。
-/

include F in
/-- ★★**[FrdI] Proposition 2.2, (i)** —— `𝒟* → 𝒟` は圏同値。

★**本質的全射性だけが中身**である —— `Definition 1.3, (i), (a)` が
「`𝒟` のどの同型類も Frobenius-trivial な対象から来る」と言い、
`Proposition 1.9, (v)`(`istrHullExists`)がそれを `𝒞^istr` へ運ぶ。 -/
theorem prop_2_2_i : (dStarToD P).IsEquivalence := by
  haveI : (dStarToD P).EssSurj := by
    constructor
    intro Y
    obtain ⟨A, -, ⟨e⟩⟩ := F.baseSurj Y
    obtain ⟨B, φ, hφ⟩ := F.isotropicHullExists A
    haveI hbi : IsIso (P.Base φ) := hφ.2.1.2
    refine ⟨(⟨B, hφ.2.2.1⟩ : Istr P), ⟨?_⟩⟩
    exact (@asIso _ _ _ _ (P.Base φ) hbi).symm ≪≫ e
  exact ⟨inferInstance, inferInstance, ‹_›⟩

/-! ## ★(ii) —— 反変関手 `𝒪^▷(−)`

原文 (FrdI p.45):
> (ii) There is a unique contravariant functor

原文 (FrdI p.45):
> is the bijection of Definition 1.3, (iii), (c). By abuse of notation, we shall

★**中身は 4 つ**:
1. 関手の存在
2. 一意性
3. (a) linear 射での値が `Proposition 1.11, (iv)` の埋め込み
4. (b) pre-step での値が `Definition 1.3, (iii), (c)` の全単射

★★**`𝒞^istr` では co-angular 性が自動である**(`Proposition 1.4, (i)`、
`isCoAngular_of_isotropic_dom`)——原文が 2 度くり返す
「necessarily co-angular」がこれである。
★**したがって (a) と (b) は 1 本の構成に統合できる** ——
pre-step は linear なので、`Proposition 1.11, (iv)` の埋め込みが
pre-step のときにちょうど `Definition 1.3, (iii), (c)` の全単射になる。
-/

section FunctorII

/-- ★★**(ii) の射の対応** —— `𝒞^istr` の linear 射 `φ : A ⟶ B` に沿う
`𝒪^▷(B) →* 𝒪^▷(A)`。

★**`Proposition 1.11, (iv)` の埋め込みそのもの**であり(＝原文の (a))、
`𝒞^istr` では co-angular 性が自動で付く。 -/
noncomputable def otriLin {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) : OTri P B →* OTri P A :=
  otriPullHom P F φ (isCoAngular_of_isotropic_dom P F hA φ) hl

/-- ★`otriLin` を特徴づける四角形。 -/
theorem otriLin_spec {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) (γ : OTri P B) :
    φ ≫ ((γ : End B) : B ⟶ B)
      = ((otriLin P F hA hl γ : End A) : A ⟶ A) ≫ φ :=
  otriPull_spec P F φ (isCoAngular_of_isotropic_dom P F hA φ) hl γ

/-- ★★**`otriLin` の一意性** —— 四角形を満たす `𝒪^▷(A)` の元は 1 つだけ。

★**これが (ii) の「unique」の中身**である。 -/
theorem otriLin_uniq {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) (γ : OTri P B) (δ : OTri P A)
    (h : φ ≫ ((γ : End B) : B ⟶ B) = ((δ : End A) : A ⟶ A) ≫ φ) :
    δ = otriLin P F hA hl γ :=
  Subtype.ext (prop_1_11_iv_unique' P F φ (isCoAngular_of_isotropic_dom P F hA φ) hl
    γ.2 δ.2 (otriLin P F hA hl γ).2 h (otriLin_spec P F hA hl γ))

/-- ★**(ii) の関手性 1** —— 恒等射では恒等準同型。 -/
theorem otriLin_id {A : C} (hA : IsIsotropic P A) :
    otriLin P F hA (isLinear_id P A) = MonoidHom.id (OTri P A) := by
  ext γ
  exact congrArg (fun x : OTri P A => ((x : End A) : A ⟶ A))
    (otriLin_uniq P F hA (isLinear_id P A) γ γ (by simp)).symm

/-- ★★**(ii) の関手性 2** —— 合成を(反変に)保つ。

★**一意性から出る**:`(φ ≫ ψ) ≫ γ_E = φ ≫ (γ_B ≫ ψ) = (γ_A ≫ φ) ≫ ψ`。 -/
theorem otriLin_comp {A B E : C} (hA : IsIsotropic P A) (hB : IsIsotropic P B)
    {φ : A ⟶ B} (hlφ : IsLinear P φ) {ψ : B ⟶ E} (hlψ : IsLinear P ψ) :
    otriLin P F hA (IsLinear.comp P hlφ hlψ)
      = (otriLin P F hA hlφ).comp (otriLin P F hB hlψ) := by
  ext γ
  exact congrArg (fun x : OTri P A => ((x : End A) : A ⟶ A))
    (otriLin_uniq P F hA (IsLinear.comp P hlφ hlψ) γ
      (otriLin P F hA hlφ (otriLin P F hB hlψ γ)) (by
        rw [Category.assoc, otriLin_spec P F hB hlψ γ, ← Category.assoc,
          otriLin_spec P F hA hlφ (otriLin P F hB hlψ γ), Category.assoc])).symm

/-- ★**(a)** `otriLin` は単射 —— `Proposition 1.11, (iv)` の「injection」。 -/
theorem otriLin_injective {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) : Function.Injective (otriLin P F hA hl) :=
  otriPullHom_injective P F φ (isCoAngular_of_isotropic_dom P F hA φ) hl

/-- ★★★**(b)** `φ` が pre-step なら `otriLin` は**全単射** ——
`Definition 1.3, (iii), (c)` の全単射に一致する。

★**全射性は `otriFwd`(順方向)から**、単射性は `Proposition 1.11, (iv)` から。
★**両者が同じ四角形で特徴づけられるので、写像として一致する。** -/
theorem otriLin_bijective_of_preStep {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hs : IsPreStep P φ) : Function.Bijective (otriLin P F hA hs.1) := by
  refine ⟨otriLin_injective P F hA hs.1, fun α => ?_⟩
  obtain ⟨β, ⟨hβm, hβe⟩, -⟩ :=
    F.otriFwd φ (isCoAngular_of_isotropic_dom P F hA φ) hs (α : End A) α.2
  exact ⟨⟨β, hβm⟩, (otriLin_uniq P F hA hs.1 ⟨β, hβm⟩ α hβe).symm⟩

/-- ★★**(b)** `Definition 1.3, (iii), (c)` の全単射との一致(逆方向)。

`otriFwd` が与える `β` は、ちょうど `otriLin` の逆像である。 -/
theorem otriLin_otriFwd {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hs : IsPreStep P φ) (α : OTri P A) {β : End B} (hβm : β ∈ OTri P B)
    (hβe : φ ≫ β = ((α : End A) : A ⟶ A) ≫ φ) :
    otriLin P F hA hs.1 ⟨β, hβm⟩ = α :=
  (otriLin_uniq P F hA hs.1 ⟨β, hβm⟩ α hβe).symm

/-! ### ★★関手として束ねる —— `(𝒞^istr)^lin → Mon`

原文 (FrdI p.45):
> the restriction of this functor on D∗to (Cistr)lin. Finally,

★**対象は `𝒞^istr` の対象、射は linear なものだけ**。 -/

/-- ★**`(𝒞^istr)^lin`** —— `𝒞^lin` の isotropic 対象がなす充満部分圏。 -/
abbrev IstrLin : Type u2 :=
  ObjectProperty.FullSubcategory (fun A : Lin P => IsIsotropic P A.obj)

/-- ★★★**[FrdI] Proposition 2.2, (ii)(の `(𝒞^istr)^lin` への制限)**
—— 反変関手 `𝒪^▷(−) : (𝒞^istr)^lin → Mon`。

★**射の対応は `Proposition 1.11, (iv)` の埋め込み**(＝原文の (a))、
★**pre-step の上ではそれが `Definition 1.3, (iii), (c)` の全単射になる**(＝原文の (b))。 -/
noncomputable def otriFunctorLin : (IstrLin P)ᵒᵖ ⥤ MonCat.{v2} where
  obj A := MonCat.of (OTri P A.unop.obj.obj)
  map {_ B} f := MonCat.ofHom (otriLin P F B.unop.property f.unop.hom.2)
  map_id A := by
    apply MonCat.hom_ext
    exact otriLin_id P F A.unop.property
  map_comp {_ B E} f g := by
    apply MonCat.hom_ext
    exact otriLin_comp P F E.unop.property B.unop.property g.unop.hom.2 f.unop.hom.2

end FunctorII

/-! ## ★(iii) —— `𝒪^×` は `𝒪^▷` の単元群、`Div` はその核を潰す

原文 (FrdI p.45):
> functor of the functor of (ii) which is equal to the subfunctor

原文 (FrdI p.45):
> the notation of §0]. Moreover, the operation “Div(−)” determines a functorial

★**中身は 3 つ**:
1. `𝒪^×(A)` は `𝒪^▷(A)` の**単元のなす部分モノイド**(＝ `𝒪^▷(A)^±`)
2. `Div : 𝒪^▷(A) → Φ(A)` は**モノイド準同型**
   (★`End` の積は `x * y = y ≫ x`、`Div` は底恒等・linear では**素直に足し算**になる)
3. ★★**その核がちょうど `𝒪^×(A)`** —— これが `𝒪^▷(A)^char ↪ Φ(A)` の中身。
   ★**`A` が isotropic であることをここで使う**(`Div = 0` の pre-step が同型になる)。
-/

/-- ★**(iii)-2** `Div` は `𝒪^▷(A)` 上で加法的。

★`End` の積は `x * y = y ≫ x` で、底恒等・linear なら
`Div (y ≫ x) = Φ(𝟙)(Div x) + 1 · Div y = Div x + Div y`。 -/
theorem otri_Div_mul {A : C} (x y : OTri P A) :
    P.Div ((x * y).1 : A ⟶ A) = P.Div (x.1 : A ⟶ A) + P.Div (y.1 : A ⟶ A) := by
  show P.Div ((y.1 : A ⟶ A) ≫ (x.1 : A ⟶ A)) = _
  rw [P.Div_comp, show P.degFr (x.1 : A ⟶ A) = 1 from x.2.2,
    show P.Base (y.1 : A ⟶ A) = 𝟙 _ from by
      have h : P.Base (y.1 : A ⟶ A) = P.Base (𝟙 A) := y.2.1
      rwa [P.Base_id] at h,
    MonoidOn.map_id]
  simp [add_comm]

/-- ★**(iii)-2** `Div` の単位元。 -/
theorem otri_Div_one {A : C} : P.Div ((1 : OTri P A).1 : A ⟶ A) = 0 := P.Div_id A

include F in
/-- ★★★**(iii)-3** `𝒪^▷(A)` の元が単元 ⟺ `Div = 0`。

★★**`𝒪^×(A)` は `Div` の核**である —— これが原文の
`𝒪^▷(A)^char = 𝒪^▷(A)/𝒪^×(A) ↪ Φ(A)` の中身。

★`⟸` に **`A` が isotropic** であることを使う(`Div = 0` の pre-step は同型)。 -/
theorem otri_isUnit_iff_Div_zero {A : C} (hiso : IsIsotropic P A) (x : OTri P A) :
    IsIso (x.1 : A ⟶ A) ↔ P.Div (x.1 : A ⟶ A) = 0 := by
  constructor
  · intro h
    haveI := h
    exact isIsometric_of_isIso P _
  · intro h
    exact hiso A (x.1 : A ⟶ A) h ⟨x.2.2, by
      have hb : P.Base (x.1 : A ⟶ A) = P.Base (𝟙 A) := x.2.1
      rw [P.Base_id] at hb
      show IsIso (P.Base (x.1 : A ⟶ A))
      rw [hb]
      infer_instance⟩

include F in
/-- ★★★**(iii)-3 の本体** —— `𝒪^▷(A)^char = 𝒪^▷(A)/𝒪^×(A) ↪ Φ(A)` の**単射性**。

★`Div x = Div y` なら `x` と `y` は**単元 1 つ分しか違わない** ——
これが「induces an inclusion」の中身である。

★★**理由は `Proposition 1.11, (vi)`(`faithfulUpToUnits`)**:
`𝒪^▷(A)` の元は base-identity かつ linear なので **pre-step** であり、
`Div` が等しければ metrically equivalent、底はどちらも `𝟙` だから base-equivalent。
★co-angular 性は `A` が isotropic であることから自動で付く。 -/
theorem otri_div_eq_iff {A : C} (hiso : IsIsotropic P A) (x y : OTri P A) :
    P.Div ((x : End A) : A ⟶ A) = P.Div ((y : End A) : A ⟶ A)
      ↔ ∃ u : OTimes P A, ((x : End A) : A ⟶ A) = ((y : End A) : A ⟶ A)
          ≫ ((u : End A) : A ⟶ A) := by
  have hbx : P.Base ((x : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((x : End A) : A ⟶ A) = P.Base (𝟙 A) := x.2.1
    rwa [P.Base_id] at h
  have hby : P.Base ((y : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((y : End A) : A ⟶ A) = P.Base (𝟙 A) := y.2.1
    rwa [P.Base_id] at h
  have hsx : IsPreStep P ((x : End A) : A ⟶ A) :=
    ⟨x.2.2, by show IsIso (P.Base ((x : End A) : A ⟶ A)); rw [hbx]; infer_instance⟩
  have hsy : IsPreStep P ((y : End A) : A ⟶ A) :=
    ⟨y.2.2, by show IsIso (P.Base ((y : End A) : A ⟶ A)); rw [hby]; infer_instance⟩
  constructor
  · intro h
    obtain ⟨u, hu, he⟩ := F.faithfulUpToUnits ((x : End A) : A ⟶ A) ((y : End A) : A ⟶ A)
      (show P.Base _ = P.Base _ from by rw [hbx, hby]) h
      (isCoAngular_of_isotropic_dom P F hiso _) hsx
      (isCoAngular_of_isotropic_dom P F hiso _) hsy
    exact ⟨⟨u, hu⟩, he⟩
  · rintro ⟨u, he⟩
    haveI : IsIso ((u : End A) : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso _).mp u.2.2
    rw [he, P.Div_comp, hby, MonoidOn.map_id,
      show P.Div ((u : End A) : A ⟶ A) = 0 from isIsometric_of_isIso P _,
      show P.degFr ((u : End A) : A ⟶ A) = 1 from u.2.1.2]
    simp

/-! ## ★(iv) —— isotropic hull が誘導する `𝒪^▷` の埋め込み

原文 (FrdI p.45):
> (iv) If φ : A →Aistr is an isotropic hull in C, then φ determines a natural

★**中身は isotropic hull の普遍性そのもの** —— `α ∈ 𝒪^▷(A)` に対し
`α ≫ φ : A ⟶ A^istr` を普遍性に通すと**一意な** `β : A^istr ⟶ A^istr` が出る。
★**単射性は `φ` が pre-step ゆえ mono**(`Definition 1.3, (v), (a)`)であることから。
-/

section IsotropicHull

variable {A B : C} (φ : A ⟶ B) (hφ : IsIsotropicHull P φ)

/-- ★普遍性が与える `β` —— `α ≫ φ = φ ≫ β`。 -/
noncomputable def hullOTriMap (α : End A) : End B :=
  (hφ.2.2.2 B hφ.2.2.1 ((α : A ⟶ A) ≫ φ)).choose

theorem hullOTriMap_sq (α : End A) :
    (α : A ⟶ A) ≫ φ = φ ≫ (hullOTriMap P φ hφ α : B ⟶ B) :=
  (hφ.2.2.2 B hφ.2.2.1 ((α : A ⟶ A) ≫ φ)).choose_spec.1

theorem hullOTriMap_uniq (α : End A) (β : End B)
    (h : (α : A ⟶ A) ≫ φ = φ ≫ (β : B ⟶ B)) : β = hullOTriMap P φ hφ α :=
  (hφ.2.2.2 B hφ.2.2.1 ((α : A ⟶ A) ≫ φ)).choose_spec.2 β h

/-- ★★**(iv)** `𝒪^▷(A) → 𝒪^▷(A^istr)` はモノイド準同型。

★`End` の積は `x * y = y ≫ x` なので、四角形を 2 つ繋げばそのまま出る。 -/
noncomputable def hullOTriHom : End A →* End B where
  toFun := hullOTriMap P φ hφ
  map_one' := (hullOTriMap_uniq P φ hφ 1 1 (by simp)).symm
  map_mul' α₁ α₂ := by
    refine (hullOTriMap_uniq P φ hφ (α₁ * α₂) (hullOTriMap P φ hφ α₁ * hullOTriMap P φ hφ α₂)
      ?_).symm
    show ((α₂ : A ⟶ A) ≫ (α₁ : A ⟶ A)) ≫ φ
      = φ ≫ ((hullOTriMap P φ hφ α₂ : B ⟶ B) ≫ (hullOTriMap P φ hφ α₁ : B ⟶ B))
    conv_lhs => rw [Category.assoc, hullOTriMap_sq P φ hφ α₁, ← Category.assoc,
      hullOTriMap_sq P φ hφ α₂, Category.assoc]

include F in
/-- ★★**(iv) の「inclusion」** —— 単射である。

★`φ` は pre-step なので `Definition 1.3, (v), (a)` により mono。
`α₁ ≫ φ = φ ≫ β = α₂ ≫ φ` から `α₁ = α₂` が出る。 -/
theorem hullOTriHom_injective : Function.Injective (hullOTriHom P φ hφ) := by
  haveI : Mono φ := F.preStepMono φ hφ.2.1
  intro α₁ α₂ h
  have h1 := hullOTriMap_sq P φ hφ α₁
  have h2 := hullOTriMap_sq P φ hφ α₂
  have : (α₁ : A ⟶ A) ≫ φ = (α₂ : A ⟶ A) ≫ φ := by
    rw [h1, h2, show hullOTriMap P φ hφ α₁ = hullOTriMap P φ hφ α₂ from h]
  exact (cancel_mono φ).mp this

/-- ★★**(iv)** 像は `𝒪^▷(A^istr)` に入る —— 底恒等かつ linear が保たれる。 -/
theorem hullOTriHom_mem (α : End A) (hα : α ∈ OTri P A) :
    hullOTriHom P φ hφ α ∈ OTri P B := by
  haveI hbi : IsIso (P.Base φ) := hφ.2.1.2
  have hsq := hullOTriMap_sq P φ hφ α
  constructor
  · -- base-identity
    have hb := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1,
      P.Base_id, Category.id_comp] at hb
    show P.Base (hullOTriMap P φ hφ α : B ⟶ B) = P.Base (𝟙 B)
    rw [P.Base_id]
    exact ((cancel_epi (P.Base φ)).mp (by rw [← hb]; simp)).symm
  · -- linear
    have hd := congrArg P.degFr hsq
    rw [P.degFr_comp, P.degFr_comp, show P.degFr (α : A ⟶ A) = 1 from hα.2, mul_one] at hd
    show P.degFr (hullOTriMap P φ hφ α : B ⟶ B) = 1
    exact (mul_right_cancel (b := P.degFr φ) (by rw [← hd, one_mul])).symm

end IsotropicHull

end ABC3.Found.FrdI
