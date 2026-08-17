import ABC3.Found.FrdI.Prop33Classes

/-!
# [FrdI] Proposition 3.3, (iv) の締め —— `𝒞^un-tr` の 2 本の圏同値

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.60。

原文 (FrdI p.60):
> which is faithful and essentially surjective; moreover, this functor determines

★`FrobenioidCore` の 21 条は `Prop33UnTr.lean` で埋まった。
★★残るのは `Frobenioid` の **`Definition 1.3, (iii), (d)` の 2 本の圏同値**である。

## ★★測定(2026-08-17)—— なぜ 3 性質が別々の道具で落ちるか

★★**行き先 `Order(Φ(A_𝒟))` は前順序圏**である。★したがって:

| 性質 | 何に帰着するか | 道具 |
|---|---|---|
| **忠実性** | 添字圏の hom が **subsingleton** | ★`Under` 側は `𝒞^un-tr` が **totally epimorphic**、`Over` 側は **pre-step が mono** |
| **充満性** | 行き先の hom は自明なので **「射が 1 本あればよい」** | ★`𝒞^istr` 側の充満性へ持ち上げて押し出す |
| **本質的全射性** | `𝒞^istr` の証人を押し出すだけ | ★`Div` が商で不変(`istrPre_Div`) |

★★**要点は「行き先が前順序圏だから充満性が存在問題に潰れる」**ことである。
★これを使わずに `map` の等式まで作ろうとすると、`Order` の hom の等式で
`Subsingleton` を引き直すことになる。

## ★★★co-angular が自動であることの効き方

★`𝒞^un-tr` は全対象が isotropic なので **すべての射が co-angular**
(`unTr_coAngular`)。★したがって `coaPreProp (unTrPre)` は
**pre-step であること**に潰れ、持ち上げのときに co-angular 性を別途示す必要がない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★0. pull-back の類型 —— 9 クラス目 -/

include P in
/-- ★★★**pull-back 射も商を渡る**(原文の 9 クラスの最後)。

原文 (FrdI p.60):
> phism; morphism of a given Frobenius degree) if and only if it arises from such

★`⟸` は `unTr_isPullBack_of_istr`、★`⟹` は `unTr_isPullBack_to`。 -/
theorem unTr_isPullBack_iff (Fc : FrobenioidCore P) {A B : Istr P} (α₀ : A ⟶ B) :
    IsPullBack (unTrPre P Fc) ((istrToUnTr P).map α₀) ↔ IsPullBack (istrPre P Fc) α₀ :=
  ⟨unTr_isPullBack_to P Fc α₀, unTr_isPullBack_of_istr P Fc α₀⟩

/-! ## ★1. 持ち上げ —— `𝒞^un-tr` の co-angular pre-step は `𝒞^istr` から来る -/

include P in
/-- ★`𝒞^un-tr` の pre-step は `𝒞^istr` の pre-step から来る。 -/
theorem unTr_isPreStep_lift (Fc : FrobenioidCore P) {A B : Istr P}
    (f : (show UnTr P from A) ⟶ (show UnTr P from B))
    (hf : IsPreStep (unTrPre P Fc) f) :
    ∃ f₀ : A ⟶ B, (istrToUnTr P).map f₀ = f ∧
      IsCoAngular (istrPre P Fc) f₀ ∧ IsPreStep (istrPre P Fc) f₀ := by
  obtain ⟨f₀, hf₀⟩ := (istrToUnTr P).map_surjective f
  refine ⟨f₀, hf₀, istr_coAngular_all P Fc f₀, ?_⟩
  refine (unTr_isPreStep_iff P Fc f₀.hom).mp ?_
  show IsPreStep (unTrPre P Fc) ((istrToUnTr P).map f₀)
  rw [hf₀]
  exact hf

/-! ## ★2. 添字圏の hom は subsingleton

★★`Under` 側は **`𝒞^un-tr` が totally epimorphic**(`unTr_totEpi`)、
`Over` 側は **pre-step が mono**(`unTr_preStepMono`)で落ちる。 -/

section Coa

variable (Fc : FrobenioidCore P)
  [MorphismProperty.IsMultiplicative (coaPreProp (unTrPre P Fc))]

include P in
/-- ★★**コスライスの hom は高々 1 本** —— `Z.hom` が epi だから。 -/
theorem unTr_coaPreUnder_subsingleton (A : UnTr P)
    (Z W : Under (⟨A⟩ : WideSubcategory (coaPreProp (unTrPre P Fc)))) :
    Subsingleton (Z ⟶ W) := by
  refine ⟨fun f g => ?_⟩
  haveI : Epi Z.hom.hom := unTr_totEpi P _ _ Z.hom.hom
  refine Under.UnderMorphism.ext (InducedWideCategory.Hom.ext ?_)
  refine (cancel_epi Z.hom.hom).mp ?_
  exact (congrArg InducedWideCategory.Hom.hom (Under.w f)).trans
    (congrArg InducedWideCategory.Hom.hom (Under.w g)).symm

include P in
/-- ★★**スライスの hom は高々 1 本** —— `W.hom` が pre-step ゆえ mono だから。 -/
theorem unTr_coaPreOver_subsingleton (A : UnTr P)
    (Z W : Over (⟨A⟩ : WideSubcategory (coaPreProp (unTrPre P Fc)))) :
    Subsingleton (Z ⟶ W) := by
  refine ⟨fun f g => ?_⟩
  haveI : Mono W.hom.hom := unTr_preStepMono P Fc W.hom.hom W.hom.property.2
  refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
  refine (cancel_mono W.hom.hom).mp ?_
  exact (congrArg InducedWideCategory.Hom.hom (Over.w f)).trans
    (congrArg InducedWideCategory.Hom.hom (Over.w g)).symm

/-! ## ★3. 忠実性 -/

include P in
theorem unTr_coaPreUnder_faithful (A : UnTr P) :
    (coaPreUnderFunctor (unTrPre P Fc) A).Faithful where
  map_injective {Z W} {f g} _ := (unTr_coaPreUnder_subsingleton P Fc A Z W).elim f g

include P in
theorem unTr_coaPreOver_faithful (A : UnTr P) :
    (coaPreOverFunctor (unTrPre P Fc) A).Faithful where
  map_injective {Z W} {f g} _ := (unTr_coaPreOver_subsingleton P Fc A Z W).elim f g

/-! ## ★4. 充満性 —— ★★**行き先が前順序圏なので「射が 1 本あればよい」**

★持ち上げて `𝒞^istr` 側の充満性を使い、押し出す。 -/

include P in
theorem unTr_isPreStep_map (Fc : FrobenioidCore P) {A B : Istr P} (f : A ⟶ B) :
    IsPreStep (unTrPre P Fc) ((istrToUnTr P).map f) ↔ IsPreStep (istrPre P Fc) f := Iff.rfl

include P in
theorem unTr_coaPreUnder_full (G : Frobenioid P) (A : UnTr P) :
    (coaPreUnderFunctor (unTrPre P Fc) A).Full := by
  letI := coaPreProp_isMultiplicative (istrPre P Fc) (istr_frobenioidCore P Fc).coAngularComp
  haveI := (istr_frobenioid P Fc G).coaPreUnderEquiv (show Istr P from A)
  constructor
  intro Z W g
  obtain ⟨-, Zr, ⟨zh, zp⟩⟩ := Z
  obtain ⟨-, Wr, ⟨wh, wp⟩⟩ := W
  obtain ⟨ζ, hζ, hζc, hζs⟩ := unTr_isPreStep_lift P Fc zh zp.2
  obtain ⟨ω, hω, hωc, hωs⟩ := unTr_isPreStep_lift P Fc wh wp.2
  subst hζ
  subst hω
  set Zι : Under (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc))) :=
    Under.mk (show (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc)))
      ⟶ ⟨show Istr P from Zr.obj⟩ from ⟨ζ, hζc, hζs⟩) with hZι
  set Wι : Under (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc))) :=
    Under.mk (show (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc)))
      ⟶ ⟨show Istr P from Wr.obj⟩ from ⟨ω, hωc, hωs⟩) with hWι
  have hle : (coaPreUnderFunctor (istrPre P Fc) (show Istr P from A)).obj Zι
      ≤ (coaPreUnderFunctor (istrPre P Fc) (show Istr P from A)).obj Wι := leOfHom g
  let h := (coaPreUnderFunctor (istrPre P Fc) (show Istr P from A)).preimage (homOfLE hle)
  have hw0 : (Zι.hom ≫ h.right).hom = Wι.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Under.w h)
  have hw : ζ ≫ h.right.hom = ω := hw0
  refine ⟨Under.homMk (show Zr ⟶ Wr from
      ⟨(istrToUnTr P).map h.right.hom, unTr_coAngular P Fc _,
        (unTr_isPreStep_map P Fc h.right.hom).mpr h.right.property.2⟩) ?_, ?_⟩
  · refine InducedWideCategory.Hom.ext ?_
    show (istrToUnTr P).map ζ ≫ (istrToUnTr P).map h.right.hom = (istrToUnTr P).map ω
    exact ((istrToUnTr P).map_comp ζ h.right.hom).symm.trans
      (congrArg (fun t => (istrToUnTr P).map t) hw)
  · exact Subsingleton.elim _ _

include P in
theorem unTr_coaPreOver_full (G : Frobenioid P) (A : UnTr P) :
    (coaPreOverFunctor (unTrPre P Fc) A).Full := by
  letI := coaPreProp_isMultiplicative (istrPre P Fc) (istr_frobenioidCore P Fc).coAngularComp
  haveI := (istr_frobenioid P Fc G).coaPreOverEquiv (show Istr P from A)
  constructor
  intro Z W g
  obtain ⟨Zl, -, ⟨zh, zp⟩⟩ := Z
  obtain ⟨Wl, -, ⟨wh, wp⟩⟩ := W
  obtain ⟨ζ, hζ, hζc, hζs⟩ := unTr_isPreStep_lift P Fc zh zp.2
  obtain ⟨ω, hω, hωc, hωs⟩ := unTr_isPreStep_lift P Fc wh wp.2
  subst hζ
  subst hω
  set Zι : Over (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc))) :=
    Over.mk (show (⟨show Istr P from Zl.obj⟩ : WideSubcategory (coaPreProp (istrPre P Fc)))
      ⟶ ⟨(show Istr P from A)⟩ from ⟨ζ, hζc, hζs⟩) with hZι
  set Wι : Over (⟨(show Istr P from A)⟩ : WideSubcategory (coaPreProp (istrPre P Fc))) :=
    Over.mk (show (⟨show Istr P from Wl.obj⟩ : WideSubcategory (coaPreProp (istrPre P Fc)))
      ⟶ ⟨(show Istr P from A)⟩ from ⟨ω, hωc, hωs⟩) with hWι
  have hle : (coaPreOverFunctor (istrPre P Fc) (show Istr P from A)).obj Zι
      ⟶ (coaPreOverFunctor (istrPre P Fc) (show Istr P from A)).obj Wι := g
  let h := (coaPreOverFunctor (istrPre P Fc) (show Istr P from A)).preimage hle
  have hw0 : (h.left ≫ Wι.hom).hom = Zι.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Over.w h)
  have hw : h.left.hom ≫ ω = ζ := hw0
  refine ⟨Over.homMk (show Zl ⟶ Wl from
      ⟨(istrToUnTr P).map h.left.hom, unTr_coAngular P Fc _,
        (unTr_isPreStep_map P Fc h.left.hom).mpr h.left.property.2⟩) ?_, ?_⟩
  · refine InducedWideCategory.Hom.ext ?_
    show (istrToUnTr P).map h.left.hom ≫ (istrToUnTr P).map ω = (istrToUnTr P).map ζ
    exact ((istrToUnTr P).map_comp h.left.hom ω).symm.trans
      (congrArg (fun t => (istrToUnTr P).map t) hw)
  · exact Subsingleton.elim _ _

/-! ## ★5. 本質的全射性 —— `𝒞^istr` の証人を押し出す -/

include P in
theorem unTr_coaPreUnder_essSurj (G : Frobenioid P) (A : UnTr P) :
    (coaPreUnderFunctor (unTrPre P Fc) A).EssSurj := by
  letI := coaPreProp_isMultiplicative (istrPre P Fc) (istr_frobenioidCore P Fc).coAngularComp
  haveI := (istr_frobenioid P Fc G).coaPreUnderEquiv (show Istr P from A)
  refine ⟨fun x => ?_⟩
  obtain ⟨Zι, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := coaPreUnderFunctor (istrPre P Fc) (show Istr P from A)) x
  refine ⟨Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (unTrPre P Fc)))
      ⟶ ⟨Zι.right.obj⟩ from ⟨(istrToUnTr P).map Zι.hom.hom, unTr_coAngular P Fc _,
        (unTr_isPreStep_map P Fc Zι.hom.hom).mpr Zι.hom.property.2⟩), ⟨e⟩⟩

include P in
theorem unTr_coaPreOver_essSurj (G : Frobenioid P) (A : UnTr P) :
    (coaPreOverFunctor (unTrPre P Fc) A).EssSurj := by
  letI := coaPreProp_isMultiplicative (istrPre P Fc) (istr_frobenioidCore P Fc).coAngularComp
  haveI := (istr_frobenioid P Fc G).coaPreOverEquiv (show Istr P from A)
  refine ⟨fun x => ?_⟩
  obtain ⟨Zι, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := coaPreOverFunctor (istrPre P Fc) (show Istr P from A)) x
  refine ⟨Over.mk (show (⟨Zι.left.obj⟩ : WideSubcategory (coaPreProp (unTrPre P Fc)))
      ⟶ ⟨A⟩ from ⟨(istrToUnTr P).map Zι.hom.hom, unTr_coAngular P Fc _,
        (unTr_isPreStep_map P Fc Zι.hom.hom).mpr Zι.hom.property.2⟩), ⟨e⟩⟩

end Coa

/-! ## ★6. ★★★★**`𝒞^un-tr` は Frobenioid である** -/

include P in
/-- ★★★★**[FrdI] Proposition 3.3, (iv)** —— `𝒞^un-tr → 𝔽_Φ` は
`𝒞^un-tr` に **Frobenioid の構造**を定める。

原文 (FrdI p.60):
> a natural structure of Frobenioid on Cun-tr, with respect to which Cun-tr is of

★`FrobenioidCore` の 21 条は `unTr_frobenioidCore`(`Prop33UnTr.lean`)。
★`Definition 1.3, (iii), (d)` の 2 本の圏同値は本ファイル。 -/
theorem unTr_frobenioid (Fc : FrobenioidCore P) (G : Frobenioid P) :
    Frobenioid (unTrPre P Fc) := by
  haveI := unTr_coaPreProp_isMultiplicative P Fc
  exact ⟨unTr_frobenioidCore P Fc,
    fun A => ⟨unTr_coaPreUnder_faithful P Fc A, unTr_coaPreUnder_full P Fc G A,
      unTr_coaPreUnder_essSurj P Fc G A⟩,
    fun A => ⟨unTr_coaPreOver_faithful P Fc A, unTr_coaPreOver_full P Fc G A,
      unTr_coaPreOver_essSurj P Fc G A⟩⟩

/-! ## ★7. `Proposition 3.3` の 5 主張の対応

| 主張 | 内容 | 実装 |
|---|---|---|
| (i) | `End(𝒞^pl-bk_A → 𝒞)^bs-iso` の 2 主張 | `prop_3_3_i_mem_oTri` / `prop_3_3_i_converse`(`Prop33i.lean`) |
| (ii) | unit-equivalent ⟺ `𝔽_Φ` で同じ | `prop_3_3_ii`(`Prop33.lean`) |
| (iii) | `𝒞^istr → 𝒞^un-tr` は full ＋ 本質的全射、忠実 ⟺ unit-trivial | `prop_3_3_iii`(`UnTr.lean`) |
| (iv) | `𝒞^un-tr → 𝔽_Φ` は忠実 ＋ 本質的全射、Frobenioid 構造、isotropic ＋ unit-trivial 型、**9 クラスの移送** | `unTrToElem_essSurj` / `unTr_frobenioid` / `unTr_isotropic` / `unTr_unitTrivial` / `UnTr.lean` ＋ `Prop33Classes.lean` ＋ `unTr_isPullBack_iff` |
| (v) | `𝒞 → 𝔽_Φ` が圏同値 ⟺ Aut-ample ＋ unit-trivial ＋ base-trivial | `prop_3_3_v`(`Prop33v.lean`) |
-/

/-- ★★★★**[FrdI] Proposition 3.3** —— 5 主張がすべて実装された。 -/
def prop_3_3.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59, item := "Proposition 3.3",
    sectionId := "frdi-prop-3-3" }

end ABC3.Found.FrdI
