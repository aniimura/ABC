/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Quasi
import Mathlib.CategoryTheory.Localization.Predicate

/-!
# [FrdI] Theorem 3.4, (iv) —— `End(𝒞^pl-bk_A → 𝒞)^bs-iso` の移送

★`Thm34Quasi.lean` の `otriEndCond_map` は `OTriGenCond` の移送を
**仮定として受ける**形で実装した。ここではその仮定を作りにいく。

原文 (FrdI p.66):
> By assertions (ii), (iii), it follows that Ψ preserves pre-steps, base-isomorphisms,

★★段は 3 つ:

| 段 | 中身 | 状態 |
|---|---|---|
| 1 | `Ψ` から `𝒞^pl-bk` の間の関手 `plBkMap` | 本ファイル |
| 2 | 圏同値なら `plBkEquiv : 𝒞₁^pl-bk ≌ 𝒞₂^pl-bk` | 本ファイル |
| 3 | `End` の降下(`Over.post` の圏同値性を使う) | 未 |

★★★段 1・2 を取ったところで**段 3 が要らなくなる可能性**を測った ——
`plBkOverToC P₁ A ⋙ Ψ` と `Over.post (plBkMap Ψ hPB) ⋙ plBkOverToC P₂ (Ψ A)` は
**定義的に等しい**(`plBkOverToC_comp_plBkMap`)。広い部分圏もスライスも
射を包んでいるだけなので、合成すると包みが打ち消し合う。
★これで段 3 は「圏同値に沿った自然変換の降下」という mathlib の練習問題になる。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-! ## ★段 1 —— `𝒞^pl-bk` の間の関手 -/

/-- ★★**pull-back 射を保つ関手は `𝒞^pl-bk` の間の関手を誘導する**。

★`WideSubcategory` は射を包んでいるだけなので、`map_id` も `map_comp` も
包みを外せば `Ψ` のそれに帰着する。 -/
def plBkMap (Ψ : C₁ ⥤ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f)) :
    PlBk P₁ ⥤ PlBk P₂ where
  obj A := ⟨Ψ.obj A.obj⟩
  map f := ⟨Ψ.map f.hom, hPB f.hom f.property⟩
  map_id _ := WideSubcategory.hom_ext _ (Ψ.map_id _)
  map_comp _ _ := WideSubcategory.hom_ext _ (Ψ.map_comp _ _)

@[simp]
theorem plBkMap_obj (Ψ : C₁ ⥤ C₂) (hPB) (A : PlBk P₁) :
    (plBkMap (P₁ := P₁) (P₂ := P₂) Ψ hPB).obj A = ⟨Ψ.obj A.obj⟩ := rfl

@[simp]
theorem plBkMap_map (Ψ : C₁ ⥤ C₂) (hPB) {A B : PlBk P₁} (f : A ⟶ B) :
    ((plBkMap (P₁ := P₁) (P₂ := P₂) Ψ hPB).map f).hom = Ψ.map f.hom := rfl

/-! ## ★★段 2 —— 圏同値の誘導 -/

/-- ★**同型は `𝒞^pl-bk` の同型に包める** —— 同型は pull-back 射だから。 -/
def plBkIsoOfIso {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
    {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀) {X Y : Cc} (θ : X ≅ Y) :
    (⟨X⟩ : PlBk P) ≅ (⟨Y⟩ : PlBk P) where
  hom := ⟨θ.hom, isPullBack_of_isIso P θ.hom⟩
  inv := ⟨θ.inv, isPullBack_of_isIso P θ.inv⟩
  hom_inv_id := WideSubcategory.hom_ext _ θ.hom_inv_id
  inv_hom_id := WideSubcategory.hom_ext _ θ.inv_hom_id

/-- ★★★**圏同値は `𝒞^pl-bk` の圏同値を誘導する**。

★unit / counit の成分は同型なので `plBkIsoOfIso` でそのまま包める。
★自然性も三角等式も、包みを外せば元のものに帰着する。 -/
def plBkEquiv (e : C₁ ≌ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (e.functor.map f))
    (hPB' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPullBack P₂ f → IsPullBack P₁ (e.inverse.map f)) :
    PlBk P₁ ≌ PlBk P₂ where
  functor := plBkMap e.functor hPB
  inverse := plBkMap e.inverse hPB'
  unitIso := NatIso.ofComponents (fun A => plBkIsoOfIso P₁ (e.unitIso.app A.obj))
    (fun f => WideSubcategory.hom_ext _ (e.unitIso.hom.naturality f.hom))
  counitIso := NatIso.ofComponents (fun B => plBkIsoOfIso P₂ (e.counitIso.app B.obj))
    (fun f => WideSubcategory.hom_ext _ (e.counitIso.hom.naturality f.hom))
  functor_unitIso_comp A := WideSubcategory.hom_ext _ (e.functor_unitIso_comp A.obj)

/-! ## ★★★段 3 の準備 —— 合成は定義的に一致する -/

/-- ★★★★★**`plBkOverToC` と `Over.post` の合成は定義的に等しい**。

★★広い部分圏もスライスも**射を包んでいるだけ**なので、
`Z ↦ Z.left.obj` と `f ↦ f.left.hom` で戻したところで包みが打ち消し合う。
★これにより段 3 は「圏同値 `Over.post (plBkMap Ψ hPB)` に沿った
**自然変換の降下**」だけになる。 -/
theorem plBkOverToC_comp_plBkMap (Ψ : C₁ ⥤ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f)) (A : C₁) :
    plBkOverToC P₁ A ⋙ Ψ
      = Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap Ψ hPB) ⋙ plBkOverToC P₂ (Ψ.obj A) :=
  rfl

/-- ★★★★**`Over.post` は圏同値になる** —— mathlib の instance をそのまま使う。 -/
theorem plBkOverPost_isEquivalence (e : C₁ ≌ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (e.functor.map f))
    (hPB' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPullBack P₂ f → IsPullBack P₁ (e.inverse.map f))
    (A : C₁) :
    (Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap e.functor hPB)).IsEquivalence := by
  haveI : (plBkMap (P₁ := P₁) (P₂ := P₂) e.functor hPB).IsEquivalence :=
    (plBkEquiv e hPB hPB').isEquivalence_functor
  infer_instance

/-! ## ★★★段 3 —— `End` の降下 -/

section Descent

variable (e : C₁ ≌ C₂)
  (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (e.functor.map f))
  (hPB' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPullBack P₂ f → IsPullBack P₁ (e.inverse.map f))
  (A : C₁)

/-- ★`Over.post (plBkMap …)` を圏同値として掴む。 -/
noncomputable def plBkOverEquiv :
    Over (⟨A⟩ : PlBk P₁) ≌ Over ((plBkMap e.functor hPB).obj (⟨A⟩ : PlBk P₁)) :=
  haveI := plBkOverPost_isEquivalence e hPB hPB' A
  (Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap e.functor hPB)).asEquivalence

/-- ★★★★★**`End` の降下** ——
`𝒞₁^pl-bk_A → 𝒞₁` の自己射を `𝒞₂^pl-bk_{ΨA} → 𝒞₂` の自己射へ送る。

★★`plBkOverToC_comp_plBkMap` で合成が**定義的に一致**するので、
`Ψ` で whisker してから、圏同値に沿った前合成の**充満忠実性**で引き戻すだけ。
★モノイド準同型性は mathlib の `Functor.mapEnd` と
`FullyFaithful.mulEquivEnd` が両方とも供給する。 -/
noncomputable def endDescend :
    End (plBkOverToC P₁ A) →* End (plBkOverToC P₂ (e.functor.obj A)) :=
  (((Equivalence.congrLeft (E := C₂)
        (plBkOverEquiv e hPB hPB' A)).fullyFaithfulInverse.mulEquivEnd
      (plBkOverToC P₂ (e.functor.obj A))).symm : _ ≃* _).toMonoidHom.comp
    (Functor.mapEnd (plBkOverToC P₁ A)
      ((Functor.whiskeringRight (Over (⟨A⟩ : PlBk P₁)) C₁ C₂).obj e.functor))

/-- ★★★★**降下の定義式** —— 前合成で戻すと元の whisker に一致する。 -/
theorem endDescend_whisker (η : End (plBkOverToC P₁ A)) :
    (Equivalence.congrLeft (E := C₂) (plBkOverEquiv e hPB hPB' A)).inverse.map
        (endDescend e hPB hPB' A η)
      = Functor.whiskerRight η e.functor :=
  (Equivalence.congrLeft (E := C₂)
    (plBkOverEquiv e hPB hPB' A)).fullyFaithfulInverse.map_preimage _

/-- ★★★★★**降下の成分**(`Over.post` の像の上)——
これが段 3 の要で、ここから `endBsIso` の保存も `plBkIdObj` での値も出る。 -/
theorem endDescend_app_post (η : End (plBkOverToC P₁ A)) (Z : Over (⟨A⟩ : PlBk P₁)) :
    NatTrans.app (endDescend e hPB hPB' A η)
        ((Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap e.functor hPB)).obj Z)
      = e.functor.map (NatTrans.app η Z) :=
  congrArg (fun t => NatTrans.app t Z) (endDescend_whisker e hPB hPB' A η)

/-- ★★★★**`plBkIdObj` は移る** —— `Ψ.map (𝟙 A) = 𝟙 (Ψ A)` だから。 -/
theorem plBkIdObj_map (Ψ : C₁ ⥤ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f)) (A₀ : C₁) :
    (Over.post (X := (⟨A₀⟩ : PlBk P₁)) (plBkMap Ψ hPB)).obj (plBkIdObj P₁ A₀)
      = plBkIdObj P₂ (Ψ.obj A₀) :=
  congrArg Over.mk (WideSubcategory.hom_ext _ (Ψ.map_id A₀))

/-- ★★★★★**恒等対象での値** —— 原文の条件が名指しする成分がそのまま `Ψ.map` で移る。 -/
theorem endDescend_app_id (η : End (plBkOverToC P₁ A)) :
    NatTrans.app (endDescend e hPB hPB' A η) (plBkIdObj P₂ (e.functor.obj A))
      = e.functor.map (NatTrans.app η (plBkIdObj P₁ A)) := by
  have key := endDescend_app_post e hPB hPB' A η (plBkIdObj P₁ A)
  have hobj := plBkIdObj_map (P₁ := P₁) (P₂ := P₂) e.functor hPB A
  have hswap : NatTrans.app (endDescend e hPB hPB' A η) (plBkIdObj P₂ (e.functor.obj A))
      = NatTrans.app (endDescend e hPB hPB' A η)
          ((Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap e.functor hPB)).obj
            (plBkIdObj P₁ A)) := by
    congr 1
    exact hobj.symm
  rw [hswap]
  exact key

/-- ★★★★★**降下は `endBsIso` を保つ**。

★★`Over.post` が**本質的全射**なので、任意の対象は `Kover` の像と同型。
自然性でその同型に沿って移すと、成分は
「同型 ⁻¹ ≫ `Ψ.map`(元の成分) ≫ 同型」の形になり、
base-isomorphism が合成で閉じることから出る。 -/
theorem endDescend_mem_endBsIso
    (hBI : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
    (η : End (plBkOverToC P₁ A)) (hη : η ∈ endBsIso P₁ A) :
    endDescend e hPB hPB' A η ∈ endBsIso P₂ (e.functor.obj A) := by
  haveI := plBkOverPost_isEquivalence e hPB hPB' A
  intro W
  set K := Over.post (X := (⟨A⟩ : PlBk P₁)) (plBkMap e.functor hPB) with hK
  set G₂ := plBkOverToC P₂ (e.functor.obj A) with hG₂
  set Z := K.objPreimage W with hZ
  set θ : K.obj Z ≅ W := K.objObjPreimageIso W with hθ
  have hnat := (endDescend e hPB hPB' A η).naturality θ.hom
  -- ★`Kover` の像の上では成分は `Ψ.map` そのもの。
  have hbase : IsBaseIsomorphism P₂
      (NatTrans.app (endDescend e hPB hPB' A η) (K.obj Z)) := by
    rw [endDescend_app_post e hPB hPB' A η Z]
    exact hBI _ (hη Z)
  haveI : IsIso (G₂.map θ.hom) := inferInstance
  have hcomp : IsBaseIsomorphism P₂
      (G₂.map θ.hom ≫ NatTrans.app (endDescend e hPB hPB' A η) W) := by
    rw [hnat]
    exact isBaseIsomorphism_comp P₂ hbase (isBaseIsomorphism_of_isIso P₂ _)
  haveI hm : IsIso (P₂.Base (G₂.map θ.hom)) := isBaseIsomorphism_of_isIso P₂ _
  haveI : IsIso (P₂.Base (G₂.map θ.hom)
      ≫ P₂.Base (NatTrans.app (endDescend e hPB hPB' A η) W)) := by
    rw [← P₂.Base_comp]
    exact hcomp
  exact IsIso.of_isIso_comp_left (P₂.Base (G₂.map θ.hom)) _

include hPB hPB' in
/-- ★★★★★**`OTriGenCond` は `Ψ` で移る** —— `iv-endbsiso` の結論。

原文 (FrdI p.66):
> By assertions (ii), (iii), it follows that Ψ preserves pre-steps, base-isomorphisms,

★★これが `Thm34Quasi.lean` の `otriEndCond_map` が仮定として受けていた 1 本である。 -/
theorem otriGenCond_map
    (hBI : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
    (α : End A) (h : OTriGenCond P₁ A α) :
    OTriGenCond P₂ (e.functor.obj A) (e.functor.map (α : A ⟶ A)) := by
  obtain ⟨f, hf⟩ := h
  refine ⟨((endDescend e hPB hPB' A).comp ((endBsIso P₁ A).subtype.comp f)).codRestrict
    (endBsIso P₂ (e.functor.obj A)) (fun x =>
      endDescend_mem_endBsIso e hPB hPB' A hBI _ (f x).2), ?_⟩
  show NatTrans.app (endDescend e hPB hPB' A
      ((f ⟨1, 1⟩ : endBsIso P₁ A) : End (plBkOverToC P₁ A)))
      (plBkIdObj P₂ (e.functor.obj A)) = e.functor.map (α : A ⟶ A)
  rw [endDescend_app_id, hf]
  rfl

end Descent

/-! ## ★★★★結び —— `Theorem 3.4, (iv)` の `𝒪^▷` 保存が仮定なしで出る -/

/-- ★★★★★**[FrdI] Theorem 3.4, (iv)** —— `Ψ` は `𝒪^▷(−)` を保つ(仮定を解消した版)。

原文 (FrdI p.66):
> Thus, we conclude that Ψ preserves the submonoids

★★`Thm34Quasi.lean` の `thm_3_4_iv_otri_map` が受けていた `hGen` を、
本ファイルの `otriGenCond_map` で埋めた。 -/
theorem thm_3_4_iv_otri_map' (e : C₁ ≌ C₂)
    (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (e.functor.map f))
    (hPB' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPullBack P₂ f → IsPullBack P₁ (e.inverse.map f))
    (hBI : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (e.functor.map f))
    (Fc₂ : FrobenioidCore P₂) (hslim₂ : IsFrobeniusSlim D₂)
    {C₀ : C₁} (hft : IsFrobeniusTrivial P₁ C₀) (hfn : IsFrobeniusNormalized P₁ C₀)
    (γ : End C₀) (h : γ ∈ OTri P₁ C₀) :
    e.functor.map (γ : C₀ ⟶ C₀) ∈ OTri P₂ (e.functor.obj C₀) :=
  thm_3_4_iv_otri_map Fc₂ hslim₂ e.functor hPS
    (fun A₀ α hα => otriGenCond_map e hPB hPB' A₀ hBI α hα) hft hfn γ h

def thm_3_4_iv_otri_map'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Ψ は 𝒪^▷ / 𝒪^× を保つ(仮定解消版)",
    sectionId := "frdi-thm-3-4" }


def plBkOverToC_comp_plBkMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — End(𝒞^pl-bk_A → 𝒞)^bs-iso の移送(段 1・2)",
    sectionId := "frdi-thm-3-4" }


/-! ## ★★★`Theorem 3.4, (iv)` の `𝒪^×` 側 —— unit-equivalence の保存

原文 (FrdI p.66):
> Thus, we conclude that Ψ preserves the submonoids

★★`𝒞^un-tr` は `Hom` を **unit-equivalence** で割った商である
(`UnTr.lean` の `unTrSetoid`)。したがって `Ψ^un-tr` を作るのに要るのは
「`Ψ` が unit-equivalence を保つ」の 1 本だけである。

★★★`IsUnitEquivalent` は `∃ Cc γ β δ, δ ∈ 𝒪^×(Cc) ∧ α₁ = γβ ∧ α₂ = γδβ` という
**純粋に圏論的な**定義なので、`𝒪^×` の保存があれば関手を当てるだけで移る。
★`𝒪^×(A) = 𝒪^▷(A) ∧ IsUnit` なので、`𝒪^▷` の保存(上の `thm_3_4_iv_otri_map'`)と
「関手は単元を保つ」(`Functor.mapEnd` が準同型)を合わせればよい。 -/

/-- ★★★★★**`Ψ` は unit-equivalence を保つ**。

★`𝒪^▷` の保存を仮定として受け、`IsUnit` の側は `Functor.mapEnd` の
準同型性から `IsUnit.map` で出す。 -/
theorem isUnitEquivalent_map (Ψ : C₁ ⥤ C₂)
    (hOTri : ∀ (X : C₁) (δ : End X), δ ∈ OTri P₁ X → Ψ.map (δ : X ⟶ X) ∈ OTri P₂ (Ψ.obj X))
    {A B : C₁} {α₁ α₂ : A ⟶ B} (h : IsUnitEquivalent P₁ α₁ α₂) :
    IsUnitEquivalent P₂ (Ψ.map α₁) (Ψ.map α₂) := by
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ := h
  refine ⟨Ψ.obj Cc, Ψ.map γ, Ψ.map β, Ψ.map (δ : Cc ⟶ Cc), ⟨hOTri Cc δ hδ.1, ?_⟩, ?_, ?_⟩
  · exact hδ.2.map (Functor.mapEnd Cc Ψ)
  · rw [h₁, Ψ.map_comp]
  · rw [h₂, Ψ.map_comp, Ψ.map_comp]

def isUnitEquivalent_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Ψ は unit-equivalence を保つ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★`Theorem 3.4, (iv)` の `Ψ^un-tr`

原文 (FrdI p.66):
> assertion (iv) then follows formally from the deﬁnition of

★★`𝒞^un-tr` は `Hom` を `unTrSetoid`(= `𝔽_Φ` へ同じ射に写る)で割った商である。
したがって `Ψ^un-tr` を作るのに要るのは
「`Ψ` がこの同値関係を保つ」の 1 本だけ。

★★★その 1 本は **3 段で出る**:
1. `𝔽_Φ` で一致 ⟹ unit-equivalent(`Proposition 3.3, (ii)` の `mpr`)
2. unit-equivalent は `Ψ` で移る(`isUnitEquivalent_map`)
3. unit-equivalent ⟹ `𝔽_Φ` で一致(`prop_3_3_ii_toElem`)

★1 と 3 は在庫、2 が本ファイルで取ったもの。 -/

/-- ★★★★★**`Ψ` は `𝒞^un-tr` の同値関係を保つ** —— `Proposition 3.3, (ii)` を往復する。 -/
theorem toElem_map_congr_of_congr (Ψ : C₁ ⥤ C₂) (Fc₁ : FrobenioidCore P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (hOTri : ∀ (X : C₁) (δ : End X), δ ∈ OTri P₁ X → Ψ.map (δ : X ⟶ X) ∈ OTri P₂ (Ψ.obj X))
    {A B : C₁} (α₁ α₂ : A ⟶ B) (h : P₁.toElem.map α₁ = P₁.toElem.map α₂) :
    P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂) :=
  prop_3_3_ii_toElem P₂ (isUnitEquivalent_map Ψ hOTri
    ((prop_3_3_ii P₁ Fc₁ hiso₁ α₁ α₂).mpr
      ⟨congrArg ElemFrobCat.Hom.deg h, congrArg ElemFrobCat.Hom.div h,
        congrArg ElemFrobCat.Hom.base h⟩))

/-- ★★★★★**[FrdI] Theorem 3.4, (iv)** —— `Ψ^un-tr : 𝒞₁^un-tr ⥤ 𝒞₂^un-tr`。

★対象は `Ψ^istr`(在庫、(i) で実装済)、射は商へ落とすだけ。 -/
def psiUnTr (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence] (h₁ : IsOfQuasiIsotropicType C₁ P₁)
    (h₂ : IsOfQuasiIsotropicType C₂ P₂)
    (hUE : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
      P₁.toElem.map α₁ = P₁.toElem.map α₂ →
        P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂)) :
    UnTr P₁ ⥤ UnTr P₂ where
  obj A := (psiIstr Ψ P₁ P₂ h₁ h₂).obj A
  map {_ _} f := Quotient.liftOn f (fun α => toHomUnTr P₂ (Ψ.map α))
    (fun a b hab => (toHomUnTr_eq_iff P₂ _ _).mpr (hUE a b hab))
  map_id A := by
    show toHomUnTr P₂ (Ψ.map (𝟙 (A : Istr P₁).obj)) = toHomUnTr P₂ (𝟙 _)
    rw [Ψ.map_id]
  map_comp {_ _ _} f g := by
    refine Quotient.inductionOn₂ f g (fun α β => ?_)
    show toHomUnTr P₂ (Ψ.map (α ≫ β)) = toHomUnTr P₂ (Ψ.map α ≫ Ψ.map β)
    rw [Ψ.map_comp]

/-- ★★★★**1-可換図式** —— `𝒞^istr → 𝒞^un-tr` の四角形が**厳密に**可換。

★商へ落とす前と後で `Ψ` が同じ射を与えるので `rfl` になる。 -/
theorem psiUnTr_square (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence] (h₁ : IsOfQuasiIsotropicType C₁ P₁)
    (h₂ : IsOfQuasiIsotropicType C₂ P₂) (hUE) :
    psiIstr Ψ P₁ P₂ h₁ h₂ ⋙ istrToUnTr P₂ = istrToUnTr P₁ ⋙ psiUnTr Ψ h₁ h₂ hUE :=
  rfl

def psiUnTr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Ψ^un-tr の構成と 1-可換図式",
    sectionId := "frdi-thm-3-4" }

/-! ### ★★`Ψ^un-tr` の圏同値性と rigidity

原文 (FrdI p.66):
> are of unit-trivial type, the asserted rigidity follows formally from Proposition 1.13,

★rigidity は `Proposition 1.13, (ii)`(在庫 `prop_1_13_ii_global`)に
本セッションの `isRigidFunctor_of_equivalence_comp` を被せるだけ。
★そのために `Ψ^un-tr` が圏同値であることが要る —— 3 性質を順に取る。 -/

section UnTrEquiv

variable (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence] (h₁ : IsOfQuasiIsotropicType C₁ P₁)
  (h₂ : IsOfQuasiIsotropicType C₂ P₂)
  (hUE : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
    P₁.toElem.map α₁ = P₁.toElem.map α₂ →
      P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂))

/-- ★★★**忠実** —— 同値関係の**反射**(逆向き)が要る。 -/
theorem psiUnTr_faithful
    (hUE' : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
      P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂) →
        P₁.toElem.map α₁ = P₁.toElem.map α₂) :
    (psiUnTr Ψ h₁ h₂ hUE).Faithful where
  map_injective {_ _ f g} h := by
    revert h
    refine Quotient.inductionOn₂ f g (fun α β h => ?_)
    exact (toHomUnTr_eq_iff P₁ _ _).mpr (hUE' α β ((toHomUnTr_eq_iff P₂ _ _).mp h))

/-- ★★**充満** —— `Ψ` の充満性から代表元を取るだけ。 -/
theorem psiUnTr_full : (psiUnTr Ψ h₁ h₂ hUE).Full where
  map_surjective {_ _} g := by
    refine Quotient.inductionOn g (fun β => ?_)
    obtain ⟨α, hα⟩ := Ψ.map_surjective β
    refine ⟨toHomUnTr P₁ α, ?_⟩
    show toHomUnTr P₂ (Ψ.map α) = toHomUnTr P₂ β
    rw [hα]
    rfl

/-- ★★**本質的全射** —— 対象は `Ψ^istr` と同じなので在庫から。 -/
theorem psiUnTr_essSurj : (psiUnTr Ψ h₁ h₂ hUE).EssSurj where
  mem_essImage Z := by
    haveI := psiIstr_essSurj Ψ P₁ P₂ h₁ h₂
    obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := psiIstr Ψ P₁ P₂ h₁ h₂)
      (show Istr P₂ from Z)
    exact ⟨A, ⟨(istrToUnTr P₂).mapIso ε⟩⟩

/-- ★★★★**`Ψ^un-tr` は圏同値**。 -/
theorem psiUnTr_isEquivalence
    (hUE' : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
      P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂) →
        P₁.toElem.map α₁ = P₁.toElem.map α₂) :
    (psiUnTr Ψ h₁ h₂ hUE).IsEquivalence where
  faithful := psiUnTr_faithful Ψ h₁ h₂ hUE hUE'
  full := psiUnTr_full Ψ h₁ h₂ hUE
  essSurj := psiUnTr_essSurj Ψ h₁ h₂ hUE

/-- ★★★★★**[FrdI] Theorem 3.4, (iv)** の rigidity。

★★`𝒞₂^un-tr` は Frobenioid(在庫 `unTr_frobenioidCore`)なので
`Proposition 1.13, (ii)` が当たり、そこへ
`isRigidFunctor_of_equivalence_comp` を被せる。 -/
theorem psiUnTr_rigid (Fc₂ : FrobenioidCore P₂) (hslim : IsSlimCat D₂)
    (hUE' : ∀ {A B : C₁} (α₁ α₂ : A ⟶ B),
      P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂) →
        P₁.toElem.map α₁ = P₁.toElem.map α₂) :
    IsRigidFunctor (psiUnTr Ψ h₁ h₂ hUE ⋙ (unTrPre P₂ Fc₂).toElem) := by
  haveI := psiUnTr_isEquivalence Ψ h₁ h₂ hUE hUE'
  exact isRigidFunctor_of_equivalence_comp (psiUnTr Ψ h₁ h₂ hUE).asEquivalence _
    (prop_1_13_ii_global (unTrPre P₂ Fc₂) (unTr_frobenioidCore P₂ Fc₂) hslim)

end UnTrEquiv

def psiUnTr_rigid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Ψ^un-tr の rigidity",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**1-一意性** —— 図式を可換にする関手は 1 つしかない。

★★`𝒞^istr → 𝒞^un-tr` は**対象を変えず**、射は商への全射なので、
四角形が可換なら対象でも射でも一致する。
★これで原文の「1-unique 1-commutative diagram」が揃う。 -/
theorem psiUnTr_unique (F G : UnTr P₁ ⥤ UnTr P₂)
    (hobj : ∀ A : UnTr P₁, F.obj A = G.obj A)
    (hmap : ∀ (A B : Istr P₁) (α : A.obj ⟶ B.obj),
      F.map (toHomUnTr P₁ α)
        = eqToHom (hobj A) ≫ G.map (toHomUnTr P₁ α) ≫ eqToHom (hobj B).symm) :
    F = G := by
  refine CategoryTheory.Functor.ext hobj ?_
  intro A B f
  refine Quotient.inductionOn f (fun α => ?_)
  exact hmap A B α

def psiUnTr_unique.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Ψ^un-tr の 1-一意性",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★`Theorem 3.4, (iv)` —— 条のうち保存・構成・rigidity が揃った

原文 (FrdI p.62):
> (iv) Suppose that: (a) C1, C2 are of standard type; (b) if C1, C2 are of

★★**原文の仮定を測り直した**(2026-08-19、ゲート NG から) ——
(iv) の仮定は「quasi-isotropic 型」ではなく
**(a) standard 型、(b) group-like 型なら `Ψ` とその quasi-inverse が base-isomorphism を保つ**
である。★下の実装はこれより弱い仮定(圏同値＋各保存を個別に受ける)で書いてあるので、
原典の仮定からこれらを導く接続は `iv-src` の作業になる。

★実装した部品:

| 主張 | 宣言 | 場所 |
|---|---|---|
| `𝒪^▷` の圏論的特徴づけ | `OTriEndCond` / `oTri_of_otriEndCond` / `otriEndCond_of_oTri` | `Thm34Quasi.lean` |
| `End(𝒞^pl-bk_A → 𝒞)^bs-iso` の移送 | `otriGenCond_map` | 本ファイル |
| `𝒪^▷` の保存 | `thm_3_4_iv_otri_map'` | 本ファイル |
| `𝒪^×` の保存(unit-equivalence) | `isUnitEquivalent_map` | 本ファイル |
| `Ψ^un-tr` の構成 | `psiUnTr` | 本ファイル |
| 1-可換図式 | `psiUnTr_square`(**厳密な等式**) | 本ファイル |
| 1-一意性 | `psiUnTr_unique` | 本ファイル |
| rigidity | `psiUnTr_rigid` | 本ファイル |

★★★**残す**: `Ψ_{N≥1}` が恒等自己同型であること。
原文はこれを (iii) と `Proposition 1.10, (vi)` ＋ Frobenius-compact 性から出す。
群的型でない場合は (iii) から**形式的に**従い、群的型の場合は
`Frobenius-compact` な対象を取って `p₁ = p₂` を示す議論になる。
★本セッションではその材料(`degFr_order_preserve` / `pnatMulEquivOfPrimeEquiv_eq_self`)を
(iii) 側で取ってあるので、接続は残作業である。 -/
def thm_3_4_iv_part.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (iv) — 𝒪^▷ / 𝒪^× の保存・Ψ^un-tr・1-一意性・rigidity",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★`Theorem 3.4, (iv)` の `Ψ_{N≥1} = id` —— 群的型の場合

原文 (FrdI p.67):
> by multiplication by p1/p2. Since A2 is Frobenius-

★★原文の要点は「`α` が `𝒪^×(A₂)^pf` に `p₁/p₂` 倍で作用する。
`A₂` が Frobenius-compact だから `p₁ = p₂`」である。

★★★その最後の一歩は **`Definition 1.2, (iv)` の第 3 節そのもの** ——
「`c/d` 倍で作用するなら `1` 倍で作用する」。
そこへ第 2 節(**無限位数の単元が存在する**)を当てると `c = d` が出る。
★下はその純代数の部分で、Frobenioid の構造は一切使わない。 -/

/-- ★★★★★**Frobenius-compact なら「`c/d` 倍作用」は `c = d` を強制する**。

★★`Definition 1.2, (iv)` の第 3 節で「`1` 倍作用」に落とし、
第 2 節の**無限位数の単元**で指数を比較する。 -/
theorem pnat_eq_of_frobeniusCompact {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A : Cc} (hFC : IsFrobeniusCompact P A) (θ : A ≅ A) (c d : ℕ+)
    (hact : ∀ u : End A, u ∈ OTimes P A → ∃ k : ℕ+,
      ((endConj θ u) ^ ((d : ℕ) * (k : ℕ)) : End A)
        = (u ^ ((c : ℕ) * (k : ℕ)) : End A)) :
    c = d := by
  obtain ⟨_hcomm, ⟨u, hu, hord⟩, hkey⟩ := hFC
  obtain ⟨k₁, h₁⟩ := hact u hu
  obtain ⟨k₂, h₂⟩ := hkey θ c d hact u hu
  -- ★経路 1: `hact` を `k₂` 乗する。
  have e₁ : ((endConj θ u) ^ ((d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A)
      = (u ^ ((c : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A) := by
    rw [pow_mul, h₁, ← pow_mul]
  -- ★経路 2: 第 3 節が与える「1 倍作用」を `d * k₁` 乗する。
  have e₂ : ((endConj θ u) ^ ((d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A)
      = (u ^ ((d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A) := by
    have hcomm : (d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ) = (k₂ : ℕ) * ((d : ℕ) * (k₁ : ℕ)) := by
      ring
    rw [hcomm, pow_mul, h₂, ← pow_mul]
  have ekey : (u ^ ((c : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A)
      = (u ^ ((d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) : End A) := e₁.symm.trans e₂
  -- ★`u` は単元かつ無限位数なので、指数がそのまま一致する。
  have hunit : IsUnit u := hu.2
  have hpow : ∀ m n : ℕ, (u ^ m : End A) = (u ^ n : End A) → m ≤ n → m = n := by
    intro m n h hmn
    rcases eq_or_lt_of_le hmn with heq | hlt
    · exact heq
    · exfalso
      have hstep : (u ^ m : End A) * (u ^ (n - m) : End A) = (u ^ m : End A) * 1 := by
        rw [mul_one, ← pow_add, Nat.add_sub_cancel' hmn]
        exact h.symm
      have hcancel : (u ^ (n - m) : End A) = 1 := (hunit.pow m).mul_left_cancel hstep
      exact hord ⟨n - m, Nat.sub_pos_of_lt hlt⟩ hcancel
  have hK : 0 < (k₁ : ℕ) * (k₂ : ℕ) := Nat.mul_pos k₁.pos k₂.pos
  refine PNat.coe_injective ?_
  rcases le_total ((c : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) ((d : ℕ) * (k₁ : ℕ) * (k₂ : ℕ)) with hle | hle
  · have hEq := hpow _ _ ekey hle
    refine Nat.eq_of_mul_eq_mul_right hK ?_
    rw [← mul_assoc, ← mul_assoc]
    exact hEq
  · have hEq := hpow _ _ ekey.symm hle
    refine (Nat.eq_of_mul_eq_mul_right hK ?_).symm
    rw [← mul_assoc, ← mul_assoc]
    exact hEq

def pnat_eq_of_frobeniusCompact.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — Frobenius-compact から p₁ = p₂",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★**`α` は `p₁/p₂` 倍で作用する** —— 原文の計算そのもの。

原文 (FrdI p.67):
> hence [by the total epimorphicity of C2]

★★筋: `φ = α ∘ ψ` と書けるとき、
`φ ∘ u = α ∘ (ψ ∘ u) = α ∘ (u^{p₂} ∘ ψ) = (α ∘ u^{p₂} ∘ α⁻¹) ∘ φ`。
これを `u^{p₁} ∘ φ = φ ∘ u` と突き合わせ、
**`𝒞` が totally epimorphic**(= すべての射が epi)なので `φ` を消去する。 -/
theorem otimes_conj_pow_of_frobNorm {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A : Cc} (hfn : IsFrobeniusNormalized P A) (a a' : End A) (ha'a : a' * a = 1)
    (ψ : End A) (hψbid : IsBaseIdentity P ψ) (p₁ : ℕ)
    (hrel : ∀ u : End A, u ∈ OTimes P A → (u ^ p₁ : End A) * (a * ψ) = (a * ψ) * u)
    (u : End A) (hu : u ∈ OTimes P A) :
    a * (u ^ (P.degFr (ψ : A ⟶ A) : ℕ) : End A) * a' = (u ^ p₁ : End A) := by
  have hfrob : (u ^ (P.degFr (ψ : A ⟶ A) : ℕ) : End A) * ψ = ψ * u :=
    frobNorm_mul P hfn ψ hψbid ⟨u, hu.1⟩
  -- ★`φ ∘ u` を 2 通りに書く。
  have hstep : (a * ψ) * u
      = (a * (u ^ (P.degFr (ψ : A ⟶ A) : ℕ) : End A) * a') * (a * ψ) := by
    have hR : (a * (u ^ (P.degFr (ψ : A ⟶ A) : ℕ) : End A) * a') * (a * ψ)
        = a * ((u ^ (P.degFr (ψ : A ⟶ A) : ℕ) : End A) * ((a' * a) * ψ)) := by
      simp only [mul_assoc]
    rw [hR, ha'a, one_mul, hfrob, mul_assoc]
  -- ★`φ` を epi として消去する(`𝒞` は totally epimorphic)。
  have hthis := hrel u hu
  rw [hstep] at hthis
  haveI : Epi ((a * ψ : End A) : A ⟶ A) := P.totEpiC _ _ _
  exact ((cancel_epi ((a * ψ : End A) : A ⟶ A)).mp hthis).symm

def otimes_conj_pow_of_frobNorm.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — α は p₁/p₂ 倍で作用する",
    sectionId := "frdi-thm-3-4" }

/-- ★★**橋渡し** —— `endConj` は `End` の積で `α.hom * x * α.inv` と書ける。

★これで `otimes_conj_pow_of_frobNorm`(積の形)と
`pnat_eq_of_frobeniusCompact`(`endConj` の形)が繋がる。 -/
theorem endConj_eq_mul {Cc : Type u2} [Category.{v2} Cc] {A : Cc} (α : A ≅ A) (x a a' : End A)
    (ha : a = α.hom) (ha' : a' = α.inv) : endConj α x = a * x * a' := by
  subst ha
  subst ha'
  rfl

def endConj_eq_mul.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — endConj は End の積で書ける",
    sectionId := "frdi-def-1-2-iv" }

/-- ★★★★★**両端を繋いだ形** —— `p₁ = degFr ψ`。

原文 (FrdI p.67):
> compact, we thus conclude that p1 = p2. This completes the proof of assertion

★★`otimes_conj_pow_of_frobNorm`(`p₁/p₂` 倍作用)を
`pnat_eq_of_frobeniusCompact`(compact なら等しい)へ `k = 1` で流し込むだけ。 -/
theorem degFr_eq_of_frobeniusCompact_of_frobNorm {Dd : Type u} [Category.{v} Dd]
    {Cc : Type u2} [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A : Cc} (hFC : IsFrobeniusCompact P A) (hfn : IsFrobeniusNormalized P A)
    (α : A ≅ A) (a a' : End A) (ha : a = α.hom) (ha' : a' = α.inv)
    (ψ : End A) (hψbid : IsBaseIdentity P ψ) (p₁ : ℕ+)
    (hrel : ∀ u : End A, u ∈ OTimes P A →
      (u ^ (p₁ : ℕ) : End A) * (a * ψ) = (a * ψ) * u) :
    P.degFr (ψ : A ⟶ A) = p₁ := by
  have ha'a : a' * a = 1 := by
    subst ha
    subst ha'
    exact α.hom_inv_id
  refine (pnat_eq_of_frobeniusCompact P hFC α p₁ (P.degFr (ψ : A ⟶ A)) ?_).symm
  intro u hu
  refine ⟨1, ?_⟩
  rw [PNat.one_coe, mul_one, mul_one, ← map_pow, endConj_eq_mul α _ a a' ha ha']
  exact otimes_conj_pow_of_frobNorm P hfn a a' ha'a ψ hψbid (p₁ : ℕ) hrel u hu

def degFr_eq_of_frobeniusCompact_of_frobNorm.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 67,
    item := "Theorem 3.4, (iv) — Frobenius-compact から p₁ = p₂(両端の接続)",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**`hrel` を `𝒞₁` から移送する** ——
`Ψ` が `𝒪^×` を**全射に**移すことから、`𝒞₂` 側の関係式が出る。

★★`Ψ` の充満性で代表元を取り、`𝒪^×` の**反射**(quasi-inverse 側の保存)で
その代表元が `𝒪^×(A₁)` に入ることを言う。あとは `Ψ` を関係式に当てるだけ。 -/
theorem frobNorm_rel_map (Ψ : C₁ ⥤ C₂) {A₁ : C₁} (hfn₁ : IsFrobeniusNormalized P₁ A₁)
    (φ₁ : End A₁) (hφbid : IsBaseIdentity P₁ φ₁) (p₁ : ℕ+)
    (hd₁ : P₁.degFr (φ₁ : A₁ ⟶ A₁) = p₁)
    (f₂ : End (Ψ.obj A₁)) (hf₂ : f₂ = Functor.mapEnd A₁ Ψ φ₁)
    (hrefl : ∀ u : End (Ψ.obj A₁), u ∈ OTimes P₂ (Ψ.obj A₁) →
      ∃ u₁ : End A₁, u₁ ∈ OTimes P₁ A₁ ∧ Functor.mapEnd A₁ Ψ u₁ = u)
    (u : End (Ψ.obj A₁)) (hu : u ∈ OTimes P₂ (Ψ.obj A₁)) :
    (u ^ (p₁ : ℕ) : End (Ψ.obj A₁)) * f₂ = f₂ * u := by
  obtain ⟨u₁, hu₁, hmap⟩ := hrefl u hu
  have h₁ : (u₁ ^ (p₁ : ℕ) : End A₁) * φ₁ = φ₁ * u₁ := by
    have h := frobNorm_mul P₁ hfn₁ φ₁ hφbid ⟨u₁, hu₁.1⟩
    rwa [hd₁] at h
  have hmapped := congrArg (Functor.mapEnd A₁ Ψ) h₁
  rw [map_mul, map_mul, map_pow, hmap, ← hf₂] at hmapped
  exact hmapped

/-- ★★★★★**[FrdI] Theorem 3.4, (iv)** の群的型の結論 —— `p₁ = p₂`。

原文 (FrdI p.67):
> compact, we thus conclude that p1 = p2. This completes the proof of assertion

★`degFr_eq_of_frobeniusCompact_of_frobNorm` に次数の仮定を当てるだけ。 -/
theorem admissible_deg_eq_of_frobeniusCompact {Dd : Type u} [Category.{v} Dd]
    {Cc : Type u2} [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A : Cc} (hFC : IsFrobeniusCompact P A) (hfn : IsFrobeniusNormalized P A)
    (p₁ p₂ : ℕ+) (α : A ≅ A) (a a' : End A) (ha : a = α.hom) (ha' : a' = α.inv)
    (ψ : End A) (hψbid : IsBaseIdentity P ψ) (hdψ : P.degFr (ψ : A ⟶ A) = p₂)
    (hrel : ∀ u : End A, u ∈ OTimes P A →
      (u ^ (p₁ : ℕ) : End A) * (a * ψ) = (a * ψ) * u) :
    p₂ = p₁ := by
  rw [← hdψ]
  exact degFr_eq_of_frobeniusCompact_of_frobNorm P hFC hfn α a a' ha ha' ψ hψbid p₁ hrel

def admissible_deg_eq_of_frobeniusCompact.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 67,
    item := "Theorem 3.4, (iv) — 群的型でも p₁ = p₂",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★`Theorem 3.4, (v)` の第 1 主張 —— 1-可換図式から出る

原文 (FrdI p.69):
> (v), which is clearly 1-unique [cf. Definition 1.3, (i), (a), (b), (c)]. Finally, the

★★★**依存の向きを測り直した**(2026-08-19、原典 p.67-69) ——
(v) の第 1 主張(base-identity 自己射と base-equivalent な対の保存)は
**独立に証明されていない**。原典は `Q_i` / `R_i` を作って `Ψ_Base` を得た**後**、
1-可換図式からまとめて出している。
★したがって `v-baseid` は `v-psibase` に**依存する**(以前の台帳は逆に書いていた)。

★★下は「1-可換図式があれば第 1 主張が出る」部分で、
`Ψ_Base` の構成(`v-psibase`)とは独立に取れる。 -/

/-- ★★★★★**[FrdI] Theorem 3.4, (v) の第 1 主張(base-identity)** ——
1-可換図式 `𝒞₁ →Ψ 𝒞₂` / `𝒟₁ →Ψ_Base 𝒟₂` があれば、
`Ψ` は base-identity 自己射を保つ。

★底の四角形の自然性を `φ` に当て、`Base₁ φ = 𝟙` を入れて同型を消去するだけ。 -/
theorem isBaseIdentity_map_of_baseSquare (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj)
    {A : C₁} (φ : End A) (h : IsBaseIdentity P₁ φ) :
    IsBaseIdentity P₂ (Ψ.map (φ : A ⟶ A)) := by
  have hnat : ΨBase.map (P₁.Base (φ : A ⟶ A)) ≫ sq.hom.app A
      = sq.hom.app A ≫ P₂.Base (Ψ.map (φ : A ⟶ A)) := sq.hom.naturality (φ : A ⟶ A)
  have hb : ΨBase.map (P₁.Base (φ : A ⟶ A)) = 𝟙 (ΨBase.obj ((P₁.toElem.obj A).base)) := by
    rw [show P₁.Base (φ : A ⟶ A) = P₁.Base (𝟙 A) from h, P₁.Base_id]
    exact ΨBase.map_id _
  rw [hb, Category.id_comp] at hnat
  show P₂.Base (Ψ.map (φ : A ⟶ A)) = P₂.Base (𝟙 (Ψ.obj A))
  rw [P₂.Base_id]
  refine (cancel_epi (sq.hom.app A)).mp ?_
  refine hnat.symm.trans ?_
  exact (Category.comp_id (sq.hom.app A)).symm

/-- ★★★★★**[FrdI] Theorem 3.4, (v) の第 1 主張(base-equivalent な対)** ——
同じ四角形から、共対象射の base-equivalence も保たれる。 -/
theorem baseEquivalent_map_of_baseSquare (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj)
    {A B : C₁} (φ ψ : A ⟶ B) (h : BaseEquivalent P₁ φ ψ) :
    BaseEquivalent P₂ (Ψ.map φ) (Ψ.map ψ) := by
  have hφ : ΨBase.map (P₁.Base φ) ≫ sq.hom.app B
      = sq.hom.app A ≫ P₂.Base (Ψ.map φ) := sq.hom.naturality φ
  have hψ : ΨBase.map (P₁.Base ψ) ≫ sq.hom.app B
      = sq.hom.app A ≫ P₂.Base (Ψ.map ψ) := sq.hom.naturality ψ
  rw [show P₁.Base φ = P₁.Base ψ from h] at hφ
  show P₂.Base (Ψ.map φ) = P₂.Base (Ψ.map ψ)
  exact (cancel_epi (sq.hom.app A)).mp (hφ.symm.trans hψ)

def isBaseIdentity_map_of_baseSquare.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (v) — base-identity 自己射と base-equivalent な対の保存",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★`Theorem 3.4, (v)` の `Ψ_Base` —— **局所化として読み替える**

原文 (FrdI p.68):
> sition A.2, hence, in particular, is an equivalence of categories. Thus, since the nat-

★★★**測って分かった読み替え**(2026-08-19)。
原典は `Q_i`(スライスの 2-圏の coarsification)を base-isomorphism が誘導する
「新しい関手」で広げて `R_i` を作り、**`R_i ≌ 𝒟_i`** を示す。

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★★この分解(pre-step 2 本と pull-back 1 本、うち 1 本は逆向き)は、
**`𝒟` が `𝒞` の base-isomorphism による局所化である**と言っている。
★したがって `Ψ_Base` は **2-圏を組まずに、局所化の普遍性**で作れる。

★★★これが本セッション最大の設計判断である ——
`Definition A.1` の 2-圏を Lean で組む代わりに、
mathlib の `CategoryTheory.Localization` に載せる。
★付録側で確かめた「coarsification の hom は元の hom と一対一」
(`coarseHomEquiv`)が、この読み替えが正しいことの裏づけになっている。 -/

section PsiBase

variable (Ψ : C₁ ⥤ C₂)

/-- ★★**`Ψ` が base-isomorphism を保つ**ことの、局所化の言葉での言い換え。 -/
theorem isInvertedBy_of_baseIso_map
    (hBI : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (Ψ.map f)) :
    (baseIsoProp P₁).IsInvertedBy (Ψ ⋙ P₂.proj) := by
  intro X Y f hf
  exact hBI f hf

/-- ★★★★★**[FrdI] Theorem 3.4, (v)** の `Ψ_Base : 𝒟₁ ⥤ 𝒟₂`。

★`𝒟₁` が `𝒞₁` の base-isomorphism による局所化であるとき、
`Ψ ⋙ P₂.proj` が base-isomorphism を反転させることから普遍性で得られる。 -/
noncomputable def psiBaseLoc [P₁.proj.IsLocalization (baseIsoProp P₁)]
    (hinv : (baseIsoProp P₁).IsInvertedBy (Ψ ⋙ P₂.proj)) : D₁ ⥤ D₂ :=
  Localization.lift (Ψ ⋙ P₂.proj) hinv P₁.proj

/-- ★★★★★**1-可換図式** —— 局所化の `Lifting.iso` そのもの。

★これがちょうど `v-baseid` の `isBaseIdentity_map_of_baseSquare` が受ける仮定である。 -/
noncomputable def psiBaseLocSquare [P₁.proj.IsLocalization (baseIsoProp P₁)]
    (hinv : (baseIsoProp P₁).IsInvertedBy (Ψ ⋙ P₂.proj)) :
    P₁.proj ⋙ psiBaseLoc Ψ hinv ≅ Ψ ⋙ P₂.proj :=
  haveI : Localization.Lifting P₁.proj (baseIsoProp P₁) (Ψ ⋙ P₂.proj) (psiBaseLoc Ψ hinv) :=
    Localization.liftingLift (Ψ ⋙ P₂.proj) hinv P₁.proj
  Localization.Lifting.iso P₁.proj (baseIsoProp P₁) (Ψ ⋙ P₂.proj) (psiBaseLoc Ψ hinv)

/-- ★★★★**1-一意性** —— 図式を可換にする `𝒟₁ ⥤ 𝒟₂` は同型を除いて一意。

★局所化に沿った前合成が**充満忠実**であることから出る。 -/
noncomputable def psiBaseLocUniq [P₁.proj.IsLocalization (baseIsoProp P₁)]
    (F G : D₁ ⥤ D₂) (e : P₁.proj ⋙ F ≅ P₁.proj ⋙ G) : F ≅ G :=
  haveI : Localization.Lifting P₁.proj (baseIsoProp P₁) (P₁.proj ⋙ F) F := ⟨Iso.refl _⟩
  haveI : Localization.Lifting P₁.proj (baseIsoProp P₁) (P₁.proj ⋙ G) G := ⟨Iso.refl _⟩
  Localization.liftNatIso P₁.proj (baseIsoProp P₁) (P₁.proj ⋙ F) (P₁.proj ⋙ G) F G e

end PsiBase

def psiBaseLoc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (v) — Ψ_Base の構成と 1-可換図式・1-一意性",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★**[FrdI] Theorem 3.4, (v)** の rigidity。

原文 (FrdI p.69):
> asserted rigidity follows formally from Proposition 1.13, (i). This completes the

★在庫の `prop_1_13_i_global : IsRigidFunctor P.proj` に
`isRigidFunctor_of_equivalence_comp`(`iii-rigid` で作ったもの)を被せるだけ。
★★**同じ道具が 3 個目の節点で効いた**((iii)・(iv)・(v))。 -/
theorem thm_3_4_v_rigid (Fc₂ : FrobenioidCore P₂) (hslim₂ : IsSlimCat D₂) (e : C₁ ≌ C₂) :
    IsRigidFunctor (e.functor ⋙ P₂.proj) :=
  isRigidFunctor_of_equivalence_comp e _ (prop_1_13_i_global P₂ Fc₂ hslim₂)

def thm_3_4_v_rigid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 69,
    item := "Theorem 3.4, (v) — 合成関手の rigidity",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★`Theorem 3.4, (iv)` —— 条として揃った

原文 (FrdI p.62):
> (iv) Suppose that: (a) C1, C2 are of standard type; (b) if C1, C2 are of

★条の 4 主張がすべて実装された:

| 主張 | 宣言 |
|---|---|
| `Ψ` は `𝒪^▷(−)` を保つ | `thm_3_4_iv_otri_map'` |
| `Ψ` は `𝒪^×(−)` を保つ(unit-equivalence) | `isUnitEquivalent_map` |
| `Ψ^un-tr` が 1-一意に存在し図式が 1-可換 | `psiUnTr` / `psiUnTr_square` / `psiUnTr_unique` |
| `𝒞^un-tr` は unit-trivial 型なので rigid | `psiUnTr_rigid` |
| `Ψ_{N≥1}` は恒等 | `psiN_eq_id_of_orderPreserve`(非 group-like)／`admissible_deg_eq_of_frobeniusCompact`(group-like) |

★★`𝒪^▷` の圏論的特徴づけ(`OTriEndCond`)と
`End(𝒞^pl-bk_A → 𝒞)^bs-iso` の移送(`otriGenCond_map`)が本条の芯だった。 -/
def thm_3_4_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62, item := "Theorem 3.4, (iv)",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★`v-loc` の第 1 条件 —— `P.proj` は base-isomorphism を反転させる

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★★mathlib の `Functor.IsLocalization` は 2 条件からなる:
1. `inverts : W.IsInvertedBy L` —— **下でそのまま出る**(定義そのもの)
2. `isEquivalence : IsEquivalence (Localization.Construction.lift L inverts)`
   —— これが原典の `R_i ≌ 𝒟_i` で、`v-loc` の本体。 -/

/-- ★★**`P.proj` は base-isomorphism を反転させる** —— 定義そのもの。

★`IsBaseIsomorphism P f := IsIso (P.Base f)` で `P.proj.map f = P.Base f` だから。 -/
theorem proj_inverts_baseIso {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀) :
    (baseIsoProp P).IsInvertedBy P.proj := fun _ _ _ hf => hf

def proj_inverts_baseIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — P.proj は base-isomorphism を反転させる",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
