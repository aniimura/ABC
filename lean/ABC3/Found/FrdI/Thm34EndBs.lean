/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Quasi

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

end ABC3.Found.FrdI
