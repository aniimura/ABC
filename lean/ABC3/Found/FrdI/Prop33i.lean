import ABC3.Found.FrdI.Remark312
import ABC3.Found.FrdI.Prop22

/-!
# [FrdI] Proposition 3.3, (i) —— `End(𝒞^pl-bk_A → 𝒞)^bs-iso`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.59–p.61。

原文 (FrdI p.59):
> Proposition 3.3. (Base-identity Pre-steps and Units)

## ★★(i) の主張は 2 つ

原文 (FrdI p.60):
> Proof. First, we consider assertion (i). Note that since the composite of the functor

| # | 内容 |
|---|---|
| 1 | `𝒟` が Frobenius-slim なら、`𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` による `1 ∈ ℤ≥0` の像の成分は**すべて base-identity pre-step** |
| 2 | 逆に `𝒞` が Frobenius-normalized 型で `A` が Frobenius-trivial なら、`A` のどの base-identity pre-step 自己射もそう現れる |

## ★★★主張 1 の骨(原文どおり、2026-08-17 に測定)

原文 (FrdI p.60):
> obtained by considering the Frobenius degree of the induced endomorphism of A

★**2 つの準同型に分けて、それぞれ `𝔽 ↠ ℕ≥1` を経由させる**:

| 成分 | 行き先 | 経由の理由 |
|---|---|---|
| **底** | `Aut(𝒟_{A_𝒟} → 𝒟)` | ★`𝒟` が **Frobenius-slim**(定義そのもの) |
| **次数** | `ℕ≥1` | ★★**`ℕ≥1` が可換かつ簡約的**(`elemFrob_factors_of_cancel`) |

★どちらも `⟨1,1⟩` の次数が `1` なので、像は `1` になる。
★**底が恒等 = base-identity、次数が 1 = linear。合わせて base-identity pre-step**
(すなわち `𝒪^▷(−)` の元)である。

★★底の側で `𝒞^pl-bk_A → 𝒟` が `𝒟_{A_𝒟} → 𝒟` を経由することは
**`Definition 1.3, (i), (c)` の圏同値**(`plBkEquiv`)である。

## ★本ファイルの範囲

★**土台(対象の定義)と、次数の側**をここで取る。
★底の側は `plBkEquiv` を通した `Aut` の移送が要る(別途)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★自然な関手 `𝒞^pl-bk_A → 𝒞` -/

/-- ★**`𝒞^pl-bk_A → 𝒞`** —— pull-back 射の広い部分圏のスライスから `𝒞` へ戻る。

★`map_id` も `map_comp` も **`rfl`** である(広い部分圏は射を包んでいるだけだから)。 -/
def plBkOverToC (A : C) : Over (⟨A⟩ : PlBk P) ⥤ C where
  obj Z := Z.left.obj
  map f := f.left.hom
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★**底へ落とすと `plBkOverFunctor` を経由する** —— `rfl`。

★★これが原文の「the composite of the functor `𝒞^pl-bk_A → 𝒞` with the natural
projection functor `𝒞 → 𝒟` factors as …」の中身である。 -/
theorem plBkOverToC_comp_proj (A : C) :
    plBkOverToC P A ⋙ P.proj = plBkOverFunctor P A ⋙ Over.forget _ := rfl

/-! ## ★`End(𝒞^pl-bk_A → 𝒞)^bs-iso` -/

/-- ★★**原文の `End(𝒞^pl-bk_A → 𝒞)^bs-iso`** ——
成分がすべて base-isomorphism である自然変換のなす部分モノイド。

★`End F` の積は `x * y = y ≫ x` なので、`mul_mem` は
**base-isomorphism が合成で閉じる**ことに帰着する。 -/
def endBsIso (A : C) : Submonoid (End (plBkOverToC P A)) where
  carrier := {η | ∀ Z, IsBaseIsomorphism P (NatTrans.app η Z)}
  one_mem' := fun Z => by
    show IsIso (P.Base (NatTrans.app (𝟙 (plBkOverToC P A)) Z))
    rw [show NatTrans.app (𝟙 (plBkOverToC P A)) Z = 𝟙 _ from rfl, P.Base_id]
    infer_instance
  mul_mem' {x y} hx hy := fun Z => by
    show IsIso (P.Base (NatTrans.app (y ≫ x) Z))
    rw [show NatTrans.app (y ≫ x) Z = NatTrans.app y Z ≫ NatTrans.app x Z from rfl]
    exact isBaseIsomorphism_comp P (hy Z) (hx Z)

/-! ## ★次数の側 —— `End(…)^bs-iso →* ℕ≥1`

★**各対象 `Z` での成分**の Frobenius 次数を取る。
★次数は乗法的(`degFr_comp`)なので準同型になる。 -/

/-- ★恒等射が定める `𝒞^pl-bk_A` の対象。 -/
def plBkIdObj (A : C) : Over (⟨A⟩ : PlBk P) :=
  Over.mk (⟨𝟙 A, isPullBack_of_isIso P (𝟙 A)⟩ :
    (⟨A⟩ : PlBk P) ⟶ (⟨A⟩ : PlBk P))

/-- ★★**次数を取る準同型** `End(𝒞^pl-bk_A → 𝒞) →* ℕ≥1`(対象 `Z` での成分)。

★`End` の積が `x * y = y ≫ x` で、`degFr` が `degFr (ψ ≫ φ) = degFr φ * degFr ψ`
なので、**両方の反転が打ち消して**準同型になる。 -/
def endDegHom (A : C) (Z : Over (⟨A⟩ : PlBk P)) : End (plBkOverToC P A) →* ℕ+ where
  toFun η := P.degFr (NatTrans.app η Z)
  map_one' := by
    show P.degFr (NatTrans.app (𝟙 (plBkOverToC P A)) Z) = 1
    rw [show NatTrans.app (𝟙 (plBkOverToC P A)) Z = 𝟙 _ from rfl]
    exact degFr_of_isIso P (𝟙 _)
  map_mul' x y := by
    show P.degFr (NatTrans.app (y ≫ x) Z) = _
    rw [show NatTrans.app (y ≫ x) Z = NatTrans.app y Z ≫ NatTrans.app x Z from rfl,
      P.degFr_comp]

/-- ★★★**主張 1 の「次数の側」** —— `𝔽 → End(𝒞^pl-bk_A → 𝒞)` による
`⟨1,1⟩` の像は **linear**(次数 1)である。

原文 (FrdI p.60):
> obtained by considering the Frobenius degree of the induced endomorphism of A

★★**`ℕ≥1` が可換かつ簡約的**なので `elemFrob_factors_of_cancel` が効き、
合成 `𝔽 → ℕ≥1` は `deg` を経由する。★`⟨1,1⟩` の次数は `1` なので像も `1`。 -/
theorem endDeg_gen_eq_one (A : C) (Z : Over (⟨A⟩ : PlBk P))
    (f : ElemFrob ℕ →* End (plBkOverToC P A)) :
    P.degFr (NatTrans.app (f ⟨1, 1⟩) Z) = 1 := by
  obtain ⟨g, hg⟩ := elemFrob_factors_of_cancel ((endDegHom P A Z).comp f)
  have h := congrArg (fun m : ElemFrob ℕ →* ℕ+ => m ⟨1, 1⟩) hg
  show (endDegHom P A Z) (f ⟨1, 1⟩) = 1
  show ((endDegHom P A Z).comp f) ⟨1, 1⟩ = 1
  rw [h]
  show g ((⟨1, 1⟩ : ElemFrob ℕ).deg) = 1
  show g 1 = 1
  exact g.map_one

/-! ## ★底の側 —— `Aut(𝒟_{A_𝒟} → 𝒟)` へ移して Frobenius-slim を当てる

原文 (FrdI p.60):
> Proof. First, we consider assertion (i). Note that since the composite of the functor

★★**`plBkOverToC ⋙ P.proj = plBkOverFunctor ⋙ Over.forget` は `rfl`**(上)であり、
`plBkOverFunctor` は `Definition 1.3, (i), (c)` により**圏同値**である。
★★★したがって **`F ↦ plBkOverFunctor ⋙ F` は関手圏の圏同値**
(`Equivalence.congrLeft`)——`Aut` はそのまま移る。
★**`Prop 3.11, (iii)` の rigidity で使ったのと同じ道具**である。
-/

variable (Fc : FrobenioidCore P)

/-- ★**`η` の底成分がなす自己同型**。★成分がすべて同型(`bs-iso`)なので
`whiskerRight` が自然同型になる。 -/
noncomputable def baseAut (A : C) (η : endBsIso P A) :
    Aut (plBkOverToC P A ⋙ P.proj) := by
  haveI : ∀ Z, IsIso (NatTrans.app
      (Functor.whiskerRight (η : End (plBkOverToC P A)) P.proj) Z) := fun Z => η.2 Z
  haveI := NatIso.isIso_of_isIso_app
    (Functor.whiskerRight (η : End (plBkOverToC P A)) P.proj)
  exact asIso (Functor.whiskerRight (η : End (plBkOverToC P A)) P.proj)

@[simp] theorem baseAut_hom_app (A : C) (η : endBsIso P A) (Z : Over (⟨A⟩ : PlBk P)) :
    NatTrans.app (baseAut P A η).hom Z
      = P.Base (NatTrans.app (η : End (plBkOverToC P A)) Z) := rfl

/-- ★★**底成分を取る準同型** `End(…)^bs-iso →* Aut(𝒞^pl-bk_A → 𝒟)`。

★`End` の積は `x * y = y ≫ x`、`Aut` の積は `x * y = y ≪≫ x` で**向きが揃う**ので、
`whiskerRight` の関手性がそのまま準同型性になる。 -/
noncomputable def baseAutHom (A : C) :
    endBsIso P A →* Aut (plBkOverToC P A ⋙ P.proj) where
  toFun := baseAut P A
  map_one' := by
    refine Aut.ext (NatTrans.ext (funext fun Z => ?_))
    show P.Base (NatTrans.app (𝟙 (plBkOverToC P A)) Z) = _
    rw [show NatTrans.app (𝟙 (plBkOverToC P A)) Z = 𝟙 _ from rfl, P.Base_id]
    rfl
  map_mul' x y := by
    refine Aut.ext (NatTrans.ext (funext fun Z => ?_))
    exact P.Base_comp _ _

/-- ★★★**`Aut` の移送** —— `Aut(𝒟_{A_𝒟} → 𝒟) ≃* Aut(𝒞^pl-bk_A → 𝒟)`。

★★**`F ↦ plBkOverFunctor ⋙ F` 自体が関手圏の圏同値**(`Equivalence.congrLeft`)
なので、`Aut` はそのまま移る。★これが原文の
「the composite of the functor `𝒞^pl-bk_A → 𝒞` with the natural projection
functor `𝒞 → 𝒟` factors as …」の使いどころである。

★`Prop 3.11, (iii)` の rigidity(`isRigidFunctor_comp_of_isEquivalence`)と
**同じ道具**で、そこでは充満忠実性を直に使い、ここでは
mathlib の `autMulEquivOfFullyFaithful` に載せている。 -/
noncomputable def baseAutEquiv (A : C) [(plBkOverFunctor P A).IsEquivalence] :
    Aut (Over.forget ((P.toElem.obj A).base)) ≃* Aut (plBkOverToC P A ⋙ P.proj) :=
  haveI : (((plBkOverFunctor P A).asEquivalence.congrLeft (E := D)).inverse).IsEquivalence :=
    inferInstance
  (Functor.FullyFaithful.ofFullyFaithful
      ((plBkOverFunctor P A).asEquivalence.congrLeft (E := D)).inverse).autMulEquivOfFullyFaithful
    (Over.forget ((P.toElem.obj A).base))

include P Fc in
/-- ★★★**主張 1 の「底の側」** —— `𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` による
`⟨1,1⟩` の像の底成分は**恒等**である。

原文 (FrdI p.60):
> Proof. First, we consider assertion (i). Note that since the composite of the functor

★★`𝒟` の **Frobenius-slim 性**を `A_𝒟` に当て、`baseAutEquiv` で戻す。
★`⟨1,1⟩` の次数は `1` なので、経由した `ℕ≥1` 側で `1` になる。 -/
theorem baseAut_gen_eq_one (hslim : IsFrobeniusSlim D) (A : C)
    (f : ElemFrob.Standard →* endBsIso P A) :
    baseAutHom P A (f ⟨1, 1⟩) = 1 := by
  haveI := Fc.plBkEquiv A
  obtain ⟨g, hg⟩ := hslim ((P.toElem.obj A).base)
    (((baseAutEquiv P A).symm.toMonoidHom).comp ((baseAutHom P A).comp f))
  have h := congrArg
    (fun m : ElemFrob.Standard →* Aut (Over.forget ((P.toElem.obj A).base)) => m ⟨1, 1⟩) hg
  have h1 : (baseAutEquiv P A).symm (baseAutHom P A (f ⟨1, 1⟩)) = 1 := by
    refine h.trans ?_
    show g ((⟨1, 1⟩ : ElemFrob ℕ).deg) = 1
    show g 1 = 1
    exact g.map_one
  have h2 := congrArg (baseAutEquiv P A) h1
  rwa [MulEquiv.apply_symm_apply, map_one] at h2

include P Fc in
/-- ★★★**[FrdI] Proposition 3.3, (i) の主張 1** ——
`𝒟` が Frobenius-slim なら、`𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` による
`1 ∈ ℤ≥0` の像の**各成分は base-identity pre-step**、すなわち `𝒪^▷(−)` の元である。

原文 (FrdI p.59):
> Proposition 3.3. (Base-identity Pre-steps and Units)

★**底が恒等**(`baseAut_gen_eq_one`、Frobenius-slim)かつ
**次数が 1**(`endDeg_gen_eq_one`、`ℕ≥1` の可換性と簡約性)。 -/
theorem prop_3_3_i_mem_oTri (hslim : IsFrobeniusSlim D) (A : C)
    (f : ElemFrob.Standard →* endBsIso P A) (Z : Over (⟨A⟩ : PlBk P)) :
    NatTrans.app (f ⟨1, 1⟩ : End (plBkOverToC P A)) Z ∈ OTri P Z.left.obj := by
  refine ⟨?_, ?_⟩
  · -- ★底の側
    have h := baseAut_gen_eq_one P Fc hslim A f
    have hZ := congrArg (fun t : Aut (plBkOverToC P A ⋙ P.proj) => NatTrans.app t.hom Z) h
    show P.Base (NatTrans.app (f ⟨1, 1⟩ : End (plBkOverToC P A)) Z) = P.Base (𝟙 _)
    rw [P.Base_id]
    exact hZ
  · -- ★次数の側
    exact endDeg_gen_eq_one P A Z ((endBsIso P A).subtype.comp f)

/-! ## ★残り —— 主張 2(逆向き)

★`𝒞` が Frobenius-normalized 型で `A` が Frobenius-trivial なら、
`A` の**どの** base-identity pre-step 自己射もこの形で現れる、という逆である。
★★これは `Frobenius-normalized` の定義そのものが与える構成なので、
そちらの語彙が要る(別途)。 -/

/-- ★★★**[FrdI] Proposition 3.3, (i) の主張 1** —— 2 主張のうち 1 つ。

★★主張 2(逆向き)は未実装なので**条つき**で記録する。 -/
def prop_3_3_i_mem_oTri.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59, item := "Proposition 3.3, (i) — 主張 1",
    sectionId := "frdi-prop-3-3" }

end ABC3.Found.FrdI
