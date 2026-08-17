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

## ★★★主張 2 の骨(原文どおり、2026-08-17 に測定)

原文 (FrdI p.59):
> if C is of Frobenius-normalized type, and A is Frobenius-trivial, then every

| 段 | 中身 | 道具 |
|---|---|---|
| A | `𝔽 → End_𝒞(A)` を作る | `A` の Frobenius-triviality の `ζ` ＋ **Frobenius-normalized の関係式** |
| B | 関手の自己射へ持ち上げる | ★**`Proposition 1.11, (iii)` の `∃!`** |

★★段 B の**自然性・準同型性・`plBkIdObj` での値**は、
★★★**すべて `∃!` の一意性の側**から出る —— ここが要点である。

## ★本ファイルの範囲

★**(i) の 2 主張がともに揃った**(`prop_3_3_i_mem_oTri` / `prop_3_3_i_converse`)。
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

/-! ## ★★主張 2(逆向き)—— どの base-identity pre-step 自己射もこの形で現れる

原文 (FrdI p.59):
> if C is of Frobenius-normalized type, and A is Frobenius-trivial, then every

★★**段は 2 つ**:

| 段 | 中身 | 道具 |
|---|---|---|
| A | `𝔽 → End_𝒞(A)` を作る | `A` の Frobenius-triviality の `ζ` ＋ **Frobenius-normalized の関係式** |
| B | それを関手の自己射へ持ち上げる | ★**`Proposition 1.11, (iii)` の `∃!`** |

★★段 B の**自然性・準同型性・`plBkIdObj` での値**は、
**すべて `∃!` の一意性の側**から出る —— これが要点である。 -/

/-! ### ★段 B —— `Prop 1.11, (iii)` による対象ごとの一意持ち上げ -/

section Lift

variable {A : C} (α : End A) (hα : IsBaseIdentity P α)

include hα in
/-- ★`Prop 1.11, (iii)` を当てるための四角形。★`α` が base-identity なので自明。 -/
theorem plBkLift_sq_aux (Z : Over (⟨A⟩ : PlBk P)) :
    P.Base Z.hom.hom ≫ P.Base (α : A ⟶ A)
      = (𝟙 ((P.toElem.obj Z.left.obj).base) : _ ⟶ _) ≫ P.Base Z.hom.hom := by
  rw [show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα, P.Base_id]
  simp

/-- ★★**対象ごとの一意持ち上げ** —— `Prop 1.11, (iii)` の `∃!` から取る。 -/
noncomputable def plBkLift (Z : Over (⟨A⟩ : PlBk P)) : End Z.left.obj :=
  (prop_1_11_iii P Z.hom.hom Z.hom.property α (𝟙 _)
    (plBkLift_sq_aux P α hα Z)).choose

theorem plBkLift_base (Z : Over (⟨A⟩ : PlBk P)) :
    P.Base (plBkLift P α hα Z : Z.left.obj ⟶ Z.left.obj) = 𝟙 _ :=
  ((prop_1_11_iii P Z.hom.hom Z.hom.property α (𝟙 _)
    (plBkLift_sq_aux P α hα Z)).choose_spec.1).1

theorem plBkLift_sq (Z : Over (⟨A⟩ : PlBk P)) :
    Z.hom.hom ≫ (α : A ⟶ A)
      = (plBkLift P α hα Z : Z.left.obj ⟶ Z.left.obj) ≫ Z.hom.hom :=
  ((prop_1_11_iii P Z.hom.hom Z.hom.property α (𝟙 _)
    (plBkLift_sq_aux P α hα Z)).choose_spec.1).2

/-- ★★**一意性** —— 「底が恒等」と「四角形」を満たす自己射は `plBkLift` だけ。
★★★**自然性も準同型性も、すべてここから出る。** -/
theorem plBkLift_uniq (Z : Over (⟨A⟩ : PlBk P)) (β : End Z.left.obj)
    (h1 : P.Base (β : Z.left.obj ⟶ Z.left.obj) = 𝟙 _)
    (h2 : Z.hom.hom ≫ (α : A ⟶ A) = (β : Z.left.obj ⟶ Z.left.obj) ≫ Z.hom.hom) :
    β = plBkLift P α hα Z :=
  (prop_1_11_iii P Z.hom.hom Z.hom.property α (𝟙 _)
    (plBkLift_sq_aux P α hα Z)).choose_spec.2 β ⟨h1, h2⟩

/-- ★`plBkLift_uniq` の向きを変えたもの(書き換えに使う)。 -/
theorem plBkLift_eq (Z : Over (⟨A⟩ : PlBk P)) (β : End Z.left.obj)
    (h1 : P.Base (β : Z.left.obj ⟶ Z.left.obj) = 𝟙 _)
    (h2 : Z.hom.hom ≫ (α : A ⟶ A) = (β : Z.left.obj ⟶ Z.left.obj) ≫ Z.hom.hom) :
    plBkLift P α hα Z = β :=
  (plBkLift_uniq P α hα Z β h1 h2).symm

/-- ★★**自然変換**になる —— 自然性は `plBkLift_uniq` を
**`g.left.hom` 自身を pull-back として使って**当てると出る。 -/
noncomputable def plBkNat : End (plBkOverToC P A) where
  app Z := plBkLift P α hα Z
  naturality {Z W} g := by
    -- ★`g.left.hom` も pull-back 射である
    have hZ : Z.hom.hom = g.left.hom ≫ W.hom.hom :=
      (congrArg InducedWideCategory.Hom.hom (Over.w g)).symm
    -- ★`g.left.hom ≫ plBkLift W` は `Z` 側の一意持ち上げの条件を満たす
    have hbase : P.Base (g.left.hom ≫ (plBkLift P α hα W : W.left.obj ⟶ W.left.obj))
        = P.Base g.left.hom := by
      rw [P.Base_comp, plBkLift_base P α hα W, Category.comp_id]
    have key : (g.left.hom ≫ (plBkLift P α hα W : W.left.obj ⟶ W.left.obj))
        = (plBkLift P α hα Z : Z.left.obj ⟶ Z.left.obj) ≫ g.left.hom := by
      -- ★`g.left.hom` を pull-back として `Prop 1.11, (iii)` を当てる
      have huniq := prop_1_11_iii P g.left.hom g.left.property
        (plBkLift P α hα W) (𝟙 _) (by
          rw [plBkLift_base P α hα W]; simp)
      obtain ⟨γ, ⟨hγ1, hγ2⟩, hγu⟩ := huniq
      have h1 : (γ : Z.left.obj ⟶ Z.left.obj) = plBkLift P α hα Z := by
        refine plBkLift_uniq P α hα Z γ hγ1 ?_
        rw [hZ, Category.assoc, plBkLift_sq P α hα W, ← Category.assoc, hγ2,
          Category.assoc]
      rw [← h1, hγ2]
    exact key

@[simp] theorem plBkNat_app (Z : Over (⟨A⟩ : PlBk P)) :
    NatTrans.app (plBkNat P α hα) Z = plBkLift P α hα Z := rfl

/-- ★`plBkNat` の成分はすべて base-isomorphism —— 底が恒等だから。 -/
theorem plBkNat_mem : plBkNat P α hα ∈ endBsIso P A := by
  intro Z
  show IsIso (P.Base (plBkLift P α hα Z : Z.left.obj ⟶ Z.left.obj))
  rw [plBkLift_base P α hα Z]
  infer_instance

end Lift

/-- ★**base-identity 自己射のなす部分モノイド**。 -/
def bsIdSubmonoid (A : C) : Submonoid (End A) where
  carrier := {φ : End A | IsBaseIdentity P φ}
  one_mem' := rfl
  mul_mem' {x y} hx hy := by
    show P.Base ((y : A ⟶ A) ≫ (x : A ⟶ A)) = P.Base (𝟙 A)
    rw [P.Base_comp, show P.Base (y : A ⟶ A) = P.Base (𝟙 A) from hy,
      show P.Base (x : A ⟶ A) = P.Base (𝟙 A) from hx, P.Base_id, Category.id_comp]

/-- ★★**持ち上げは準同型** —— `map_one'` も `map_mul'` も `plBkLift_uniq` から。 -/
noncomputable def plBkNatHom (A : C) : bsIdSubmonoid P A →* endBsIso P A where
  toFun α := ⟨plBkNat P (α : End A) α.2, plBkNat_mem P (α : End A) α.2⟩
  map_one' := by
    refine Subtype.ext (NatTrans.ext (funext fun Z => ?_))
    refine plBkLift_eq P (1 : End A) (Submonoid.one_mem (bsIdSubmonoid P A)) Z (1 : End _)
      (P.Base_id _) ?_
    show Z.hom.hom ≫ 𝟙 A = 𝟙 Z.left.obj ≫ Z.hom.hom
    simp
  map_mul' x y := by
    refine Subtype.ext (NatTrans.ext (funext fun Z => ?_))
    refine plBkLift_eq P ((x : End A) * (y : End A))
      (Submonoid.mul_mem (bsIdSubmonoid P A) x.2 y.2) Z
      ((plBkLift P (x : End A) x.2 Z) * (plBkLift P (y : End A) y.2 Z)) ?_ ?_
    · show P.Base ((plBkLift P (y : End A) y.2 Z : Z.left.obj ⟶ Z.left.obj)
          ≫ (plBkLift P (x : End A) x.2 Z : Z.left.obj ⟶ Z.left.obj)) = 𝟙 _
      rw [P.Base_comp, plBkLift_base P (x : End A) x.2 Z,
        plBkLift_base P (y : End A) y.2 Z, Category.id_comp]
    · show Z.hom.hom ≫ (((y : End A) : A ⟶ A) ≫ ((x : End A) : A ⟶ A))
        = ((plBkLift P (y : End A) y.2 Z : Z.left.obj ⟶ Z.left.obj)
            ≫ (plBkLift P (x : End A) x.2 Z : Z.left.obj ⟶ Z.left.obj)) ≫ Z.hom.hom
      rw [← Category.assoc, plBkLift_sq P (y : End A) y.2 Z, Category.assoc,
        plBkLift_sq P (x : End A) x.2 Z, ← Category.assoc]

/-- ★**恒等対象での値は `α` 自身** —— `𝟙 A` を pull-back として一意性を当てる。 -/
theorem plBkNatHom_at_id (A : C) (α : bsIdSubmonoid P A) :
    NatTrans.app ((plBkNatHom P A α : endBsIso P A) : End (plBkOverToC P A))
        (plBkIdObj P A)
      = (α : End A) := by
  refine plBkLift_eq P (α : End A) α.2 (plBkIdObj P A) (α : End A) ?_ ?_
  · exact (show P.Base ((α : End A) : A ⟶ A) = P.Base (𝟙 A) from α.2).trans (P.Base_id A)
  · exact (Category.id_comp ((α : End A) : A ⟶ A)).trans
      (Category.comp_id ((α : End A) : A ⟶ A)).symm

/-! ### ★段 A —— `𝔽 → End_𝒞(A)` を Frobenius-normalized から作る -/

section StepA

variable {A : C} (hfn : IsFrobeniusNormalized P A)

include hfn in
/-- ★★**Frobenius-normalized の関係式を `End A` の積で書き直したもの**。

★原文は `α^d ◦ φ = φ ◦ α`。★`End A` の積は `x * y = y ≫ x` なので
**`α^d * φ = φ * α`** になる。 -/
theorem frobNorm_mul (φ : End A) (hφ : IsBaseIdentity P φ) (α : OTri P A) :
    ((α : End A) ^ (P.degFr (φ : A ⟶ A) : ℕ)) * φ = φ * (α : End A) := by
  have h := hfn φ hφ (α : End A) α.2
  exact h

include hfn in
/-- ★★**関係式の `b` 乗版** —— `α^(d·b) * φ = φ * α^b`。★`b` についての帰納。 -/
theorem frobNorm_mul_pow (φ : End A) (hφ : IsBaseIdentity P φ) (α : OTri P A) (b : ℕ) :
    ((α : End A) ^ ((P.degFr (φ : A ⟶ A) : ℕ) * b)) * φ = φ * ((α : End A) ^ b) := by
  induction b with
  | zero => simp
  | succ n ih =>
    have hstep : (P.degFr (φ : A ⟶ A) : ℕ) * (n + 1)
        = (P.degFr (φ : A ⟶ A) : ℕ) * n + (P.degFr (φ : A ⟶ A) : ℕ) := by ring
    rw [hstep, pow_add, mul_assoc, frobNorm_mul P hfn φ hφ α, ← mul_assoc, ih,
      mul_assoc, ← pow_succ]

end StepA

/-- ★★★**主張 2 の段 A** —— `𝔽 → End_𝒞(A)`。

原文 (FrdI p.61):
> as in the definition of the term “Frobenius-trivial” [cf. the homomorphism “ζ” of

★`⟨a, n⟩ ↦ θ^a * ζ n`。★**乗法性がちょうど Frobenius-normalized の関係式**である。 -/
noncomputable def frobHomOfNormalized {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : ℕ+ →* End A) (hζd : ∀ n : ℕ+, P.degFr (ζ n) = n)
    (hζb : ∀ n : ℕ+, IsBaseIdentity P (ζ n)) (θ : OTri P A) :
    ElemFrob.Standard →* End A where
  toFun x := ((θ : End A) ^ x.div) * ζ x.deg
  map_one' := by
    show ((θ : End A) ^ (0 : ℕ)) * ζ 1 = 1
    rw [pow_zero, one_mul, ζ.map_one]
  map_mul' x y := by
    show ((θ : End A) ^ (x.div + (x.deg : ℕ) • y.div)) * ζ (x.deg * y.deg)
      = (((θ : End A) ^ x.div) * ζ x.deg) * (((θ : End A) ^ y.div) * ζ y.deg)
    have hkey : ((θ : End A) ^ ((x.deg : ℕ) * y.div)) * ζ x.deg
        = ζ x.deg * ((θ : End A) ^ y.div) := by
      have := frobNorm_mul_pow P hfn (ζ x.deg) (hζb x.deg) θ y.div
      rwa [hζd x.deg] at this
    rw [ζ.map_mul, smul_eq_mul, pow_add]
    calc ((θ : End A) ^ x.div * (θ : End A) ^ ((x.deg : ℕ) * y.div)) * (ζ x.deg * ζ y.deg)
        = (θ : End A) ^ x.div
            * ((((θ : End A) ^ ((x.deg : ℕ) * y.div)) * ζ x.deg) * ζ y.deg) := by
          simp only [mul_assoc]
      _ = (θ : End A) ^ x.div * ((ζ x.deg * ((θ : End A) ^ y.div)) * ζ y.deg) := by
          rw [hkey]
      _ = ((θ : End A) ^ x.div * ζ x.deg) * (((θ : End A) ^ y.div) * ζ y.deg) := by
          simp only [mul_assoc]

/-- ★★★**[FrdI] Proposition 3.3, (i) の主張 2** ——
`A` が Frobenius-trivial かつ Frobenius-normalized なら、
`A` の**どの** base-identity pre-step 自己射(`𝒪^▷(A)` の元)も
`𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` による `1 ∈ ℤ≥0` の像として現れる。

原文 (FrdI p.59):
> if C is of Frobenius-normalized type, and A is Frobenius-trivial, then every

★★段 A で `𝔽 → End_𝒞(A)` を作り、段 B(`Proposition 1.11, (iii)`)で
**関手の自己射へ一意に持ち上げる**。★`plBkIdObj` での値が `θ` そのものになる。 -/
theorem prop_3_3_i_converse {A : C} (hft : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (θ : OTri P A) :
    ∃ f : ElemFrob.Standard →* endBsIso P A,
      NatTrans.app ((f ⟨1, 1⟩ : endBsIso P A) : End (plBkOverToC P A))
          (plBkIdObj P A)
        = (θ : End A) := by
  obtain ⟨ζ, hζd, hζb⟩ := hft
  set F := frobHomOfNormalized P hfn ζ hζd (fun n => (hζb n).1) θ with hF
  have hmem : ∀ x : ElemFrob.Standard, F x ∈ bsIdSubmonoid P A := by
    intro x
    refine Submonoid.mul_mem _ (pow_mem ?_ _) ?_
    · exact θ.2.1
    · exact (hζb x.deg).1
  refine ⟨(plBkNatHom P A).comp (F.codRestrict (bsIdSubmonoid P A) hmem), ?_⟩
  rw [MonoidHom.comp_apply, plBkNatHom_at_id]
  show ((θ : End A) ^ (1 : ℕ)) * ζ 1 = (θ : End A)
  rw [pow_one, ζ.map_one, mul_one]

/-- ★★★**[FrdI] Proposition 3.3, (i)** —— 2 主張がともに実装された。

★主張 1 は `prop_3_3_i_mem_oTri`(`𝒟` が Frobenius-slim)、
★主張 2 は `prop_3_3_i_converse`(`A` が Frobenius-trivial かつ Frobenius-normalized)。 -/
def prop_3_3_i_mem_oTri.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59, item := "Proposition 3.3, (i)",
    sectionId := "frdi-prop-3-3" }

end ABC3.Found.FrdI
