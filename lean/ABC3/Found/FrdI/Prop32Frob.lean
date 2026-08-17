import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle

/-!
# [FrdI] Proposition 3.2, (iii) —— `𝒞^pf` が Frobenioid であること(道具立て)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.59。

原文 (FrdI p.59):
> (iii) The category Cpf, equipped with the functor Cpf →FΦpf of the diagram of

★`Prop32.lean` までで **(i)**(1-可換図式と pre-Frobenioid 構造 `pfRootPre`)と
**(ii) の 5 項**(linear / base-isomorphism / isometry / pre-step / 次数)が出ている。

## ★本ファイルの役割 —— 「代表元で計算する」道具

★★`Hom^pf` は **Frobenius 型の根**についての filtered colimit なので、
`𝒞^birat` のとき(`Prop44Core.lean`)と違い、
★**`𝒪^▷` の元も `𝒞` の射そのもので代表される**(比ではない)。
★したがって `Definition 1.3` の 21 条は原則として
「代表元へ降ろす → `𝒞` の条件を使う → 押し出す」で渡る。

★本ファイルではまずその**土台**を置く:

| 宣言 | 内容 |
|---|---|
| `idxMk3` | 3 つ組の添字対象を 3 本の Frobenius 型射から作る |
| `compPf_mk_pair` | ★**真ん中を共有する代表元での合成則**(`compPf_mk` の使いやすい形) |
| `mk_diag_id` | 対角の添字での `𝟙` は `toHomPf (𝟙)` |
| `homPf_inv_pair` | ★★`𝒞` の同型は `Hom^pf` でも両側逆射を与える |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★1. 3 つ組の添字対象 -/

variable {P F} in
/-- ★3 つ組の添字対象を、次数の等しい 3 本の Frobenius 型射から作る。 -/
def idxMk3 {A B E A' B' E' : C} (a : A ⟶ A') (b : B ⟶ B') (e : E ⟶ E')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (he : IsFrobeniusType P e)
    (hab : P.degFr a = P.degFr b) (hbe : P.degFr b = P.degFr e) : IdxPf3 P F A B E :=
  Under.mk (Y := (⟨(A', B', E')⟩ : TriFr P F))
    (show triFrObj P F A B E ⟶ (⟨(A', B', E')⟩ : TriFr P F) from
      ⟨(a, b, e), ha, hb, he, hab, hbe⟩)

/-! ## ★2. 真ん中を共有する代表元での合成則 -/

variable {P F} in
/-- ★★**合成則(使いやすい形)** —— 3 本の根を明示して `compPf_mk` を当てる。 -/
theorem compPf_mk_pair {A B E A' B' E' : C} (a : A ⟶ A') (b : B ⟶ B') (e : E ⟶ E')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (he : IsFrobeniusType P e)
    (hab : P.degFr a = P.degFr b) (hbe : P.degFr b = P.degFr e)
    (φ : A' ⟶ B') (ψ : B' ⟶ E') :
    compPf P F (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hab) φ)
        (HomPf.mk (idxMk (P := P) (F := F) b e hb he hbe) ψ)
      = HomPf.mk (idxMk (P := P) (F := F) a e ha he (hab.trans hbe)) (φ ≫ ψ) :=
  compPf_mk (idxMk3 (F := F) a b e ha hb he hab hbe) φ ψ

/-! ## ★3. 対角の添字での恒等射 -/

variable {P F} in
/-- ★添字の始対象から対角の添字への遷移射。 -/
def idxOneToDiag {A A' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (hd : P.degFr a = P.degFr a) :
    idxOne P F A A ⟶ idxMk (P := P) (F := F) a a ha ha hd :=
  Under.homMk (show (⟨(A, A)⟩ : BiFr P F) ⟶ (⟨(A', A')⟩ : BiFr P F) from ⟨(a, a), ha, ha, rfl⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp a) (Category.id_comp a)))

variable {P F} in
/-- ★**対角の添字での `𝟙` は `𝒞 → 𝒞^pf` の像の `𝟙`**。 -/
theorem mk_diag_id {A A' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (hd : P.degFr a = P.degFr a) :
    HomPf.mk (idxMk (P := P) (F := F) a a ha ha hd) (𝟙 A') = toHomPf (F := F) (𝟙 A) := by
  have ht : idxTransport P F (idxOneToDiag (F := F) a ha hd) (𝟙 A) = 𝟙 A' := by
    refine frobTransport_eq (F := F) _ _ _ _ _ _ _ ?_
    show 𝟙 A ≫ a = a ≫ 𝟙 A'
    rw [Category.id_comp, Category.comp_id]
  rw [← ht, HomPf.mk_map]
  rfl

/-! ## ★4. `𝒞` の同型は `Hom^pf` でも両側逆射を与える -/

variable {P F} in
/-- ★★**`𝒞` の同型を `Hom^pf` へ送ると両側逆射になる**。

★★`𝒞^birat` と違い、添字は **Frobenius 型の根**なので、
`(a, b)` を `(b, a)` に**入れ替えるだけ**で逆側の添字が作れる。 -/
theorem homPf_inv_pair {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A' ⟶ B') (ψ : B' ⟶ A') (h1 : φ ≫ ψ = 𝟙 A') (h2 : ψ ≫ φ = 𝟙 B') :
    compPf P F (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) φ)
          (HomPf.mk (idxMk (P := P) (F := F) b a hb ha hd.symm) ψ)
        = toHomPf (F := F) (𝟙 A)
      ∧ compPf P F (HomPf.mk (idxMk (P := P) (F := F) b a hb ha hd.symm) ψ)
          (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) φ)
        = toHomPf (F := F) (𝟙 B) := by
  refine ⟨(compPf_mk_pair a b a ha hb ha hd hd.symm φ ψ).trans ?_,
    (compPf_mk_pair b a b hb ha hb hd.symm hd ψ φ).trans ?_⟩
  · exact (congrArg (HomPf.mk (idxMk (P := P) (F := F) a a ha ha (hd.trans hd.symm))) h1).trans
      (mk_diag_id a ha _)
  · exact (congrArg (HomPf.mk (idxMk (P := P) (F := F) b b hb hb (hd.symm.trans hd))) h2).trans
      (mk_diag_id b hb _)

variable {P F} in
/-- ★★**入れ替えた添字での合成**(片側だけ)。 -/
theorem compPf_mk_flip {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A' ⟶ B') (ψ : B' ⟶ A') (h1 : φ ≫ ψ = 𝟙 A') :
    compPf P F (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) φ)
        (HomPf.mk (idxMk (P := P) (F := F) b a hb ha hd.symm) ψ)
      = toHomPf (F := F) (𝟙 A) :=
  (compPf_mk_pair a b a ha hb ha hd hd.symm φ ψ).trans
    ((congrArg (HomPf.mk (idxMk (P := P) (F := F) a a ha ha (hd.trans hd.symm))) h1).trans
      (mk_diag_id a ha _))

/-! ## ★5. 全射性から**片側逆射だけで**同型が出る

★★`𝒞^pf` は totally epimorphic(`pfRoot_totEpi`)なので、
`f ≫ g = 𝟙` の 1 本から `g ≫ f = 𝟙` が出る。★**逆側の根の指数合わせを省ける**。 -/

theorem isIso_of_comp_eq_id {K : Type u3} [Category.{v3} K] (h : IsTotallyEpimorphic K)
    {X Y : K} (f : X ⟶ Y) (g : Y ⟶ X) (hfg : f ≫ g = 𝟙 X) : IsIso f := by
  haveI : Epi f := h _ _ f
  refine ⟨g, hfg, ?_⟩
  have hh : f ≫ (g ≫ f) = f ≫ 𝟙 Y := by
    rw [← Category.assoc, hfg, Category.id_comp, Category.comp_id]
  exact (cancel_epi f).mp hh

/-! ## ★6. ★★★代表元が `𝒞` の同型なら、`𝒞^pf` でも同型

★★根の指数は `compRoot_eq_lift` で
`PA = PE = X.root * Y.root`、`PB = X.root * X.root` と**対称に**選ぶ。
★これで逆側の添字(`(a,b)` を入れ替えたもの)が**型どおりに収まる**。 -/

variable {P F} in
/-- ★★★★**代表元が同型なら `𝒞^pf` の射も同型**。 -/
theorem pfRoot_isIso_of_rep {X Y : PfRootObj P F} (f : X ⟶ Y)
    (Z : IdxPf P F (rtObj P F X.obj (X.root * Y.root))
      (rtObj P F Y.obj (X.root * X.root)))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) [hφ : IsIso φ]
    (hf : (rtRootIso P F X.obj Y.obj
        (show X.root * Y.root = X.root * Y.root from rfl)
        (show X.root * X.root = X.root * X.root from rfl)).inv f = HomPf.mk Z φ) :
    IsIso f := by
  obtain ⟨ha, hb, hd⟩ := Z.hom.property
  refine isIso_of_comp_eq_id (pfRoot_totEpi P F) f
    ((rtRootIso P F Y.obj X.obj
        (show X.root * X.root = X.root * X.root from rfl)
        (show X.root * Y.root = X.root * Y.root from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) Z.hom.hom.2 Z.hom.hom.1 hb ha hd.symm) (@inv _ _ _ _ φ hφ)))
    ?_
  have hcomp := compRoot_eq_lift (P := P) (F := F) f
    ((rtRootIso P F Y.obj X.obj
        (show X.root * X.root = X.root * X.root from rfl)
        (show X.root * Y.root = X.root * Y.root from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) Z.hom.hom.2 Z.hom.hom.1 hb ha hd.symm) (@inv _ _ _ _ φ hφ)))
    (c := 1) (PA := X.root * Y.root) (PB := X.root * X.root) (PE := X.root * Y.root)
    (hcA := by rw [one_mul]) (hcB := by rw [one_mul])
    (hcE := by rw [one_mul]; exact mul_comm _ _)
    (ef := X.root) (eg := X.root) (er := Y.root)
    (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl)
    (hrA := mul_comm _ _) (hrE := mul_comm _ _)
  show compRoot P F f _ = idRoot P F X
  have hflip : compPf P F (HomPf.mk Z φ)
      (HomPf.mk (idxMk (P := P) (F := F) Z.hom.hom.2 Z.hom.hom.1 hb ha hd.symm)
        (@inv _ _ _ _ φ hφ))
      = toHomPf (F := F) (𝟙 (rtObj P F X.obj (X.root * Y.root))) :=
    compPf_mk_flip Z.hom.hom.1 Z.hom.hom.2 ha hb hd φ (@inv _ _ _ _ φ hφ)
      (@IsIso.hom_inv_id _ _ _ _ φ hφ)
  refine hcomp.trans ?_
  rw [hf, Iso.hom_inv_id_apply]
  exact (congrArg (ConcreteCategory.hom (rtRootIso P F X.obj X.obj
      (mul_comm X.root Y.root) (mul_comm X.root Y.root)).hom) hflip).trans
    (rtRootIso_hom_id P F X.obj (mul_comm X.root Y.root))

/-! ## ★7. ★★★`𝒞^pf` は isotropic 型

原文 (FrdI p.59):
> (iii) The category Cpf, equipped with the functor Cpf →FΦpf of the diagram of

★★`𝒞` が **Frobenius-isotropic 型**なら、添字を押し上げて
**始域が isotropic になる代表元**が取れる。★そこでは isometric pre-step は同型であり、
`pfRoot_isIso_of_rep` で `𝒞^pf` の同型に上がる。 -/

variable {P F} in
/-- ★添字を「始域が isotropic になる所」まで押し上げる。 -/
theorem exists_idx_isotropic (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (V : IdxPf P F A B) :
    ∃ (W : IdxPf P F A B) (u : V ⟶ W), IsIsotropic P W.right.obj.1 := by
  obtain ⟨Dd, α, hα, hDd⟩ := hfi V.right.obj.1
  obtain ⟨B₂, β, hβ, hβd⟩ := F.frobDegSurj V.right.obj.2 (P.degFr α)
  obtain ⟨hva, hvb, hvd⟩ := V.hom.property
  refine ⟨idxMk (P := P) (F := F) (V.hom.hom.1 ≫ α) (V.hom.hom.2 ≫ β)
      (IsFrobeniusType.comp P F hva hα) (IsFrobeniusType.comp P F hvb hβ) ?_,
    Under.homMk (show V.right ⟶ (⟨(Dd, B₂)⟩ : BiFr P F) from
      ⟨(α, β), hα, hβ, hβd.symm⟩) (WideSubcategory.hom_ext _ rfl), hDd⟩
  rw [P.degFr_comp, P.degFr_comp, hβd]
  exact congrArg (fun t => P.degFr α * t) hvd

variable {P F} in
/-- ★★★★**[FrdI] Proposition 3.2, (iii) の「isotropic 型」の条**。

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural
-/
theorem pfRoot_isOfIsotropicType (hfi : IsOfFrobeniusIsotropicType P) :
    IsOfIsotropicType (pfRootPre P F) := by
  intro X Y f hisom hstep
  -- ★段 1: 「持ち上げた」代表元を取る
  obtain ⟨V, ψ, hV⟩ := HomPf.exists_rep (P := P) (F := F)
    ((rtRootIso P F X.obj Y.obj
      (show X.root * Y.root = X.root * Y.root from rfl)
      (show X.root * X.root = X.root * X.root from rfl)).inv f)
  -- ★段 2: 添字を isotropic な所まで押し上げる
  obtain ⟨W, u, hW⟩ := exists_idx_isotropic (F := F) hfi V
  obtain ⟨ψ', hψ'⟩ : ∃ t : W.right.obj.1 ⟶ W.right.obj.2,
      t = idxTransport P F u ψ := ⟨_, rfl⟩
  have hfW : (rtRootIso P F X.obj Y.obj
      (show X.root * Y.root = X.root * Y.root from rfl)
      (show X.root * X.root = X.root * X.root from rfl)).inv f
      = HomPf.mk W ψ' := by
    rw [hψ', HomPf.mk_map (P := P) (F := F) (Z := V) (W := W) u ψ, hV]
  -- ★段 3: `f` を押し出した添字の代表元として書き直す
  have hpush := rtRootIso_hom_mk (F := F) X.obj Y.obj
    (show X.root * Y.root = X.root * Y.root from rfl)
    (show X.root * X.root = X.root * X.root from rfl) W ψ'
  have hf2 : f = (rtRootIso P F X.obj Y.obj
      (show X.root * Y.root = X.root * Y.root from rfl)
      (show X.root * X.root = X.root * X.root from rfl)).hom (HomPf.mk W ψ') := by
    rw [← hfW, Iso.inv_hom_id_apply]
  have hfmk := hf2.trans hpush
  -- ★段 4: 代表元は `𝒞` の isometric pre-step
  have hps : IsPreStep P ψ' :=
    (isPreStep_mk_iff (X := X) (Y := Y) _ _).mp (hfmk ▸ hstep)
  have his : IsIsometric P ψ' :=
    (isIsometric_mk_iff (X := X) (Y := Y) _ _).mp (hfmk ▸ hisom)
  -- ★段 5: isotropic な始域から出る isometric pre-step は同型
  haveI : IsIso ψ' := hW _ _ his hps
  exact pfRoot_isIso_of_rep f W ψ' hfW

/-! ## ★8. isotropic 型から**ただで出る** 7 条

★★`Proposition 1.4, (i)` により **isotropic 型の Frobenioid ではすべての射が co-angular**
であり、★**isometric pre-step はすべて同型**である。
★これで `Definition 1.3` の (iii)(a)(b)・(v)(b)(c)・(vii)(a)(b) が落ちる。 -/

variable {P F} in
/-- ★★**`𝒞^pf` のすべての射は co-angular**(`Proposition 1.4, (i)`)。 -/
theorem pfRoot_isCoAngular (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) : IsCoAngular (pfRootPre P F) f :=
  fun _ _ _ β _ _ _ hβi hβs _ =>
    pfRoot_isOfIsotropicType (F := F) hfi _ _ β hβi hβs

variable {P F} in
/-- ★**(vii)(b)** —— 全対象が isotropic なので自明。 -/
theorem pfRoot_isotropicClosed (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (_f : X ⟶ Y) (_hX : IsIsotropic (pfRootPre P F) X) : IsIsotropic (pfRootPre P F) Y :=
  pfRoot_isOfIsotropicType (F := F) hfi Y

variable {P F} in
/-- ★**(vii)(a)** —— `𝟙` が isotropic hull。 -/
theorem pfRoot_isotropicHullExists (hfi : IsOfFrobeniusIsotropicType P) (X : PfRootObj P F) :
    ∃ (Y : PfRootObj P F) (φ : X ⟶ Y), IsIsotropicHull (pfRootPre P F) φ := by
  refine ⟨X, 𝟙 X, isIsometric_id _ X, isPreStep_id _ X,
    pfRoot_isOfIsotropicType (F := F) hfi X, ?_⟩
  intro Cc _ γ
  exact ⟨γ, (Category.id_comp γ).symm, fun y hy => by simpa using hy.symm⟩

variable {P F} in
/-- ★**(iii)(a)** —— 全射が co-angular なので自明。 -/
theorem pfRoot_coAngularComp (hfi : IsOfFrobeniusIsotropicType P)
    {X Y Z : PfRootObj P F} (ψ : X ⟶ Y) (φ : Y ⟶ Z)
    (_hψ : IsCoAngular (pfRootPre P F) ψ) (_hφ : IsCoAngular (pfRootPre P F) φ) :
    IsCoAngular (pfRootPre P F) (ψ ≫ φ) :=
  pfRoot_isCoAngular hfi _

variable {P F} in
/-- ★**(iii)(b)** —— 同上。 -/
theorem pfRoot_coAngularOfPreStep (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (α : X ⟶ Y) (_hc : IsCoAngular (pfRootPre P F) α)
    (_hs : IsPreStep (pfRootPre P F) α) (φ : X ⟶ Y) : IsCoAngular (pfRootPre P F) φ :=
  pfRoot_isCoAngular hfi φ

variable {P F} in
/-- ★**(v)(b)** —— `φ = φ ≫ 𝟙`。 -/
theorem pfRoot_preStepFactor (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (φ : X ⟶ Y) (hφ : IsPreStep (pfRootPre P F) φ) :
    ∃ (M : PfRootObj P F) (β : X ⟶ M) (α : M ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (pfRootPre P F) β ∧ IsPreStep (pfRootPre P F) β ∧
        IsIsometric (pfRootPre P F) α ∧ IsPreStep (pfRootPre P F) α :=
  ⟨Y, φ, 𝟙 Y, (Category.comp_id φ).symm, pfRoot_isCoAngular hfi φ, hφ,
    isIsometric_id _ Y, isPreStep_id _ Y⟩

variable {P F} in
/-- ★**(v)(c)** —— `φ = 𝟙 ≫ φ`。 -/
theorem pfRoot_preStepFactor' (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (φ : X ⟶ Y) (hφ : IsPreStep (pfRootPre P F) φ) :
    ∃ (M : PfRootObj P F) (β : X ⟶ M) (α : M ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (pfRootPre P F) β ∧ IsPreStep (pfRootPre P F) β ∧
        IsCoAngular (pfRootPre P F) α ∧ IsPreStep (pfRootPre P F) α :=
  ⟨X, 𝟙 X, φ, (Category.id_comp φ).symm, isIsometric_id _ X, isPreStep_id _ X,
    pfRoot_isCoAngular hfi φ, hφ⟩

variable {P F} in
/-- ★**(v)(b) の一意性** —— `α`・`α'` は isometric pre-step ゆえ同型。 -/
theorem pfRoot_preStepFactorUniq (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (M M' : PfRootObj P F) (β : X ⟶ M) (α : M ⟶ Y)
    (β' : X ⟶ M') (α' : M' ⟶ Y) (heq : β ≫ α = β' ≫ α')
    (_hβc : IsCoAngular (pfRootPre P F) β) (_hβs : IsPreStep (pfRootPre P F) β)
    (hαi : IsIsometric (pfRootPre P F) α) (hαs : IsPreStep (pfRootPre P F) α)
    (_hβ'c : IsCoAngular (pfRootPre P F) β') (_hβ's : IsPreStep (pfRootPre P F) β')
    (hα'i : IsIsometric (pfRootPre P F) α') (hα's : IsPreStep (pfRootPre P F) α') :
    ∃ γ : M ≅ M', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI : IsIso α := pfRoot_isOfIsotropicType (F := F) hfi M _ α hαi hαs
  haveI : IsIso α' := pfRoot_isOfIsotropicType (F := F) hfi M' _ α' hα'i hα's
  refine ⟨⟨α ≫ inv α', α' ≫ inv α, by simp, by simp⟩, by simp, ?_⟩
  show β' = β ≫ (α ≫ inv α')
  rw [← Category.assoc, heq, Category.assoc]
  simp

variable {P F} in
/-- ★**(v)(c) の一意性** —— `β`・`β'` は isometric pre-step ゆえ同型。 -/
theorem pfRoot_preStepFactorUniq' (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (M M' : PfRootObj P F) (β : X ⟶ M) (α : M ⟶ Y)
    (β' : X ⟶ M') (α' : M' ⟶ Y) (heq : β ≫ α = β' ≫ α')
    (hβi : IsIsometric (pfRootPre P F) β) (hβs : IsPreStep (pfRootPre P F) β)
    (_hαc : IsCoAngular (pfRootPre P F) α) (_hαs : IsPreStep (pfRootPre P F) α)
    (hβ'i : IsIsometric (pfRootPre P F) β') (hβ's : IsPreStep (pfRootPre P F) β')
    (_hα'c : IsCoAngular (pfRootPre P F) α') (_hα's : IsPreStep (pfRootPre P F) α') :
    ∃ γ : M ≅ M', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI : IsIso β := pfRoot_isOfIsotropicType (F := F) hfi X _ β hβi hβs
  haveI : IsIso β' := pfRoot_isOfIsotropicType (F := F) hfi X _ β' hβ'i hβ's
  refine ⟨⟨inv β ≫ β', inv β' ≫ β, by simp, by simp⟩, ?_, by simp⟩
  show α' = (inv β' ≫ β) ≫ α
  rw [Category.assoc, heq, ← Category.assoc]
  simp

/-! ## ★9. ★★★★合成の代表元 —— これが残り 13 条の共通の足場

★★`compRoot` の定義そのものが使う根の取り方
(`rfl` / `rfl`、`mul_comm` / `mul_comm`、`mul_comm` / `rfl`)で
`exists_rep3` を当てるだけでよい。★`compRoot_eq_lift` は要らない。 -/

variable {P F} in
/-- ★★★★**`𝒞^pf` の合成を `𝒞` の合成へ落とす**。 -/
theorem compRoot_rep {X Y Z : PfRootObj P F} (f : X ⟶ Y) (g : Y ⟶ Z) :
    ∃ (V : IdxPf3 P F (rtObj P F X.obj (Z.root * Y.root))
        (rtObj P F Y.obj (Z.root * X.root)) (rtObj P F Z.obj (Y.root * X.root)))
      (φ : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      (rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
          (show Z.root * X.root = Z.root * X.root from rfl)).inv f
        = HomPf.mk ((idx12 P F _ _ _).obj V) φ ∧
      (rtRootIso P F Y.obj Z.obj (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
          (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g
        = HomPf.mk ((idx23 P F _ _ _).obj V) ψ ∧
      compRoot P F f g
        = (rtRootIso P F X.obj Z.obj
            (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
            (show Y.root * X.root = Y.root * X.root from rfl)).hom
          (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)) := by
  obtain ⟨V, φ, ψ, hφ, hψ⟩ := exists_rep3 (P := P) (F := F)
    ((rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv f)
    ((rtRootIso P F Y.obj Z.obj (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).inv g)
  refine ⟨V, φ, ψ, hφ, hψ, ?_⟩
  show compRoot P F f g = _
  unfold compRoot
  rw [hφ, hψ, compPf_mk]

/-! ## ★10. (ii) —— 各次数の Frobenius 型射

★★`𝒞^pf` では `(A, m)` から `(A^n, m)` への次数 `n` の Frobenius 型射を取る。
★その代表元は `𝒞` の中の `rtObj A m ⟶ rtObj (rtObj A n) m` であり、
`frobDegSurj` ＋ `frobDegUniq` で作れる。 -/

variable {P F} in
/-- ★`rtObj A d` から `rtObj (rtObj A n) d` への次数 `n` の Frobenius 型射。 -/
theorem exists_rt_frob (A : C) (d n : ℕ+) :
    ∃ θ : rtObj P F A d ⟶ rtObj P F (rtObj P F A n) d,
      IsFrobeniusType P θ ∧ P.degFr θ = n := by
  obtain ⟨B', θ₀, hθ₀, hθ₀d⟩ := F.frobDegSurj (rtObj P F A d) n
  have h1 : IsFrobeniusType P (rtExt P F A d ≫ θ₀) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A d) hθ₀
  have h2 : IsFrobeniusType P (rtExt P F A n ≫ rtExt P F (rtObj P F A n) d) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A n) (rtExt_frobType P F (rtObj P F A n) d)
  have hdeg : P.degFr (rtExt P F A d ≫ θ₀)
      = P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) d) := by
    rw [P.degFr_comp, P.degFr_comp, hθ₀d, rtExt_degFr, rtExt_degFr, rtExt_degFr, mul_comm]
  obtain ⟨β, hβ, hβe⟩ := F.frobDegUniq A B' _ _ _ h1 h2 hdeg
  haveI : IsIso β := hβ
  refine ⟨θ₀ ≫ β, IsFrobeniusType.comp P F hθ₀ (isFrobeniusType_of_isIso P β), ?_⟩
  rw [P.degFr_comp, show P.degFr β = 1 from isLinear_of_isIso P β, one_mul, hθ₀d]

set_option maxHeartbeats 1600000 in
variable {P F} in
/-- ★★**(ii) の存在** —— `𝒞^pf` にも各次数の Frobenius 型射がある。 -/
theorem pfRoot_frobDegSurj (hfi : IsOfFrobeniusIsotropicType P) (X : PfRootObj P F) (n : ℕ+) :
    ∃ (Y : PfRootObj P F) (φ : X ⟶ Y),
      IsFrobeniusType (pfRootPre P F) φ ∧ (pfRootPre P F).degFr φ = n := by
  obtain ⟨θ, hθ, hθd⟩ := exists_rt_frob (F := F) X.obj X.root n
  refine ⟨⟨rtObj P F X.obj n, X.root⟩, ?_⟩
  refine ⟨HomPf.mk (idxOne P F (rtObj P F X.obj X.root)
    (rtObj P F (rtObj P F X.obj n) X.root)) θ, ?_⟩
  refine ⟨⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · exact pfRoot_isCoAngular hfi _
  · exact (isIsometric_mk_iff (X := X) (Y := ⟨rtObj P F X.obj n, X.root⟩)
      (idxOne P F (rtObj P F X.obj X.root)
        (rtObj P F (rtObj P F X.obj n) X.root)) θ).mpr hθ.1.2
  · exact (isBaseIsomorphism_mk_iff (X := X) (Y := ⟨rtObj P F X.obj n, X.root⟩)
      (idxOne P F (rtObj P F X.obj X.root)
        (rtObj P F (rtObj P F X.obj n) X.root)) θ).mpr hθ.2
  · exact (degFr_mk_iff (X := X) (Y := ⟨rtObj P F X.obj n, X.root⟩)
      (idxOne P F (rtObj P F X.obj X.root)
        (rtObj P F (rtObj P F X.obj n) X.root)) θ n).mpr hθd

/-! ## ★11. 遷移は pre-step を保つ -/

variable {P F} in
/-- ★**添字の遷移は pre-step を保つ**。 -/
theorem idxTransport_isPreStep {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsPreStep P φ) :
    IsPreStep P (idxTransport P F u φ) := by
  have hsq := idxTransport_spec (F := F) u φ
  obtain ⟨hua, hub, hud⟩ := u.right.property
  refine ⟨?_, ?_⟩
  · show P.degFr (idxTransport P F u φ) = 1
    exact (repDeg_map (F := F) u φ).trans h.1
  · haveI hia : IsIso (P.Base u.right.hom.1) := hua.2
    haveI hib : IsIso (P.Base u.right.hom.2) := hub.2
    haveI hif : IsIso (P.Base φ) := h.2
    have h1 : P.Base (φ ≫ u.right.hom.2)
        = P.Base (u.right.hom.1 ≫ idxTransport P F u φ) := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp] at h1
    haveI : IsIso (P.Base u.right.hom.1 ≫ P.Base (idxTransport P F u φ)) := by
      rw [← h1]; infer_instance
    exact IsIso.of_isIso_comp_left (P.Base u.right.hom.1) _

/-! ## ★12. ★★★`Hom^pf` の中で pre-step を消す

★★添字の等式は「ある上界で一致」としか言わないので、
その上界を **3 脚の像へ持ち上げ**(`idx13` の cofinal 性)、
`IdxPf3` の filtered 性で `V` との上界を取る。
★★`idx_hom_ext`(添字圏は細い)で 2 本の遷移射を同一視すれば、
`idxTransport_comp_pair` が使えて `𝒞` の `preStepMono` に落ちる。 -/

variable {P F} in
/-- ★★★**`Hom^pf` の中で pre-step は消せる**。 -/
theorem homPf_cancel_preStep {A B E : C} (V : IdxPf3 P F A B E)
    (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hψ : IsPreStep P ψ)
    (h : HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ' ≫ ψ)) :
    HomPf.mk ((idx12 P F A B E).obj V) φ = HomPf.mk ((idx12 P F A B E).obj V) φ' := by
  obtain ⟨V', t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  obtain ⟨V'', ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) V'
  set s : V ⟶ IsFiltered.max V V'' := IsFiltered.leftToMax V V'' with hs
  set r : V'' ⟶ IsFiltered.max V V'' := IsFiltered.rightToMax V V'' with hr
  have hm : t ≫ k ≫ (idx13 P F A B E).map r = (idx13 P F A B E).map s :=
    idx_hom_ext _ _
  have hA : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have hB : idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ' ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have key : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ) := by
    rw [hA, hB, ht]
  have e1 := idxTransport_comp_pair (F := F) s φ ψ
  have e2 := idxTransport_comp_pair (F := F) s φ' ψ
  have hcan : idxTransport P F ((idx12 P F A B E).map s) φ
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ
      = idxTransport P F ((idx12 P F A B E).map s) φ'
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ := by
    rw [e1, e2]
    exact key
  haveI : Mono (idxTransport P F ((idx23 P F A B E).map s) ψ) :=
    F.preStepMono _ (idxTransport_isPreStep (F := F) ((idx23 P F A B E).map s) ψ hψ)
  have hφφ' : idxTransport P F ((idx12 P F A B E).map s) φ
      = idxTransport P F ((idx12 P F A B E).map s) φ' :=
    (cancel_mono (idxTransport P F ((idx23 P F A B E).map s) ψ)).mp hcan
  rw [← HomPf.mk_map ((idx12 P F A B E).map s) φ, ← HomPf.mk_map ((idx12 P F A B E).map s) φ',
    hφφ']

/-! ## ★13. 平行 2 射 ＋ 1 射を共通の 3 脚添字へ -/

variable {P F} in
/-- ★**左側が平行**な版(`Prop32.lean` の `exists_rep3_pair` は右側が平行)。 -/
theorem exists_rep3_pairL {A B E : C} (u v : HomPf P F A B) (g : HomPf P F B E) :
    ∃ (V : IdxPf3 P F A B E) (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      u = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      v = HomPf.mk ((idx12 P F A B E).obj V) φ' ∧
      g = HomPf.mk ((idx23 P F A B E).obj V) ψ := by
  obtain ⟨V₁, φ₁, ψ₁, hu, hg₁⟩ := exists_rep3 (P := P) (F := F) u g
  obtain ⟨V₂, φ₂, ψ₂, hv, hg₂⟩ := exists_rep3 (P := P) (F := F) v g
  refine ⟨IsFiltered.max V₁ V₂,
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.leftToMax V₁ V₂)) φ₁,
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.rightToMax V₁ V₂)) φ₂,
    idxTransport P F ((idx23 P F A B E).map (IsFiltered.leftToMax V₁ V₂)) ψ₁,
    ?_, ?_, ?_⟩
  · rw [HomPf.mk_map]; exact hu
  · rw [HomPf.mk_map]; exact hv
  · rw [HomPf.mk_map]; exact hg₁

variable {P F} in
/-- ★★**`𝒞^pf` の合成(平行 2 射を揃えた版)**。 -/
theorem compRoot_rep_pairL {W X Y : PfRootObj P F} (u v : W ⟶ X) (f : X ⟶ Y) :
    ∃ (V : IdxPf3 P F (rtObj P F W.obj (Y.root * X.root))
        (rtObj P F X.obj (Y.root * W.root)) (rtObj P F Y.obj (X.root * W.root)))
      (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      u = (rtRootIso P F W.obj X.obj
          (show Y.root * X.root = Y.root * X.root from rfl)
          (show Y.root * W.root = Y.root * W.root from rfl)).hom
        (HomPf.mk ((idx12 P F _ _ _).obj V) φ) ∧
      v = (rtRootIso P F W.obj X.obj
          (show Y.root * X.root = Y.root * X.root from rfl)
          (show Y.root * W.root = Y.root * W.root from rfl)).hom
        (HomPf.mk ((idx12 P F _ _ _).obj V) φ') ∧
      f = (rtRootIso P F X.obj Y.obj
          (show Y.root * W.root = W.root * Y.root from mul_comm _ _)
          (show X.root * W.root = W.root * X.root from mul_comm _ _)).hom
        (HomPf.mk ((idx23 P F _ _ _).obj V) ψ) ∧
      compRoot P F u f
        = (rtRootIso P F W.obj Y.obj
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
            (show X.root * W.root = X.root * W.root from rfl)).hom
          (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)) ∧
      compRoot P F v f
        = (rtRootIso P F W.obj Y.obj
            (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
            (show X.root * W.root = X.root * W.root from rfl)).hom
          (HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ)) := by
  obtain ⟨V, φ, φ', ψ, hu, hv, hf⟩ := exists_rep3_pairL (P := P) (F := F)
    ((rtRootIso P F W.obj X.obj (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * W.root = Y.root * W.root from rfl)).inv u)
    ((rtRootIso P F W.obj X.obj (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * W.root = Y.root * W.root from rfl)).inv v)
    ((rtRootIso P F X.obj Y.obj (show Y.root * W.root = W.root * Y.root from mul_comm _ _)
      (show X.root * W.root = W.root * X.root from mul_comm _ _)).inv f)
  refine ⟨V, φ, φ', ψ, ?_, ?_, ?_, ?_, ?_⟩
  · rw [← hu, Iso.inv_hom_id_apply]
  · rw [← hv, Iso.inv_hom_id_apply]
  · rw [← hf, Iso.inv_hom_id_apply]
  · show compRoot P F u f = _
    unfold compRoot
    rw [hu, hf, compPf_mk]
  · show compRoot P F v f = _
    unfold compRoot
    rw [hv, hf, compPf_mk]

/-! ## ★14. (v)(a) —— `𝒞^pf` の pre-step は mono -/

variable {P F} in
/-- ★★★**(v)(a)** —— `𝒞^pf` の pre-step は monomorphism。 -/
theorem pfRoot_preStepMono {X Y : PfRootObj P F} (f : X ⟶ Y)
    (hf : IsPreStep (pfRootPre P F) f) : Mono f := by
  refine ⟨fun {W} u v huv => ?_⟩
  obtain ⟨V, φ, φ', ψ, hurep, hvrep, hfrep, hu, hv⟩ := compRoot_rep_pairL (F := F) u v f
  have hψ : IsPreStep P ψ := by
    refine (isPreStep_mk_iff (X := X) (Y := Y)
      ((pushIdx (F := F)
        (rtLift P F X.obj (show Y.root * W.root = W.root * Y.root from mul_comm _ _))
        (rtLift_frobType P F X.obj _)
        (rtLift P F Y.obj (show X.root * W.root = W.root * X.root from mul_comm _ _))
        (rtLift_frobType P F Y.obj _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx23 P F _ _ _).obj V)) ψ).mp ?_
    rw [← rtRootIso_hom_mk (F := F) X.obj Y.obj _ _ ((idx23 P F _ _ _).obj V) ψ, ← hfrep]
    exact hf
  have h0 : (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ))
      = (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ)) := by
    rw [← hu, ← hv]; exact huv
  have h1 : HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ) :=
    ((Iso.hom_inv_id_apply _ _).symm.trans (congrArg _ h0)).trans (Iso.hom_inv_id_apply _ _)
  rw [hurep, hvrep, homPf_cancel_preStep (F := F) V φ φ' ψ hψ h1]

/-! ## ★15. (i)(a) —— 底の同型類は Frobenius-trivial 対象の像

★★`𝒞^pf` の対象 `(A, n)` の底は `(P.toElem.obj A).base` そのもの
(`pfRootToElem` の `obj`)なので、★`𝒞` の `baseSurj` の証人 `A` を
`(A, 1)` として送るだけでよい。 -/

variable {P F} in
/-- ★`𝒞 → 𝒞^pf` は Frobenius 型射を Frobenius 型射に送る。 -/
theorem toPfRoot_isFrobeniusType (hfi : IsOfFrobeniusIsotropicType P) {A B : C} (φ : A ⟶ B)
    (hφ : IsFrobeniusType P φ) : IsFrobeniusType (pfRootPre P F) ((toPfRoot P F).map φ) := by
  refine ⟨⟨pfRoot_isCoAngular hfi _, ?_⟩, ?_⟩
  · show rootDiv (toRootHom (F := F) φ) = 0
    rw [rootDiv_toRootHom]
    show Pf.mk (P.Div φ) 1 = 0
    rw [show P.Div φ = 0 from hφ.1.2, Pf.mk_eq_zero_iff]
    exact ⟨1, by simp⟩
  · show IsIso (rootBase (toRootHom (F := F) φ))
    rw [rootBase_toRootHom]
    exact hφ.2

variable {P F} in
/-- ★`𝒞 → 𝒞^pf` は Frobenius 次数を保つ。 -/
theorem toPfRoot_degFr {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).degFr ((toPfRoot P F).map φ) = P.degFr φ :=
  rootDeg_toRootHom (F := F) φ

variable {P F} in
/-- ★`𝒞 → 𝒞^pf` は base-identity を保つ。 -/
theorem toPfRoot_isBaseIdentity {A : C} (φ : End A) (hφ : IsBaseIdentity P φ) :
    IsBaseIdentity (pfRootPre P F) ((toPfRoot P F).map (φ : A ⟶ A)) := by
  show rootBase (toRootHom (F := F) (φ : A ⟶ A))
    = (pfRootPre P F).Base (𝟙 ((toPfRoot P F).obj A))
  rw [rootBase_toRootHom, ← (toPfRoot P F).map_id A]
  show P.Base (φ : A ⟶ A) = rootBase (toRootHom (F := F) (𝟙 A))
  rw [rootBase_toRootHom]
  exact hφ

variable {P F} in
/-- ★★**(i)(a)** —— `𝒞^pf` の底の同型類も Frobenius-trivial 対象の像。 -/
theorem pfRoot_baseSurj (hfi : IsOfFrobeniusIsotropicType P) (Y : D) :
    ∃ X : PfRootObj P F, IsFrobeniusTrivial (pfRootPre P F) X ∧
      Nonempty (((pfRootPre P F).toElem.obj X).base ≅ Y) := by
  obtain ⟨A, ⟨ζ, hdeg, hbf⟩, ⟨e⟩⟩ := F.baseSurj Y
  refine ⟨(toPfRoot P F).obj A, ⟨{
      toFun := fun n => (toPfRoot P F).map ((ζ n : A ⟶ A))
      map_one' := by
        show (toPfRoot P F).map ((ζ 1 : End A) : A ⟶ A) = 𝟙 _
        rw [show ((ζ 1 : End A) : A ⟶ A) = 𝟙 A from
          congrArg (fun t : End A => (t : A ⟶ A)) ζ.map_one]
        exact (toPfRoot P F).map_id A
      map_mul' := fun x y => by
        show (toPfRoot P F).map ((ζ (x * y) : End A) : A ⟶ A)
          = (toPfRoot P F).map ((ζ y : End A) : A ⟶ A)
            ≫ (toPfRoot P F).map ((ζ x : End A) : A ⟶ A)
        rw [show ((ζ (x * y) : End A) : A ⟶ A)
          = ((ζ y : End A) : A ⟶ A) ≫ ((ζ x : End A) : A ⟶ A) from
            congrArg (fun t : End A => (t : A ⟶ A)) (ζ.map_mul x y)]
        exact (toPfRoot P F).map_comp _ _ }, ?_, ?_⟩, ⟨e⟩⟩
  · intro n
    exact (toPfRoot_degFr (F := F) ((ζ n : End A) : A ⟶ A)).trans (hdeg n)
  · intro n
    exact ⟨toPfRoot_isBaseIdentity (F := F) (ζ n) (hbf n).1,
      toPfRoot_isFrobeniusType hfi _ (hbf n).2⟩

/-! ## ★16. 根を上げる道具

★★`𝒞` の射を「両側の `k` 乗根の間」へ運ぶ。★添字圏の遷移として書けば
`idxTransport_isPreStep`(と `repDeg_map`)がそのまま使える。 -/

variable {P F} in
/-- ★`k` 乗根へ上げる添字。 -/
noncomputable def idxPow (U V : C) (k : ℕ+) : IdxPf P F U V :=
  idxMk (P := P) (F := F) (rtExt P F U k) (rtExt P F V k)
    (rtExt_frobType P F U k) (rtExt_frobType P F V k)
    (by rw [rtExt_degFr, rtExt_degFr])

variable {P F} in
/-- ★始対象からの遷移射。 -/
noncomputable def idxPowHom (U V : C) (k : ℕ+) : idxOne P F U V ⟶ idxPow (F := F) U V k :=
  Under.homMk (show (⟨(U, V)⟩ : BiFr P F) ⟶ (⟨(rtObj P F U k, rtObj P F V k)⟩ : BiFr P F) from
    ⟨(rtExt P F U k, rtExt P F V k), rtExt_frobType P F U k, rtExt_frobType P F V k,
      by rw [rtExt_degFr, rtExt_degFr]⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))

variable {P F} in
/-- ★★**`k` 乗へ運んだ射**。 -/
noncomputable def liftPow {U V : C} (p : U ⟶ V) (k : ℕ+) :
    rtObj P F U k ⟶ rtObj P F V k :=
  idxTransport P F (idxPowHom (F := F) U V k) p

variable {P F} in
/-- ★四角形。 -/
theorem liftPow_spec {U V : C} (p : U ⟶ V) (k : ℕ+) :
    p ≫ rtExt P F V k = rtExt P F U k ≫ liftPow (F := F) p k :=
  idxTransport_spec (F := F) (idxPowHom (F := F) U V k) p

variable {P F} in
/-- ★pre-step は保たれる。 -/
theorem liftPow_isPreStep {U V : C} (p : U ⟶ V) (k : ℕ+) (hp : IsPreStep P p) :
    IsPreStep P (liftPow (F := F) p k) :=
  idxTransport_isPreStep (F := F) (idxPowHom (F := F) U V k) p hp

variable {P F} in
/-- ★★**`rtObj (rtObj A m) n` と `rtObj A (n * m)` を同一視する**。 -/
theorem exists_rtObj_assoc (A : C) (m n t : ℕ+) (ht : t = n * m) :
    ∃ β : rtObj P F (rtObj P F A m) n ⟶ rtObj P F A t,
      IsIso β ∧ rtExt P F A m ≫ rtExt P F (rtObj P F A m) n ≫ β = rtExt P F A t := by
  have h1 : IsFrobeniusType P (rtExt P F A m ≫ rtExt P F (rtObj P F A m) n) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A m) (rtExt_frobType P F (rtObj P F A m) n)
  have h2 : IsFrobeniusType P (rtExt P F A t) := rtExt_frobType P F A t
  have hdeg : P.degFr (rtExt P F A m ≫ rtExt P F (rtObj P F A m) n)
      = P.degFr (rtExt P F A t) := by
    rw [P.degFr_comp, rtExt_degFr, rtExt_degFr, rtExt_degFr, ht]
  obtain ⟨β, hβ, hβe⟩ := F.frobDegUniq A _ _ _ _ h1 h2 hdeg
  exact ⟨β, hβ, (Category.assoc _ _ _).symm.trans hβe⟩

/-! ## ★17. (i)(b) —— 底の同型は pre-step の span で持ち上がる

★★根の取り方が要点である。`X = (A, n)`、`Y = (B, m)` に対し
**`𝒞` の span を `A^m` と `B^n` の間で取り**、その頂点 `V₀` を
`W := (V₀, n * m)` として使う。
★モデルで言えば `W = V₀/(nm) ≤ A/n ⟺ V₀ ≤ mA` であり、
これがちょうど `V₀ ⟶ A^m` という pre-step である。 -/

variable {P F} in
/-- ★始対象の添字での `repBase` は `Base` そのもの。 -/
theorem repBase_idxOne {A B : C} (φ : A ⟶ B) :
    repBase (idxOne P F A B) φ = P.Base φ := by
  have h := repBase_spec (F := F) (idxOne P F A B) φ
  show repBase (idxOne P F A B) φ = P.Base φ
  have h2 : repBase (idxOne P F A B) φ ≫ P.Base (𝟙 B) = P.Base (𝟙 A) ≫ P.Base φ := h
  rw [P.Base_id, P.Base_id, Category.comp_id, Category.id_comp] at h2
  exact h2

/-- ★`f ≫ α = g` から `α = f⁻¹ ≫ g`。 -/
theorem eq_inv_comp_of {W X Y : D} (f : W ⟶ X) (hf : IsIso f) (g : W ⟶ Y) (α : X ⟶ Y)
    (h : f ≫ α = g) : α = @inv _ _ _ _ f hf ≫ g := by
  haveI := hf
  rw [← h, IsIso.inv_hom_id_assoc]

variable {P F} in
/-- ★★**(i)(b)** —— `𝒞^pf` 版。 -/
theorem pfRoot_preStepSpan (X Y : PfRootObj P F)
    (α : ((pfRootPre P F).toElem.obj X).base ⟶ ((pfRootPre P F).toElem.obj Y).base)
    (hα : IsIso α) :
    ∃ (W : PfRootObj P F) (φ : W ⟶ X) (ψ : W ⟶ Y) (hφ : IsPreStep (pfRootPre P F) φ),
      IsPreStep (pfRootPre P F) ψ ∧
        α = @inv _ _ _ _ ((pfRootPre P F).Base φ) hφ.2 ≫ (pfRootPre P F).Base ψ := by
  haveI := hα
  haveI hieA : IsIso (P.Base (rtExt P F X.obj Y.root)) := (rtExt_frobType P F X.obj Y.root).2
  haveI hieB : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  haveI hinv : IsIso (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
      ≫ P.Base (rtExt P F Y.obj X.root)) :=
    IsIso.comp_isIso' (IsIso.inv_isIso) (IsIso.comp_isIso' hα hieB)
  obtain ⟨V₀, p, q, hp, hq, hspan⟩ := F.preStepSpan (rtObj P F X.obj Y.root)
    (rtObj P F Y.obj X.root)
    (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
      ≫ P.Base (rtExt P F Y.obj X.root)) hinv
  haveI hip : IsIso (P.Base p) := hp.2
  haveI hiq : IsIso (P.Base q) := hq.2
  obtain ⟨βA, hβA, hβAe⟩ := exists_rtObj_assoc (F := F) X.obj Y.root X.root (X.root * Y.root) rfl
  obtain ⟨βB, hβB, hβBe⟩ := exists_rtObj_assoc (F := F) Y.obj X.root Y.root (X.root * Y.root)
    (mul_comm _ _)
  haveI := hβA
  haveI := hβB
  have hφs : IsPreStep P (liftPow (F := F) p X.root ≫ βA) :=
    IsPreStep.comp P (liftPow_isPreStep (F := F) p X.root hp) (isPreStep_of_isIso P βA)
  have hψs : IsPreStep P (liftPow (F := F) q Y.root ≫ βB) :=
    IsPreStep.comp P (liftPow_isPreStep (F := F) q Y.root hq) (isPreStep_of_isIso P βB)
  have hbase : ∀ {U V : C} (r : V₀ ⟶ U) (k t : ℕ+)
      (L : rtObj P F V₀ k ⟶ rtObj P F U k)
      (_hL : r ≫ rtExt P F U k = rtExt P F V₀ k ≫ L)
      (γ : rtObj P F U k ⟶ rtObj P F V t) (eU : V ⟶ U) (heU : IsIso (P.Base eU))
      (_hcomp : eU ≫ rtExt P F U k ≫ γ = rtExt P F V t),
      rootBase (show HomRoot P F ⟨V₀, t⟩ ⟨V, k⟩ from
          HomPf.mk (idxOne P F (rtObj P F V₀ k) (rtObj P F V t)) (L ≫ γ))
        = P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU := by
    intro U V r k t L hL γ eU heU hcomp
    haveI := heU
    have hmor : rtExt P F V₀ k ≫ (L ≫ γ) = r ≫ (rtExt P F U k ≫ γ) := by
      rw [← Category.assoc, ← hL, Category.assoc]
    refine (rootBase_uniq _ _ ?_).symm
    show (P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU) ≫ P.Base (rtExt P F V t)
      = P.Base (rtExt P F V₀ k) ≫ pfBase (HomPf.mk (idxOne P F _ _) _)
    rw [pfBase_mk, repBase_idxOne, ← hcomp]
    have hRHS : P.Base (rtExt P F V₀ k) ≫ P.Base (L ≫ γ)
        = P.Base r ≫ P.Base (rtExt P F U k ≫ γ) :=
      (P.Base_comp _ _).symm.trans ((congrArg P.Base hmor).trans (P.Base_comp _ _))
    have hLHS : (P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU)
          ≫ P.Base (eU ≫ (rtExt P F U k ≫ γ))
        = P.Base r ≫ P.Base (rtExt P F U k ≫ γ) := by
      rw [P.Base_comp eU (rtExt P F U k ≫ γ)]
      simp
    exact hLHS.trans hRHS.symm
  have hbφ := hbase p X.root (X.root * Y.root) (liftPow (F := F) p X.root)
    (liftPow_spec (F := F) p X.root) βA (rtExt P F X.obj Y.root) hieA hβAe
  have hbψ := hbase q Y.root (X.root * Y.root) (liftPow (F := F) q Y.root)
    (liftPow_spec (F := F) q Y.root) βB (rtExt P F Y.obj X.root) hieB hβBe
  refine ⟨⟨V₀, X.root * Y.root⟩,
    HomPf.mk (idxOne P F _ _) (liftPow (F := F) p X.root ≫ βA),
    HomPf.mk (idxOne P F _ _) (liftPow (F := F) q Y.root ≫ βB),
    (isPreStep_mk_iff (X := (⟨V₀, X.root * Y.root⟩ : PfRootObj P F)) (Y := X) _ _).mpr hφs,
    (isPreStep_mk_iff (X := (⟨V₀, X.root * Y.root⟩ : PfRootObj P F)) (Y := Y) _ _).mpr hψs, ?_⟩
  have hkey : rootBase (show HomRoot P F ⟨V₀, X.root * Y.root⟩ X from
        HomPf.mk (idxOne P F _ _) (liftPow (F := F) p X.root ≫ βA)) ≫ α
      = rootBase (show HomRoot P F ⟨V₀, X.root * Y.root⟩ Y from
        HomPf.mk (idxOne P F _ _) (liftPow (F := F) q Y.root ≫ βB)) := by
    have h1 : P.Base p ≫ (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
        ≫ P.Base (rtExt P F Y.obj X.root)) = P.Base q := by
      refine (congrArg (fun t => P.Base p ≫ t) hspan).trans ?_
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    have hcancel : P.Base (rtExt P F Y.obj X.root)
        ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) hieB = 𝟙 _ :=
      @IsIso.hom_inv_id _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) hieB
    rw [hbφ, hbψ, ← h1]
    simp only [Category.assoc]
    refine congrArg (fun t => P.Base p
      ≫ (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ t)) ?_
    exact (Category.comp_id α).symm.trans (congrArg (fun t => α ≫ t) hcancel.symm)
  exact eq_inv_comp_of _ _ _ _ hkey

/-! ## ★18. 遷移は同型を保つ

★★逆射の遷移が逆射になる。★示すのは `a ≫ (T₁ ≫ T₂) = a ≫ 𝟙` で、
`a` が epi(`𝒞` は totally epimorphic)なので消せる。 -/

variable {P F} in
/-- ★★**遷移は同型を保つ**。 -/
theorem frobTransport_isIso {A' B' A'' B'' : C}
    (a : A' ⟶ A'') (ha : IsFrobeniusType P a) (b : B' ⟶ B'') (hb : IsFrobeniusType P b)
    (hd : P.degFr a = P.degFr b) (θ : A' ⟶ B') (hθ : IsIso θ) :
    IsIso (frobTransport (F := F) a ha b hb hd θ) := by
  haveI := hθ
  haveI hea : Epi a := P.totEpiC _ _ _
  haveI heb : Epi b := P.totEpiC _ _ _
  have s1 : θ ≫ b = a ≫ frobTransport (F := F) a ha b hb hd θ :=
    frobTransport_spec _ _ _ _ _ θ
  have s2 : inv θ ≫ a = b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) :=
    frobTransport_spec _ _ _ _ _ (inv θ)
  refine ⟨frobTransport (F := F) b hb a ha hd.symm (inv θ), ?_, ?_⟩
  · refine (cancel_epi a).mp ?_
    have e1 : a ≫ (frobTransport (F := F) a ha b hb hd θ
        ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
        = (a ≫ frobTransport (F := F) a ha b hb hd θ)
          ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) := (Category.assoc _ _ _).symm
    have e2 : (a ≫ frobTransport (F := F) a ha b hb hd θ)
          ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)
        = (θ ≫ b) ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) :=
      congrArg (fun t => t ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)) s1.symm
    have e3 : (θ ≫ b) ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)
        = θ ≫ (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)) := Category.assoc _ _ _
    have e4 : θ ≫ (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
        = θ ≫ (inv θ ≫ a) := congrArg (fun t => θ ≫ t) s2.symm
    have e5 : θ ≫ (inv θ ≫ a) = a := by
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    exact ((((e1.trans e2).trans e3).trans e4).trans e5).trans (Category.comp_id a).symm
  · refine (cancel_epi b).mp ?_
    have e1 : b ≫ (frobTransport (F := F) b hb a ha hd.symm (inv θ)
        ≫ frobTransport (F := F) a ha b hb hd θ)
        = (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
          ≫ frobTransport (F := F) a ha b hb hd θ := (Category.assoc _ _ _).symm
    have e2 : (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
          ≫ frobTransport (F := F) a ha b hb hd θ
        = (inv θ ≫ a) ≫ frobTransport (F := F) a ha b hb hd θ :=
      congrArg (fun t => t ≫ frobTransport (F := F) a ha b hb hd θ) s2.symm
    have e3 : (inv θ ≫ a) ≫ frobTransport (F := F) a ha b hb hd θ
        = inv θ ≫ (a ≫ frobTransport (F := F) a ha b hb hd θ) := Category.assoc _ _ _
    have e4 : inv θ ≫ (a ≫ frobTransport (F := F) a ha b hb hd θ)
        = inv θ ≫ (θ ≫ b) := congrArg (fun t => inv θ ≫ t) s1.symm
    have e5 : inv θ ≫ (θ ≫ b) = b := by
      rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    exact ((((e1.trans e2).trans e3).trans e4).trans e5).trans (Category.comp_id b).symm

variable {P F} in
/-- ★添字の遷移は同型を保つ。 -/
theorem idxTransport_isIso {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (θ : Z.right.obj.1 ⟶ Z.right.obj.2) (hθ : IsIso θ) :
    IsIso (idxTransport P F u θ) :=
  frobTransport_isIso _ _ _ _ _ θ hθ

/-! ## ★19. ★★代表元が同型なら同型(添字を持ち上げない形)

★★`pfRoot_isIso_of_rep` は `compRoot` の定義に合わせた「持ち上げた添字」の形だが、
★実際に使うのは**素の添字**の形である。`pushIdx` の cofinal 性で移す。 -/

variable {P F} in
/-- ★★★**素の添字で代表元が同型なら `𝒞^pf` でも同型**。 -/
theorem pfRoot_isIso_mk {Y Z : PfRootObj P F}
    (W : IdxPf P F (rtObj P F Y.obj Z.root) (rtObj P F Z.obj Y.root))
    (β₀ : W.right.obj.1 ⟶ W.right.obj.2) (hβ₀ : IsIso β₀) :
    IsIso (show Y ⟶ Z from HomPf.mk W β₀) := by
  obtain ⟨V, ⟨u⟩⟩ := exists_hom_of_final
    (pushIdx (F := F)
      (rtLift P F Y.obj (show Y.root * Z.root = Y.root * Z.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (rtLift P F Z.obj (show Y.root * Y.root = Y.root * Y.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])) W
  set β₁ := idxTransport P F u β₀ with hβ₁
  have hiso : IsIso β₁ := by rw [hβ₁]; exact idxTransport_isIso (F := F) u β₀ hβ₀
  have h1 : (show Y ⟶ Z from HomPf.mk W β₀)
      = (rtRootIso P F Y.obj Z.obj (show Y.root * Z.root = Y.root * Z.root from rfl)
          (show Y.root * Y.root = Y.root * Y.root from rfl)).hom
        (HomPf.mk V β₁) := by
    rw [rtRootIso_hom_mk, hβ₁]
    exact (HomPf.mk_map (P := P) (F := F) u β₀).symm
  refine pfRoot_isIso_of_rep (show Y ⟶ Z from HomPf.mk W β₀) V β₁ (hφ := hiso) ?_
  rw [h1, Iso.hom_inv_id_apply]

/-! ## ★20. 同じ始域から出る 2 射を揃える(co-span)

★★`exists_rep3` は「合成できる 2 射」を揃えるが、
`frobDegUniq` に要るのは**同じ始域から出る 2 射**である。
★3 脚添字の `idx12` と `idx13` を使えば同じ形で作れる。 -/

variable {P F} in
/-- ★**co-span を共通の 3 脚添字へ**。 -/
theorem exists_rep_cospan {A B E : C} (f : HomPf P F A B) (g : HomPf P F A E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (χ : V.right.obj.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx13 P F A B E).obj V) χ := by
  obtain ⟨Zf, φ₀, hf⟩ := HomPf.exists_rep (P := P) (F := F) f
  obtain ⟨Zg, χ₀, hg⟩ := HomPf.exists_rep (P := P) (F := F) g
  obtain ⟨hfa, hfb, hfab⟩ := Zf.hom.property
  obtain ⟨hga, hge, hgae⟩ := Zg.hom.property
  obtain ⟨E₁, e₁, he₁, he₁d⟩ := F.frobDegSurj E (P.degFr Zf.hom.hom.1)
  obtain ⟨B₂, b₂, hb₂, hb₂d⟩ := F.frobDegSurj B (P.degFr Zg.hom.hom.1)
  refine ⟨IsFiltered.max
      (Under.mk (Y := (⟨(Zf.right.obj.1, Zf.right.obj.2, E₁)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(Zf.hom.hom.1, Zf.hom.hom.2, e₁), hfa, hfb, he₁, hfab,
            hfab.symm.trans he₁d.symm⟩))
      (Under.mk (Y := (⟨(Zg.right.obj.1, B₂, Zg.right.obj.2)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(Zg.hom.hom.1, b₂, Zg.hom.hom.2), hga, hb₂, hge, hb₂d.symm,
            hb₂d.trans hgae⟩)),
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.leftToMax _ _)) φ₀,
    idxTransport P F ((idx13 P F A B E).map (IsFiltered.rightToMax _ _)) χ₀, ?_, ?_⟩
  · rw [HomPf.mk_map]; exact hf.symm
  · rw [HomPf.mk_map]; exact hg.symm

variable {P F} in
/-- ★3 脚添字も「第 1 脚が isotropic」な所まで押し上げられる。 -/
theorem exists_idx3_isotropic (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E) :
    ∃ (W : IdxPf3 P F A B E) (u : V ⟶ W), IsIsotropic P W.right.obj.1 := by
  obtain ⟨Dd, a, ha, hDd⟩ := hfi V.right.obj.1
  obtain ⟨B₂, b, hb, hbd⟩ := F.frobDegSurj V.right.obj.2.1 (P.degFr a)
  obtain ⟨E₂, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr a)
  obtain ⟨hva, hvb, hve, hvab, hvbe⟩ := V.hom.property
  refine ⟨Under.mk (Y := (⟨(Dd, B₂, E₂)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(V.hom.hom.1 ≫ a, V.hom.hom.2.1 ≫ b, V.hom.hom.2.2 ≫ e),
          IsFrobeniusType.comp P F hva ha, IsFrobeniusType.comp P F hvb hb,
          IsFrobeniusType.comp P F hve he, ?_, ?_⟩),
    Under.homMk (show V.right ⟶ (⟨(Dd, B₂, E₂)⟩ : TriFr P F) from
      ⟨(a, b, e), ha, hb, he, hbd.symm, hbd.trans hed.symm⟩)
      (WideSubcategory.hom_ext _ rfl), hDd⟩
  · rw [P.degFr_comp, P.degFr_comp, hbd, hvab]
  · rw [P.degFr_comp, P.degFr_comp, hbd, hed, hvbe]

/-! ## ★21. (ii) の一意性 —— `frobDegUniq`

★★`φ`・`ψ` を **`compRoot` が使う根の高さ**へ持ち上げて co-span を揃え、
添字を isotropic まで押し上げてから `𝒞` の `frobDegUniq` を当てる。
★得た同型 `β₀` を `idx23` の添字で戻せば、合成則は `compPf_mk` そのもの。 -/

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★**(ii) の一意性** —— `𝒞^pf` 版。 -/
theorem pfRoot_frobDegUniq (hfi : IsOfFrobeniusIsotropicType P)
    (X Y Z : PfRootObj P F) (φ : X ⟶ Y) (ψ : X ⟶ Z)
    (hφ : IsFrobeniusType (pfRootPre P F) φ) (hψ : IsFrobeniusType (pfRootPre P F) ψ)
    (hdeg : (pfRootPre P F).degFr φ = (pfRootPre P F).degFr ψ) :
    ∃ β : Y ⟶ Z, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨V, φ₀, χ₀, hφ0, hχ0⟩ := exists_rep_cospan (P := P) (F := F)
    ((rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
      (show Z.root * X.root = Z.root * X.root from rfl)).inv φ)
    ((rtRootIso P F X.obj Z.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = Y.root * X.root from rfl)).inv ψ)
  obtain ⟨W, u, hW⟩ := exists_idx3_isotropic (F := F) hfi V
  set φ₁ := idxTransport P F ((idx12 P F _ _ _).map u) φ₀ with hφ₁
  set χ₁ := idxTransport P F ((idx13 P F _ _ _).map u) χ₀ with hχ₁
  have hφW : (rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv φ
      = HomPf.mk ((idx12 P F _ _ _).obj W) φ₁ := by
    rw [hφ₁, HomPf.mk_map]; exact hφ0
  have hχW : (rtRootIso P F X.obj Z.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = Y.root * X.root from rfl)).inv ψ
      = HomPf.mk ((idx13 P F _ _ _).obj W) χ₁ := by
    rw [hχ₁, HomPf.mk_map]; exact hχ0
  -- ★`φ`・`ψ` を押し出した添字の代表元として書き直す
  have hφmk : φ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Z.root * Y.root from rfl))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Y.obj (show Z.root * X.root = Z.root * X.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx12 P F _ _ _).obj W)) φ₁ := by
    rw [← rtRootIso_hom_mk (F := F) X.obj Y.obj _ _ ((idx12 P F _ _ _).obj W) φ₁,
      ← hφW, Iso.inv_hom_id_apply]
  have hψmk : ψ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Z.obj (show Y.root * X.root = Y.root * X.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx13 P F _ _ _).obj W)) χ₁ := by
    rw [← rtRootIso_hom_mk (F := F) X.obj Z.obj _ _ ((idx13 P F _ _ _).obj W) χ₁,
      ← hχW, Iso.inv_hom_id_apply]
  set Wφ : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root) :=
    (pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Z.root * Y.root from rfl))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Y.obj (show Z.root * X.root = Z.root * X.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx12 P F _ _ _).obj W) with hWφ
  set Wψ : IdxPf P F (rtObj P F X.obj Z.root) (rtObj P F Z.obj X.root) :=
    (pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Z.obj (show Y.root * X.root = Y.root * X.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx13 P F _ _ _).obj W) with hWψ
  -- ★代表元は `𝒞` の Frobenius 型
  have hco : ∀ {U : C} (t : W.right.obj.1 ⟶ U), IsCoAngular P t :=
    fun t => prop_1_4_i P t (fun _ g => F.isotropicClosed g hW)
  have hφ₁F : IsFrobeniusType P φ₁ := by
    refine ⟨⟨hco φ₁, ?_⟩, ?_⟩
    · refine (isIsometric_mk_iff (X := X) (Y := Y) Wφ φ₁).mp ?_
      rw [← hφmk]; exact hφ.1.2
    · refine (isBaseIsomorphism_mk_iff (X := X) (Y := Y) Wφ φ₁).mp ?_
      rw [← hφmk]; exact hφ.2
  have hχ₁F : IsFrobeniusType P χ₁ := by
    refine ⟨⟨hco χ₁, ?_⟩, ?_⟩
    · refine (isIsometric_mk_iff (X := X) (Y := Z) Wψ χ₁).mp ?_
      rw [← hψmk]; exact hψ.1.2
    · refine (isBaseIsomorphism_mk_iff (X := X) (Y := Z) Wψ χ₁).mp ?_
      rw [← hψmk]; exact hψ.2
  have hd1 : P.degFr φ₁ = P.degFr χ₁ := by
    have e1 : (pfRootPre P F).degFr φ = P.degFr φ₁ := by
      rw [hφmk]
      exact (degFr_mk_iff (X := X) (Y := Y) Wφ φ₁ (P.degFr φ₁)).mpr rfl
    have e2 : (pfRootPre P F).degFr ψ = P.degFr χ₁ := by
      rw [hψmk]
      exact (degFr_mk_iff (X := X) (Y := Z) Wψ χ₁ (P.degFr χ₁)).mpr rfl
    rw [← e1, ← e2, hdeg]
  obtain ⟨β₀, hβ₀, hβ₀e⟩ := F.frobDegUniq _ _ _ φ₁ χ₁ hφ₁F hχ₁F hd1
  -- ★戻す
  refine ⟨(rtRootIso P F Y.obj Z.obj
      (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).hom
    (HomPf.mk ((idx23 P F _ _ _).obj W) β₀), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact pfRoot_isIso_mk _ _ hβ₀
  · show compRoot P F φ _ = ψ
    unfold compRoot
    rw [hφW, Iso.hom_inv_id_apply, compPf_mk]
    refine Eq.trans (congrArg (ConcreteCategory.hom (rtRootIso P F X.obj Z.obj
        (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = Y.root * X.root from rfl)).hom)
      (congrArg (HomPf.mk ((idx13 P F _ _ _).obj W)) hβ₀e)) ?_
    rw [← hχW, Iso.inv_hom_id_apply]

/-! ## ★22. ★★★根の指数は上げられる —— `(A, n) ≅ (A^m, n*m)`

★★これで **2 つの対象を同じ根の指数に揃えられる**。
`X = (A, n)`、`Y = (B, m)` なら `X ≅ (A^m, n*m)`、`Y ≅ (B^n, n*m)` である。

★★★根が等しいと `compRoot` の 3 つの `rtRootIso` が**すべて同じ命題の証明**になり
(`r*r = r*r`)、証明無関係で**同一の射**になる。★これが (iii)(c) や (vi) の鍵。 -/

variable {P F} in
/-- ★★**根の指数を上げる同型** —— `(A, n) ≅ (A^m, t)`(`t = n * m`)。 -/
theorem pfRoot_exists_iso_root (A : C) (n m t : ℕ+) (ht : t = n * m) :
    ∃ e : (⟨A, n⟩ : PfRootObj P F) ⟶ ⟨rtObj P F A m, t⟩, IsIso e := by
  obtain ⟨β, hβ, -⟩ := exists_rtObj_assoc (F := F) A m n t ht
  haveI := hβ
  exact ⟨HomPf.mk (idxOne P F _ _) (@inv _ _ _ _ β hβ),
    pfRoot_isIso_mk _ _ (@IsIso.inv_isIso _ _ _ _ β hβ)⟩

/-! ## ★23. `𝒪^▷` と (iii)(c) は同型で移せる(圏の言葉だけ) -/

variable {P F} in
/-- ★同型で挟むと `𝒪^▷` の元は `𝒪^▷` の元に移る。 -/
theorem oTri_conj {X X' : PfRootObj P F} (e : X ⟶ X') (e' : X' ⟶ X)
    (he1 : e ≫ e' = 𝟙 X) (he2 : e' ≫ e = 𝟙 X') (α : End X)
    (hα : α ∈ OTri (pfRootPre P F) X) :
    (show End X' from e' ≫ (α : X ⟶ X) ≫ e) ∈ OTri (pfRootPre P F) X' := by
  haveI : IsIso e := ⟨⟨e', he1, he2⟩⟩
  haveI : IsIso e' := ⟨⟨e, he2, he1⟩⟩
  constructor
  · show (pfRootPre P F).Base (e' ≫ (α : X ⟶ X) ≫ e)
      = (pfRootPre P F).Base (𝟙 X')
    rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp,
      show (pfRootPre P F).Base (α : X ⟶ X) = (pfRootPre P F).Base (𝟙 X) from hα.1,
      (pfRootPre P F).Base_id, Category.id_comp, (pfRootPre P F).Base_id,
      ← (pfRootPre P F).Base_comp, he2, (pfRootPre P F).Base_id]
  · show (pfRootPre P F).degFr (e' ≫ (α : X ⟶ X) ≫ e) = 1
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (α : X ⟶ X) = 1 from hα.2,
      show (pfRootPre P F).degFr e = 1 from isLinear_of_isIso (pfRootPre P F) e,
      show (pfRootPre P F).degFr e' = 1 from isLinear_of_isIso (pfRootPre P F) e']
    simp

/-! ## ★24. 3 脚添字を「対角かつ isotropic」へ

★★`𝒪^▷` の元を `𝒞` の `𝒪^▷` の元として読むには、
**第 1・第 2 脚が同じ射**でなければならない(そうでないと
`Base` が恒等にならず、共役が残る)。
★`frob_common_upper` で 2 脚を揃え、そのうえで `hfi` で isotropic まで押す。 -/

variable {P F} in
/-- ★★**3 脚添字を「第 1・第 2 脚が一致し、終域が isotropic」な所まで押し上げる**。 -/
theorem exists_idx3_diag (hfi : IsOfFrobeniusIsotropicType P) {A E : C}
    (V : IdxPf3 P F A A E) :
    ∃ (Dd E₃ : C) (l : A ⟶ Dd) (hl : IsFrobeniusType P l) (m : E ⟶ E₃)
      (hm : IsFrobeniusType P m) (hd : P.degFr l = P.degFr m),
      IsIsotropic P Dd ∧ Nonempty (V ⟶ idxMk3 (F := F) l l m hl hl hm rfl hd) := by
  obtain ⟨hv1, hv2, hv3, h12, h23⟩ := V.hom.property
  obtain ⟨X₃, a, c, ha, hc, had, hcd, hac⟩ :=
    frob_common_upper P F V.hom.hom.1 hv1 V.hom.hom.2.1 hv2
  obtain ⟨Dd, δ, hδ, hDd⟩ := hfi X₃
  obtain ⟨E₃, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr (a ≫ δ))
  have hacδ : V.hom.hom.1 ≫ (a ≫ δ) = V.hom.hom.2.1 ≫ (c ≫ δ) := by
    rw [← Category.assoc, ← Category.assoc, hac]
    rfl
  have hdaδ : P.degFr (a ≫ δ) = P.degFr (c ≫ δ) := by
    rw [P.degFr_comp a δ, P.degFr_comp c δ, had, hcd, h12]
    rfl
  have hlm : P.degFr (V.hom.hom.1 ≫ (a ≫ δ)) = P.degFr (V.hom.hom.2.2 ≫ e) := by
    rw [P.degFr_comp V.hom.hom.1 (a ≫ δ), P.degFr_comp V.hom.hom.2.2 e, hed, h12, h23]
  refine ⟨Dd, E₃, V.hom.hom.1 ≫ (a ≫ δ),
    IsFrobeniusType.comp P F hv1 (IsFrobeniusType.comp P F ha hδ),
    V.hom.hom.2.2 ≫ e, IsFrobeniusType.comp P F hv3 he, hlm, hDd,
    ⟨Under.homMk (show V.right ⟶ (⟨(Dd, Dd, E₃)⟩ : TriFr P F) from
      ⟨(a ≫ δ, c ≫ δ, e), IsFrobeniusType.comp P F ha hδ,
        IsFrobeniusType.comp P F hc hδ, he, hdaδ, ?_⟩)
      (WideSubcategory.hom_ext _ (Prod.ext rfl (Prod.ext hacδ.symm rfl)))⟩⟩
  rw [← hdaδ, hed]

/-! ## ★25. 対角の添字での `𝒪^▷` の判定

★★第 1・第 2 脚が**同じ射** `l` なら、`repBase = Base l ≫ Base a ≫ (Base l)⁻¹` が
共役なので、★`Base a = 𝟙` と `rootBase = 𝟙` は同値になる。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★**対角の添字での `𝒪^▷` の判定**。 -/
theorem oTri_mk_diag {X : PfRootObj P F} {D : C} (l : rtObj P F X.obj X.root ⟶ D)
    (hl : IsFrobeniusType P l) (a : D ⟶ D) :
    (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      ∈ OTri (pfRootPre P F) X ↔ (P.Base a = P.Base (𝟙 D) ∧ P.degFr a = 1) := by
  haveI hil : IsIso (P.Base l) := hl.2
  haveI hie : IsIso (P.Base (rtExt P F X.obj X.root)) := (rtExt_frobType P F X.obj X.root).2
  have hrep : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
      = P.Base l ≫ P.Base a :=
    repBase_spec (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hroot : rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root)
        ≫ pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a) :=
    rootBase_spec (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
  have hpf : pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
    pfBase_mk (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hidB : 𝟙 ((pfRootPre P F).toElem.obj X).base ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root) := Category.id_comp _
  have hidL : 𝟙 ((P.toElem.obj (rtObj P F X.obj X.root)).base) ≫ P.Base l
      = P.Base l := Category.id_comp _
  constructor
  · rintro ⟨hb, hd⟩
    refine ⟨?_, (degFr_mk_iff (X := X) (Y := X)
      (idxMk (P := P) (F := F) l l hl hl rfl) a 1).mp hd⟩
    have h0 : rootBase (show End X from
        HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        = (pfRootPre P F).Base (𝟙 X) := hb
    rw [(pfRootPre P F).Base_id] at h0
    rw [h0, hpf] at hroot
    -- hroot : 𝟙 ≫ Base (rtExt) = Base (rtExt) ≫ repBase
    have h1 : P.Base (rtExt P F X.obj X.root) ≫ 𝟙 _
        = P.Base (rtExt P F X.obj X.root)
          ≫ repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
      (Category.comp_id _).trans (hidB.symm.trans hroot)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      ((cancel_epi (P.Base (rtExt P F X.obj X.root))).mp h1).symm
    rw [h2] at hrep
    -- hrep : 𝟙 ≫ Base l = Base l ≫ Base a
    have h3 : P.Base l ≫ 𝟙 ((P.toElem.obj D).base) = P.Base l ≫ P.Base a :=
      (Category.comp_id _).trans (hidL.symm.trans hrep)
    rw [P.Base_id]
    exact ((cancel_epi (P.Base l)).mp h3).symm
  · rintro ⟨hb, hd⟩
    refine ⟨?_, (degFr_mk_iff (X := X) (Y := X)
      (idxMk (P := P) (F := F) l l hl hl rfl) a 1).mpr hd⟩
    rw [P.Base_id] at hb
    rw [hb] at hrep
    -- hrep : repBase ≫ Base l = Base l ≫ 𝟙
    have h1 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
        = 𝟙 _ ≫ P.Base l := hrep.trans ((Category.comp_id _).trans hidL.symm)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      (cancel_mono (P.Base l)).mp h1
    rw [hpf, h2] at hroot
    -- hroot : rootBase ≫ Base (rtExt) = Base (rtExt) ≫ 𝟙
    show rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = (pfRootPre P F).Base (𝟙 X)
    rw [(pfRootPre P F).Base_id]
    refine (cancel_mono (P.Base (rtExt P F X.obj X.root))).mp ?_
    exact hroot.trans ((Category.comp_id _).trans hidB.symm)

/-! ## ★26. (iii)(c) の順方向(根が等しい場合)

★★根が等しいと `compRoot α φ` と `compRoot φ β` が**同じ `rtRootIso`** を使う
(どちらも `r*r = r*r` の証明で、証明無関係により同一)。
★したがって両者は `compPf_mk` の 1 本で比較でき、
`𝒞` の `otriFwd` がそのまま渡る。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**(iii)(c) の順方向**(根が等しい場合)。 -/
theorem pfRoot_otriFwd_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ)
    (α : End (⟨A, r⟩ : PfRootObj P F)) (hα : α ∈ OTri (pfRootPre P F) ⟨A, r⟩) :
    ∃ β : End (⟨B, r⟩ : PfRootObj P F), β ∈ OTri (pfRootPre P F) ⟨B, r⟩ ∧
      (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
        = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
  obtain ⟨V, a₀, f₀, ha₀, hf₀⟩ := exists_rep3 (P := P) (F := F)
    ((rtRootIso P F A A (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv
      (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩))
    ((rtRootIso P F A B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).inv φ)
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag (F := F) hfi V
  obtain ⟨a₁, ha₁⟩ : ∃ t : Dd ⟶ Dd,
      t = idxTransport P F ((idx12 P F _ _ _).map u) a₀ := ⟨_, rfl⟩
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  have haW : (rtRootIso P F A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁ := by
    rw [ha₁]
    exact ha₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) a₀).symm
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀).symm
  have hαmk : (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) rfl) a₁ := by
    exact ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A A
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) haW)).trans
      (rtRootIso_hom_mk (F := F) A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ := by
    exact ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have ha₁O : P.Base a₁ = P.Base (𝟙 Dd) ∧ P.degFr a₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) a₁).mp ?_
    rw [← hαmk]; exact hα
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨b₁, ⟨hb₁O, hb₁e⟩, -⟩ :=
    F.otriFwd f₁ hf₁c hf₁s (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
  refine ⟨(rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    refine (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) _).mpr ⟨hb₁O.1, hb₁O.2⟩
  · show compRoot P F φ _ = compRoot P F (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) φ
    have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
      compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ _
    have e2 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
        (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
      compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁
    have e3 : HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃))
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
      congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e
    unfold compRoot
    rw [hfW, Iso.hom_inv_id_apply, haW, e1, e2, e3]

/-! ## ★27. 3 脚添字を「第 2・第 3 脚が一致し、第 1 脚の終域が isotropic」へ

★(iii)(c) の**逆方向**では `β` が第 2・第 3 脚にまたがる自己射なので、
今度はそちら側を対角にする。 -/

variable {P F} in
/-- ★★**3 脚添字を「第 2・第 3 脚が同じ射で、第 1 脚の終域が isotropic」へ**。 -/
theorem exists_idx3_diag23 (hfi : IsOfFrobeniusIsotropicType P) {A E : C}
    (V : IdxPf3 P F A E E) :
    ∃ (Dd E₃ : C) (l : A ⟶ Dd) (hl : IsFrobeniusType P l) (m : E ⟶ E₃)
      (hm : IsFrobeniusType P m) (hd : P.degFr l = P.degFr m),
      IsIsotropic P Dd ∧ Nonempty (V ⟶ idxMk3 (F := F) l m m hl hm hm hd rfl) := by
  obtain ⟨hv1, hv2, hv3, h12, h23⟩ := V.hom.property
  obtain ⟨X₃, s, t, hs, ht, hsd, htd, hst⟩ :=
    frob_common_upper P F V.hom.hom.2.1 hv2 V.hom.hom.2.2 hv3
  obtain ⟨A₁, a, ha, had⟩ := F.frobDegSurj V.right.obj.1 (P.degFr s)
  obtain ⟨Dd, δ, hδ, hDd⟩ := hfi A₁
  obtain ⟨E₃, ε, hε, hεd⟩ := F.frobDegSurj X₃ (P.degFr δ)
  have hstε : V.hom.hom.2.1 ≫ (s ≫ ε) = V.hom.hom.2.2 ≫ (t ≫ ε) := by
    rw [← Category.assoc, ← Category.assoc, hst]
    rfl
  have hsdt : P.degFr s = P.degFr t := by
    rw [hsd, htd, h23]
    rfl
  have hdeg1 : P.degFr (a ≫ δ) = P.degFr (s ≫ ε) := by
    rw [P.degFr_comp a δ, P.degFr_comp s ε, had, hεd]
  have hdeg2 : P.degFr (s ≫ ε) = P.degFr (t ≫ ε) := by
    rw [P.degFr_comp s ε, P.degFr_comp t ε, hsdt]
  have hlm : P.degFr (V.hom.hom.1 ≫ (a ≫ δ)) = P.degFr (V.hom.hom.2.1 ≫ (s ≫ ε)) := by
    rw [P.degFr_comp V.hom.hom.1 (a ≫ δ), P.degFr_comp V.hom.hom.2.1 (s ≫ ε),
      hdeg1, h12]
  refine ⟨Dd, E₃, V.hom.hom.1 ≫ (a ≫ δ),
    IsFrobeniusType.comp P F hv1 (IsFrobeniusType.comp P F ha hδ),
    V.hom.hom.2.1 ≫ (s ≫ ε), IsFrobeniusType.comp P F hv2
      (IsFrobeniusType.comp P F hs hε), hlm, hDd,
    ⟨Under.homMk (show V.right ⟶ (⟨(Dd, E₃, E₃)⟩ : TriFr P F) from
      ⟨(a ≫ δ, s ≫ ε, t ≫ ε), IsFrobeniusType.comp P F ha hδ,
        IsFrobeniusType.comp P F hs hε, IsFrobeniusType.comp P F ht hε,
        hdeg1, hdeg2⟩)
      (WideSubcategory.hom_ext _ (Prod.ext rfl (Prod.ext rfl hstε.symm)))⟩⟩

/-! ## ★28. (iii)(c) の逆方向(根が等しい場合) -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**(iii)(c) の逆方向**(根が等しい場合)。 -/
theorem pfRoot_otriBwd_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ)
    (β : End (⟨B, r⟩ : PfRootObj P F)) (hβ : β ∈ OTri (pfRootPre P F) ⟨B, r⟩) :
    ∃ α : End (⟨A, r⟩ : PfRootObj P F), α ∈ OTri (pfRootPre P F) ⟨A, r⟩ ∧
      (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
        = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
  obtain ⟨V, f₀, b₀, hf₀, hb₀⟩ := exists_rep3 (P := P) (F := F)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).inv (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag23 (F := F) hfi V
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx12 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  obtain ⟨b₁, hb₁⟩ : ∃ t : E₃ ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) b₀ := ⟨_, rfl⟩
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) f₀).symm
  have hbW : (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) b₁ := by
    rw [hb₁]
    exact hb₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) b₀).symm
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have hβmk : (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) rfl) b₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F B B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hbW)).trans
      (rtRootIso_hom_mk (F := F) B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) m m hm hm rfl) b₁)
  have hb₁O : P.Base b₁ = P.Base (𝟙 E₃) ∧ P.degFr b₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) b₁).mp ?_
    rw [← hβmk]; exact hβ
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨a₁, ⟨ha₁O, ha₁e⟩, -⟩ :=
    F.otriBwd f₁ hf₁c hf₁s (b₁ : End E₃) ⟨hb₁O.1, hb₁O.2⟩
  refine ⟨(rtRootIso P F A A (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) ((a₁ : End Dd) : Dd ⟶ Dd)), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) _).mpr ⟨ha₁O.1, ha₁O.2⟩
  · have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) b₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (f₁ ≫ b₁) :=
      compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ b₁
    have e2 : compPf P F
        (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) ((a₁ : End Dd) : Dd ⟶ Dd))
        (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁) :=
      compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) _ f₁
    have e3 : HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (f₁ ≫ b₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁) :=
      congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) ha₁e
    show compRoot P F φ _ = compRoot P F _ φ
    unfold compRoot
    rw [hfW, hbW, Iso.hom_inv_id_apply, e1, e2, e3]

/-! ## ★29. 同じ添字の 2 射の底が一致すれば、代表元の底も一致する -/

variable {P F} in
/-- ★**同じ添字なら `rootBase` の一致から `Base` の一致が出る**。 -/
theorem base_eq_of_rootBase_eq {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (f f' : Z.right.obj.1 ⟶ Z.right.obj.2)
    (h : rootBase (show HomRoot P F X Y from HomPf.mk Z f)
      = rootBase (show HomRoot P F X Y from HomPf.mk Z f')) :
    P.Base f = P.Base f' := by
  haveI h1 : IsIso (P.Base Z.hom.hom.1) := Z.hom.property.1.2
  haveI he : IsIso (P.Base (rtExt P F X.obj Y.root)) :=
    (rtExt_frobType P F X.obj Y.root).2
  have s1 := rootBase_spec (show HomRoot P F X Y from HomPf.mk Z f)
  have s2 := rootBase_spec (show HomRoot P F X Y from HomPf.mk Z f')
  rw [pfBase_mk] at s1 s2
  have e1 : P.Base (rtExt P F X.obj Y.root) ≫ repBase Z f
      = P.Base (rtExt P F X.obj Y.root) ≫ repBase Z f' := by
    rw [← s1, ← s2, h]
  have e2 : repBase Z f = repBase Z f' :=
    (cancel_epi (P.Base (rtExt P F X.obj Y.root))).mp e1
  have r1 := repBase_spec Z f
  have r2 := repBase_spec Z f'
  rw [e2] at r1
  exact (cancel_epi (P.Base Z.hom.hom.1)).mp (r1.symm.trans r2)

/-! ## ★30. (iii)(c) の「底にしか依らない」条(根が等しい場合) -/

set_option maxHeartbeats 4000000 in
variable {P F} in
/-- ★★★**(iii)(c) の第 3 条**(根が等しい場合)。 -/
theorem pfRoot_otriBase_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ φ' : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ) (hφ's : IsPreStep (pfRootPre P F) φ')
    (hbase : (pfRootPre P F).Base φ = (pfRootPre P F).Base φ')
    (α : End (⟨A, r⟩ : PfRootObj P F)) (hα : α ∈ OTri (pfRootPre P F) ⟨A, r⟩)
    (β : End (⟨B, r⟩ : PfRootObj P F)) (_hβ : β ∈ OTri (pfRootPre P F) ⟨B, r⟩)
    (heq : (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ)) :
    (φ' ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ') := by
  obtain ⟨V, a₀, f₀, f₀', ha₀, hf₀, hf₀'⟩ := exists_rep3_pair (P := P) (F := F)
    ((rtRootIso P F A A (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv
      (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩))
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ')
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag (F := F) hfi V
  obtain ⟨a₁, ha₁⟩ : ∃ t : Dd ⟶ Dd,
      t = idxTransport P F ((idx12 P F _ _ _).map u) a₀ := ⟨_, rfl⟩
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  obtain ⟨f₁', hf₁'⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀' := ⟨_, rfl⟩
  have haW : (rtRootIso P F A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁ := by
    rw [ha₁]
    exact ha₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) a₀).symm
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀).symm
  have hfW' : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ'
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁' := by
    rw [hf₁']
    exact hf₀'.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀').symm
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have hφ'mk : φ' = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁' :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW')).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
  have hαmk : (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) rfl) a₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A A
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) haW)).trans
      (rtRootIso_hom_mk (F := F) A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
  rw [hφmk, hφ'mk] at hbase
  have hbf : P.Base f₁ = P.Base f₁' :=
    base_eq_of_rootBase_eq (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ f₁' hbase
  have ha₁O : P.Base a₁ = P.Base (𝟙 Dd) ∧ P.degFr a₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) a₁).mp ?_
    rw [← hαmk]; exact hα
  have hidx : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) :=
    idxMk (P := P) (F := F) l m hl hm hdlm
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁'s : IsPreStep P f₁' := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁').mp ?_
    rw [← hφ'mk]; exact hφ's
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  have hf₁'c : IsCoAngular P f₁' :=
    prop_1_4_i P f₁' (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨b₁, ⟨hb₁O, hb₁e⟩, -⟩ :=
    F.otriFwd f₁ hf₁c hf₁s (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
  have hb₁e' : (f₁' ≫ ((b₁ : End E₃) : E₃ ⟶ E₃))
      = ((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁' :=
    F.otriBase f₁ f₁' hf₁c hf₁s hf₁'c hf₁'s hbf (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
      (b₁ : End E₃) hb₁O hb₁e
  have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ _
  have e1' : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁' ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁' _
  have e2 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
      (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
    compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁
  have e2' : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
      (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁') :=
    compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁'
  have heq0 : (φ ≫ (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
    show compRoot P F φ _ = compRoot P F _ φ
    unfold compRoot
    rw [hfW, Iso.hom_inv_id_apply, haW, e1, e2]
    exact congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom)
      (congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e)
  have heq0' : (φ' ≫ (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ') := by
    show compRoot P F φ' _ = compRoot P F _ φ'
    unfold compRoot
    rw [hfW', Iso.hom_inv_id_apply, haW, e1', e2']
    exact congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom)
      (congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e')
  haveI : Epi φ := pfRoot_totEpi P F _ _ φ
  have hββ : (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    (cancel_epi φ).mp (heq.trans heq0.symm)
  rw [hββ]
  exact heq0'

/-! ## ★31. 根を揃えて (iii)(c) を一般の場合へ

★★`X = (A,n)`、`Y = (B,m)` を `X' = (A^m, n*m)`、`Y' = (B^n, n*m)` に取り替える。
★どちらも根が `n*m` なので `sameRoot` 版が当たり、結果は共役で戻す。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の順方向、`𝒞^pf` 版。 -/
theorem pfRoot_otriFwd (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (α : End X) (hα : α ∈ OTri (pfRootPre P F) X) :
    ∃! β : End Y, β ∈ OTri (pfRootPre P F) Y ∧
      ((φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hα' := oTri_conj (F := F) eX (inv eX) (IsIso.hom_inv_id eX) (IsIso.inv_hom_id eX)
    α hα
  obtain ⟨β', hβ'O, hβ'e⟩ := pfRoot_otriFwd_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY) hφ'
    (inv eX ≫ (α : X ⟶ X) ≫ eX) hα'
  have hβO := oTri_conj (F := F) (inv eY) eY (IsIso.inv_hom_id eY) (IsIso.hom_inv_id eY)
    β' hβ'O
  have h1 : inv eX ≫ (φ ≫ eY ≫ (β' : _ ⟶ _))
      = inv eX ≫ ((α : X ⟶ X) ≫ φ ≫ eY) := by
    have h := hβ'e
    simp only [Category.assoc] at h ⊢
    rw [h]
    simp only [Category.assoc, IsIso.hom_inv_id_assoc]
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  have h2 : φ ≫ eY ≫ (β' : _ ⟶ _) = (α : X ⟶ X) ≫ φ ≫ eY :=
    (cancel_epi (inv eX)).mp h1
  have hβe : (φ ≫ (eY ≫ (β' : _ ⟶ _) ≫ inv eY)) = ((α : X ⟶ X) ≫ φ) := by
    calc φ ≫ (eY ≫ (β' : _ ⟶ _) ≫ inv eY)
        = (φ ≫ eY ≫ (β' : _ ⟶ _)) ≫ inv eY := by simp only [Category.assoc]
      _ = ((α : X ⟶ X) ≫ φ ≫ eY) ≫ inv eY := by rw [h2]
      _ = (α : X ⟶ X) ≫ φ := by
          simp only [Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  refine ⟨eY ≫ (β' : _ ⟶ _) ≫ inv eY, ⟨hβO, hβe⟩, ?_⟩
  rintro y ⟨-, hy⟩
  haveI : Epi φ := pfRoot_totEpi P F _ _ φ
  exact (cancel_epi φ).mp (hy.trans hβe.symm)

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の逆方向、`𝒞^pf` 版。 -/
theorem pfRoot_otriBwd (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (β : End Y) (hβ : β ∈ OTri (pfRootPre P F) Y) :
    ∃! α : End X, α ∈ OTri (pfRootPre P F) X ∧
      ((φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hβ' := oTri_conj (F := F) eY (inv eY) (IsIso.hom_inv_id eY) (IsIso.inv_hom_id eY)
    β hβ
  obtain ⟨α', hα'O, hα'e⟩ := pfRoot_otriBwd_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY) hφ'
    (inv eY ≫ (β : Y ⟶ Y) ≫ eY) hβ'
  have hαO := oTri_conj (F := F) (inv eX) eX (IsIso.inv_hom_id eX) (IsIso.hom_inv_id eX)
    α' hα'O
  have h1 : inv eX ≫ (φ ≫ (β : Y ⟶ Y) ≫ eY)
      = inv eX ≫ ((eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ ≫ eY) := by
    have h := hα'e
    simp only [Category.assoc, IsIso.hom_inv_id_assoc] at h
    simp only [Category.assoc, IsIso.inv_hom_id_assoc]
    exact h
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  have h2 : φ ≫ (β : Y ⟶ Y) ≫ eY = (eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ ≫ eY :=
    (cancel_epi (inv eX)).mp h1
  haveI : Mono eY := IsIso.mono_of_iso _
  have hαe : (φ ≫ (β : Y ⟶ Y)) = ((eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ) := by
    refine (cancel_mono eY).mp ?_
    simp only [Category.assoc]
    simp only [Category.assoc] at h2
    exact h2
  refine ⟨eX ≫ (α' : _ ⟶ _) ≫ inv eX, ⟨hαO, hαe⟩, ?_⟩
  rintro y ⟨-, hy⟩
  haveI : Mono φ := pfRoot_preStepMono (F := F) φ hφs
  exact (cancel_mono φ).mp (hy.symm.trans hαe)

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の第 3 条、`𝒞^pf` 版。 -/
theorem pfRoot_otriBase (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ φ' : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (_hφ'c : IsCoAngular (pfRootPre P F) φ')
    (hφ's : IsPreStep (pfRootPre P F) φ')
    (hbase : (pfRootPre P F).Base φ = (pfRootPre P F).Base φ')
    (α : End X) (hα : α ∈ OTri (pfRootPre P F) X)
    (β : End Y) (hβ : β ∈ OTri (pfRootPre P F) Y)
    (heq : (φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) :
    (φ' ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ') := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hψ : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hψ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ' ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφ's (isPreStep_of_isIso (pfRootPre P F) eY))
  have hbase' : (pfRootPre P F).Base (inv eX ≫ φ ≫ eY)
      = (pfRootPre P F).Base (inv eX ≫ φ' ≫ eY) := by
    rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp,
      (pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp, hbase]
  have hα' := oTri_conj (F := F) eX (inv eX) (IsIso.hom_inv_id eX) (IsIso.inv_hom_id eX)
    α hα
  have hβ' := oTri_conj (F := F) eY (inv eY) (IsIso.hom_inv_id eY) (IsIso.inv_hom_id eY)
    β hβ
  have heq' : ((inv eX ≫ φ ≫ eY) ≫ (inv eY ≫ (β : Y ⟶ Y) ≫ eY))
      = ((inv eX ≫ (α : X ⟶ X) ≫ eX) ≫ (inv eX ≫ φ ≫ eY)) := by
    simp only [Category.assoc, IsIso.hom_inv_id_assoc, IsIso.inv_hom_id_assoc]
    have h3 : (φ ≫ (β : Y ⟶ Y)) ≫ eY = ((α : X ⟶ X) ≫ φ) ≫ eY :=
      congrArg (fun t => t ≫ eY) heq
    simp only [Category.assoc] at h3
    exact congrArg (fun t => inv eX ≫ t) h3
  have hres := pfRoot_otriBase_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY)
    (inv eX ≫ φ' ≫ eY) hψ hψ' hbase' (inv eX ≫ (α : X ⟶ X) ≫ eX) hα'
    (inv eY ≫ (β : Y ⟶ Y) ≫ eY) hβ' heq'
  simp only [Category.assoc, IsIso.hom_inv_id_assoc, IsIso.inv_hom_id_assoc] at hres
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  haveI : Mono eY := IsIso.mono_of_iso _
  refine (cancel_mono eY).mp ?_
  refine (cancel_epi (inv eX)).mp ?_
  simp only [Category.assoc]
  exact hres

/-! ## ★33. 単系の補題 —— divisorial なら `n • a = n • b → a = b`

★★`Φ(A)` は divisorial(integral ＋ saturated ＋ sharp)。
`n • (a − b) = 0` は像に入るので **saturated** で `a − b` 自身が像に入り、
`b ≼ a` が出る。対称に `a ≼ b`。**integral の簡約**と **sharp** で `a = b`。

★これが `MetricallyEquivalent` を `Pf` から代表元へ降ろす鍵である。 -/

theorem eq_of_nsmul_eq_of_divisorial {M : Type w} [AddCommMonoid M]
    (hdiv : IsDivisorial M) {a b : M} {n : ℕ} (hn : 0 < n) (hab : n • a = n • b) : a = b := by
  obtain ⟨⟨hint, hsat, -⟩, hsharp⟩ := hdiv
  haveI : IsCancelAdd M := isCancelAdd_of_isIntegralMonoid (M := M) hint
  have key : ∀ x y : M, n • x = n • y → ∃ c : M, y + c = x := by
    intro x y hxy
    have h1 : n • (toGp M x - toGp M y) = 0 := by
      rw [smul_sub, ← toGp_nsmul, ← toGp_nsmul, hxy, sub_self]
    have h2 : n • (toGp M x - toGp M y) ∈ Set.range (toGp M) := by
      rw [h1]; exact ⟨0, toGp_zero M⟩
    obtain ⟨c, hc⟩ := hsat _ n hn h2
    refine ⟨c, hint ?_⟩
    rw [toGp_add, hc]
    abel
  obtain ⟨c, hc⟩ := key a b hab
  obtain ⟨d, hd⟩ := key b a hab.symm
  have h3 : b + (c + d) = b + 0 := by
    rw [← add_assoc, hc, hd, add_zero]
  have h4 : c + d = 0 := add_left_cancel h3
  have h5 : c = 0 := hsharp c ⟨⟨c, d, h4, by rw [add_comm]; exact h4⟩, rfl⟩
  rw [← hc, h5, add_zero]

/-- ★**`Pf` の同じ分母どうしは分子で決まる**(divisorial なら)。 -/
theorem Pf.mk_inj_of_divisorial {M : Type w} [AddCommMonoid M] (hdiv : IsDivisorial M)
    {a b : M} {n : ℕ+} (h : Pf.mk a n = Pf.mk b n) : a = b := by
  obtain ⟨k, hk⟩ := Quotient.exact h
  exact eq_of_nsmul_eq_of_divisorial hdiv (Nat.mul_pos k.pos n.pos) hk

variable {P F} in
/-- ★★**同じ添字なら `rootDiv` の一致から `Div` の一致が出る**。 -/
theorem div_eq_of_rootDiv_eq {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (f f' : Z.right.obj.1 ⟶ Z.right.obj.2)
    (h : rootDiv (show HomRoot P F X Y from HomPf.mk Z f)
      = rootDiv (show HomRoot P F X Y from HomPf.mk Z f')) :
    P.Div f = P.Div f' := by
  replace h : Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv (HomPf.mk Z f)))
    = Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv (HomPf.mk Z f'))) := h
  rw [pfDiv_mk, pfDiv_mk] at h
  replace h : Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root)))
        (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div f)) (repRoot Z)))
    = Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root)))
        (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div f')) (repRoot Z))) := h
  rw [Pf.map_mk, Pf.map_mk, Pf.divBy_mk, Pf.divBy_mk] at h
  have h1 := Pf.mk_inj_of_divisorial (P.divisorial _) h
  exact Φ.map_injective (P.Base Z.hom.hom.1)
    (Φ.map_injective (P.Base (rtExt P F X.obj Y.root)) h1)

/-! ## ★34. (vi) 単元を除く忠実性(根が等しい場合) -/

variable {P F} in
/-- ★平行 2 射を共通の添字へ。 -/
theorem exists_rep_par {A B : C} (f g : HomPf P F A B) :
    ∃ (Z : IdxPf P F A B) (φ ψ : Z.right.obj.1 ⟶ Z.right.obj.2),
      f = HomPf.mk Z φ ∧ g = HomPf.mk Z ψ := by
  obtain ⟨Zf, φ₀, hf⟩ := HomPf.exists_rep (P := P) (F := F) f
  obtain ⟨Zg, ψ₀, hg⟩ := HomPf.exists_rep (P := P) (F := F) g
  exact ⟨IsFiltered.max Zf Zg,
    idxTransport P F (IsFiltered.leftToMax Zf Zg) φ₀,
    idxTransport P F (IsFiltered.rightToMax Zf Zg) ψ₀,
    (by rw [HomPf.mk_map]; exact hf.symm),
    (by rw [HomPf.mk_map]; exact hg.symm)⟩

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★**(vi)**(根が等しい場合)。 -/
theorem pfRoot_faithfulUpToUnits_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hbe : BaseEquivalent (pfRootPre P F) φ ψ)
    (hme : MetricallyEquivalent (pfRootPre P F) φ ψ)
    (hφs : IsPreStep (pfRootPre P F) φ) (hψs : IsPreStep (pfRootPre P F) ψ) :
    ∃ α : End (⟨B, r⟩ : PfRootObj P F), α ∈ OTimes (pfRootPre P F) ⟨B, r⟩ ∧
      φ = (ψ ≫ ((α : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))) := by
  obtain ⟨Z₀, φ₀, ψ₀, hφ0, hψ0⟩ := exists_rep_par (P := P) (F := F)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv ψ)
  obtain ⟨W, u, hW⟩ := exists_idx_isotropic (F := F) hfi Z₀
  obtain ⟨φ₁, hφ₁⟩ : ∃ t : W.right.obj.1 ⟶ W.right.obj.2,
      t = idxTransport P F u φ₀ := ⟨_, rfl⟩
  obtain ⟨ψ₁, hψ₁⟩ : ∃ t : W.right.obj.1 ⟶ W.right.obj.2,
      t = idxTransport P F u ψ₀ := ⟨_, rfl⟩
  obtain ⟨hw1, hw2, hwd⟩ := W.hom.property
  have hφW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ = HomPf.mk W φ₁ := by
    rw [hφ₁, HomPf.mk_map]; exact hφ0
  have hψW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv ψ = HomPf.mk W ψ₁ := by
    rw [hψ₁, HomPf.mk_map]; exact hψ0
  have hφmk : φ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
      (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hφW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) W φ₁)
  have hψmk : ψ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
      (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj W) ψ₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hψW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) W ψ₁)
  rw [hφmk, hψmk] at hbe hme
  have hbe' : BaseEquivalent P φ₁ ψ₁ :=
    base_eq_of_rootBase_eq (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ ψ₁ hbe
  have hme' : MetricallyEquivalent P φ₁ ψ₁ :=
    div_eq_of_rootDiv_eq (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ ψ₁ hme
  have hφ₁s : IsPreStep P φ₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁).mp ?_
    rw [← hφmk]; exact hφs
  have hψ₁s : IsPreStep P ψ₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) ψ₁).mp ?_
    rw [← hψmk]; exact hψs
  have hφ₁c : IsCoAngular P φ₁ :=
    prop_1_4_i P φ₁ (fun _ g => F.isotropicClosed g hW)
  have hψ₁c : IsCoAngular P ψ₁ :=
    prop_1_4_i P ψ₁ (fun _ g => F.isotropicClosed g hW)
  obtain ⟨α₀, hα₀U, hα₀e⟩ :=
    F.faithfulUpToUnits φ₁ ψ₁ hbe' hme' hφ₁c hφ₁s hψ₁c hψ₁s
  refine ⟨(rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) W.hom.hom.2 W.hom.hom.2 hw2 hw2 rfl)
      ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2)), ?_, ?_⟩
  · refine ⟨?_, ?_⟩
    · rw [rtRootIso_hom_mk]
      exact (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
        (rtLift P F B (show r * r = r * r from rfl) ≫ W.hom.hom.2)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hw2) _).mpr
        ⟨hα₀U.1.1, hα₀U.1.2⟩
    · refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
      rw [rtRootIso_hom_mk]
      exact pfRoot_isIso_mk _ _ ((CategoryTheory.isUnit_iff_isIso _).mp hα₀U.2)
  · have e1 : compPf P F (HomPf.mk W ψ₁)
        (HomPf.mk (idxMk (P := P) (F := F) W.hom.hom.2 W.hom.hom.2 hw2 hw2 rfl)
          ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2))
        = HomPf.mk W (ψ₁ ≫ ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2)) :=
      compPf_mk (idxMk3 (F := F) W.hom.hom.1 W.hom.hom.2 W.hom.hom.2 hw1 hw2 hw2 hwd rfl)
        ψ₁ _
    show φ = compRoot P F ψ _
    unfold compRoot
    rw [hψW, Iso.hom_inv_id_apply]
    refine Eq.trans ?_ (congrArg (ConcreteCategory.hom (rtRootIso P F A B
      (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) e1).symm
    refine Eq.trans ?_ (congrArg (fun t => ConcreteCategory.hom (rtRootIso P F A B
      (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom
      (HomPf.mk W t)) hα₀e)
    rw [← hφW, Iso.inv_hom_id_apply]

/-! ## ★35. (vi) を一般の根へ -/

/-- ★同型で挟んでも `Div` の一致は保たれる。 -/
theorem metricallyEquivalent_conj {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3) {X X' Y Y' : C3}
    (a : X' ⟶ X) (b : Y ⟶ Y') (ha : IsIso a) (hb : IsIso b) {f g : X ⟶ Y}
    (h : MetricallyEquivalent Q f g) :
    MetricallyEquivalent Q (a ≫ f ≫ b) (a ≫ g ≫ b) := by
  haveI := ha
  haveI := hb
  show Q.Div (a ≫ f ≫ b) = Q.Div (a ≫ g ≫ b)
  rw [Q.Div_comp, Q.Div_comp, Q.Div_comp, Q.Div_comp,
    show Q.Div a = 0 from isIsometric_of_isIso Q a,
    show Q.Div b = 0 from isIsometric_of_isIso Q b,
    show Q.Div f = Q.Div g from h]
  simp

/-- ★同型で挟んでも `Base` の一致は保たれる。 -/
theorem baseEquivalent_conj {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3) {X X' Y Y' : C3}
    (a : X' ⟶ X) (b : Y ⟶ Y') {f g : X ⟶ Y} (h : BaseEquivalent Q f g) :
    BaseEquivalent Q (a ≫ f ≫ b) (a ≫ g ≫ b) := by
  show Q.Base (a ≫ f ≫ b) = Q.Base (a ≫ g ≫ b)
  rw [Q.Base_comp, Q.Base_comp, Q.Base_comp, Q.Base_comp,
    show Q.Base f = Q.Base g from h]

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (vi)** —— `𝒞^pf` 版。 -/
theorem pfRoot_faithfulUpToUnits (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (φ ψ : X ⟶ Y)
    (hbe : BaseEquivalent (pfRootPre P F) φ ψ)
    (hme : MetricallyEquivalent (pfRootPre P F) φ ψ)
    (_hφc : IsCoAngular (pfRootPre P F) φ) (hφs : IsPreStep (pfRootPre P F) φ)
    (_hψc : IsCoAngular (pfRootPre P F) ψ) (hψs : IsPreStep (pfRootPre P F) ψ) :
    ∃ α : End Y, α ∈ OTimes (pfRootPre P F) Y ∧ φ = (ψ ≫ ((α : Y ⟶ Y))) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hψ' : IsPreStep (pfRootPre P F) (inv eX ≫ ψ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hψs (isPreStep_of_isIso (pfRootPre P F) eY))
  obtain ⟨α', hα'U, hα'e⟩ := pfRoot_faithfulUpToUnits_sameRoot (F := F) hfi
    (inv eX ≫ φ ≫ eY) (inv eX ≫ ψ ≫ eY)
    (baseEquivalent_conj (pfRootPre P F) (inv eX) eY hbe)
    (metricallyEquivalent_conj (pfRootPre P F) (inv eX) eY inferInstance inferInstance hme)
    hφ' hψ'
  refine ⟨eY ≫ ((α' : _ ⟶ _) ≫ inv eY), ⟨?_, ?_⟩, ?_⟩
  · exact (oTri_conj (F := F) (inv eY) eY (IsIso.inv_hom_id eY) (IsIso.hom_inv_id eY)
      α' hα'U.1)
  · refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
    haveI : IsIso ((α' : _ ⟶ _)) := (CategoryTheory.isUnit_iff_isIso _).mp hα'U.2
    show IsIso (eY ≫ ((α' : _ ⟶ _) ≫ inv eY))
    infer_instance
  · have h := hα'e
    simp only [Category.assoc] at h ⊢
    have h2 : eX ≫ inv eX ≫ φ ≫ eY = eX ≫ inv eX ≫ ψ ≫ eY ≫ (α' : _ ⟶ _) := by
      rw [h]
    simp only [IsIso.hom_inv_id_assoc] at h2
    have h3 : (φ ≫ eY) ≫ inv eY = (ψ ≫ eY ≫ (α' : _ ⟶ _)) ≫ inv eY :=
      congrArg (fun t => t ≫ inv eY) h2
    simpa only [Category.assoc, IsIso.hom_inv_id, Category.comp_id] using h3

/-! ## ★36. ★★★★`k` 乗根の高さへの持ち上げ `Λ_k`

★★**動機**(2026-08-17): 残る 4 条(`pullBackLB` / `arbFactor` / `arbFactorUniq` /
`plBkEquiv`)はどれも「`𝒞` の射を `𝒞^pf` へ**根つきで**送る」道具を要る。
★`toRootHom` は根 `1` にしか送れないが、必要なのは**任意の根 `k`** である。

★★**構成**: `rtExt A k : A ⟶ A^{(k)}` と `rtExt B k : B ⟶ B^{(k)}` は
同じ次数 `k` の Frobenius 型射なので、`Proposition 1.10, (i)` の遷移
(`frobTransport`)が `φ : A ⟶ B` を `A^{(k)} ⟶ B^{(k)}` へ一意に持ち上げる。 -/

variable {P F} in
/-- ★★`φ : A ⟶ B` を **`k` 乗根の高さへ持ち上げた射** `A^{(k)} ⟶ B^{(k)}`。 -/
noncomputable def rtMap (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    rtObj P F A k ⟶ rtObj P F B k :=
  frobTransport (F := F) (rtExt P F A k) (rtExt_frobType P F A k)
    (rtExt P F B k) (rtExt_frobType P F B k)
    (by rw [rtExt_degFr, rtExt_degFr]) φ

variable {P F} in
/-- ★持ち上げの特徴づけ(四角形)。 -/
theorem rtMap_spec (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    φ ≫ rtExt P F B k = rtExt P F A k ≫ rtMap (F := F) k φ :=
  frobTransport_spec _ _ _ _ _ φ

variable {P F} in
/-- ★持ち上げの一意性。 -/
theorem rtMap_eq (k : ℕ+) {A B : C} (φ : A ⟶ B) (ψ : rtObj P F A k ⟶ rtObj P F B k)
    (h : φ ≫ rtExt P F B k = rtExt P F A k ≫ ψ) : rtMap (F := F) k φ = ψ :=
  frobTransport_eq _ _ _ _ _ φ ψ h

variable {P F} in
/-- ★★**次数は変わらない**。 -/
theorem rtMap_degFr (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.degFr (rtMap (F := F) k φ) = P.degFr φ := by
  have h := congrArg P.degFr (rtMap_spec (F := F) k φ)
  rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr] at h
  have h2 : P.degFr φ * k = P.degFr (rtMap (F := F) k φ) * k := (mul_comm _ _).trans h
  exact (mul_right_cancel h2).symm

variable {P F} in
/-- ★★**零因子は `k` 倍で写る** —— `(Base ext)^*(Div (Λφ)) = k • Div φ`。 -/
theorem rtMap_Div (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    Φ.map (P.Base (rtExt P F A k)) (P.Div (rtMap (F := F) k φ))
      = ((k : ℕ+) : ℕ) • P.Div φ := by
  have h := congrArg P.Div (rtMap_spec (F := F) k φ)
  rw [P.Div_comp, P.Div_comp, show P.Div (rtExt P F B k) = 0 from (rtExt_frobType P F B k).1.2,
    show P.Div (rtExt P F A k) = 0 from (rtExt_frobType P F A k).1.2,
    rtExt_degFr, map_zero, zero_add, smul_zero, add_zero] at h
  exact h.symm

variable {P F} in
/-- ★★**等長性は行き来する**。 -/
theorem rtMap_isIsometric_iff (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    IsIsometric P (rtMap (F := F) k φ) ↔ IsIsometric P φ := by
  constructor
  · intro h
    have h1 : ((k : ℕ+) : ℕ) • P.Div φ = 0 := by
      rw [← rtMap_Div (F := F) k φ, show P.Div (rtMap (F := F) k φ) = 0 from h, map_zero]
    exact nsmul_eq_zero_of_isSharp (P.divisorial _).2 h1
  · intro h
    refine Φ.map_injective (P.Base (rtExt P F A k)) ?_
    rw [rtMap_Div, show P.Div φ = 0 from h, smul_zero, map_zero]

variable {P F} in
/-- ★★**底の同型性も行き来する**。 -/
theorem rtMap_isBaseIso_iff (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    IsBaseIsomorphism P (rtMap (F := F) k φ) ↔ IsBaseIsomorphism P φ := by
  haveI hA : IsIso (P.Base (rtExt P F A k)) := (rtExt_frobType P F A k).2
  haveI hB : IsIso (P.Base (rtExt P F B k)) := (rtExt_frobType P F B k).2
  refine (isIso_iff_of_sq (P.Base φ) (P.Base (rtExt P F B k))
    (P.Base (rtExt P F A k)) (P.Base (rtMap (F := F) k φ)) ?_).symm
  rw [← P.Base_comp, ← P.Base_comp, rtMap_spec]

variable {P F} in
/-- ★恒等射の持ち上げは恒等射。 -/
@[simp] theorem rtMap_id (k : ℕ+) (A : C) :
    rtMap (F := F) k (𝟙 A) = 𝟙 (rtObj P F A k) :=
  rtMap_eq k _ _ (by rw [Category.id_comp, Category.comp_id])

variable {P F} in
/-- ★★**持ち上げは合成を保つ**。 -/
theorem rtMap_comp (k : ℕ+) {A B E : C} (φ : A ⟶ B) (ψ : B ⟶ E) :
    rtMap (F := F) k (φ ≫ ψ) = rtMap (F := F) k φ ≫ rtMap (F := F) k ψ := by
  refine rtMap_eq k _ _ ?_
  rw [Category.assoc, rtMap_spec, ← Category.assoc, rtMap_spec, Category.assoc]

variable {P F} in
/-- ★同型の持ち上げは同型。 -/
theorem rtMap_isIso (k : ℕ+) {A B : C} (u : A ⟶ B) [IsIso u] :
    IsIso (rtMap (F := F) k u) :=
  ⟨rtMap (F := F) k (inv u),
    by rw [← rtMap_comp, IsIso.hom_inv_id, rtMap_id],
    by rw [← rtMap_comp, IsIso.inv_hom_id, rtMap_id]⟩

/-! ### ★`Λ_k` —— `𝒞` の射から `𝒞^pf` の根 `k` の射へ -/

variable {P F} in
/-- ★★★**`φ : A ⟶ B` の根 `k` での像** `(A,k) ⟶ (B,k)`。 -/
noncomputable def lamHom (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (⟨A, k⟩ : PfRootObj P F) ⟶ (⟨B, k⟩ : PfRootObj P F) :=
  HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)

variable {P F} in
@[simp] theorem lamHom_degFr (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).degFr (lamHom (F := F) k φ) = P.degFr φ := by
  show rootDeg (show HomRoot P F ⟨A, k⟩ ⟨B, k⟩ from
    HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)) = _
  rw [rootDeg_mk]
  exact rtMap_degFr (F := F) k φ

variable {P F} in
theorem lamHom_isLinear (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsLinear P φ) :
    IsLinear (pfRootPre P F) (lamHom (F := F) k φ) := by
  show (pfRootPre P F).degFr (lamHom (F := F) k φ) = 1
  rw [lamHom_degFr]; exact h

variable {P F} in
theorem lamHom_isIsometric (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsIsometric P φ) :
    IsIsometric (pfRootPre P F) (lamHom (F := F) k φ) :=
  (isIsometric_mk_iff (X := (⟨A, k⟩ : PfRootObj P F)) (Y := (⟨B, k⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)).mpr
    ((rtMap_isIsometric_iff (F := F) k φ).mpr h)

variable {P F} in
theorem lamHom_isBaseIso (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism (pfRootPre P F) (lamHom (F := F) k φ) :=
  (isBaseIsomorphism_mk_iff (X := (⟨A, k⟩ : PfRootObj P F)) (Y := (⟨B, k⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)).mpr
    ((rtMap_isBaseIso_iff (F := F) k φ).mpr h)

variable {P F} in
theorem lamHom_isPreStep (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsPreStep P φ) :
    IsPreStep (pfRootPre P F) (lamHom (F := F) k φ) :=
  ⟨lamHom_isLinear k φ h.1, lamHom_isBaseIso k φ h.2⟩

variable {P F} in
/-- ★★同型の像は同型。 -/
theorem lamHom_isIso (k : ℕ+) {A B : C} (u : A ⟶ B) (hu : IsIso u) :
    IsIso (lamHom (F := F) k u) := by
  haveI := hu
  exact pfRoot_isIso_mk _ _ (rtMap_isIso (F := F) k u)

variable {P F} in
/-- ★★**Frobenius 型の像は Frobenius 型**(`𝒞^pf` が isotropic 型のとき)。 -/
theorem lamHom_isFrobeniusType (hfi : IsOfFrobeniusIsotropicType P) (k : ℕ+)
    {A B : C} (φ : A ⟶ B) (h : IsFrobeniusType P φ) :
    IsFrobeniusType (pfRootPre P F) (lamHom (F := F) k φ) :=
  ⟨⟨pfRoot_isCoAngular hfi _, lamHom_isIsometric k φ h.1.2⟩, lamHom_isBaseIso k φ h.2⟩

/-! ## ★37. 根を上げる同型の一般形と、根が等しいときの合成則 -/

variable {P F} in
/-- ★★**任意の Frobenius 型射で根を上げられる** —— `pfRoot_exists_iso_root` の一般形。

★`z : A ⟶ U` が次数 `d` の Frobenius 型射なら `(A, k) ≅ (U, k*d)`。
★`pfRoot_exists_iso_root` は `z = rtExt A d` の場合であった。 -/
theorem pfRoot_iso_of_frobType {A U : C} (z : A ⟶ U) (hz : IsFrobeniusType P z)
    {d k t : ℕ+} (hzd : P.degFr z = d) (ht : t = k * d) :
    ∃ e : (⟨A, k⟩ : PfRootObj P F) ⟶ ⟨U, t⟩, IsIso e := by
  obtain ⟨e₀, he₀⟩ := pfRoot_exists_iso_root (F := F) A k d t ht
  obtain ⟨θ, hθ, -⟩ := F.frobDegUniq A U (rtObj P F A d) z (rtExt P F A d) hz
    (rtExt_frobType P F A d) (by rw [hzd, rtExt_degFr])
  haveI := hθ
  haveI := he₀
  haveI : IsIso (lamHom (F := F) t (@inv _ _ _ _ θ hθ)) :=
    lamHom_isIso t _ (@IsIso.inv_isIso _ _ _ _ θ hθ)
  exact ⟨e₀ ≫ lamHom (F := F) t (@inv _ _ _ _ θ hθ), inferInstance⟩

/-- ★**根 `r` の対象どうしの `Hom` を `r·r` 乗の高さで書くときの添字の押し出し**。 -/
noncomputable def pushRt (A B : C) (r : ℕ+) :
    IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) ⥤
      IdxPf P F (rtObj P F A r) (rtObj P F B r) :=
  pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
    (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
    (by rw [rtLift_degFr, rtLift_degFr])

/-- ★**根が等しいときの、`Hom` の同一視** ——
`Hom^pf(A^{(r·r)}, B^{(r·r)}) ≅ ((A,r) ⟶ (B,r))`。 -/
noncomputable def rtJ (A B : C) (r : ℕ+) :
    HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r))
      ≅ HomPf P F (rtObj P F A r) (rtObj P F B r) :=
  rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)

variable {P F} in
/-- ★★**根が等しいときの合成則** —— `compRoot` は `rtJ` で共役した `compPf` である。

★★3 つの `rtRootIso` の添字はすべて `(dA, dB, e, tA, tB) = (r, r, r, r·r, r·r)` で、
命題の証明部分は証明無関係で潰れる。★これが (iii)(c) や (vi) で効いた仕掛けと同じもの。 -/
theorem compRoot_sameRoot {A B E : C} {r : ℕ+}
    (f : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) (g : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) :
    f ≫ g = (rtJ P F A E r).hom
      (compPf P F ((rtJ P F A B r).inv f) ((rtJ P F B E r).inv g)) := rfl

variable {P F} in
/-- ★★**`rtJ` は合成を移す**。 -/
theorem rtJ_comp {A B E : C} {r : ℕ+}
    (x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (y : HomPf P F (rtObj P F B (r * r)) (rtObj P F E (r * r))) :
    (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from (rtJ P F A B r).hom x)
        ≫ (show ((⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) from (rtJ P F B E r).hom y)
      = (rtJ P F A E r).hom (compPf P F x y) := by
  rw [compRoot_sameRoot, Iso.hom_inv_id_apply, Iso.hom_inv_id_apply]

variable {P F} in
/-- ★`rtJ` の代表元での計算則。 -/
theorem rtJ_hom_mk {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (rtJ P F A B r).hom (HomPf.mk V φ)
      = HomPf.mk ((pushRt P F A B r).obj V) φ :=
  rtRootIso_hom_mk (F := F) A B _ _ V φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 次数。 -/
theorem rtJ_degFr {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (pfRootPre P F).degFr
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      = P.degFr φ := by
  rw [rtJ_hom_mk]
  exact rootDeg_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
    ((pushRt P F A B r).obj V) φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 等長性。 -/
theorem rtJ_isIsometric_iff {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    IsIsometric (pfRootPre P F)
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      ↔ IsIsometric P φ := by
  rw [rtJ_hom_mk]
  exact isIsometric_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
    ((pushRt P F A B r).obj V) φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 底の同型性。 -/
theorem rtJ_isBaseIso_iff {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    IsBaseIsomorphism (pfRootPre P F)
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      ↔ IsBaseIsomorphism P φ := by
  rw [rtJ_hom_mk]
  exact isBaseIsomorphism_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
    (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushRt P F A B r).obj V) φ

/-! ## ★38. ★★★★`Hom^pf` の Frobenius 分解

★★**これが `pullBackLB` / `arbFactor` の共通の足場**である。
`𝒞` の `arbFactor` を代表元に当て、**Frobenius 型の部分を添字の脚へ吸収する**。
★吸収先は **`A₀` の次数 `n` の標準拡大 `A₀^{(n)}`** に取る —— そうすると
`𝒞^pf` 側で中間対象を `(A^{(n)}, r)` の形に書ける。 -/

variable {P F} in
/-- ★添字の始対象から `idxMk` への遷移射。 -/
noncomputable def idxOneHom {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    idxOne P F A B ⟶ idxMk (P := P) (F := F) a b ha hb hd :=
  Under.homMk (show (⟨(A, B)⟩ : BiFr P F) ⟶ (⟨(A', B')⟩ : BiFr P F) from
    ⟨(a, b), ha, hb, hd⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))

variable {P F} in
/-- ★★**`𝒞` の射の像を、脚を伸ばした代表元で書く**。 -/
theorem toHomPf_eq_mk {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A ⟶ B) (ψ : A' ⟶ B') (h : φ ≫ b = a ≫ ψ) :
    toHomPf (F := F) φ = HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) ψ := by
  have h1 : idxTransport P F (idxOneHom (F := F) a b ha hb hd) φ = ψ :=
    frobTransport_eq a ha b hb hd φ ψ h
  rw [← h1]
  exact (HomPf.mk_map (idxOneHom (F := F) a b ha hb hd) φ).symm

variable {P F} in
/-- ★★★★**`Hom^pf` の Frobenius 分解** ——
`x : A₀ ⟶ B₀`(`Hom^pf` の元、次数 `n`)は
**標準の Frobenius 拡大 `A₀ → A₀^{(n)}` の像**と**次数 `1` の射**の合成になる。 -/
theorem homPf_frobSplit {A₀ B₀ : C} (x : HomPf P F A₀ B₀) {n : ℕ+} (hn : pfDeg x = n) :
    ∃ z : HomPf P F (rtObj P F A₀ n) B₀,
      pfDeg z = 1 ∧ compPf P F (toHomPf (F := F) (rtExt P F A₀ n)) z = x := by
  obtain ⟨V, χ, hx⟩ := HomPf.exists_rep (P := P) (F := F) x
  have hχ : P.degFr χ = n := by rw [← hn, ← hx, pfDeg_mk]; rfl
  obtain ⟨hv1, hv2, hvd⟩ := V.hom.property
  obtain ⟨G, H, γ, β, a, hfac, hγ, hβ, ha⟩ := F.arbFactor χ
  have haLin : P.degFr a = 1 := (F.pullBackLB a ha).2
  have hρ : P.degFr (β ≫ a) = 1 := by
    rw [P.degFr_comp, haLin, show P.degFr β = 1 from hβ.1, mul_one]
  have hγn : P.degFr γ = n := by
    have h0 := congrArg P.degFr hfac
    rw [P.degFr_comp, hρ, one_mul, hχ] at h0
    exact h0.symm
  -- ★`A₀^{(n)}` の `d` 次拡大と `G` を合わせる
  have h1 : IsFrobeniusType P (rtExt P F A₀ n ≫
      rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1)) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A₀ n)
      (rtExt_frobType P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
  have h2 : IsFrobeniusType P (V.hom.hom.1 ≫ γ) := IsFrobeniusType.comp P F hv1 hγ
  have hdeg : P.degFr (rtExt P F A₀ n ≫
        rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
      = P.degFr (V.hom.hom.1 ≫ γ) := by
    rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr, hγn, mul_comm]
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A₀ _ G _ _ h1 h2 hdeg
  haveI := hθ
  obtain ⟨a₂, ha₂eq⟩ : ∃ t : rtObj P F A₀ n ⟶ G,
      t = rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1) ≫ θ := ⟨_, rfl⟩
  have ha₂F : IsFrobeniusType P a₂ := by
    rw [ha₂eq]
    exact IsFrobeniusType.comp P F
      (rtExt_frobType P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
      (isFrobeniusType_of_isIso P θ)
  have ha₂d : P.degFr a₂ = P.degFr V.hom.hom.1 := by
    rw [ha₂eq, P.degFr_comp, show P.degFr θ = 1 from isLinear_of_isIso P θ, rtExt_degFr,
      one_mul]
  have hsq : rtExt P F A₀ n ≫ a₂ = V.hom.hom.1 ≫ γ := by
    rw [ha₂eq, ← Category.assoc]
    exact hθe
  -- ★2 つの代表元
  have hy : toHomPf (F := F) (rtExt P F A₀ n)
      = HomPf.mk (idxMk (P := P) (F := F) V.hom.hom.1 a₂ hv1 ha₂F ha₂d.symm) γ :=
    toHomPf_eq_mk V.hom.hom.1 a₂ hv1 ha₂F ha₂d.symm _ γ hsq
  refine ⟨HomPf.mk (idxMk (P := P) (F := F) a₂ V.hom.hom.2 ha₂F hv2 (ha₂d.trans hvd))
      (β ≫ a), ?_, ?_⟩
  · rw [pfDeg_mk]; exact hρ
  · rw [hy]
    exact (compPf_mk_pair (P := P) (F := F) V.hom.hom.1 a₂ V.hom.hom.2
      hv1 ha₂F hv2 ha₂d.symm (ha₂d.trans hvd) γ (β ≫ a)).trans
      ((congrArg (HomPf.mk V) hfac.symm).trans hx)

/-! ## ★39. ★★★★`Λ_k` は関手である

★★`lamHom` を `rtJ` で書き直すと、合成則が `rtJ_comp` ＋ `toHomPf_comp` ＋
`rtMap_comp` の 3 本から出る。 -/

variable {P F} in
/-- ★添字の始対象から「押し出した始対象」への遷移射。 -/
noncomputable def idxPushOneHom {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    idxOne P F A B ⟶ (pushIdx (F := F) a ha b hb hd).obj (idxOne P F A' B') :=
  Under.homMk (show (⟨(A, B)⟩ : BiFr P F) ⟶ (⟨(A', B')⟩ : BiFr P F) from
    ⟨(a, b), ha, hb, hd⟩)
    (WideSubcategory.hom_ext _ (Prod.ext
      ((Category.id_comp _).trans (Category.comp_id _).symm)
      ((Category.id_comp _).trans (Category.comp_id _).symm)))

variable {P F} in
/-- ★★**根の不変性で `𝒞` の射の像がどう写るか**。 -/
theorem toHomPf_rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A ⟶ B) (ψ : A' ⟶ B') (h : φ ≫ b = a ≫ ψ) :
    (rootIso (F := F) a ha b hb hd).hom (toHomPf (F := F) ψ) = toHomPf (F := F) φ := by
  have h1 : idxTransport P F (idxPushOneHom (F := F) a ha b hb hd) φ = ψ :=
    frobTransport_eq a ha b hb hd φ ψ h
  refine Eq.trans (rootIso_hom_mk a ha b hb hd (idxOne P F A' B') ψ) ?_
  rw [← h1]
  exact HomPf.mk_map (idxPushOneHom (F := F) a ha b hb hd) φ

variable {P F} in
/-- ★★**`rtMap` は 1 段上げても四角形で繋がる**。 -/
theorem rtMap_rtLift (r : ℕ+) {A B : C} (φ : A ⟶ B) :
    rtMap (F := F) r φ ≫ rtLift P F B (show r * r = r * r from rfl)
      = rtLift P F A (show r * r = r * r from rfl) ≫ rtMap (F := F) (r * r) φ := by
  haveI : Epi (rtExt P F A r) := P.totEpiC _ _ _
  refine (cancel_epi (rtExt P F A r)).mp ?_
  rw [← Category.assoc, ← rtMap_spec, Category.assoc, rtLift_ext, rtMap_spec,
    ← Category.assoc, rtLift_ext]

variable {P F} in
/-- ★★★**`Λ_r` を `rtJ` で書く**。 -/
theorem lamHom_eq_rtJ (r : ℕ+) {A B : C} (φ : A ⟶ B) :
    lamHom (F := F) r φ
      = (rtJ P F A B r).hom (toHomPf (F := F) (rtMap (F := F) (r * r) φ)) :=
  (toHomPf_rootIso (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
    (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
    (by rw [rtLift_degFr, rtLift_degFr]) (rtMap (F := F) r φ)
    (rtMap (F := F) (r * r) φ) (rtMap_rtLift (F := F) r φ)).symm

variable {P F} in
/-- ★★★**`Λ_k` は合成を保つ**。 -/
theorem lamHom_comp (k : ℕ+) {A B E : C} (φ : A ⟶ B) (ψ : B ⟶ E) :
    lamHom (F := F) k φ ≫ lamHom (F := F) k ψ = lamHom (F := F) k (φ ≫ ψ) := by
  rw [lamHom_eq_rtJ, lamHom_eq_rtJ, lamHom_eq_rtJ, rtJ_comp, ← toHomPf_comp, rtMap_comp]

variable {P F} in
/-- ★★★**`Λ_k` は恒等射を保つ**。 -/
theorem lamHom_id (k : ℕ+) (A : C) :
    lamHom (F := F) k (𝟙 A) = 𝟙 (⟨A, k⟩ : PfRootObj P F) := by
  show HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F A k)) (rtMap (F := F) k (𝟙 A))
    = idRoot P F ⟨A, k⟩
  rw [rtMap_id]
  rfl

/-! ## ★40. (iv)(b) の第 1 歩 —— pull-back は linear(根が等しい場合)

★★手は `𝒞^birat` の `birat_pullBack_repr` と同じ:
**`α` を「Frobenius 型 ≫ 次数 1」に分解し、次数 1 の側を `α` に沿って持ち上げる**と
`1 = degFr α · degFr g` になる。
★`𝒞^pf` では分解の中間対象を **`(A^{(n)}, r)`** に取れる(`homPf_frobSplit` の
吸収先が標準拡大 `A₀^{(n)}` だから)——ここが根つきでも閉じる理由である。 -/

variable {P F} in
/-- ★`rtJ` の像の次数(代表元を経由しない形)。 -/
theorem rtJ_degFr' {A B : C} (r : ℕ+)
    (x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r))) :
    (pfRootPre P F).degFr
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from (rtJ P F A B r).hom x)
      = pfDeg x := by
  obtain ⟨V, φ, hφ⟩ := HomPf.exists_rep (P := P) (F := F) x
  rw [← hφ, rtJ_degFr, pfDeg_mk]
  rfl

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**pull-back は linear**(根が等しい場合)。 -/
theorem pfRoot_isLinear_of_pullBack_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hα : IsPullBack (pfRootPre P F) α) : IsLinear (pfRootPre P F) α := by
  obtain ⟨x, hx⟩ : ∃ x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)),
      (rtJ P F A B r).hom x = α := ⟨(rtJ P F A B r).inv α, Iso.inv_hom_id_apply _ _⟩
  obtain ⟨n, hn⟩ : ∃ n : ℕ+, pfDeg x = n := ⟨_, rfl⟩
  obtain ⟨z, hz1, hz2⟩ := homPf_frobSplit (P := P) (F := F) x hn
  -- ★`A^{(n)}` を `r·r` 乗したものと、`A^{(r·r)}` を `n` 乗したものを合わせる
  have h1 : IsFrobeniusType P
      (rtExt P F A n ≫ rtExt P F (rtObj P F A n) (r * r)) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A n)
      (rtExt_frobType P F (rtObj P F A n) (r * r))
  have h2 : IsFrobeniusType P
      (rtExt P F A (r * r) ≫ rtExt P F (rtObj P F A (r * r)) n) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A (r * r))
      (rtExt_frobType P F (rtObj P F A (r * r)) n)
  have hdeg : P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) (r * r))
      = P.degFr (rtExt P F A (r * r) ≫ rtExt P F (rtObj P F A (r * r)) n) := by
    rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr, rtExt_degFr, rtExt_degFr,
      mul_comm]
  obtain ⟨w, hw, -⟩ := F.frobDegUniq A _ _ _ _ h1 h2 hdeg
  haveI := hw
  -- ★分解の 2 本
  obtain ⟨gam, hgam⟩ : ∃ t : rtObj P F A (r * r) ⟶ rtObj P F (rtObj P F A n) (r * r),
      t = rtExt P F (rtObj P F A (r * r)) n ≫ @inv _ _ _ _ w hw := ⟨_, rfl⟩
  obtain ⟨Gam, hGam⟩ : ∃ t : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨rtObj P F A n, r⟩,
      t = (rtJ P F A (rtObj P F A n) r).hom (toHomPf (F := F) gam) := ⟨_, rfl⟩
  obtain ⟨R, hR⟩ : ∃ t : (⟨rtObj P F A n, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      t = (rtJ P F (rtObj P F A n) B r).hom
        (compPf P F (toHomPf (F := F) w) z) := ⟨_, rfl⟩
  -- ★合成すると `α`
  have hcomp : Gam ≫ R = α := by
    rw [hGam, hR, rtJ_comp, ← compPf_assoc, ← toHomPf_comp, hgam, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id, hz2, hx]
  -- ★`R` は次数 1
  have hRdeg : (pfRootPre P F).degFr R = 1 := by
    rw [hR, rtJ_degFr', pfDeg_comp, hz1, pfDeg_toHomPf,
      show P.degFr w = 1 from isLinear_of_isIso P w, mul_one]
  -- ★`Gam` の底は同型
  have hGamBase : IsIso ((pfRootPre P F).Base Gam) := by
    have hb : IsBaseIsomorphism P gam := by
      haveI : IsIso (@inv _ _ _ _ w hw) := @IsIso.inv_isIso _ _ _ _ w hw
      rw [hgam]
      exact (IsFrobeniusType.comp P F (rtExt_frobType P F (rtObj P F A (r * r)) n)
        (isFrobeniusType_of_isIso P (@inv _ _ _ _ w hw))).2
    have := (rtJ_isBaseIso_iff (F := F) (A := A) (B := rtObj P F A n) r
      (idxOne P F (rtObj P F A (r * r)) (rtObj P F (rtObj P F A n) (r * r))) gam).mpr hb
    rw [hGam]
    exact this
  -- ★持ち上げ
  have hbase : (pfRootPre P F).Base R
      = @inv _ _ _ _ ((pfRootPre P F).Base Gam) hGamBase ≫ (pfRootPre P F).Base α := by
    refine eq_inv_comp_of ((pfRootPre P F).Base Gam) hGamBase _ _ ?_
    rw [← (pfRootPre P F).Base_comp, hcomp]
  obtain ⟨g, ⟨hg1, -⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα
    (⟨rtObj P F A n, r⟩ : PfRootObj P F) R _ hbase
  have hmul : (pfRootPre P F).degFr α * (pfRootPre P F).degFr g = 1 := by
    rw [← hRdeg, ← hg1, (pfRootPre P F).degFr_comp]
  have hcoe : (((pfRootPre P F).degFr α : ℕ+) : ℕ)
      * (((pfRootPre P F).degFr g : ℕ+) : ℕ) = 1 := by
    rw [← PNat.mul_coe, hmul]; rfl
  show (pfRootPre P F).degFr α = 1
  exact PNat.coe_injective (Nat.dvd_one.mp ⟨_, hcoe.symm⟩)

/-! ## ★41. ★★★★根が違ってもよい合成の**組み立て**規則

★`compRoot_rep` は「合成を分解する」向きだった。★ここでは**逆向き** ——
3 脚の添字と 2 本の `𝒞` の射から、`𝒞^pf` の 2 射とその合成を**作る**。
★★これがあれば `arbFactor` の中間対象を根の違う対象に取れる。 -/

variable {P F} in
/-- ★★★★**3 脚の添字から `𝒞^pf` の合成を組み立てる**。 -/
theorem compRoot_mk3 {X Y Z : PfRootObj P F} {A' B' E' : C}
    (a : rtObj P F X.obj (Z.root * Y.root) ⟶ A')
    (b : rtObj P F Y.obj (Z.root * X.root) ⟶ B')
    (e : rtObj P F Z.obj (Y.root * X.root) ⟶ E')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (he : IsFrobeniusType P e)
    (hab : P.degFr a = P.degFr b) (hbe : P.degFr b = P.degFr e)
    (φ : A' ⟶ B') (ψ : B' ⟶ E') :
    (show (X ⟶ Y) from (rtRootIso P F X.obj Y.obj
        (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hab) φ))
      ≫ (show (Y ⟶ Z) from (rtRootIso P F Y.obj Z.obj
        (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).hom
      (HomPf.mk (idxMk (P := P) (F := F) b e hb he hbe) ψ))
      = (rtRootIso P F X.obj Z.obj
          (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
          (show Y.root * X.root = Y.root * X.root from rfl)).hom
        (HomPf.mk (idxMk (P := P) (F := F) a e ha he (hab.trans hbe)) (φ ≫ ψ)) := by
  show compRoot P F _ _ = _
  unfold compRoot
  rw [Iso.hom_inv_id_apply, Iso.hom_inv_id_apply, compPf_mk_pair]

/-! ## ★42. 代表元の脚を `k` 乗ぶん伸ばす -/

variable {P F} in
/-- ★添字の脚を `k` 乗ぶん伸ばす遷移射。 -/
noncomputable def idxExtHom {A B : C} (V : IdxPf P F A B) (k : ℕ+) :
    V ⟶ idxMk (P := P) (F := F)
      (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k)
      (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k)
      (IsFrobeniusType.comp P F V.hom.property.1 (rtExt_frobType P F _ k))
      (IsFrobeniusType.comp P F V.hom.property.2.1 (rtExt_frobType P F _ k))
      (by rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr]
          exact congrArg (fun m => k * m) V.hom.property.2.2) :=
  Under.homMk (show V.right ⟶ (⟨(rtObj P F V.right.obj.1 k,
      rtObj P F V.right.obj.2 k)⟩ : BiFr P F) from
    ⟨(rtExt P F V.right.obj.1 k, rtExt P F V.right.obj.2 k),
      rtExt_frobType P F _ k, rtExt_frobType P F _ k,
      by rw [rtExt_degFr, rtExt_degFr]⟩)
    (WideSubcategory.hom_ext _ (Prod.ext rfl rfl))

variable {P F} in
/-- ★★**脚を伸ばした代表元** —— 射は `rtMap` で持ち上がる。 -/
theorem mk_rtExt {A B : C} (V : IdxPf P F A B) (k : ℕ+)
    (χ : V.right.obj.1 ⟶ V.right.obj.2) :
    HomPf.mk (idxMk (P := P) (F := F)
      (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k)
      (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k)
      (IsFrobeniusType.comp P F V.hom.property.1 (rtExt_frobType P F _ k))
      (IsFrobeniusType.comp P F V.hom.property.2.1 (rtExt_frobType P F _ k))
      (by rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr]
          exact congrArg (fun m => k * m) V.hom.property.2.2))
      (rtMap (F := F) k χ)
      = HomPf.mk V χ := by
  have h1 : idxTransport P F (idxExtHom (F := F) V k) χ = rtMap (F := F) k χ :=
    frobTransport_eq _ _ _ _ _ χ _ (rtMap_spec (F := F) k χ)
  rw [← h1]
  exact HomPf.mk_map (idxExtHom (F := F) V k) χ

variable {P F} in
/-- ★★**脚を伸ばした代表元(伸ばし先を与える形)**。 -/
theorem mk_rtExt_gen {A B : C} (V : IdxPf P F A B) (k : ℕ+)
    (a : A ⟶ rtObj P F V.right.obj.1 k) (b : B ⟶ rtObj P F V.right.obj.2 k)
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (hva : V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k = a)
    (hvb : V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k = b)
    (χ : V.right.obj.1 ⟶ V.right.obj.2) :
    HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) (rtMap (F := F) k χ)
      = HomPf.mk V χ := by
  refine Eq.trans ?_ (HomPf.mk_map
    (Under.homMk (show V.right ⟶ (⟨(rtObj P F V.right.obj.1 k,
        rtObj P F V.right.obj.2 k)⟩ : BiFr P F) from
      ⟨(rtExt P F V.right.obj.1 k, rtExt P F V.right.obj.2 k),
        rtExt_frobType P F _ k, rtExt_frobType P F _ k,
        by rw [rtExt_degFr, rtExt_degFr]⟩)
      (WideSubcategory.hom_ext _ (Prod.ext hva hvb)) :
        V ⟶ idxMk (P := P) (F := F) a b ha hb hd) χ)
  exact congrArg (HomPf.mk _)
    (frobTransport_eq _ _ _ _ _ χ _ (rtMap_spec (F := F) k χ)).symm

/-! ## ★43. ★★★★linear な射は「pre-step ≫ 等長」に分解する

★★**中間対象の根が変わる**のが要点である。`X = (A,r)`、`Z = (B,r)` で
代表元の脚の次数が `c` のとき、中間対象は **`(Cc, r·r·c)`** に取る。
★`𝒞` の `Proposition 1.7, (iii)` の分解 `χ = α₁ ≫ γ₁` を
**`rtMap (r·r)` で丸ごと持ち上げ**、3 脚の添字の脚を**同型(次数 1)**に取る。 -/

variable {P F} in
/-- ★`idxMk` の脚が等しければ代表元も等しい。 -/
theorem mk_idxMk_congr {A B A' B' : C} {a a' : A ⟶ A'} {b b' : B ⟶ B'}
    (ha : a = a') (hb : b = b')
    {ha1 : IsFrobeniusType P a} {hb1 : IsFrobeniusType P b} {hd1 : P.degFr a = P.degFr b}
    {ha2 : IsFrobeniusType P a'} {hb2 : IsFrobeniusType P b'}
    {hd2 : P.degFr a' = P.degFr b'} (φ : A' ⟶ B') :
    HomPf.mk (idxMk (P := P) (F := F) a b ha1 hb1 hd1) φ
      = HomPf.mk (idxMk (P := P) (F := F) a' b' ha2 hb2 hd2) φ := by
  subst ha
  subst hb
  rfl

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**linear ⟹ 「pre-step ≫ 等長」**(根が等しい場合)。 -/
theorem pfRoot_linearSplit_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) (hlin : IsLinear (pfRootPre P F) α) :
    ∃ (T : PfRootObj P F) (pp : (⟨A, r⟩ : PfRootObj P F) ⟶ T) (ll : T ⟶ ⟨B, r⟩),
      IsPreStep (pfRootPre P F) pp ∧ IsIsometric (pfRootPre P F) ll ∧ α = pp ≫ ll := by
  obtain ⟨V, χ, hV⟩ := HomPf.exists_rep (P := P) (F := F)
    (show HomPf P F (rtObj P F A r) (rtObj P F B r) from α)
  have hχlin : IsLinear P χ :=
    (rootDeg_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F)) V χ).symm.trans
      (by rw [hV]; exact hlin)
  obtain ⟨Cc, α₁, γ₁, hfac, hα₁, hγ₁⟩ := (prop_1_7_iii_linear_factor P F χ).mp hχlin
  obtain ⟨hv1, hv2, hvd⟩ := V.hom.property
  obtain ⟨t, ht⟩ : ∃ t : ℕ+, t = r * r * P.degFr V.hom.hom.1 := ⟨_, rfl⟩
  -- ★2 本の脚(次数 1、すなわち同型)
  have hd1 : P.degFr (rtLift P F A (show r * t = t * r from mul_comm _ _))
      = P.degFr (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 (r * r)) := by
    rw [rtLift_degFr, P.degFr_comp, rtExt_degFr, ht]
  have hd3 : P.degFr (rtLift P F B (show t * r = t * r from rfl))
      = P.degFr (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 (r * r)) := by
    rw [rtLift_degFr, P.degFr_comp, rtExt_degFr, ht, hvd]
  obtain ⟨l₁, hl₁iso, hl₁⟩ := F.frobDegUniq (rtObj P F A r) _ _
    (rtLift P F A (show r * t = t * r from mul_comm _ _))
    (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 (r * r))
    (rtLift_frobType P F A _)
    (IsFrobeniusType.comp P F hv1 (rtExt_frobType P F _ (r * r))) hd1
  obtain ⟨l₃, hl₃iso, hl₃⟩ := F.frobDegUniq (rtObj P F B r) _ _
    (rtLift P F B (show t * r = t * r from rfl))
    (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 (r * r))
    (rtLift_frobType P F B _)
    (IsFrobeniusType.comp P F hv2 (rtExt_frobType P F _ (r * r))) hd3
  haveI := hl₁iso
  haveI := hl₃iso
  have hl₁F : IsFrobeniusType P l₁ := isFrobeniusType_of_isIso P l₁
  have hl₃F : IsFrobeniusType P l₃ := isFrobeniusType_of_isIso P l₃
  have hl₂F : IsFrobeniusType P (𝟙 (rtObj P F Cc (r * r))) :=
    isFrobeniusType_of_isIso P (𝟙 _)
  have hab : P.degFr l₁ = P.degFr (𝟙 (rtObj P F Cc (r * r))) := by
    rw [show P.degFr l₁ = 1 from isLinear_of_isIso P l₁, P.degFr_id]
  have hbe : P.degFr (𝟙 (rtObj P F Cc (r * r))) = P.degFr l₃ := by
    rw [show P.degFr l₃ = 1 from isLinear_of_isIso P l₃, P.degFr_id]
  refine ⟨⟨Cc, t⟩, ?_, ?_, ?_, ?_, ?_⟩
  · exact (rtRootIso P F A Cc (show r * t = r * t from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) l₁ (𝟙 (rtObj P F Cc (r * r))) hl₁F hl₂F hab)
        (rtMap (F := F) (r * r) α₁))
  · exact (rtRootIso P F Cc B (show r * r = r * r from mul_comm _ _)
      (show t * r = r * t from mul_comm _ _)).hom
      (HomPf.mk (idxMk (P := P) (F := F) (𝟙 (rtObj P F Cc (r * r))) l₃ hl₂F hl₃F hbe)
        (rtMap (F := F) (r * r) γ₁))
  · rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨Cc, t⟩ : PfRootObj P F)) _ _).mpr
      ⟨(rtMap_degFr (F := F) (r * r) α₁).trans hα₁.1,
        (rtMap_isBaseIso_iff (F := F) (r * r) α₁).mpr hα₁.2⟩
  · rw [rtRootIso_hom_mk]
    exact (isIsometric_mk_iff (X := (⟨Cc, t⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) _ _).mpr
      ((rtMap_isIsometric_iff (F := F) (r * r) γ₁).mpr (F.pullBackLB γ₁ hγ₁).1.2)
  · have hcc : rtMap (F := F) (r * r) α₁ ≫ rtMap (F := F) (r * r) γ₁
        = rtMap (F := F) (r * r) χ := by rw [← rtMap_comp, ← hfac]
    refine Eq.trans ?_ (compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨Cc, t⟩ : PfRootObj P F)) (Z := (⟨B, r⟩ : PfRootObj P F))
      l₁ (𝟙 (rtObj P F Cc (r * r))) l₃ hl₁F hl₂F hl₃F hab hbe
      (rtMap (F := F) (r * r) α₁) (rtMap (F := F) (r * r) γ₁)).symm
    refine Eq.trans ?_ (congrArg
      (fun u : rtObj P F V.right.obj.1 (r * r) ⟶ rtObj P F V.right.obj.2 (r * r) =>
        (rtRootIso P F A B (show r * t = t * r from mul_comm _ _)
          (show t * r = t * r from rfl)).hom
          (HomPf.mk (idxMk (P := P) (F := F) l₁ l₃ hl₁F hl₃F (hab.trans hbe)) u))
      hcc.symm)
    rw [rtRootIso_hom_mk]
    symm
    exact (mk_rtExt_gen (F := F) V (r * r)
      (rtLift P F A (show r * t = t * r from mul_comm _ _) ≫ l₁)
      (rtLift P F B (show t * r = t * r from rfl) ≫ l₃)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl₁F)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hl₃F)
      (by rw [P.degFr_comp, P.degFr_comp,
          show P.degFr l₁ = 1 from isLinear_of_isIso P l₁,
          show P.degFr l₃ = 1 from isLinear_of_isIso P l₃,
          rtLift_degFr, rtLift_degFr])
      hl₁.symm hl₃.symm χ).trans hV

/-! ## ★44. (iv)(b) —— pull-back は LB-invertible かつ linear

★★等長性の手も `𝒞^birat` と同じ: `α = pp ≫ ll`(pre-step ≫ 等長)に分解し、
`ll` を `α` に沿って持ち上げると **`pp` が同型**になる。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**pull-back は等長**(根が等しい場合)。 -/
theorem pfRoot_isIsometric_of_pullBack_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hα : IsPullBack (pfRootPre P F) α) : IsIsometric (pfRootPre P F) α := by
  obtain ⟨T, pp, ll, hpp, hll, heq⟩ := pfRoot_linearSplit_sameRoot α
    (pfRoot_isLinear_of_pullBack_sameRoot α hα)
  haveI hppb : IsIso ((pfRootPre P F).Base pp) := hpp.2
  have hbase : (pfRootPre P F).Base ll
      = @inv _ _ _ _ ((pfRootPre P F).Base pp) hppb ≫ (pfRootPre P F).Base α := by
    refine eq_inv_comp_of _ hppb _ _ ?_
    rw [← (pfRootPre P F).Base_comp, ← heq]
  obtain ⟨u, ⟨hu1, hu2⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα T ll _ hbase
  have hkey : pp ≫ u = 𝟙 (⟨A, r⟩ : PfRootObj P F) := by
    refine IsPullBack.hom_ext (pfRootPre P F) hα _ _ ?_ ?_
    · rw [Category.assoc, hu1, ← heq, Category.id_comp]
    · rw [(pfRootPre P F).Base_comp, hu2, IsIso.hom_inv_id, (pfRootPre P F).Base_id]
  haveI : IsIso pp := isIso_of_comp_eq_id (pfRoot_totEpi P F) pp u hkey
  show (pfRootPre P F).Div α = 0
  rw [heq, (pfRootPre P F).Div_comp, show (pfRootPre P F).Div ll = 0 from hll,
    show (pfRootPre P F).Div pp = 0 from isIsometric_of_isIso (pfRootPre P F) pp,
    map_zero, smul_zero, add_zero]

variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iv)(b)** —— `𝒞^pf` 版。 -/
theorem pfRoot_pullBackLB (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (α : X ⟶ Y) (hα : IsPullBack (pfRootPre P F) α) :
    IsLBInvertible (pfRootPre P F) α ∧ IsLinear (pfRootPre P F) α := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hpb' : IsPullBack (pfRootPre P F) (inv eX ≫ α ≫ eY) :=
    IsPullBack.comp (pfRootPre P F) (isPullBack_of_isIso (pfRootPre P F) (inv eX))
      (IsPullBack.comp (pfRootPre P F) hα (isPullBack_of_isIso (pfRootPre P F) eY))
  have hlin' := pfRoot_isLinear_of_pullBack_sameRoot (inv eX ≫ α ≫ eY) hpb'
  have hiso' := pfRoot_isIsometric_of_pullBack_sameRoot (inv eX ≫ α ≫ eY) hpb'
  have hlin : IsLinear (pfRootPre P F) α := by
    have h0 : (pfRootPre P F).degFr (inv eX ≫ α ≫ eY) = 1 := hlin'
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eX) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eX),
      show (pfRootPre P F).degFr eY = 1 from isLinear_of_isIso (pfRootPre P F) eY,
      one_mul, mul_one] at h0
    exact h0
  refine ⟨⟨pfRoot_isCoAngular hfi α, ?_⟩, hlin⟩
  have h1 : (pfRootPre P F).Div (inv eX ≫ α ≫ eY) = 0 := hiso'
  rw [(pfRootPre P F).Div_comp, (pfRootPre P F).Div_comp,
    show (pfRootPre P F).Div (inv eX) = 0 from
      isIsometric_of_isIso (pfRootPre P F) (inv eX),
    show (pfRootPre P F).Div eY = 0 from isIsometric_of_isIso (pfRootPre P F) eY,
    smul_zero, add_zero, map_zero, zero_add,
    show (pfRootPre P F).degFr eY = 1 from isLinear_of_isIso (pfRootPre P F) eY] at h1
  simp only [PNat.one_coe, one_smul] at h1
  show (pfRootPre P F).Div α = 0
  exact (Φ.pfOn (phiSharp P)).map_injective _ (h1.trans (map_zero _).symm)

/-! ## ★45. pull-back を `Hom^pf` の中で消す

★★残る 3 条(`arbFactor` / `arbFactorUniq` / `plBkEquiv`)はすべて
**「`𝒞` の pull-back は `𝒞^pf` でも pull-back」**に帰する。
★その普遍性のうち**一意性の側**をここで用意する ——
`homPf_cancel_preStep` と同じ形で、`Mono` の代わりに
`IsPullBack.hom_ext`(底の一致も要る)を使う。 -/

variable {P F} in
/-- ★★**pull-back 性は Frobenius 遷移で戻せる**(co-angular が要る)。

★`prop_1_10_i_pullBack_of` の逆向き。★`Proposition 1.4, (ii)` を両向きに使う:
linear と等長は四角形から**両向きに**出るが、co-angular だけは仮定に置く
(`𝒞^pf` の側では始域を isotropic に取れるので、そこはただで手に入る)。 -/
theorem pullBack_of_transport_back (Fc : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {a : A ⟶ A'} {b : B ⟶ B'} {φ' : A' ⟶ B'}
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b)
    (hd : P.degFr a = P.degFr b) (hsq : φ ≫ b = a ≫ φ')
    (hφ' : IsPullBack P φ') (hco : IsCoAngular P φ) : IsPullBack P φ := by
  obtain ⟨hlb', hlin'⟩ := Fc.pullBackLB φ' hφ'
  have hlin : IsLinear P φ := by
    have h0 := congrArg P.degFr hsq
    rw [P.degFr_comp, P.degFr_comp, show P.degFr φ' = 1 from hlin', one_mul, hd] at h0
    show P.degFr φ = 1
    exact mul_left_cancel (h0.trans (mul_one _).symm)
  refine prop_1_4_ii_mpr P Fc φ ⟨hco, ?_⟩ hlin
  have h0 := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, show P.Div b = 0 from hb.1.2,
    show P.Div a = 0 from ha.1.2, show P.Div φ' = 0 from hlb'.2,
    map_zero, map_zero, zero_add, smul_zero, add_zero] at h0
  exact nsmul_eq_zero_of_isSharp (P.divisorial _).2 h0

variable {P F} in
/-- ★遷移は底の一致を保つ。 -/
theorem idxTransport_Base_eq {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ φ' : Z.right.obj.1 ⟶ Z.right.obj.2) (h : P.Base φ = P.Base φ') :
    P.Base (idxTransport P F u φ) = P.Base (idxTransport P F u φ') := by
  haveI : IsIso (P.Base u.right.hom.1) := u.right.property.1.2
  refine (cancel_epi (P.Base u.right.hom.1)).mp ?_
  rw [← P.Base_comp, ← P.Base_comp, ← idxTransport_spec, ← idxTransport_spec,
    P.Base_comp, P.Base_comp, h]

variable {P F} in
/-- ★3 脚添字を「**第 2 脚**が isotropic」な所まで押し上げる。 -/
theorem exists_idx3_isotropic2 (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E) :
    ∃ (W : IdxPf3 P F A B E) (u : V ⟶ W), IsIsotropic P W.right.obj.2.1 := by
  obtain ⟨Dd, a, ha, hDd⟩ := hfi V.right.obj.2.1
  obtain ⟨A₂, p, hp, hpd⟩ := F.frobDegSurj V.right.obj.1 (P.degFr a)
  obtain ⟨E₂, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr a)
  obtain ⟨hva, hvb, hve, hvab, hvbe⟩ := V.hom.property
  refine ⟨Under.mk (Y := (⟨(A₂, Dd, E₂)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(V.hom.hom.1 ≫ p, V.hom.hom.2.1 ≫ a, V.hom.hom.2.2 ≫ e),
          IsFrobeniusType.comp P F hva hp, IsFrobeniusType.comp P F hvb ha,
          IsFrobeniusType.comp P F hve he, ?_, ?_⟩),
    Under.homMk (show V.right ⟶ (⟨(A₂, Dd, E₂)⟩ : TriFr P F) from
      ⟨(p, a, e), hp, ha, he, hpd, hed.symm⟩)
      (WideSubcategory.hom_ext _ rfl), hDd⟩
  · rw [P.degFr_comp, P.degFr_comp, hpd, hvab]
  · rw [P.degFr_comp, P.degFr_comp, hed, hvbe]

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**`Hom^pf` の中で pull-back は消せる**(底の一致つき)。 -/
theorem homPf_cancel_pullBack {A B E : C} (V : IdxPf3 P F A B E)
    (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hψ : IsPullBack P ψ) (hb : P.Base φ = P.Base φ')
    (h : HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ' ≫ ψ)) :
    HomPf.mk ((idx12 P F A B E).obj V) φ = HomPf.mk ((idx12 P F A B E).obj V) φ' := by
  obtain ⟨V', t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  obtain ⟨V'', ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) V'
  set s : V ⟶ IsFiltered.max V V'' := IsFiltered.leftToMax V V'' with hs
  set r : V'' ⟶ IsFiltered.max V V'' := IsFiltered.rightToMax V V'' with hr
  have hm : t ≫ k ≫ (idx13 P F A B E).map r = (idx13 P F A B E).map s :=
    idx_hom_ext _ _
  have hA : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have hB : idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ' ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have key : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ) := by
    rw [hA, hB, ht]
  have e1 := idxTransport_comp_pair (F := F) s φ ψ
  have e2 := idxTransport_comp_pair (F := F) s φ' ψ
  have hcan : idxTransport P F ((idx12 P F A B E).map s) φ
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ
      = idxTransport P F ((idx12 P F A B E).map s) φ'
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ := by
    rw [e1, e2]
    exact key
  have hψ' : IsPullBack P (idxTransport P F ((idx23 P F A B E).map s) ψ) :=
    prop_1_10_i_pullBack_of P F
      ((idx23 P F A B E).map s).right.property.1
      ((idx23 P F A B E).map s).right.property.2.1
      ((idx23 P F A B E).map s).right.property.2.2
      (idxTransport_spec ((idx23 P F A B E).map s) ψ) hψ
  have hb' : P.Base (idxTransport P F ((idx12 P F A B E).map s) φ)
      = P.Base (idxTransport P F ((idx12 P F A B E).map s) φ') :=
    idxTransport_Base_eq ((idx12 P F A B E).map s) φ φ' hb
  have hφφ' : idxTransport P F ((idx12 P F A B E).map s) φ
      = idxTransport P F ((idx12 P F A B E).map s) φ' :=
    IsPullBack.hom_ext P hψ' _ _ hcan hb'
  rw [← HomPf.mk_map ((idx12 P F A B E).map s) φ, ← HomPf.mk_map ((idx12 P F A B E).map s) φ',
    hφφ']

/-! ## ★32. ★★★★現在地と残り(2026-08-17 の測定)

### 埋まった 18 条(`FrobenioidCore (pfRootPre P F)` の 21 条のうち)

| 条 | 宣言 | 手 |
|---|---|---|
| (i)(a) `baseSurj` | `pfRoot_baseSurj` | `(A,1)` として送るだけ |
| (i)(b) `preStepSpan` | `pfRoot_preStepSpan` | ★span を `A^m`・`B^n` の間で取り `W := (V₀, n*m)` |
| (ii) `frobDegSurj` | `pfRoot_frobDegSurj` | `exists_rt_frob` |
| (ii) `frobDegUniq` | `pfRoot_frobDegUniq` | co-span を揃え、isotropic に押し上げる |
| (iii)(a) `coAngularComp` | `pfRoot_coAngularComp` | ★isotropic 型ゆえ全射が co-angular |
| (iii)(b) `coAngularOfPreStep` | `pfRoot_coAngularOfPreStep` | 同上 |
| (iii)(c) `otriFwd` | `pfRoot_otriFwd` | 根を揃える ＋ 対角 ＋ `𝒞` の `otriFwd` |
| (iii)(c) `otriBwd` | `pfRoot_otriBwd` | 同上 |
| (iii)(c) `otriBase` | `pfRoot_otriBase` | 同上 |
| (v)(a) `preStepMono` | `pfRoot_preStepMono` | `homPf_cancel_preStep` |
| (v)(b) `preStepFactor` | `pfRoot_preStepFactor` | `φ = φ ≫ 𝟙` |
| (v)(b) `preStepFactorUniq` | `pfRoot_preStepFactorUniq` | isometric pre-step は同型 |
| (v)(c) `preStepFactor'` | `pfRoot_preStepFactor'` | `φ = 𝟙 ≫ φ` |
| (v)(c) `preStepFactorUniq'` | `pfRoot_preStepFactorUniq'` | 同上 |
| (vi) `faithfulUpToUnits` | `pfRoot_faithfulUpToUnits` | ★根を揃える ＋ `rootDiv` を `Div` に降ろす(divisorial) |
| (iv)(b) `pullBackLB` | `pfRoot_pullBackLB` | ★★Frobenius 分解 ＋ linear 分解の 2 段で持ち上げる |
| (vii)(a) `isotropicHullExists` | `pfRoot_isotropicHullExists` | `𝟙` |
| (vii)(b) `isotropicClosed` | `pfRoot_isotropicClosed` | 全対象 isotropic |

★土台は `pfRoot_isOfIsotropicType`(`𝒞` が Frobenius-isotropic 型 ⟹ `𝒞^pf` は isotropic 型)。

### 残る 3 条と、それぞれの**具体的な**残り作業

★**(iv)(b) `pullBackLB` は 2026-08-17 に埋まった**。手は 2 段:

1. ★**linear** —— `homPf_frobSplit` で `α = Γ ≫ R`(`Γ` は Frobenius 型、`R` は次数 1)。
   中間対象は **`(A^{(n)}, r)`** に取れる(`homPf_frobSplit` の吸収先が標準拡大だから)。
   `R` を `α` に沿って持ち上げると `1 = degFr α · degFr g` になる。
2. ★★**等長** —— `pfRoot_linearSplit_sameRoot` で `α = pp ≫ ll`
   (`pp` は pre-step、`ll` は等長)。★ここは**中間対象の根が変わる**:
   代表元の脚の次数が `c` なら中間対象は **`(Cc, r·r·c)`** で、
   `𝒞` の `Proposition 1.7, (iii)` の分解を **`rtMap (r·r)` で丸ごと持ち上げ**、
   3 脚の添字の脚を**同型(次数 1)**に取る(`compRoot_mk3`)。
   `ll` を `α` に沿って持ち上げると `pp ≫ u = 𝟙`、全射性で `pp` が同型になり
   `Div α = 0` が出る。

**(iv)(a) `arbFactor` / (iv)(a) `arbFactorUniq` / (i)(c) `plBkEquiv`**
—— 残る 3 条はどれも ★★**「`𝒞` の pull-back は `𝒞^pf` でも pull-back」**
(`pfRoot_isPullBack_map`)を要る。
★`𝒞^birat` では `birat_isPullBack_map` として `isPullBack_of_lift`
(`Prop44Core.lean` の汎用道具、pre-Frobenioid 一般)＋ 底の公式で出してある。
★`𝒞^pf` 版では **`rtMap` が pull-back を保つこと**(`PlBkShuffle.lean` の
基底変換)と、`compRoot_mk3` による持ち上げの構成が要る。

★`arbFactor` そのものは上の 2 段(`homPf_frobSplit` ＋ `pfRoot_linearSplit_sameRoot`)を
つなぐだけで形は出ており、**足りないのは `ll` が pull-back だと言う所だけ**である。
★`arbFactorUniq` は `𝒞^birat` と同じく
`frobDegUniq` ＋ 全射性 ＋ `IsPullBack.lift` の 3 点で出る。

### `Proposition 3.2` 全体として残るもの

1. 上の 3 条(＋ `Frobenioid` の 2 本の圏同値 `coaPreUnderEquiv` / `coaPreOverEquiv`)
2. **(ii) の辞書の残り** —— Frobenius 型・pull-back・co-angular・
   base-identity 自己射・同型(`isPreStep` など 5 項は `Prop32.lean` で済み)
3. **(iii) の後半** —— `𝒞^pf ≃ (𝒞^pf)^pf`
-/

end ABC3.Found.FrdI
