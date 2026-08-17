import ABC3.Found.FrdI.Prop32

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

end ABC3.Found.FrdI
