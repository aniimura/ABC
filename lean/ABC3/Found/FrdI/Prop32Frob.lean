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

end ABC3.Found.FrdI
