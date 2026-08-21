/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Equiv
import ABC3.Found.FrdI.Prop55Birat
import ABC3.Found.FrdI.Prop55PfKappa
import ABC3.Found.FrdI.Thm52Rem272

/-!
# [FrdI] Proposition 5.5, (ii) —— 左辺 `(𝒞^pf)^birat` の側

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★測って分かった正体

模型 Frobenioid で計算すると、`(𝒞^pf)^birat` と `(𝒞^birat)^pf` は
**どちらも同じ 1 本の余極限** `colim_k Hom^birat(A^{(k)}, B^{(k)})` である。

* 右辺は `idxToBirat_final`(`Prop55Birat.lean`)で**すでにその形**になっている ——
  `colim_{W ∈ IdxPf(𝒞)(A,B)} Hom^birat(W₁, W₂)`。
* 左辺は、外の添字(`𝒞^pf` の co-angular pre-step `⟨A″,n⟩ → ⟨A,1⟩`)と
  内の `Hom^pf` を**共通の根に揃える**と、同じスパン `A^{(k)} ← X → B^{(k)}` になる。

## ★★★本ファイルの一歩 —— 右辺の代表元から左辺の元を作る

★★鍵は **`⟨A,1⟩ ≅ ⟨A₁,k⟩`** である ——
`α : A ⟶ A₁` が次数 `k` の Frobenius 型なら、`𝒞^pf` では
**`A₁` の `k` 乗根が `A`** である。
★★★**その同型は方程式 `e ≫ κ = [α]` で特徴づける**(`rootLift`、`Prop55PfKappa.lean`)——
`.choose` で取ると `W` を大きくしたときの整合性が言えないからである。

| 定理 | 中身 |
|---|---|
| `compBirat_mk_of_sq` | ★★合成は**どの引き戻しデータで計算してもよい** |
| `biratPfIsoA` / `biratPfIsoB` | `⟨A,1⟩ ≅ ⟨W₁,k⟩`(`rootLift` で作る) |
| `biratPfIdx` / `biratPfMk` | ★★★★右辺の代表元から左辺の元 |
| `biratPfMk_map` | ★★`IdxBirat` の遷移で不変(well-definedness の半分) |
| `biratPf` | ★★★★内側の余極限から降ろした写像 |

## ★★実務メモ

★**`Frobenioid (pfRootPre P F)` は仮引数で受けること。**
`pfRoot_frobenioid hfi G` を型の中に直接書くと、`whnf` が巨大な構造体を展開しようとして
**200000 heartbeats を超える**(実測 73 秒で timeout)。★`F'` を仮引数で受ける
`Prop55Birat.lean` と同じ流儀にすると 0.2 秒で通る。

★`Category.assoc` の `rw` は `𝒞^pf` の側で「`instances` 透明度で型が合わない」と言って落ちる。
**`Eq.trans` と `congrArg` で繋ぐこと**(本ファイルでも 2 箇所)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★0. `𝒞^birat` の合成は引き戻しデータの取り方に依らない -/

section CompSq

variable {G : Frobenioid P}

/-- ★★★**合成は「どの引き戻しデータで計算してもよい」**。

★★`Proposition 4.4` の合成は `birat_pull_exists` の `.choose` で定義されているが、
**`W` の構造射が mono** なので `birat_lift_unique` により
どの引き戻しデータを使っても同じ元になる。 -/
theorem compBirat_mk_of_sq (Fc : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B)
    (ψ : W.unop.left.obj ⟶ E)
    {Dd : C} (γ : Dd ⟶ Z.unop.left.obj) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (α₀ : Dd ⟶ W.unop.left.obj) (hsq : γ ≫ φ = α₀ ≫ W.unop.hom.hom) :
    compBirat P G Fc (HomBirat.mk Z φ) (HomBirat.mk W ψ)
      = HomBirat.mk (idxBiratMk P G (γ ≫ Z.unop.hom.hom)
          (G.core.coAngularComp _ _ hγc Z.unop.hom.property.1)
          (IsPreStep.comp P hγs Z.unop.hom.property.2)) (α₀ ≫ ψ) := by
  rw [compBirat_mk]
  haveI hb : Mono W.unop.hom.hom := G.core.preStepMono _ W.unop.hom.property.2
  set Z₂ : IdxBirat P G A := idxBiratMk P G (γ ≫ Z.unop.hom.hom)
    (G.core.coAngularComp _ _ hγc Z.unop.hom.property.1)
    (IsPreStep.comp P hγs Z.unop.hom.property.2) with hZ₂
  set V := IsFiltered.max (biratPullIdx Fc Z φ W) Z₂ with hV
  set cc := IsFiltered.leftToMax (biratPullIdx Fc Z φ W) Z₂ with hcc
  set cc' := IsFiltered.rightToMax (biratPullIdx Fc Z φ W) Z₂ with hcc'
  refine HomBirat.sound V cc cc' ?_
  have hraw := idxBirat_left_ext ((idxBiratHomMk γ hγc hγs rfl) ≫ cc')
    (biratPullHom Fc Z φ W ≫ cc)
  have hc : cc'.unop.left.hom ≫ γ = cc.unop.left.hom ≫ biratPullGamma Fc Z φ W := hraw
  have key := birat_lift_unique φ hb hsq (biratPull_sq Fc Z φ W) hc
  simp only [← Category.assoc] at key ⊢
  exact congrArg (fun t => t ≫ ψ) key.symm

end CompSq

/-! ## ★1. 右辺の代表元から左辺の元へ -/

/-- ★右辺の代表元 `W` の Frobenius 次数 `k`。 -/
abbrev biratPfDeg {A B : C} (W : IdxPf P F A B) : ℕ+ := P.degFr W.hom.hom.1

/-- ★★`⟨A,1⟩ ≅ ⟨W₁,k⟩`(**方程式つき**の `rootLift` で作る)。 -/
noncomputable def biratPfIsoA (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) : (⟨A, 1⟩ : PfRootObj P F) ≅ ⟨W.right.obj.1, biratPfDeg W⟩ :=
  @asIso _ _ _ _ (rootLift (F := F) hfi W.hom.hom.1 (biratPfDeg W) rfl)
    (rootLift_isIso hfi W.hom.hom.1 W.hom.property.1 _ rfl)

@[simp] theorem biratPfIsoA_hom (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) :
    (biratPfIsoA hfi W).hom = rootLift (F := F) hfi W.hom.hom.1 (biratPfDeg W) rfl := rfl

/-- ★★`⟨B,1⟩ ≅ ⟨W₂,k⟩`(次数が揃うのは `IdxPf` の条件そのもの)。 -/
noncomputable def biratPfIsoB (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) : (⟨B, 1⟩ : PfRootObj P F) ≅ ⟨W.right.obj.2, biratPfDeg W⟩ :=
  @asIso _ _ _ _
    (rootLift (F := F) hfi W.hom.hom.2 (biratPfDeg W) W.hom.property.2.2.symm)
    (rootLift_isIso hfi W.hom.hom.2 W.hom.property.2.1 _ W.hom.property.2.2.symm)

@[simp] theorem biratPfIsoB_hom (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    (W : IdxPf P F A B) :
    (biratPfIsoB hfi W).hom
      = rootLift (F := F) hfi W.hom.hom.2 (biratPfDeg W) W.hom.property.2.2.symm := rfl

/-- ★★**左辺の添字対象** —— `⟨X,k⟩ → ⟨A,1⟩` の co-angular pre-step。 -/
noncomputable def biratPfIdx (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) (Z : IdxBirat P G W.right.obj.1) :
    IdxBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) :=
  idxBiratMk (pfRootPre P F) Gpf
    (rootMap (F := F) hfi Z.unop.hom.hom (biratPfDeg W) ≫ (biratPfIsoA hfi W).inv)
    (pfRoot_isCoAngular hfi _)
    (IsPreStep.comp (pfRootPre P F)
      (rootMap_preStep hfi _ _ Z.unop.hom.property.2)
      (isPreStep_of_isIso (pfRootPre P F) _))

/-- ★★★★**右辺の代表元 `(W, Z, ψ)` から左辺 `(𝒞^pf)^birat` の元を作る**。 -/
noncomputable def biratPfMk (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) (Z : IdxBirat P G W.right.obj.1)
    (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ :=
  HomBirat.mk (biratPfIdx hfi Gpf W Z)
    (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)

/-- ★★**`IdxBirat` の遷移で不変** —— well-definedness の半分。

★`rootMap` が関手的であること(`rootMap_comp`)だけで出る。 -/
theorem biratPfMk_map (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) {Z Z' : IdxBirat P G W.right.obj.1} (u : Z ⟶ Z')
    (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    biratPfMk hfi Gpf W Z' (u.unop.left.hom ≫ ψ) = biratPfMk hfi Gpf W Z ψ := by
  have htri : u.unop.left.hom ≫ Z.unop.hom.hom = Z'.unop.hom.hom :=
    congrArg (fun t : Z'.unop.left ⟶ (coaPreObj P G W.right.obj.1) => t.hom) (Over.w u.unop)
  have hw : rootMap (F := F) hfi u.unop.left.hom (biratPfDeg W)
        ≫ (biratPfIdx hfi Gpf W Z).unop.hom.hom
      = (biratPfIdx hfi Gpf W Z').unop.hom.hom :=
    (Category.assoc _ _ _).symm.trans
      (congrArg (fun t : (⟨Z'.unop.left.obj, biratPfDeg W⟩ : PfRootObj P F)
          ⟶ ⟨W.right.obj.1, biratPfDeg W⟩ => t ≫ (biratPfIsoA hfi W).inv)
        ((rootMap_comp hfi _ _ _).symm.trans
          (congrArg (fun t => rootMap (F := F) hfi t (biratPfDeg W)) htri)))
  have hmap := HomBirat.mk_map (P := pfRootPre P F) (G := Gpf)
    (idxBiratHomMk (Z := biratPfIdx hfi Gpf W Z) (W := biratPfIdx hfi Gpf W Z')
      (rootMap (F := F) hfi u.unop.left.hom (biratPfDeg W)) (pfRoot_isCoAngular hfi _)
      (rootMap_preStep hfi _ _ u.unop.left.property.2) hw)
    (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)
  refine Eq.trans ?_ hmap
  refine congrArg (HomBirat.mk (biratPfIdx hfi Gpf W Z')) ?_
  exact (congrArg (fun t : (⟨Z'.unop.left.obj, biratPfDeg W⟩ : PfRootObj P F)
      ⟶ ⟨W.right.obj.2, biratPfDeg W⟩ => t ≫ (biratPfIsoB hfi W).inv)
    (rootMap_comp (F := F) hfi u.unop.left.hom ψ (biratPfDeg W))).trans
    (Category.assoc _ _ _)

/-! ## ★2. 右辺の内側の余極限からの写像 -/

/-- ★★右辺の内側の余極限 `Hom^birat(W₁,W₂)` からの余錐。 -/
noncomputable def biratPfCocone (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C} (W : IdxPf P F A B) :
    CategoryTheory.Limits.Cocone (homFunctorBirat P G W.right.obj.1 W.right.obj.2) :=
  CategoryTheory.Limits.Cocone.mk
    (HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩)
    { app := fun Z => TypeCat.ofHom fun ψ => biratPfMk hfi Gpf W Z ψ.down
      naturality := fun Z Z' u => by
        ext ψ
        exact biratPfMk_map hfi Gpf W u ψ.down }

/-- ★★★★**`Hom^birat(W₁,W₂) → Hom^birat_{𝒞^pf}(⟨A,1⟩,⟨B,1⟩)`**。 -/
noncomputable def biratPf (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C} (W : IdxPf P F A B) :
    HomBirat P G W.right.obj.1 W.right.obj.2
      → HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ :=
  fun z => CategoryTheory.Limits.colimit.desc _ (biratPfCocone hfi Gpf W) z

@[simp] theorem biratPf_mk (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F)) {A B : C} (W : IdxPf P F A B)
    (Z : IdxBirat P G W.right.obj.1) (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    biratPf hfi Gpf W (HomBirat.mk Z ψ) = biratPfMk hfi Gpf W Z ψ := by
  show CategoryTheory.Limits.colimit.desc _ (biratPfCocone hfi Gpf W)
    (CategoryTheory.Limits.colimit.ι
      (homFunctorBirat P G W.right.obj.1 W.right.obj.2) Z (ULift.up ψ)) = _
  rw [← types_comp_apply (CategoryTheory.Limits.colimit.ι
      (homFunctorBirat P G W.right.obj.1 W.right.obj.2) Z)
    (CategoryTheory.Limits.colimit.desc _ (biratPfCocone hfi Gpf W)),
    CategoryTheory.Limits.colimit.ι_desc]
  rfl

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (ii)` の左辺の側(右辺の代表元からの写像)。 -/
def biratPfMk.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 右辺の代表元から (𝒞^pf)^birat の元を作る",
    sectionId := "frdi-prop-5-5" }

/-! ## ★3. `IdxPf` の遷移で不変 —— well-definedness のもう半分

★★ここが「選択された同型」の問題が効くところである。
`rootLift` を**方程式で特徴づけた**(`Prop55PfKappa.lean`)おかげで、
`W` を `W′` へ大きくしたときの比較 `iA_{W′} = iA_W ≫ rootStep a` が**言える**
(`biratPfIsoA_step`)。★あとは `rootMap` の自然性の四角形を当てるだけ。 -/

theorem rootLift_congr (hfi : IsOfFrobeniusIsotropicType P) {A A₁ : C} {α α' : A ⟶ A₁}
    (h : α = α') (k : ℕ+) (hk : P.degFr α = k) (hk' : P.degFr α' = k) :
    rootLift (F := F) hfi α k hk = rootLift (F := F) hfi α' k hk' := by
  subst h; rfl

/-- ★遷移射の次数の関係 `k′ = k · d`。 -/
theorem biratPf_hk' {A B : C} {W W' : IdxPf P F A B} (u : W ⟶ W') :
    biratPfDeg W' = biratPfDeg W * P.degFr u.right.hom.1 := by
  have htriA : W.hom.hom.1 ≫ u.right.hom.1 = W'.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A B ⟶ W'.right => t.hom.1) (Under.w u)
  show P.degFr W'.hom.hom.1 = P.degFr W.hom.hom.1 * P.degFr u.right.hom.1
  rw [← htriA, P.degFr_comp, mul_comm]

/-- ★★★★**`iA_{W′} = iA_W ≫ rootStep a`** —— ここが整合性の要。 -/
theorem biratPfIsoA_step (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    {W W' : IdxPf P F A B} (u : W ⟶ W') :
    (biratPfIsoA hfi W').hom
      = (biratPfIsoA hfi W).hom ≫ rootStep (F := F) hfi u.right.hom.1 u.right.property.1
          (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W') rfl (biratPf_hk' u) := by
  have htriA : W.hom.hom.1 ≫ u.right.hom.1 = W'.hom.hom.1 :=
    congrArg (fun t : biFrObj P F A B ⟶ W'.right => t.hom.1) (Under.w u)
  have hka : P.degFr (W.hom.hom.1 ≫ u.right.hom.1) = biratPfDeg W' := by rw [htriA]
  show rootLift (F := F) hfi W'.hom.hom.1 (biratPfDeg W') rfl = _
  rw [rootLift_congr hfi htriA.symm (biratPfDeg W') rfl hka]
  exact rootLift_comp hfi W.hom.hom.1 u.right.hom.1 u.right.property.1
    (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W') rfl rfl (biratPf_hk' u) hka

/-- ★★★★**`iB_{W′} = iB_W ≫ rootStep b`**。 -/
theorem biratPfIsoB_step (hfi : IsOfFrobeniusIsotropicType P) {A B : C}
    {W W' : IdxPf P F A B} (u : W ⟶ W') :
    (biratPfIsoB hfi W').hom
      = (biratPfIsoB hfi W).hom ≫ rootStep (F := F) hfi u.right.hom.2 u.right.property.2.1
          (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W')
          u.right.property.2.2.symm (biratPf_hk' u) := by
  have htriB : W.hom.hom.2 ≫ u.right.hom.2 = W'.hom.hom.2 :=
    congrArg (fun t : biFrObj P F A B ⟶ W'.right => t.hom.2) (Under.w u)
  have hkb : P.degFr (W.hom.hom.2 ≫ u.right.hom.2) = biratPfDeg W' := by
    rw [htriB]; exact W'.hom.property.2.2.symm
  show rootLift (F := F) hfi W'.hom.hom.2 (biratPfDeg W') W'.hom.property.2.2.symm = _
  rw [rootLift_congr hfi htriB.symm (biratPfDeg W') W'.hom.property.2.2.symm hkb]
  exact rootLift_comp hfi W.hom.hom.2 u.right.hom.2 u.right.property.2.1
    (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W')
    W.hom.property.2.2.symm u.right.property.2.2.symm (biratPf_hk' u) hkb

/-- ★★★★★★**`IdxPf` の遷移で不変**(四角形の形)。

★★`z ≫ a = β ≫ α`(底の側)と `ψ ≫ b = β ≫ ψ′`(有理関数の側)の 2 本の四角形から、
`(𝒞^pf)^birat` の元が一致することを言う。
★中身は `rootMap_rootStep_sq` を 2 回と `HomBirat.mk_map` 1 回。 -/
theorem biratPfMk_step (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} {W W' : IdxPf P F A B} (u : W ⟶ W')
    {X : C} (z : X ⟶ W.right.obj.1) (hzc : IsCoAngular P z) (hzs : IsPreStep P z)
    (ψ : X ⟶ W.right.obj.2)
    {Y : C} (β : X ⟶ Y) (hβ : IsFrobeniusType P β)
    (hdβ : P.degFr β = P.degFr u.right.hom.1)
    (α : Y ⟶ W'.right.obj.1) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (hsqA : z ≫ u.right.hom.1 = β ≫ α)
    (ψ' : Y ⟶ W'.right.obj.2) (hsqB : ψ ≫ u.right.hom.2 = β ≫ ψ') :
    biratPfMk hfi Gpf W' (idxBiratMk P G α hαc hαs) ψ'
      = biratPfMk hfi Gpf W (idxBiratMk P G z hzc hzs) ψ := by
  have hk' := biratPf_hk' (F := F) u
  set Sa := rootStep (F := F) hfi u.right.hom.1 u.right.property.1
    (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W') rfl hk' with hSa
  set Sb := rootStep (F := F) hfi u.right.hom.2 u.right.property.2.1
    (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W')
    u.right.property.2.2.symm hk' with hSb
  set Sβ := rootStep (F := F) hfi β hβ
    (biratPfDeg W) (P.degFr u.right.hom.1) (biratPfDeg W') hdβ hk' with hSβ
  haveI hia : IsIso Sa := rootStep_isIso hfi _ _ _ _ _ _ _
  haveI hib : IsIso Sb := rootStep_isIso hfi _ _ _ _ _ _ _
  haveI hiβ : IsIso Sβ := rootStep_isIso hfi _ _ _ _ _ _ _
  have sqA : rootMap (F := F) hfi z (biratPfDeg W) ≫ Sa
      = Sβ ≫ rootMap (F := F) hfi α (biratPfDeg W') :=
    rootMap_rootStep_sq hfi z u.right.hom.1 u.right.property.1 β hβ α
      _ _ _ rfl hdβ hk' hsqA
  have sqB : rootMap (F := F) hfi ψ (biratPfDeg W) ≫ Sb
      = Sβ ≫ rootMap (F := F) hfi ψ' (biratPfDeg W') :=
    rootMap_rootStep_sq hfi ψ u.right.hom.2 u.right.property.2.1 β hβ ψ'
      _ _ _ u.right.property.2.2.symm hdβ hk' hsqB
  have hIA : biratPfIsoA hfi W' = biratPfIsoA hfi W ≪≫ asIso Sa :=
    Iso.ext (biratPfIsoA_step hfi u)
  have hIB : biratPfIsoB hfi W' = biratPfIsoB hfi W ≪≫ asIso Sb :=
    Iso.ext (biratPfIsoB_step hfi u)
  have keyA : (@inv _ _ _ _ Sβ hiβ) ≫ rootMap (F := F) hfi z (biratPfDeg W)
      = rootMap (F := F) hfi α (biratPfDeg W') ≫ (@inv _ _ _ _ Sa hia) := by
    rw [IsIso.inv_comp_eq, ← Category.assoc, ← sqA, Category.assoc, IsIso.hom_inv_id,
      Category.comp_id]
  have keyB : (@inv _ _ _ _ Sβ hiβ) ≫ rootMap (F := F) hfi ψ (biratPfDeg W)
      = rootMap (F := F) hfi ψ' (biratPfDeg W') ≫ (@inv _ _ _ _ Sb hib) := by
    rw [IsIso.inv_comp_eq, ← Category.assoc, ← sqB, Category.assoc, IsIso.hom_inv_id,
      Category.comp_id]
  have hwA : (@inv _ _ _ _ Sβ hiβ)
        ≫ (biratPfIdx hfi Gpf W (idxBiratMk P G z hzc hzs)).unop.hom.hom
      = (biratPfIdx hfi Gpf W' (idxBiratMk P G α hαc hαs)).unop.hom.hom := by
    show (@inv _ _ _ _ Sβ hiβ)
        ≫ (rootMap (F := F) hfi z (biratPfDeg W) ≫ (biratPfIsoA hfi W).inv)
      = rootMap (F := F) hfi α (biratPfDeg W') ≫ (biratPfIsoA hfi W').inv
    rw [hIA, Iso.trans_inv, asIso_inv, ← Category.assoc, ← Category.assoc, keyA]
  have hmap := HomBirat.mk_map (P := pfRootPre P F) (G := Gpf)
    (idxBiratHomMk (Z := biratPfIdx hfi Gpf W (idxBiratMk P G z hzc hzs))
      (W := biratPfIdx hfi Gpf W' (idxBiratMk P G α hαc hαs))
      (@inv _ _ _ _ Sβ hiβ) (pfRoot_isCoAngular hfi _)
      (isPreStep_of_isIso (pfRootPre P F) _) hwA)
    (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)
  refine Eq.trans ?_ hmap
  refine congrArg (HomBirat.mk (biratPfIdx hfi Gpf W' (idxBiratMk P G α hαc hαs))) ?_
  show rootMap (F := F) hfi ψ' (biratPfDeg W') ≫ (biratPfIsoB hfi W').inv
    = (@inv _ _ _ _ Sβ hiβ)
      ≫ (rootMap (F := F) hfi ψ (biratPfDeg W) ≫ (biratPfIsoB hfi W).inv)
  rw [hIB, Iso.trans_inv, asIso_inv, ← Category.assoc, ← Category.assoc, keyB]

/-- ★★★★★locator —— `Proposition 5.5, (ii)` の `IdxPf` 遷移での不変性。 -/
def biratPfMk_step.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — IdxPf の遷移で不変(well-definedness のもう半分)",
    sectionId := "frdi-prop-5-5" }

/-! ## ★4. 四角形の存在

★★`u : W ⟶ W′` と代表元 `(z, ψ)` から、
`biratPfMk_step` が要求する 2 本の四角形を作る。
★道具は `Proposition 1.10, (ii)`(在庫、`pre-step ≫ Frobenius` を
`Frobenius ≫ pre-step` へ組み替える)と `frobTransport`(在庫)だけである。 -/

/-- ★★★★**四角形の存在** —— `Proposition 1.10, (ii)` と `frobTransport` で作る。 -/
theorem exists_biratPf_step {A B : C} {W W' : IdxPf P F A B} (u : W ⟶ W')
    {X : C} (z : X ⟶ W.right.obj.1) (hzc : IsCoAngular P z) (hzs : IsPreStep P z)
    (ψ : X ⟶ W.right.obj.2) :
    ∃ (Y : C) (β : X ⟶ Y) (_ : IsFrobeniusType P β), P.degFr β = P.degFr u.right.hom.1 ∧
      ∃ (α : Y ⟶ W'.right.obj.1), IsCoAngular P α ∧ IsPreStep P α ∧
        ∃ ψ' : Y ⟶ W'.right.obj.2,
          z ≫ u.right.hom.1 = β ≫ α ∧ ψ ≫ u.right.hom.2 = β ≫ ψ' := by
  obtain ⟨Y, β, α, hβ, hdβ, hαs, hsq⟩ :=
    prop_1_10_ii P F z hzs u.right.hom.1 u.right.property.1
  refine ⟨Y, β, hβ, hdβ, α,
    prop_1_10_i_coAngular_of P F hβ u.right.property.1 hsq.symm hzc, hαs,
    frobTransport (F := F) β hβ u.right.hom.2 u.right.property.2.1
      (hdβ.trans u.right.property.2.2) ψ, hsq.symm, ?_⟩
  exact frobTransport_spec (F := F) β hβ u.right.hom.2 u.right.property.2.1
    (hdβ.trans u.right.property.2.2) ψ

/-- ★★★locator —— `Proposition 5.5, (ii)` の四角形の存在。 -/
def exists_biratPf_step.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — IdxPf の遷移が与える四角形の存在",
    sectionId := "frdi-prop-5-5" }

/-! ## ★5. 余錐であること —— well-definedness の完成

★★★ここで 3 つが繋がる:
* `exists_biratPf_step` —— `𝒞` の側で四角形を作る(`Proposition 1.10, (ii)`)
* `biratPf_sq` —— それが `𝒞^birat` の四角形になる(`compBirat_mk_of_sq′`)
* `biratPfMk_step` —— そこから左辺の元が一致する(`rootMap_rootStep_sq`)

★`frobTransport` の**一意性**(`frobTransport_eq`、在庫)が
「作った四角形」と「余極限の遷移写像」を結ぶ。 -/

section Cocone

variable [IsConnected D] {G : Frobenioid P}

theorem idxBiratMk_congr {A A' : C} {a a' : A' ⟶ A} (h : a = a')
    (hac : IsCoAngular P a) (has : IsPreStep P a)
    (hac' : IsCoAngular P a') (has' : IsPreStep P a') :
    idxBiratMk P G a hac has = idxBiratMk P G a' hac' has' := by subst h; rfl

/-- ★★合成の計算則(添字の構造射を好きな形で受け取る版)。 -/
theorem compBirat_mk_of_sq' (Fc : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B)
    (ψ : W.unop.left.obj ⟶ E)
    {Dd : C} (γ : Dd ⟶ Z.unop.left.obj) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (α₀ : Dd ⟶ W.unop.left.obj) (hsq : γ ≫ φ = α₀ ≫ W.unop.hom.hom)
    {δ : Dd ⟶ A} (hδ : δ = γ ≫ Z.unop.hom.hom)
    (hδc : IsCoAngular P δ) (hδs : IsPreStep P δ) :
    compBirat P G Fc (HomBirat.mk Z φ) (HomBirat.mk W ψ)
      = HomBirat.mk (idxBiratMk P G δ hδc hδs) (α₀ ≫ ψ) := by
  subst hδ
  exact compBirat_mk_of_sq Fc Z φ W ψ γ hγc hγs α₀ hsq

/-- ★★★★**`𝒞` の 2 本の四角形は `𝒞^birat` の四角形になる**。 -/
theorem biratPf_sq {A B : C} {W W' : IdxPf P F A B} (u : W ⟶ W')
    {X : C} (z : X ⟶ W.right.obj.1) (hzc : IsCoAngular P z) (hzs : IsPreStep P z)
    (ψ : X ⟶ W.right.obj.2)
    {Y : C} (β : X ⟶ Y) (α : Y ⟶ W'.right.obj.1)
    (hαc : IsCoAngular P α) (hαs : IsPreStep P α) (ψ' : Y ⟶ W'.right.obj.2)
    (hsqA : z ≫ u.right.hom.1 = β ≫ α) (hsqB : ψ ≫ u.right.hom.2 = β ≫ ψ') :
    compBirat P G F (HomBirat.mk (idxBiratMk P G z hzc hzs) ψ)
        (toHomBirat (P := P) (G := G) u.right.hom.2)
      = compBirat P G F (toHomBirat (P := P) (G := G) u.right.hom.1)
          (HomBirat.mk (idxBiratMk P G α hαc hαs) ψ') := by
  have hone : z = z ≫ (idxBiratOne P G W.right.obj.1).unop.hom.hom :=
    (Category.comp_id z).symm
  rw [compBirat_mk_toHomBirat]
  rw [show (toHomBirat (P := P) (G := G) u.right.hom.1)
      = HomBirat.mk (idxBiratOne P G W.right.obj.1) u.right.hom.1 from rfl]
  rw [compBirat_mk_of_sq' (G := G) F (idxBiratOne P G W.right.obj.1) u.right.hom.1
    (idxBiratMk P G α hαc hαs) ψ' z hzc hzs β hsqA hone hzc hzs]
  exact congrArg (HomBirat.mk (idxBiratMk P G z hzc hzs)) hsqB

/-- ★★★★★★★**`biratPf` は `IdxPf` 上の余錐** —— well-definedness の完成。

★★★これで「右辺の余極限から左辺 `(𝒞^pf)^birat` への写像」が
**代表元の取り方に依らずに定まる**。 -/
theorem biratPf_map (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    {A B : C} {W W' : IdxPf P F A B} (u : W ⟶ W')
    (zz : HomBirat P G W.right.obj.1 W.right.obj.2) :
    biratPf hfi Gpf W'
        (frobTransport (F := F') (toHomBirat (P := P) (G := G) u.right.hom.1)
          ((birat_isFrobeniusType_iff P G _).mpr ⟨u.right.property.1.1.1, u.right.property.1.2⟩)
          (toHomBirat (P := P) (G := G) u.right.hom.2)
          ((birat_isFrobeniusType_iff P G _).mpr
            ⟨u.right.property.2.1.1.1, u.right.property.2.1.2⟩)
          ((biratDeg_toHomBirat (P := P) (G := G) u.right.hom.1).trans
            (u.right.property.2.2.trans
              (biratDeg_toHomBirat (P := P) (G := G) u.right.hom.2).symm)) zz)
      = biratPf hfi Gpf W zz := by
  obtain ⟨Z, x, rfl⟩ :=
    HomColim.exists_rep (homFunctorBirat P G W.right.obj.1 W.right.obj.2) zz
  obtain ⟨ψ⟩ := x
  obtain ⟨Y, β, hβ, hdβ, α, hαc, hαs, ψ', hsqA, hsqB⟩ :=
    exists_biratPf_step u Z.unop.hom.hom Z.unop.hom.property.1 Z.unop.hom.property.2 ψ
  have hsq := biratPf_sq (F := G.core) (G := G) u Z.unop.hom.hom Z.unop.hom.property.1
    Z.unop.hom.property.2 ψ β α hαc hαs ψ' hsqA hsqB
  have htr := frobTransport_eq (F := F') (toHomBirat (P := P) (G := G) u.right.hom.1)
    ((birat_isFrobeniusType_iff P G _).mpr ⟨u.right.property.1.1.1, u.right.property.1.2⟩)
    (toHomBirat (P := P) (G := G) u.right.hom.2)
    ((birat_isFrobeniusType_iff P G _).mpr ⟨u.right.property.2.1.1.1, u.right.property.2.1.2⟩)
    ((biratDeg_toHomBirat (P := P) (G := G) u.right.hom.1).trans
      (u.right.property.2.2.trans
        (biratDeg_toHomBirat (P := P) (G := G) u.right.hom.2).symm))
    (HomBirat.mk Z ψ) (HomBirat.mk (idxBiratMk P G α hαc hαs) ψ') hsq
  have hgoal : HomColim.mk (homFunctorBirat P G W.right.obj.1 W.right.obj.2) Z (ULift.up ψ)
      = HomBirat.mk Z ψ := rfl
  rw [hgoal]
  refine ((congrArg (biratPf hfi Gpf W') htr).trans
    (biratPf_mk hfi Gpf W' (idxBiratMk P G α hαc hαs) ψ')).trans ?_
  refine Eq.trans ?_ (biratPf_mk hfi Gpf W Z ψ).symm
  exact biratPfMk_step hfi Gpf u Z.unop.hom.hom Z.unop.hom.property.1
    Z.unop.hom.property.2 ψ β hβ hdβ α hαc hαs hsqA ψ' hsqB

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の余錐性(well-definedness の完成)。 -/
def biratPf_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 右辺の余極限から (𝒞^pf)^birat への写像(well-defined)",
    sectionId := "frdi-prop-5-5" }

/-! ## ★6. 余極限から降ろす —— 右辺全体からの写像

★★`idxToBirat` が**終尾**(`idxToBirat_final`、`Prop55Birat.lean`)なので、
`colim_{W ∈ IdxPf(𝒞)(A,B)} Hom^birat(W₁,W₂) ≅ Hom^pf_{𝒞^birat}(A,B)` である。
★`biratPf_map`(余錐性)を `colimit.desc` に渡し、その終尾同型で戻せば
**`(𝒞^birat)^pf` の射から `(𝒞^pf)^birat` の射への写像**が得られる。 -/

/-- ★★★★★★左辺への写像を与える余錐(右辺の外側の余極限から)。 -/
noncomputable def biratPfOuterCocone (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C) :
    Limits.Cocone (idxToBirat P F G F' A B ⋙
      homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) :=
  Limits.Cocone.mk (HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩)
    { app := fun W => TypeCat.ofHom fun z => biratPf hfi Gpf W z.down
      naturality := fun W W' u => by
        ext z
        exact biratPf_map hfi Gpf F' u z.down }

/-- ★★★★★★★**`(𝒞^birat)^pf` の射から `(𝒞^pf)^birat` の射への写像**。 -/
noncomputable def biratPfHom (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C) :
    HomPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)
      → HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ :=
  fun z => Limits.colimit.desc _ (biratPfOuterCocone hfi Gpf F' A B)
    ((Functor.Final.colimitIso (idxToBirat P F G F' A B)
      (homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B))).inv z)

/-- ★★**代表元での計算則**。 -/
@[simp] theorem biratPfHom_mk (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A B : C)
    (W : IdxPf P F A B) (z : HomBirat P G W.right.obj.1 W.right.obj.2) :
    biratPfHom hfi Gpf F' A B
        (HomPf.mk ((idxToBirat P F G F' A B).obj W) z) = biratPf hfi Gpf W z := by
  have h2 : Limits.colimit.pre
        (homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B))
        (idxToBirat P F G F' A B)
        (Limits.colimit.ι (idxToBirat P F G F' A B ⋙
          homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) W (ULift.up z))
      = HomPf.mk ((idxToBirat P F G F' A B).obj W) z := by
    rw [← types_comp_apply (Limits.colimit.ι (idxToBirat P F G F' A B ⋙
        homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) W)
      (Limits.colimit.pre (homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B))
        (idxToBirat P F G F' A B)), Limits.colimit.ι_pre]
    rfl
  have hinv : (Functor.Final.colimitIso (idxToBirat P F G F' A B)
        (homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B))).inv
        (HomPf.mk ((idxToBirat P F G F' A B).obj W) z)
      = Limits.colimit.ι (idxToBirat P F G F' A B ⋙
          homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) W (ULift.up z) := by
    rw [← h2]
    exact Iso.hom_inv_id_apply _ _
  show Limits.colimit.desc _ (biratPfOuterCocone hfi Gpf F' A B)
      ((Functor.Final.colimitIso (idxToBirat P F G F' A B)
        (homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B))).inv
        (HomPf.mk ((idxToBirat P F G F' A B).obj W) z)) = _
  rw [hinv, ← types_comp_apply (Limits.colimit.ι (idxToBirat P F G F' A B ⋙
      homFunctorPf (biratPre P G) F' (biratUp P G A) (biratUp P G B)) W)
    (Limits.colimit.desc _ (biratPfOuterCocone hfi Gpf F' A B)), Limits.colimit.ι_desc]
  rfl

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の右辺から左辺への写像。 -/
def biratPfHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — (𝒞^birat)^pf の射から (𝒞^pf)^birat の射への写像",
    sectionId := "frdi-prop-5-5" }

end Cocone


section Surj

variable {G : Frobenioid P}

/-! ## ★全射性 —— 右辺の代表元は必ず左辺から来る

★★★右辺 `(𝒞^birat)^pf` の任意の元は、**添字を粗くすれば**
`⟨A″,n⟩ → ⟨A,1⟩`(co-angular pre-step)と
`ψ₁ : (A″)^{(k)} ⟶ A^{(n,k)}` の対で書ける。
★このとき **`biratPfMk` がその元をちょうど返す** ——
これが `Proposition 5.5, (ii)` の全射性の中身である。

★鍵は `Prop55PfKappa.lean` の**全射性の三角形** `surj_triangle'`:
同じ 1 本の `e_β = kappaLift` が `A` の側(構造射)にも `B` の側(値)にも当たる。
`e_β` は `A″`・`n`・`k` にしか依らないからである。 -/

/-- ★★全射性で使う添字 —— `θ_A = rtExt A n ≫ rtExt A^{(n)} k` と `θ_B` の対。 -/
noncomputable def surjW (A B : C) (n k : ℕ+) : IdxPf P F A B :=
  idxMk (P := P) (F := F) (rtExt P F A n ≫ rtExt P F (rtObj P F A n) k)
    (rtExt P F B n ≫ rtExt P F (rtObj P F B n) k)
    (IsFrobeniusType.comp P F (rtExt_frobType P F A n) (rtExt_frobType P F _ k))
    (IsFrobeniusType.comp P F (rtExt_frobType P F B n) (rtExt_frobType P F _ k))
    (by rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr, rtExt_degFr, rtExt_degFr])

/-- ★`surjW` の次数は `k · n`。 -/
theorem surjW_deg (A B : C) (n k : ℕ+) :
    biratPfDeg (P := P) (F := F) (surjW (P := P) (F := F) A B n k) = k * n :=
  surj_degA (P := P) (F := F) A n k

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**全射性** —— `biratPfMk` は右辺の代表元をそのまま返す。

★★`ε = [ψ₁]`(co-angular pre-step)と値 `[ψ₂]` を与えると、
`surjW A B n k` を添字にした `biratPfMk` が
**`HomBirat.mk ⟨ε⟩ [ψ₂]` に一致する**。

★★証明は 3 段:
1. `e_β := kappaLift`(`κ ≫ [β]` に沿った持ち上げ)は同型 —— `κ ≫ [β]` が Frobenius 型だから。
2. `surj_triangle'` を `A` の側と `B` の側に 1 回ずつ当てる。
   出てくる `e_β` は**同じ 1 本**である。
3. `e_β⁻¹` を `IdxBirat` の射とみなし(`idxBiratHomMk`)、`HomBirat.mk_map` に流し込む。 -/
theorem biratPfMk_surj (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (A B A'' : C) (n k : ℕ+)
    (ψ₁ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F A n) k)
    (hψ₁c : IsCoAngular P ψ₁) (hψ₁s : IsPreStep P ψ₁)
    (ψ₂ : rtObj P F (rtObj P F A'' 1) k ⟶ rtObj P F (rtObj P F B n) k)
    (hεc : IsCoAngular (pfRootPre P F)
      (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁))
    (hεs : IsPreStep (pfRootPre P F)
      (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)) :
    biratPfMk hfi Gpf (surjW (P := P) (F := F) A B n k)
        (idxBiratMk P G ψ₁ hψ₁c hψ₁s) ψ₂
      = HomBirat.mk (idxBiratMk (pfRootPre P F) Gpf
          (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
            HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁) hεc hεs)
        (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨B, 1⟩ from
          HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F B n) k) ψ₂) := by
  have hBdeg : (pfRootPre P F).degFr (pfKappa (F := F) A'' n
      ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k))
      = biratPfDeg (P := P) (F := F) (surjW (P := P) (F := F) A B n k) :=
    (surj_degB (P := P) (F := F) A'' n k).trans
      (surj_degA (P := P) (F := F) A n k).symm
  set eβ := kappaLift (F := F) hfi (biratPfDeg (P := P) (F := F)
      (surjW (P := P) (F := F) A B n k)) (pfKappa (F := F) A'' n
    ≫ toRootHom (F := F) (rtExt P F A'' 1 ≫ rtExt P F (rtObj P F A'' 1) k)) hBdeg with heβ
  haveI hiso : IsIso eβ := kappaLift_isIso hfi _ _
    (IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi) (pfKappa_frobType hfi A'' n)
      (toRootHom_frobType hfi _
        (IsFrobeniusType.comp P F (rtExt_frobType P F A'' 1) (rtExt_frobType P F _ k)))) hBdeg
  have htriA : (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)
        ≫ (biratPfIsoA hfi (surjW (P := P) (F := F) A B n k)).hom
      = eβ ≫ rootMap (F := F) hfi ψ₁ (biratPfDeg (P := P) (F := F)
        (surjW (P := P) (F := F) A B n k)) :=
    surj_triangle' (F := F) hfi A'' A n k _ ψ₁ rfl hBdeg
  have htriB : (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨B, 1⟩ from
        HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F B n) k) ψ₂)
        ≫ (biratPfIsoB hfi (surjW (P := P) (F := F) A B n k)).hom
      = eβ ≫ rootMap (F := F) hfi ψ₂ (biratPfDeg (P := P) (F := F)
        (surjW (P := P) (F := F) A B n k)) :=
    surj_triangle' (F := F) hfi A'' B n k _ ψ₂
      (surjW (P := P) (F := F) A B n k).hom.property.2.2.symm hBdeg
  have hw : (@inv _ _ _ _ eβ hiso)
        ≫ (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨A, 1⟩ from
            HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F A n) k) ψ₁)
      = rootMap (F := F) hfi ψ₁ (biratPfDeg (P := P) (F := F)
          (surjW (P := P) (F := F) A B n k))
        ≫ (biratPfIsoA hfi (surjW (P := P) (F := F) A B n k)).inv :=
    (Iso.eq_comp_inv (biratPfIsoA hfi (surjW (P := P) (F := F) A B n k))).mpr
      ((Category.assoc _ _ _).trans
        ((congrArg (fun t => (@inv _ _ _ _ eβ hiso) ≫ t) htriA).trans
          (IsIso.inv_hom_id_assoc eβ _)))
  have hval : (@inv _ _ _ _ eβ hiso)
        ≫ (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨B, 1⟩ from
            HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F B n) k) ψ₂)
      = rootMap (F := F) hfi ψ₂ (biratPfDeg (P := P) (F := F)
          (surjW (P := P) (F := F) A B n k))
        ≫ (biratPfIsoB hfi (surjW (P := P) (F := F) A B n k)).inv :=
    (Iso.eq_comp_inv (biratPfIsoB hfi (surjW (P := P) (F := F) A B n k))).mpr
      ((Category.assoc _ _ _).trans
        ((congrArg (fun t => (@inv _ _ _ _ eβ hiso) ≫ t) htriB).trans
          (IsIso.inv_hom_id_assoc eβ _)))
  refine Eq.trans ?_ (HomBirat.mk_map (P := pfRootPre P F) (G := Gpf)
    (idxBiratHomMk (Z := idxBiratMk (pfRootPre P F) Gpf _ hεc hεs)
      (W := biratPfIdx hfi Gpf (surjW (P := P) (F := F) A B n k)
        (idxBiratMk P G ψ₁ hψ₁c hψ₁s))
      (@inv _ _ _ _ eβ hiso) (pfRoot_isCoAngular hfi _)
      (isPreStep_of_isIso (pfRootPre P F) _) hw)
    (show HomRoot P F (⟨A'', n⟩ : PfRootObj P F) ⟨B, 1⟩ from
      HomPf.mk (idxPow (F := F) (rtObj P F A'' 1) (rtObj P F B n) k) ψ₂))
  exact congrArg (HomBirat.mk (biratPfIdx hfi Gpf (surjW (P := P) (F := F) A B n k)
    (idxBiratMk P G ψ₁ hψ₁c hψ₁s))) hval.symm

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の全射性。 -/
def biratPfMk_surj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 全射性(右辺の代表元は biratPfMk の像)",
    sectionId := "frdi-prop-5-5" }

end Surj

end ABC3.Found.FrdI
