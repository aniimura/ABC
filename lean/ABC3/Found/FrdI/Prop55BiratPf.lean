/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Equiv
import ABC3.Found.FrdI.Prop55Birat

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

★★鍵は **`⟨A,1⟩ ≅ ⟨A₁,k⟩`**(`pfRootIsoOfFrob`)である ——
`α : A ⟶ A₁` が次数 `k` の Frobenius 型なら、`𝒞^pf` では
**`A₁` の `k` 乗根が `A`** である。★これは
`pfRoot_exists_iso_root`(`⟨A,1⟩ ≅ ⟨A^{(k)}, k⟩`)と
`Definition 1.3, (ii)` の本質的一意性(`frobDegUniq`、`A^{(k)} ≅ A₁`)の合成である。

★これで右辺の代表元 `(W, Z, ψ)` から
「`⟨X,k⟩ → ⟨A,1⟩` の co-angular pre-step ＋ `⟨X,k⟩ → ⟨B,1⟩`」が作れる。

| 定理 | 中身 |
|---|---|
| `pfRootIsoOfFrob` | ★★`⟨A,1⟩ ≅ ⟨A₁,k⟩` |
| `lamHom_coaPre` | `lamHom` は co-angular pre-step を保つ |
| `coaPre_comp_iso` | co-angular pre-step に同型を後合成しても co-angular pre-step |
| `biratPfIdx` / `biratPfMk` | ★★★★右辺の代表元から左辺の元 |
| `biratPfMk_map` | ★★`IdxBirat` の遷移で不変(well-definedness の半分) |

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

/-! ## ★1. `⟨A,1⟩ ≅ ⟨A₁,k⟩` -/

/-- ★★★**`α : A ⟶ A₁` が次数 `k` の Frobenius 型なら `⟨A,1⟩ ≅ ⟨A₁,k⟩`**。

★`pfRoot_exists_iso_root`(`⟨A,1⟩ ≅ ⟨A^{(k)}, k⟩`)と
`Definition 1.3, (ii)` の本質的一意性(`frobDegUniq`、`A^{(k)} ≅ A₁`)の合成。 -/
noncomputable def pfRootIsoOfFrob {A A₁ : C} (α : A ⟶ A₁) (hα : IsFrobeniusType P α)
    (k : ℕ+) (hk : P.degFr α = k) :
    (⟨A, 1⟩ : PfRootObj P F) ≅ ⟨A₁, k⟩ :=
  haveI h₀ := (pfRoot_exists_iso_root (P := P) (F := F) A 1 k k (one_mul k).symm).choose_spec
  haveI hθ := (F.frobDegUniq A (rtObj P F A k) A₁ (rtExt P F A k) α
    (rtExt_frobType P F A k) hα (by rw [rtExt_degFr, hk])).choose_spec.1
  haveI : IsIso (lamHom (F := F) k
      (F.frobDegUniq A (rtObj P F A k) A₁ (rtExt P F A k) α
        (rtExt_frobType P F A k) hα (by rw [rtExt_degFr, hk])).choose) :=
    lamHom_isIso k _ hθ
  asIso ((pfRoot_exists_iso_root (P := P) (F := F) A 1 k k (one_mul k).symm).choose
    ≫ lamHom (F := F) k
      (F.frobDegUniq A (rtObj P F A k) A₁ (rtExt P F A k) α
        (rtExt_frobType P F A k) hα (by rw [rtExt_degFr, hk])).choose)

/-- ★`lamHom` は co-angular pre-step を co-angular pre-step へ移す。 -/
theorem lamHom_coaPre (hfi : IsOfFrobeniusIsotropicType P) (k : ℕ+) {A B : C} (φ : A ⟶ B)
    (h : IsPreStep P φ) :
    IsCoAngular (pfRootPre P F) (lamHom (F := F) k φ)
      ∧ IsPreStep (pfRootPre P F) (lamHom (F := F) k φ) :=
  ⟨pfRoot_isCoAngular hfi _, lamHom_isPreStep k φ h⟩

/-- ★co-angular pre-step に同型を後合成しても co-angular pre-step。 -/
theorem coaPre_comp_iso (hfi : IsOfFrobeniusIsotropicType P) {X Y Z : PfRootObj P F}
    (f : X ⟶ Y) (hf : IsPreStep (pfRootPre P F) f) (g : Y ⟶ Z) (hg : IsIso g) :
    IsCoAngular (pfRootPre P F) (f ≫ g) ∧ IsPreStep (pfRootPre P F) (f ≫ g) :=
  haveI := hg
  ⟨pfRoot_isCoAngular hfi _,
    IsPreStep.comp (pfRootPre P F) hf (isPreStep_of_isIso (pfRootPre P F) g)⟩

/-! ## ★2. 右辺の代表元から左辺の元へ -/

/-- ★右辺の代表元 `W` の Frobenius 次数 `k`。 -/
abbrev biratPfDeg {A B : C} (W : IdxPf P F A B) : ℕ+ := P.degFr W.hom.hom.1

/-- ★`⟨A,1⟩ ≅ ⟨W₁,k⟩`。 -/
noncomputable abbrev biratPfIsoA {A B : C} (W : IdxPf P F A B) :
    (⟨A, 1⟩ : PfRootObj P F) ≅ ⟨W.right.obj.1, biratPfDeg W⟩ :=
  pfRootIsoOfFrob (F := F) W.hom.hom.1 W.hom.property.1 _ rfl

/-- ★`⟨B,1⟩ ≅ ⟨W₂,k⟩`(次数が揃うのは `IdxPf` の条件そのもの)。 -/
noncomputable abbrev biratPfIsoB {A B : C} (W : IdxPf P F A B) :
    (⟨B, 1⟩ : PfRootObj P F) ≅ ⟨W.right.obj.2, biratPfDeg W⟩ :=
  pfRootIsoOfFrob (F := F) W.hom.hom.2 W.hom.property.2.1 _ W.hom.property.2.2.symm

/-- ★★**左辺の添字対象** —— `⟨X,k⟩ → ⟨A,1⟩` の co-angular pre-step。 -/
noncomputable def biratPfIdx (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) (Z : IdxBirat P G W.right.obj.1) :
    IdxBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) :=
  idxBiratMk (pfRootPre P F) Gpf
    (lamHom (F := F) (biratPfDeg W) Z.unop.hom.hom ≫ (biratPfIsoA W).inv)
    (pfRoot_isCoAngular hfi _)
    (IsPreStep.comp (pfRootPre P F)
      (lamHom_isPreStep _ _ Z.unop.hom.property.2)
      (isPreStep_of_isIso (pfRootPre P F) _))

/-- ★★★★**右辺の代表元 `(W, Z, ψ)` から左辺 `(𝒞^pf)^birat` の元を作る**。 -/
noncomputable def biratPfMk (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) (Z : IdxBirat P G W.right.obj.1)
    (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    HomBirat (pfRootPre P F) Gpf (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩ :=
  HomBirat.mk (biratPfIdx hfi Gpf W Z)
    (lamHom (F := F) (biratPfDeg W) ψ ≫ (biratPfIsoB W).inv)

/-- ★★**`IdxBirat` の遷移で不変** —— well-definedness の半分。

★`lamHom` が関手的であること(`lamHom_comp`)だけで出る。 -/
theorem biratPfMk_map (hfi : IsOfFrobeniusIsotropicType P) {G : Frobenioid P}
    (Gpf : Frobenioid (pfRootPre P F))
    {A B : C} (W : IdxPf P F A B) {Z Z' : IdxBirat P G W.right.obj.1} (u : Z ⟶ Z')
    (ψ : Z.unop.left.obj ⟶ W.right.obj.2) :
    biratPfMk hfi Gpf W Z' (u.unop.left.hom ≫ ψ) = biratPfMk hfi Gpf W Z ψ := by
  have htri : u.unop.left.hom ≫ Z.unop.hom.hom = Z'.unop.hom.hom :=
    congrArg (fun t : Z'.unop.left ⟶ (coaPreObj P G W.right.obj.1) => t.hom) (Over.w u.unop)
  have hw : lamHom (F := F) (biratPfDeg W) u.unop.left.hom
        ≫ (biratPfIdx hfi Gpf W Z).unop.hom.hom
      = (biratPfIdx hfi Gpf W Z').unop.hom.hom :=
    (Category.assoc _ _ _).symm.trans
      (congrArg (fun t : (⟨Z'.unop.left.obj, biratPfDeg W⟩ : PfRootObj P F)
          ⟶ ⟨W.right.obj.1, biratPfDeg W⟩ => t ≫ (biratPfIsoA W).inv)
        ((lamHom_comp _ _ _).trans
          (congrArg (lamHom (F := F) (biratPfDeg W)) htri)))
  have hmap := HomBirat.mk_map (P := pfRootPre P F) (G := Gpf)
    (idxBiratHomMk (Z := biratPfIdx hfi Gpf W Z) (W := biratPfIdx hfi Gpf W Z')
      (lamHom (F := F) (biratPfDeg W) u.unop.left.hom) (pfRoot_isCoAngular hfi _)
      (lamHom_isPreStep _ _ u.unop.left.property.2) hw)
    (lamHom (F := F) (biratPfDeg W) ψ ≫ (biratPfIsoB W).inv)
  refine Eq.trans ?_ hmap
  refine congrArg (HomBirat.mk (biratPfIdx hfi Gpf W Z')) ?_
  exact (congrArg (fun t : (⟨Z'.unop.left.obj, biratPfDeg W⟩ : PfRootObj P F)
      ⟶ ⟨W.right.obj.2, biratPfDeg W⟩ => t ≫ (biratPfIsoB W).inv)
    (lamHom_comp (F := F) (biratPfDeg W) u.unop.left.hom ψ).symm).trans
    (Category.assoc _ _ _)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (ii)` の左辺の側(右辺の代表元からの写像)。 -/
def biratPfMk.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 右辺の代表元から (𝒞^pf)^birat の元を作る",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
