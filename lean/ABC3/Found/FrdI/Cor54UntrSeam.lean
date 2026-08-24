/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54SeamSec
import ABC3.Found.FrdI.Cor54Compat

/-!
# [FrdI] Corollary 5.4 の縦の矢印の継ぎ目 —— `𝒞^un-tr` の段(帳尻の仮定なし)

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★`Cor54SeamSec.lean` の `cSeamIsoOfMap` を **`𝒞^un-tr` に特殊化**したもの。
`Cor54Compat.lean` が「本ファイルの範囲外」と書いていた継ぎ目そのものである。

## ★★入れるもの

| `cSeamIsoOfMap` の引数 | `𝒞^un-tr` での中身 |
|---|---|
| `P₁` / `P₂` | `unTrPre P₁ Fc₁` / `unTrPre P₂ Fc₂` |
| `R₁` / `R₂` | `unTr_ratFnData …`(有理関数の単系は `Φ^birat`) |
| `hiso` | `unTr_isotropic`(`Proposition 3.3`) |
| `hfn` | `unTr_isOfModelType` の第 2 成分 |
| `Fo` | `untrModelHomOver`(`Corollary 4.11, (iii)` の `η` から) |
| `hphi` | ★**`rfl`** —— `untrModelHomOver` の `phiHom` は `phiIsoApp` そのもの |
| `hinj` | ★**無料** —— `Div_B` は `Φ^birat ⊆ Φ^gp` の包含 |

★★★これで `Cor54Compat.lean` の「まだ実装していない条」は**なくなった**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UntrSeam

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  [IsConnected D₁] [IsConnected D₂]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目(`𝒞^un-tr` の段)** ——

```
𝒞₁^un-tr ≌ Model(Φ₁, Φ₁^birat) --Ψ^model--> Model(Φ₂, Φ₂^birat)
   |                                              ‖
 Ψ^un-tr                                          ‖
   v                                              ‖
𝒞₂^un-tr ≌ Model(Φ₂, Φ₂^birat) ═══════════════════╝
```

★★**帳尻の仮定は 1 つも残っていない。**
入力は `Corollary 4.11, (iii)(iv)` と `Corollary 5.7, (i)` のものだけである。 -/
noncomputable def untrSeamIsoSec (ΨB : D₁ ⥤ D₂)
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (Φ₁.val A)) (hfsmD₁ : IsOfFSMType D₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (Φ₂.val A)) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)),
      y ∈ phiBiratOn (unTr_frobenioid P₁ Fc₁ G₁) d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn (unTr_frobenioid P₂ Fc₂ G₂) (ΨB.obj d))
    (Ψu : UnTr P₁ ⥤ UnTr P₂)
    [Ψu.Full] [Ψu.Faithful] [ΨB.Full] [ΨB.Faithful] [ΨB.EssSurj]
    (sq : (unTrPre P₁ Fc₁).proj ⋙ ΨB ≅ Ψu ⋙ (unTrPre P₂ Fc₂).proj)
    (hdivc : ∀ {A B : UnTr P₁} (φ : A ⟶ B),
      phiIsoApp ΨB η (((unTrPre P₁ Fc₁).toElem.obj A).base) ((unTrPre P₁ Fc₁).Div φ)
        = Φ₂.map (sq.hom.app A) ((unTrPre P₂ Fc₂).Div (Ψu.map φ)))
    (hdegc : ∀ {A B : UnTr P₁} (φ : A ⟶ B),
      (unTrPre P₁ Fc₁).degFr φ = (unTrPre P₂ Fc₂).degFr (Ψu.map φ))
    (hPS : ∀ {A B : UnTr P₁} (f : A ⟶ B),
      IsPreStep (unTrPre P₁ Fc₁) f → IsPreStep (unTrPre P₂ Fc₂) (Ψu.map f))
    (S₁ : BaseSection (unTrPre P₁ Fc₁))
    (Fs₁ : ℕ+ →* SectionEnd S₁) (hFs₁ : IsFrobeniusSection S₁ Fs₁)
    (hpb : ∀ {A B : UnTr P₁} (f : A ⟶ B),
      IsPullBack (unTrPre P₁ Fc₁) f → IsPullBack (unTrPre P₂ Fc₂) (Ψu.map f))
    (hft : ∀ {A : UnTr P₁}, S₁.objP A →
      IsFrobeniusTrivial (unTrPre P₂ Fc₂) (Ψu.obj A))
    (Fs₂ : ℕ+ →* SectionEnd (S₁.map (BaseSquare.ofNatIso Ψu ΨB sq) hpb hft))
    (hFs₂ : IsFrobeniusSection
      (S₁.map (BaseSquare.ofNatIso Ψu ΨB sq) hpb hft) Fs₂) :
    (thm_5_2_iv (unTr_ratFnData Fc₁ G₁ hint₁ hfsmD₁) (unTr_isotropic P₁ Fc₁)
          (fun Z => (unTr_isOfModelType Fc₁ G₁).2 Z) S₁ Fs₁ hFs₁).functor
        ⋙ (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat).functor
      ≅ Ψu ⋙ (thm_5_2_iv (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂) (unTr_isotropic P₂ Fc₂)
          (fun Z => (unTr_isOfModelType Fc₂ G₂).2 Z)
          (S₁.map (BaseSquare.ofNatIso Ψu ΨB sq) hpb hft) Fs₂ hFs₂).functor :=
  cSeamIsoOfMap Ψu ΨB η sq hdivc hdegc hPS
    (unTr_ratFnData Fc₁ G₁ hint₁ hfsmD₁) (unTr_isotropic P₁ Fc₁)
    (fun Z => (unTr_isOfModelType Fc₁ G₁).2 Z) S₁ Fs₁ hFs₁ hpb hft
    (unTr_ratFnData Fc₂ G₂ hint₂ hfsmD₂) (unTr_isotropic P₂ Fc₂)
    (fun Z => (unTr_isOfModelType Fc₂ G₂).2 Z) Fs₂ hFs₂
    (untrModelHomOver ΨB η Fc₁ G₁ hint₁ hfsmD₁ Fc₂ G₂ hint₂ hfsmD₂ hbirat)
    (fun _ => rfl) (fun _ => Subtype.val_injective)

end UntrSeam

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目(`𝒞^un-tr` の段)。

★`Cor54Compat.lean` が「本ファイルの範囲外」と書いていた条がこれで閉じた。 -/
def untrSeamIsoSec.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(𝒞^un-tr の段、帳尻の仮定なし)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
