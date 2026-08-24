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

/-- ★★★★★★★**[FrdI] Corollary 5.4(系全体)** の locator。

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : C1rlf → C2rlf that fits into a 1-commutative

★★**原文の 4 つの主張と実装の対応**:

| 原文 | 実装 |
|---|---|
| `Ψ^rlf : 𝒞₁^rlf ⥤ 𝒞₂^rlf` の**存在** | `psiSc`(`Cor54.lean`) |
| 図式の下の行(係数の拡大)との 1-可換性 | `scBaseFunctor_comp_psiSc`(`Cor54Compat.lean`) |
| 縦の矢印 `𝒞^un-tr ⥤ 𝒞^rlf` との 1-可換性 | `untrToSc_comp_psiSc`(同上) |
| 縦の矢印 `𝒞 ⥤ 𝒞^un-tr` との 1-可換性 | `cToUnTr_comp_psiUnTr`(同上) |
| ★縦の矢印の**継ぎ目**(`𝒞^un-tr ≌ Model`) | `untrSeamIsoSec`(本ファイル) |
| **1-一意性**(等式版・同型版) | `psiSc_congr` / `phiIso_ext`(`Cor54Uniq.lean`)、`psiSc_iso_of_baseIso`(`Cor54UniqIso.lean`) |
| 「Moreover」の **rigidity** | `cor54_rigid_both`(`Cor54Rigid.lean`) |
| 「Finally」(完全化・実化との両立) | `scBaseHom_compOver_scModelHomOver`(`Cor54Compat.lean`) |

★★**1-一意性の読み方(明記)**: 本実装の 1-一意性は**原文の論証どおり**、
`Corollary 4.11, (ii)(iii)` の 1-一意性から降ろした形である ——
`Ψ^rlf` は `(Ψ_𝒟, η)` の**同型類**で決まり(`psiSc_iso_of_baseIso`)、
`(Ψ_𝒟, η)` は `Corollary 4.11, (ii)(iii)` で 1-一意(`cor_4_11_ii_uniq` / `phiIso_ext`)。
★**「図式に入る任意の関手が `Ψ^rlf` と同型」というさらに強い読みは本実装には無い** ——
`𝒞^rlf` の対象は実係数を持ち、縦の矢印 `cToSc` は本質的全射でないので
`projPrecompIsoGen` の道は通らない(`Cor54Rigid.lean` の ★4 の測定)。

★★★**継ぎ目についての測定(2026-08-25)**: 継ぎ目は
`modelType_equiv`(`Classical.choice` が基準切断を選ぶ)のままでは**偽**であり、
`thm_5_2_iv` で切断を明示して `S₂ := Ψ_*(S₁)`(`BaseSection.map`、
`Corollary 5.7, (i)`)と揃えて初めて成り立つ。 -/
def cor_5_4.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
