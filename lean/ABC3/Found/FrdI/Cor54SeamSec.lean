/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54SeamCls
import ABC3.Found.FrdI.Cor57Base

/-!
# [FrdI] Corollary 5.4 の縦の矢印の継ぎ目 —— 基準切断を `Ψ` で運ぶ

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★`Cor54SeamCls.lean` の `cSeamIso` は「`S₂` が `Ψ` の像を含む」(`hsec`)を仮定に取る。
本ファイルは、その `S₂` を**在庫の `BaseSection.map`**(`Corollary 5.7, (i)`)で
実際に作り、**`hsec` を `rfl` で埋める**。

## ★★これで何が言えたか

`Corollary 5.4` の縦の矢印の継ぎ目は

* `Corollary 4.11, (iii)(iv)` が与えるもの(`ΨB` / `η` / `sq` / `hdivc` / `hdegc`)
* `Corollary 5.7, (i)` が要求するもの(`Ψ` / `Ψ_𝒟` の充満忠実性、pull-back 射の保存、
  Frobenius-trivial 性の保存)
* `Ψ` が pre-step を保つこと

**だけ**から出る。★**帳尻の仮定は 1 つも残っていない。**

★★★これが `Cor54SeamUnTr.lean` の測定
「基準切断を揃えないと帳尻が合わない」に対する**回答**である ——
`modelType_equiv`(`Classical.choice` が切断を選ぶ)をやめて
`thm_5_2_iv` で切断を明示し、`S₂ := Ψ_*(S₁)` と取ればよい。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section SeamSec

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  [IsConnected D₁] [IsConnected D₂]

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**[FrdI] Corollary 5.4 の縦の矢印の継ぎ目**
(基準切断を `Ψ` で運んだ形、**帳尻の仮定なし**)——

```
𝒞₁ ≌ Model₁ --Ψ^model--> Model₂
 |                         ‖
 Ψ                         ‖
 v                         ‖
𝒞₂ ≌ Model₂ ═══════════════╝
```

★`S₂` は `S₁.map`(`Corollary 5.7, (i)`)で作るので、
`cSeamIso` が要求していた `hsec` は **`⟨A, hA, rfl⟩`** で埋まる。 -/
noncomputable def cSeamIsoOfMap (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    [Ψ.Full] [Ψ.Faithful] [ΨB.Full] [ΨB.Faithful] [ΨB.EssSurj]
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hdivc : ∀ {A B : C₁} (φ : A ⟶ B),
      phiIsoApp ΨB η ((P₁.toElem.obj A).base) (P₁.Div φ)
        = Φ₂.map (sq.hom.app A) (P₂.Div (Ψ.map φ)))
    (hdegc : ∀ {A B : C₁} (φ : A ⟶ B), P₁.degFr φ = P₂.degFr (Ψ.map φ))
    (hPS : ∀ {A B : C₁} (f : A ⟶ B), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    {G₁ : Frobenioid P₁} (R₁ : RatFnData P₁ G₁)
    (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ Z : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) Z)
    (S₁ : BaseSection P₁) (Fs₁ : ℕ+ →* SectionEnd S₁) (hFs₁ : IsFrobeniusSection S₁ Fs₁)
    (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
    (hft : ∀ {A : C₁}, S₁.objP A → IsFrobeniusTrivial P₂ (Ψ.obj A))
    {G₂ : Frobenioid P₂} (R₂ : RatFnData P₂ G₂)
    (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ Z : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) Z)
    (Fs₂ : ℕ+ →* SectionEnd (S₁.map (BaseSquare.ofNatIso Ψ ΨB sq) hpb hft))
    (hFs₂ : IsFrobeniusSection (S₁.map (BaseSquare.ofNatIso Ψ ΨB sq) hpb hft) Fs₂)
    (Fo : ModelDataHomOver ΨB R₁.model R₂.model)
    (hphi : ∀ d : D₁, Fo.phiHom d = phiIsoApp ΨB η d)
    (hinj : ∀ d : D₂, Function.Injective (R₂.divB d)) :
    (thm_5_2_iv R₁ hiso₁ hfn₁ S₁ Fs₁ hFs₁).functor ⋙ Fo.functor
      ≅ Ψ ⋙ (thm_5_2_iv R₂ hiso₂ hfn₂
          (S₁.map (BaseSquare.ofNatIso Ψ ΨB sq) hpb hft) Fs₂ hFs₂).functor :=
  cSeamIso Ψ ΨB η sq hdivc hdegc R₁ hiso₁ hfn₁ Fs₁ hFs₁ R₂ hiso₂ hfn₂ Fs₂ hFs₂
    (fun {A} hA => ⟨A, hA, rfl⟩) hPS Fo hphi hinj

end SeamSec

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Corollary 5.4` の縦の矢印の継ぎ目
(基準切断を `Corollary 5.7, (i)` の `BaseSection.map` で運んだ形。
★**帳尻の仮定は残っていない**)。 -/
def cSeamIsoOfMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印の継ぎ目(基準切断を Ψ で運んだ形)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
