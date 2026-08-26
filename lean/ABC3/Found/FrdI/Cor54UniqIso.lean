/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54Seam
import ABC3.Found.FrdI.Cor54Uniq

/-!
# [FrdI] Corollary 5.4 の 1-一意性 —— `Ψ^rlf` は `(Ψ_𝒟, η)` の**同型類**で決まる

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : C1rlf → C2rlf that fits into a 1-commutative

★★`Cor54Uniq.lean` の `psiSc_congr` は `(Ψ_𝒟, η)` が**等しい**場合しか扱えない。
しかし `Corollary 4.11, (ii)` が与える `Ψ_𝒟` の 1-一意性は**同型**であって等式ではない
(`cor_4_11_ii_uniq` の型を見よ)。★したがって「同型で取り替えてよい」が要る。

## ★本ファイルが閉じること

| 定理 | 中身 |
|---|---|
| `psiSc_iso_of_baseIso` | ★★★★★★**`Ψ^rlf` は `Ψ_𝒟` を同型で取り替えても同型** |

★中身は `ModelData.modelDataHomOverIsoOfBase`(`Cor54Seam.lean`)そのもの ——
`u` の帳尻は `Div_B` の単射性(`𝒞^rlf` では `S·Φ^birat ⊆ Φ^rlf,gp` の包含なので無料)で
自動になり、対象の類の帳尻は `η` の両立性を `Gp` へ持ち上げるだけで `u = 0` で足りる。

★★★これで `Corollary 5.4` の 1-一意性は

* `Corollary 4.11, (ii)` …… `Ψ_𝒟` が同型を除いて一意(`cor_4_11_ii_uniq`)
* `Corollary 4.11, (iii)` …… `η` が `𝒞₁` の上の成分で決まる(`phiIso_ext`)
* **本ファイル** …… その同型を `Ψ^rlf` の同型へ運ぶ

の 3 本で閉じる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UniqIso

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  {S : Type} [CommSemiring S]

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**[FrdI] Corollary 5.4 の 1-一意性(後半)** ——
`Ψ^rlf` は **`Ψ_𝒟` を同型で取り替えても同型**である。

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : C1rlf → C2rlf that fits into a 1-commutative

★仮定 `hcompat` は「`η` と `η'` が `θ` を通して両立する」——
`Corollary 4.11, (iii)` の `η` の作り方(`psiPhiOnBase`)から出るものである。 -/
noncomputable def psiSc_iso_of_baseIso
    {ΨB ΨB' : D₁ ⥤ D₂} (θ : ΨB ≅ ΨB')
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (η' : Φ₁.functor ≅ ΨB'.op ⋙ Φ₂.functor)
    (hcompat : ∀ (d : D₁) (x : Φ₁.val d),
      phiIsoApp ΨB η d x = Φ₂.map (θ.hom.app d) (phiIsoApp ΨB' η' d x))
    (G₁ : Frobenioid P₁) (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A)))
    (G₂ : Frobenioid P₂) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A)))
    (hfsmD₁ : IsOfFSMType D₁) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d))
    (hbirat' : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB' η' d) y ∈ phiBiratOn G₂ (ΨB'.obj d)) :
    psiSc S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
        hfsmD₁ hfsmD₂ hbirat
      ≅ psiSc S ΨB' η' G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
        hfsmD₁ hfsmD₂ hbirat' :=
  by
  refine modelDataHomOverIsoOfBase ?_ ?_ ?_ θ
    (scModelHomOver S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
      hfsmD₁ hfsmD₂ hbirat)
    (scModelHomOver S ΨB' η' G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
      hfsmD₁ hfsmD₂ hbirat') ?_
  · exact fun _ => Subtype.val_injective
  · exact fun d (x : ↥(sPhiBiratOn S G₂ d)) => -x
  · exact fun d (x : ↥(sPhiBiratOn S G₂ d)) => neg_add_cancel x
  · intro d x
    have hhom : phiIsoApp ΨB η d
        = (Φ₂.map (θ.hom.app d)).comp (phiIsoApp ΨB' η' d) :=
      AddMonoidHom.ext (hcompat d)
    show scMap (S := S) (phiIsoApp ΨB η d) x
      = scMap (Φ₂.map (θ.hom.app d)) (scMap (phiIsoApp ΨB' η' d) x)
    rw [hhom, scMap_comp]
    rfl

end UniqIso

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Corollary 5.4` の 1-一意性(後半)
(`Ψ^rlf` は `Ψ_𝒟` を**同型**で取り替えても同型)。 -/
def psiSc_iso_of_baseIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ^rlf は Ψ_𝒟 の同型類で決まる(1-一意性の後半)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
