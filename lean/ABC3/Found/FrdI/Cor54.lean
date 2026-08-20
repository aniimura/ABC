/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Rlf

/-!
# [FrdI] Corollary 5.4 —— 実化どうしを結ぶ関手 `Ψ^rlf`

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : Crlf

★`Corollary 5.4` は「`Ψ : 𝒞₁ ≅ 𝒞₂` があれば `Ψ^rlf : 𝒞₁^rlf → 𝒞₂^rlf` が
1-一意に定まり、`Proposition 5.3` の縦の矢印と 1-可換な図式をなす」と言う。
原文の証明は「`Corollary 4.10`; `4.11, (iii), (iv)` から直ちに」である。

★★**本ファイルは「`Corollary 4.11, (iii)` が与えるもの」を入力として受け取り、
`Ψ^rlf` を作る**ところまでを実装する。すなわち

* 底の関手 `Ψ_𝒟 : 𝒟₁ ⥤ 𝒟₂`
* 因子の単系の同型 `Φ₁ ≅ Ψ_𝒟^* Φ₂`(`psiPhiOnBase` がこれを与える)
* それが `Φ^birat` を `Φ^birat` へ移すこと(`hbirat`)

から `Ψ^rlf` を組み立てる。★係数 `S` は引数なので、`S = ℝ≥0` が実化、
`S = ℚ≥0` が完全化の側になる。

## ★★まだ実装していない条(記録)

1. `hbirat` を `Corollary 4.11` から**導く**こと(いまは仮定に置いている)。
2. 原文の **1-一意性**(`Ψ^rlf` が同型を除いて一意)。
3. 原文の **rigidity**(図式の合成関手がすべて rigid)。
4. `Proposition 5.3` の縦の矢印との 1-可換性。

★これらは鎖 `rlf` の節点 `cor54` に残っている。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

variable {S : Type} [CommSemiring S]

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`Φ^birat` を移す写像は `S·Φ^birat` も移す。

★生成元 `r • (1 ⊗ y)` の形で確かめれば済む(`scMap` はスカラー倍とも
`toSc` とも可換だから)。 -/
theorem sPhiBiratOn_mapOver (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂) {d : D₁} {e : D₂}
    (θ : Φ₁.val d →+ Φ₂.val e)
    (hb : ∀ y ∈ phiBiratOn G₁ d, gpMap _ θ y ∈ phiBiratOn G₂ e)
    {x : Gp (ScT S (Φ₁.val d))} (hx : x ∈ sPhiBiratOn S G₁ d) :
    gpMap _ (scMap θ) x ∈ sPhiBiratOn S G₂ e := by
  have hle : sPhiBiratOn S G₁ d ≤
      AddSubgroup.comap (gpMap _ (scMap θ)) (sPhiBiratOn S G₂ e) := by
    refine AddSubgroup.closure_le _ |>.mpr ?_
    rintro _ ⟨r, y, hy, rfl⟩
    show gpMap _ (scMap θ) (sSmulGp r (toScGp y)) ∈ sPhiBiratOn S G₂ e
    rw [gpMap_scMap_sSmulGp, gpMap_scMap_toScGp]
    exact mem_sPhiBiratOn_gen G₂ r (hb y hy)
  exact hle hx

/-- ★`𝒟` の射に沿った `Φ` の同型の成分。 -/
noncomputable abbrev phiIsoApp (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (d : D₁) : Φ₁.val d →+ Φ₂.val (ΨB.obj d) :=
  (η.hom.app (Opposite.op d)).hom

/-- ★`η` の自然性を `MonoidOn.map` の言葉で。 -/
theorem phiIsoApp_nat (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    {A B : D₁} (f : B ⟶ A) (x : Φ₁.val A) :
    phiIsoApp ΨB η B (Φ₁.map f x) = Φ₂.map (ΨB.map f) (phiIsoApp ΨB η A x) := by
  have h := η.hom.naturality f.op
  exact congrArg (fun t : (Φ₁.functor).obj (Opposite.op A) ⟶
      (ΨB.op ⋙ Φ₂.functor).obj (Opposite.op B) => (AddCommMonCat.Hom.hom t) x) h

variable (S) in
/-- ★★★★★**`Ψ^rlf` を与える `ModelData` の射**。

★入力は `Corollary 4.11, (iii)` が与えるもの:
底の関手 `Ψ_𝒟`、因子の単系の同型 `η : Φ₁ ≅ Ψ_𝒟^* Φ₂`、
そして `η` が `Φ^birat` を `Φ^birat` へ移すこと。 -/
noncomputable def scModelHomOver (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
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
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    ModelDataHomOver ΨB
      (scModel S G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hfsmD₁)
      (scModel S G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₂) where
  phiHom d := scMap (phiIsoApp ΨB η d)
  phiNat := by
    intro A B f x
    have hc : (scMap (S := S) (phiIsoApp ΨB η B)).comp (scMap (Φ₁.map f))
        = (scMap (Φ₂.map (ΨB.map f))).comp (scMap (phiIsoApp ΨB η A)) := by
      rw [← scMap_comp, ← scMap_comp]
      congr 1
      ext y
      exact phiIsoApp_nat ΨB η f y
    exact congrArg (fun t : ScT S (Φ₁.val A) →+ ScT S (Φ₂.val (ΨB.obj B)) => t x) hc
  bmonHom d :=
    AddMonoidHom.codRestrict
      ((gpMap _ (scMap (phiIsoApp ΨB η d))).comp (sPhiBiratOn S G₁ d).subtype) _
      (fun x => sPhiBiratOn_mapOver G₁ G₂ (phiIsoApp ΨB η d) (hbirat d) x.2)
  bmonNat := by
    intro A B f x
    refine Subtype.ext ?_
    have hc : (scMap (S := S) (phiIsoApp ΨB η B)).comp (scMap (Φ₁.map f))
        = (scMap (Φ₂.map (ΨB.map f))).comp (scMap (phiIsoApp ΨB η A)) := by
      rw [← scMap_comp, ← scMap_comp]
      congr 1
      ext y
      exact phiIsoApp_nat ΨB η f y
    have key : ∀ z : Gp (ScT S (Φ₁.val A)),
        gpMap _ (scMap (phiIsoApp ΨB η B)) (gpMap _ (scMap (Φ₁.map f)) z)
          = gpMap _ (scMap (Φ₂.map (ΨB.map f))) (gpMap _ (scMap (phiIsoApp ΨB η A)) z) := by
      intro z
      rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]
    exact key _
  divCompat _ _ := rfl

variable (S) in
/-- ★★★★★★**`Ψ^rlf`** —— 実化どうしを結ぶ関手(`Corollary 5.4` の主張の第 1 半)。

★`S = ℝ≥0` なら原文の `Ψ^rlf : 𝒞₁^rlf → 𝒞₂^rlf`。 -/
noncomputable def psiSc (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
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
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    ScModelObj S G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hfsmD₁ ⥤
      ScModelObj S G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₂ :=
  (scModelHomOver S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
    hfsmD₁ hfsmD₂ hbirat).functor

/-! ### ★出典の紐付け -/

/-- ★locator —— `Corollary 5.4` の `Ψ^rlf`(★**条つき**:
`Corollary 4.11, (iii)` の出力を仮定として受け取っており、
1-一意性・rigidity・1-可換性は未実装)。 -/
def psiSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ^rlf(実化どうしを結ぶ関手)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
