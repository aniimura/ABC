/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54
import ABC3.Found.FrdI.Cor411Phi

/-!
# [FrdI] Corollary 5.4 の 1-一意性 —— `Ψ^rlf` は `(Ψ_𝒟, η)` で決まる

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : C1rlf → C2rlf that fits into a 1-commutative

★★`Cor54Rigid.lean` の ★4 で測ったとおり、`Corollary 4.11, (ii)` で効いた
`projPrecompIsoGen`(「`P.proj` との前合成は関手の同型を反射する」)は
**`Corollary 5.4` の縦の矢印には効かない** —— `cToSc : 𝒞 ⥤ 𝒞^rlf` は
本質的全射ではないからである。

★★★そこで原文どおり**`Corollary 4.11, (iii)` の側から降ろす**:
`Ψ^rlf`(`psiSc`)は **`(Ψ_𝒟, η)` だけで決まる**(`scModelHomOver` の定義)。
したがって `(Ψ_𝒟, η)` の 1-一意性がそのまま `Ψ^rlf` の 1-一意性になる。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `natTrans_ext_of_essSurj` | ★本質的全射な関手の像で一致する自然変換は等しい |
| `phiIso_ext` | ★★`η : Φ₁ ≅ Ψ_𝒟^* Φ₂` は `P₁.proj.op` の像での成分で決まる |
| `psiSc_congr` | ★★★`Ψ^rlf` は `(Ψ_𝒟, η)` の関数 |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 uA vA uB vB uE vE

/-! ## ★1. 本質的全射な関手の像で決まる -/

section EssSurjExt

variable {A : Type uA} [Category.{vA} A] {B : Type uB} [Category.{vB} B]
  {E : Type uE} [Category.{vE} E]

/-- ★★★**本質的全射な関手の像で一致する自然変換は等しい**。

★`Corollary 4.11` の `projPrecompIsoGen` / `projOpPrecompIsoGen` が使っている
「像で決まる」の部分を、**一般の本質的全射な関手**について取り出したもの。 -/
theorem natTrans_ext_of_essSurj (K : A ⥤ B) [K.EssSurj] {G G' : B ⥤ E}
    (α β : G ⟶ G') (h : ∀ a : A, α.app (K.obj a) = β.app (K.obj a)) : α = β := by
  ext b
  set p := K.objPreimage b with hp
  set e := K.objObjPreimageIso b with he
  haveI : IsIso (G.map e.hom) := by
    refine ⟨G.map e.inv, ?_, ?_⟩
    · rw [← G.map_comp, e.hom_inv_id, G.map_id]
    · rw [← G.map_comp, e.inv_hom_id, G.map_id]
  refine (cancel_epi (G.map e.hom)).mp ?_
  rw [α.naturality e.hom, β.naturality e.hom, h p]

/-- ★★同上、**関手の同型**の版。 -/
theorem iso_ext_of_essSurj (K : A ⥤ B) [K.EssSurj] {G G' : B ⥤ E}
    (α β : G ≅ G') (h : ∀ a : A, α.hom.app (K.obj a) = β.hom.app (K.obj a)) : α = β :=
  Iso.ext (natTrans_ext_of_essSurj K α.hom β.hom h)

end EssSurjExt

/-! ## ★2. `η` は `P₁.proj.op` の像での成分で決まる -/

section PhiIso

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {Φ₂ : MonoidOn.{v, u, w} D₂}

/-- ★★★★★**`η : Φ₁ ≅ Ψ_𝒟^* Φ₂` は `P₁.proj.op` の像での成分で決まる**。

★`P₁.proj.op` は本質的全射(`projOp_essSurj`)なので、
`Corollary 4.11, (iii)` の `η` は `𝒞₁` の対象の上の成分だけで一意に定まる。 -/
theorem phiIso_ext (F₁ : FrobenioidCore P₁) (ΨB : D₁ ⥤ D₂)
    (η η' : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (h : ∀ Z : C₁ᵒᵖ, η.hom.app (P₁.proj.op.obj Z) = η'.hom.app (P₁.proj.op.obj Z)) :
    η = η' :=
  haveI := projOp_essSurj P₁ F₁
  iso_ext_of_essSurj P₁.proj.op η η' h

end PhiIso

/-! ## ★3. `Ψ^rlf` は `(Ψ_𝒟, η)` の関数 -/

section PsiScCongr

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  {S : Type} [CommSemiring S]

/-- ★★★★★★**[FrdI] Corollary 5.4 の 1-一意性(前半)** ——
`Ψ^rlf` は `(Ψ_𝒟, η)` の**関数**である。

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : C1rlf → C2rlf that fits into a 1-commutative

★★したがって `Corollary 4.11, (ii)(iii)` が `(Ψ_𝒟, η)` の 1-一意性を与えれば
(`cor_4_11_ii_uniq` と `phiIso_ext`)、`Ψ^rlf` の 1-一意性がそのまま出る。
★★★**`cToSc` が本質的全射でないので `projPrecompIsoGen` の道は通らない**
(`Cor54Rigid.lean` の ★4 の測定)。原文が `Corollary 4.11` を名指しするのは
このためである。 -/
theorem psiSc_congr (ΨB : D₁ ⥤ D₂) (η η' : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (hη : η = η')
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
      gpMap _ (phiIsoApp ΨB η' d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    psiSc S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
        hfsmD₁ hfsmD₂ hbirat
      = psiSc S ΨB η' G₁ hiso₁ hfn₁ hcharInj₁ hint₁ G₂ hiso₂ hfn₂ hcharInj₂ hint₂
        hfsmD₁ hfsmD₂ hbirat' := by
  subst hη
  rfl

end PsiScCongr

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Corollary 5.4` の 1-一意性の前半
(`Ψ^rlf` は `(Ψ_𝒟, η)` の関数、`η` は `𝒞₁` の上の成分で決まる)。 -/
def psiSc_congr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ^rlf は (Ψ_𝒟, η) で決まる(1-一意性の前半)",
    sectionId := "frdi-cor-5-4" }

/-- ★★★★locator —— `Corollary 4.11, (iii)` の `η` の 1-一意性
(`P₁.proj.op` の像での成分で決まる)。 -/
def phiIso_ext.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iii) — η は 𝒞₁ の上の成分で一意に決まる",
    sectionId := "frdi-cor-4-11" }

end ABC3.Found.FrdI
