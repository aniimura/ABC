import ABC3.Found.GenEll.HeightBaseChange

/-!
# [GenEll] Definition 1.2, (i) —— **有限素点側の仮定を discharge する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★`HeightBaseChange.lean` の仮定 `hfin` を消す

`htArith_baseChange` は 2 つの仮定を受けていた:

- `hfin` —— 引き戻したイデアルが拡大イデアルになる
- `harch` —— アルキメデス側が `baseChangeArc` になる

★★★**`hfin` は仮定でなく定理である**——本ファイルで証明する。

## ★★機構は `ΓSpecIso` の自然性 1 本

`baseRingHom F K (Spec.map φ) = φ` である:

    `ΓF⁻¹ ≫ (Spec.map φ).appTop ≫ ΓK = ΓF⁻¹ ≫ ΓF ≫ φ = φ`

★`Scheme.ΓSpecIso_naturality` と `Iso.inv_hom_id` だけ。
★★§9-13 の規則どおり、**圏の言葉のまま**終わる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]

/-! ## ★★`Spec.map` の底変換環準同型はもとの射である -/

/-- ★★**`baseRingHom F K (Spec.map φ) = φ`**。

★機構は `ΓSpecIso` の自然性 1 本。
★★`Spec` を取って `Γ` を取ると元に戻る、という当たり前のことだが、
`baseRingHom` の定義が `ΓSpecIso` を 2 回通すので**書いておく必要がある**。 -/
theorem baseRingHom_specMap (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K)) :
    baseRingHom F K (Spec.map φ) = φ := by
  show (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv ≫
    (Spec.map φ).appTop ≫ (Scheme.ΓSpecIso (CommRingCat.of (𝓞 K))).hom = φ
  rw [Scheme.ΓSpecIso_naturality, ← Category.assoc, Iso.inv_hom_id, Category.id_comp]

/-! ## ★★★有限素点側の仮定は定理である -/

variable {X : Scheme.{0}}

/-- ★★★**引き戻したイデアルは、自然な射に沿って拡大イデアルになる**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `x_K^* D = (x_F^* D) · 𝓞_K`

★★★**`HeightBaseChange.lean` の仮定 `hfin` が、これで定理になる。**

★機構は `pullbackIdeal_comp`(`PullbackBase.lean`)+ `baseRingHom_specMap`。 -/
theorem pullbackIdeal_specMap (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K)) :
    pullbackIdeal K D (Spec.map φ ≫ xF)
      = (pullbackIdeal F D xF).map φ.hom := by
  rw [pullbackIdeal_comp F K D xF (Spec.map φ), baseRingHom_specMap]

/-- ★★★**自然な包含 `𝓞_F → 𝓞_K` に沿った形**(`htArith_baseChange` がそのまま使える形)。 -/
theorem pullbackIdeal_specMap_algebraMap [Algebra F K] [Algebra (𝓞 F) (𝓞 K)]
    (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X) :
    pullbackIdeal K D (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = (pullbackIdeal F D xF).map (algebraMap (𝓞 F) (𝓞 K)) :=
  pullbackIdeal_specMap F K D xF _

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
アルキメデス側の仮定(埋め込みの両立)と `X(ℚ̄)` の型は `AlgPointClass.lean`(§9-744)で入った。 -/

def pullbackIdeal_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(有限素点側の底変換——仮定でなく定理として)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
