/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55BiratOmega
import ABC3.Found.FrdI.Prop55UntrFun
import ABC3.Found.FrdI.Prop55Std
import ABC3.Found.FrdI.Prop55PfArbFull

/-!
# [FrdI] Proposition 5.5 —— 命題全体の紐付け

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

本ファイルは `Proposition 5.5` の **4 つの主張がすべて実装されたこと**を
1 か所に集めて記録する(新しい数学は無い)。

## ★原文の主張と実装の対応

| 原文 | 実装 |
|---|---|
| **(i)** `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)` | `otriPfMap` / `otri_pfRoot_exists_rep`(`Prop55Pf.lean`) |
| **(ii)** `(𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf` | `unTrPfFunctor` ＋ `unTrPfFunctor_full` / `_faithful`(`Prop55UntrFun.lean`) |
| **(ii)** `(𝒞^pf)^birat ≌ (𝒞^birat)^pf` | ★`thetaFunctor_isEquivalence`(`Prop55BiratOmega.lean`) |
| **(iii)** standard 型は `𝒞^pf` / `𝒞^un-tr` / `𝒞^rlf` へ移る | `pfRoot_standardType` / `unTr_standardType` / `scModel_standardType_of_perfFactorial`(`Prop55Std.lean` / `Prop55RlfRefl.lean`) |
| **(iii)** model 型 | `pfRoot_isOfModelType` / `unTr_isOfModelType` / `scModel_isOfModelType` |
| **(iv)** 単系の同定 `(Φ^pf)^birat = ℚ·Φ^birat` ほか | `phiBiratOn_pf_eq_qPhiBiratOn_full`(`Prop55PfArbFull.lean`) |

## ★★(ii) の birat の側について(本増分で閉じた)

```
Ω : 𝒞^pf ⥤ (𝒞^birat)^pf        (omegaFunctor)
Θ := biratDescFunctor Ω hΩ : (𝒞^pf)^birat ⥤ (𝒞^birat)^pf   (thetaFunctor)
```

* **本質的全射** …… 対象の写像が `⟨A,n⟩ ↦ ⟨A^birat, n⟩` なので**無料**
  (`thetaFunctor_obj_surjective`)。
* **充満忠実** …… 根 1 では橋渡し `Θ.map ∘ biratPfHom = pfCatToRoot.map`
  (`theta_biratPfHom`)と在庫の `biratPfHom_bijective` から。
  根 `k` へは `Θ` が `Σ_k` と可換(`theta_scaleRootBirat`)なことで運び、
  一般の根 `⟨A,n⟩,⟨B,m⟩` へは `pfRoot_exists_iso_root` で**両方を根 `n·m` に揃えて**
  `map_bijective_of_iso` で運ぶ。

## ★逸脱の記録

`Prop55BiratOmega.lean` は `{C : Type u2} [Category.{u2} C]`, `[Category.{u2} D]` を置く
(`homPfMap` が `C₁` と `C₂` に同じ hom 宇宙を要求するため)。
★数学的な仮定ではなく**宇宙パラメータの特殊化**であり、原典の主張には何も足していない。
-/

namespace ABC3.Found.FrdI

/-- ★★★★★★★★locator —— **[FrdI] Proposition 5.5 の全体**。

★(i)(ii)(iii)(iv) の 4 主張がすべて実装された(上の対応表)。
★本増分で閉じたのは (ii) の birat の側
(`thetaFunctor_isEquivalence : (𝒞^pf)^birat ≌ (𝒞^birat)^pf`)。 -/
def prop_5_5.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
