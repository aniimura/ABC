/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.RatTowerColimit
import Mathlib.AlgebraicGeometry.AffineTransitionLimit
import Mathlib.CategoryTheory.Limits.Preserves.Opposites

/-!
# [GenEll] Remark 1.5.1 —— **`Spec ℚ = lim Spec ℤ[1/n!]` と点の spreading out**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★ここで mathlib と噛み合う

`RatTowerColimit.lean` で `ℚ = colim ℤ[1/n!]`（`ratTowerIsColimit`）を取った。
★`Spec` は `ΓSpec.adjunction : Γ ⊣ Spec` の**右随伴**なので極限を保つ
——`PreservesLimits Scheme.Spec` は**インスタンスで通る**（2026-08-27 実測）。

したがって

    Spec ℚ = lim_n Spec ℤ[1/n!]   （`specRatTowerIsLimit`）

★★これで mathlib の `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`
（`AffineTransitionLimit.lean`、EGA IV §8）が要求する形が揃う。

## ★★★★要求されたインスタンスは全部通った

| 要求 | 通し方 |
|---|---|
| `IsCofiltered ℕᵒᵖ` | ★そのまま `infer_instance` |
| `IsAffineHom (D.map f)` | `Spec.map` はアフィン射 |
| `CompactSpace (D.obj i)` | `Spec R` はアフィンだから |
| `QuasiSeparatedSpace (D.obj i)` | 同上 |
| `t : D ⟶ const (Spec ℤ)` | ★★**`Spec ℤ` は終対象**なので一意——自然性も自動 |

★★★**`Spec ℤ` が終対象であることが効いた**——底への自然変換を手で作らずに済む。

## ★★★★★★到達点

`exists_factor_specRatTower` —— **`Spec ℚ ⟶ X` は `Spec ℤ[1/n!] ⟶ X` を経由する**
（`X` が `Spec ℤ` 上有限表示なら）。★これが**点の spreading out** である。

## ★残っている段

原文 `Remark 1.5.1` の後半には、点ではなく**モデルの同型**の spreading out が要る:

    X_ℚ ≅ X'_ℚ  ⟹  ∃Σ, X ×_ℤ ℤ[Σ⁻¹] ≅ X' ×_ℤ ℤ[Σ⁻¹]

★★そのためには `X_ℚ = lim (X ×_ℤ ℤ[1/n!])` の図式が要る
——`X ×_ℤ (−)`（`Over.mapPullbackAdj` の右随伴）で本ファイルの極限を押し出す。
★★★「同型であることの降下」は**両向きの射を降ろして単射側
（`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`）に流せば出る**
——mathlib の欠落は無い。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★`Spec` で極限に移す -/

/-- ★★**`n ↦ Spec ℤ[1/n!]`** の余フィルター図式。 -/
noncomputable def specRatTowerDiagram : ℕᵒᵖ ⥤ Scheme.{0} := ratTowerDiagram.op ⋙ Scheme.Spec

instance specRatTowerDiagram_isAffine (i : ℕᵒᵖ) : IsAffine (specRatTowerDiagram.obj i) := by
  show IsAffine (Scheme.Spec.obj (Opposite.op (ratTowerDiagram.obj i.unop)))
  infer_instance

instance specRatTowerDiagram_compact (i : ℕᵒᵖ) : CompactSpace (specRatTowerDiagram.obj i) :=
  inferInstance

instance specRatTowerDiagram_qs (i : ℕᵒᵖ) : QuasiSeparatedSpace (specRatTowerDiagram.obj i) :=
  inferInstance

instance specRatTowerDiagram_affineHom {i j : ℕᵒᵖ} (f : i ⟶ j) :
    IsAffineHom (specRatTowerDiagram.map f) := by
  show IsAffineHom (Scheme.Spec.map (ratTowerDiagram.map f.unop).op)
  infer_instance

/-- ★★`Spec ℚ` を頂点とする錐。 -/
noncomputable def specRatTowerCone : Cone specRatTowerDiagram :=
  Scheme.Spec.mapCone ratTowerCocone.op

/-- ★★★★★★★★**`Spec ℚ = lim_n Spec ℤ[1/n!]`**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★`Spec` は `ΓSpec.adjunction` の右随伴なので極限を保つ——
`PreservesLimits Scheme.Spec` は**インスタンスで通る**。
★★あとは `ratTowerIsColimit` を反対圏の極限に読み替えて流すだけである。 -/
noncomputable def specRatTowerIsLimit : IsLimit specRatTowerCone :=
  isLimitOfPreserves Scheme.Spec ratTowerIsColimit.op

/-! ## ★★★★★★底への自然変換 —— `Spec ℤ` は終対象 -/

/-- ★**底 `Spec ℤ` への自然変換**。

★★`Spec ℤ` は終対象（`specZIsTerminal`）なので、射も自然性も**一意に決まる**。 -/
noncomputable def specRatTowerToZ :
    specRatTowerDiagram ⟶
      (Functor.const ℕᵒᵖ).obj (Spec (CommRingCat.of ℤ)) where
  app _ := specZIsTerminal.from _
  naturality _ _ _ := specZIsTerminal.hom_ext _ _

/-! ## ★★★★★★★★点の spreading out -/

/-- ★★★★★★★★**点の spreading out** ——
`Spec ℚ ⟶ X` は `Spec ℤ[1/n!] ⟶ X` を経由する。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★★これが mathlib の `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`
（`AffineTransitionLimit.lean`、EGA IV §8）に `specRatTowerIsLimit` を流したものである。

★★仮定は `X` が `Spec ℤ` 上**有限表示**であること——原文の `X` は `ℤ`-固有・`ℤ`-平坦な
有限型スキームなので満たす（★ただしその接続は本ファイルには入っていない）。

★★★★**モデルの同型**の spreading out（原文が実際に使う形）には、
`X_ℚ = lim (X ×_ℤ ℤ[1/n!])` の図式が要る——`X ×_ℤ (−)` で本補題の極限を押し出す段である。 -/
theorem exists_factor_specRatTower {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFinitePresentation f]
    (a : specRatTowerCone.pt ⟶ X) :
    ∃ (i : ℕᵒᵖ) (g : specRatTowerDiagram.obj i ⟶ X),
      specRatTowerCone.π.app i ≫ g = a ∧ g ≫ f = specRatTowerToZ.app i :=
  Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation
    specRatTowerDiagram specRatTowerToZ f specRatTowerCone specRatTowerIsLimit a
    (by
      apply NatTrans.ext
      funext i
      exact specZIsTerminal.hom_ext _ _)

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` の後半には
モデルの同型の spreading out（`X_ℚ ≅ X'_ℚ ⟹ ℤ[Σ⁻¹] 上で同型`）と、
因子 `D` の降下、`Σ` の外での conductor の一致が残っている。 -/

def specRatTowerIsLimit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(Spec ℚ = lim_n Spec ℤ[1/n!])",
    sectionId := "genell-rem-1-5-1" }

def exists_factor_specRatTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(点の spreading out——モデルの同型の降下は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_factor_specRatTower.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "specRatTowerIsLimit(Spec ℚ = lim_n Spec ℤ[1/n!])"
      (.inProject "ABC3" "ABC3.Found.GenEll.specRatTowerIsLimit") 9,
    .citation "[mathlib]" "Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation(EGA IV §8)"
      (.inMathlib "AlgebraicGeometry.Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation") 9,
    .citation "[ABC3]" "ratTowerIsColimit(ℚ = colim ℤ[1/n!])"
      (.inProject "ABC3" "ABC3.Found.GenEll.ratTowerIsColimit") 9,
    .implicitStep
      ("★★Spec ℤ が終対象(specZIsTerminal)であることが効いた——底への自然変換を" ++
       "手で作らずに済む。射も自然性も一意に決まる") 9,
    .implicitStep
      ("★仮定は X が Spec ℤ 上有限表示であること。原文の X は ℤ-固有・ℤ-平坦な" ++
       "有限型スキームなので満たすが、その接続は本ファイルには入っていない") 9,
    .implicitStep
      ("★★★★残る段: モデルの同型の spreading out には X_ℚ = lim (X ×_ℤ ℤ[1/n!]) の" ++
       "図式が要る(X ×_ℤ (−) は Over.mapPullbackAdj の右随伴なので極限を保つ)。" ++
       "「同型であることの降下」は両向きの射を降ろして単射側" ++
       "(Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType)に流せば出る" ++
       "——mathlib の欠落は無い") 9 ]

end ABC3.Found.GenEll
