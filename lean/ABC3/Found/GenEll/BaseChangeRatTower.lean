/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SpecRatTower
import Mathlib.CategoryTheory.Limits.Constructions.Over.Connected
import Mathlib.CategoryTheory.Comma.Over.Pullback

/-!
# [GenEll] Remark 1.5.1 —— **`X_ℚ = lim (X ×_ℤ ℤ[1/n!])` とモデルの spreading out**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★到達点 —— 原文が使う形の spreading out

原文 `Remark 1.5.1` の証明はこう言う（p.9 の本文）——生成ファイバーの同型が
有限個の素数 `Σ` を除いて `ℤ[Σ⁻¹]` へ延びる。★その**射の側**が本ファイルで通った:

    X_ℚ ⟶ X′  は  X ×_ℤ ℤ[1/n!] ⟶ X′  を経由する   （`exists_factor_baseChangeRatTower`）

## ★★★★経路（4 段、すべて随伴で押す）

| 段 | 使ったもの |
|---|---|
| `ℚ = colim ℤ[1/n!]` | ★`ratTowerIsColimit`（自前、`RatTowerColimit.lean`） |
| `Spec ℚ = lim Spec ℤ[1/n!]` | `Spec` は `ΓSpec.adjunction` の右随伴 |
| `Over (Spec ℤ)` でも極限 | ★`Over.forget` が**連結形の極限を創出する**（`ℕᵒᵖ` は余フィルターゆえ連結） |
| `X_ℚ = lim (X ×_ℤ ℤ[1/n!])` | ★★`Over.pullback f` は `Over.mapPullbackAdj` の右随伴 |

★★★**`Spec ℤ` が終対象であることが 2 回効いた** —— 底への自然変換と、
mathlib の補題の最後の仮定がどちらも `specZIsTerminal.hom_ext` 1 行で出る。

## ★★★★★★インスタンスはすべて既存の道具で通った

| 要求 | 通し方 |
|---|---|
| `IsAffineHom (D.map g)` | ★`MorphismProperty.overPullbackMap`（底変換で保たれる） |
| `CompactSpace (D.obj i)` | `pullback.snd` がアフィン射 ⟹ 準コンパクト ⟹ `X` から降りる |
| `QuasiSeparatedSpace (D.obj i)` | 同上（`quasiSeparatedSpace_of_quasiSeparated`） |
| `IsCofiltered ℕᵒᵖ` / `IsConnected ℕᵒᵖ` | ★そのまま `infer_instance` |

## ★残っている段

★**「同型であることの降下」** —— 本ファイルは**射**を降ろす。同型を降ろすには
両向きの射を降ろし、合成が恒等と生成ファイバーで一致することを
mathlib の `Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`（単射側）に流す。
★★**mathlib の欠落は無い。**

★★★そのあと **因子 `D` の降下**と、**`Σ` の外での conductor の一致**を経て、
`LogCondSigma.lean` の `abs_logCond_sub_le_sum_log`（既存）で `Remark 1.5.1` が閉じる。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★`Over (Spec ℤ)` へ持ち上げる -/

/-- ★★図式を `Over (Spec ℤ)` へ持ち上げる（`Spec ℤ` は終対象なので構造射は一意）。 -/
noncomputable def overRatTowerDiagram :
    ℕᵒᵖ ⥤ Over (Spec (CommRingCat.of ℤ)) :=
  specRatTowerDiagram.toOver _ (fun _ => specZIsTerminal.from _)
    (fun _ => specZIsTerminal.hom_ext _ _)

instance overRatTower_obj_hom_affine (i : ℕᵒᵖ) : IsAffineHom (overRatTowerDiagram.obj i).hom := by
  haveI : IsAffine ((overRatTowerDiagram.obj i).left) := specRatTowerDiagram_isAffine i
  exact isAffineHom_of_isAffine _

/-- ★★`Over (Spec ℤ)` へ持ち上げた錐。 -/
noncomputable def overRatTowerCone : Cone overRatTowerDiagram where
  pt := Over.mk (specZIsTerminal.from specRatTowerCone.pt)
  π := { app := fun i => Over.homMk (specRatTowerCone.π.app i) (specZIsTerminal.hom_ext _ _)
         naturality := fun i j f => by
           apply Over.OverMorphism.ext
           exact specRatTowerCone.π.naturality f }

/-- ★★★★★★★**`Over (Spec ℤ)` でも極限である**。

★`Over.forget` は**連結形の極限を創出する**
（`Mathlib/CategoryTheory/Limits/Constructions/Over/Connected.lean`）。
★★`ℕᵒᵖ` は余フィルターゆえ連結なので、`isLimitOfReflects` がそのまま効く。 -/
noncomputable def overRatTowerIsLimit : IsLimit overRatTowerCone :=
  isLimitOfReflects (Over.forget _) specRatTowerIsLimit

/-! ## ★★★★★★★★底変換 `X ×_ℤ (−)` -/

/-- ★★★★**底変換した図式** `n ↦ X ×_ℤ Spec ℤ[1/n!]`。 -/
noncomputable def baseChangeRatTowerDiagram {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) : ℕᵒᵖ ⥤ Scheme.{0} :=
  overRatTowerDiagram ⋙ Over.pullback f ⋙ Over.forget X

/-- ★★★★**その錐** —— 頂点は `X ×_ℤ Spec ℚ`（＝ `X_ℚ`）。 -/
noncomputable def baseChangeRatTowerCone {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) :
    Cone (baseChangeRatTowerDiagram f) :=
  (Over.forget X).mapCone ((Over.pullback f).mapCone overRatTowerCone)

/-- ★★★★★★★★**`X_ℚ = lim_n (X ×_ℤ ℤ[1/n!])`**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★`Over.pullback f` は `Over.mapPullbackAdj f` の**右随伴**なので極限を保つ。
★★`Over.forget X` も（創出するので）保つ。★★★あとは `overRatTowerIsLimit` を流すだけ。 -/
noncomputable def baseChangeRatTowerIsLimit {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) :
    IsLimit (baseChangeRatTowerCone f) := by
  haveI : PreservesLimits (Over.pullback f) := (Over.mapPullbackAdj f).rightAdjoint_preservesLimits
  exact isLimitOfPreserves (Over.forget X)
    (isLimitOfPreserves (Over.pullback f) overRatTowerIsLimit)

/-! ## ★★★★★★mathlib の降下補題が要求するインスタンス -/

instance baseChangeRatTower_affineHom {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    {i j : ℕᵒᵖ} (g : i ⟶ j) : IsAffineHom ((baseChangeRatTowerDiagram f).map g) := by
  show IsAffineHom (((Over.pullback f).map (overRatTowerDiagram.map g)).left)
  exact MorphismProperty.overPullbackMap f (overRatTowerDiagram.map g)
    (specRatTowerDiagram_affineHom g)

instance baseChangeRatTower_compact {X : Scheme.{0}} [CompactSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (i : ℕᵒᵖ) :
    CompactSpace ((baseChangeRatTowerDiagram f).obj i) := by
  have hA : IsAffineHom (pullback.snd (overRatTowerDiagram.obj i).hom f) :=
    MorphismProperty.pullback_snd _ _ (overRatTower_obj_hom_affine i)
  have hQ : QuasiCompact (pullback.snd (overRatTowerDiagram.obj i).hom f) :=
    @instQuasiCompactOfIsAffineHom _ _ _ hA
  show CompactSpace ((Over.pullback f).obj (overRatTowerDiagram.obj i)).left
  exact @QuasiCompact.compactSpace_of_compactSpace _ _ _ hQ _

instance baseChangeRatTower_qs {X : Scheme.{0}} [QuasiSeparatedSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (i : ℕᵒᵖ) :
    QuasiSeparatedSpace ((baseChangeRatTowerDiagram f).obj i) := by
  have hA : IsAffineHom (pullback.snd (overRatTowerDiagram.obj i).hom f) :=
    MorphismProperty.pullback_snd _ _ (overRatTower_obj_hom_affine i)
  have hQ : QuasiSeparated (pullback.snd (overRatTowerDiagram.obj i).hom f) := by
    haveI := hA
    infer_instance
  show QuasiSeparatedSpace ((Over.pullback f).obj (overRatTowerDiagram.obj i)).left
  exact @quasiSeparatedSpace_of_quasiSeparated _ _ _ _ hQ

/-- ★底 `Spec ℤ` への自然変換（終対象なので一意）。 -/
noncomputable def baseChangeRatTowerToZ {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) :
    baseChangeRatTowerDiagram f ⟶
      (Functor.const ℕᵒᵖ).obj (Spec (CommRingCat.of ℤ)) where
  app _ := specZIsTerminal.from _
  naturality _ _ _ := specZIsTerminal.hom_ext _ _

/-! ## ★★★★★★★★★モデルの spreading out -/

/-- ★★★★★★★★★**モデルの spreading out** ——
`X_ℚ ⟶ X′` は `X ×_ℤ ℤ[1/n!] ⟶ X′` を経由する。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★★**これが原文 `Remark 1.5.1` の証明が使う spreading out の射の側である。**
mathlib の `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation`
（`AffineTransitionLimit.lean`、EGA IV §8）に `baseChangeRatTowerIsLimit` を流した。

★★仮定: `X` は準コンパクト・準分離（原文の `ℤ`-固有な `X` は満たす）、
`X′` は `Spec ℤ` 上**有限表示**（原文の `ℤ`-固有・`ℤ`-平坦な有限型スキームは満たす）。
★その接続は本ファイルには入っていない。

★★★★**同型の降下**は本定理の両向きを取り、合成が恒等と生成ファイバーで一致することを
`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`（単射側）に流せば出る。 -/
theorem exists_factor_baseChangeRatTower {X X' : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    [LocallyOfFinitePresentation f']
    (a : (baseChangeRatTowerCone f).pt ⟶ X') :
    ∃ (i : ℕᵒᵖ) (g : (baseChangeRatTowerDiagram f).obj i ⟶ X'),
      (baseChangeRatTowerCone f).π.app i ≫ g = a ∧ g ≫ f' = (baseChangeRatTowerToZ f).app i :=
  Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation
    (baseChangeRatTowerDiagram f) (baseChangeRatTowerToZ f) f'
    (baseChangeRatTowerCone f) (baseChangeRatTowerIsLimit f) a
    (by
      apply NatTrans.ext
      funext i
      exact specZIsTerminal.hom_ext _ _)

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` にはまだ
(1) 同型であることの降下、(2) 因子 `D` の降下、(3) `Σ` の外での conductor の一致、
が残っている。★ただし (1) は mathlib の単射側から出るので**欠落は無い**。 -/

def baseChangeRatTowerIsLimit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(X_ℚ = lim_n (X ×_ℤ ℤ[1/n!]))",
    sectionId := "genell-rem-1-5-1" }

def exists_factor_baseChangeRatTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(モデルの spreading out——射の側。同型の降下は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_factor_baseChangeRatTower.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "baseChangeRatTowerIsLimit(X_ℚ = lim_n (X ×_ℤ ℤ[1/n!]))"
      (.inProject "ABC3" "ABC3.Found.GenEll.baseChangeRatTowerIsLimit") 9,
    .citation "[mathlib]" "Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation(EGA IV §8)"
      (.inMathlib "AlgebraicGeometry.Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation") 9,
    .citation "[mathlib]" "Over.mapPullbackAdj(底変換は右随伴——極限を保つ)"
      (.inMathlib "CategoryTheory.Over.mapPullbackAdj") 9,
    .citation "[mathlib]" "MorphismProperty.overPullbackMap(アフィン射は底変換で保たれる)"
      (.inMathlib "CategoryTheory.MorphismProperty.overPullbackMap") 9,
    .implicitStep
      ("★★★経路は 4 段、すべて随伴で押す: ℚ = colim ℤ[1/n!](自前) → " ++
       "Spec ℚ = lim Spec ℤ[1/n!](Spec は ΓSpec.adjunction の右随伴) → " ++
       "Over (Spec ℤ) でも極限(Over.forget が連結形の極限を創出、ℕᵒᵖ は連結) → " ++
       "X_ℚ = lim (X ×_ℤ ℤ[1/n!])(Over.pullback は Over.mapPullbackAdj の右随伴)") 9,
    .implicitStep
      ("★★Spec ℤ が終対象であることが 2 回効いた——底への自然変換と、" ++
       "mathlib の補題の最後の仮定がどちらも specZIsTerminal.hom_ext 1 行で出る") 9,
    .implicitStep
      ("★仮定: X は準コンパクト・準分離(原文の ℤ-固有な X は満たす)、" ++
       "X′ は Spec ℤ 上有限表示(原文の ℤ-固有・ℤ-平坦な有限型スキームは満たす)。" ++
       "その接続は本ファイルには入っていない") 9,
    .implicitStep
      ("★★★★同型の降下は本定理の両向きを取り、合成が恒等と生成ファイバーで" ++
       "一致することを Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType" ++
       "(単射側)に流せば出る——mathlib の欠落は無い") 9 ]

end ABC3.Found.GenEll
