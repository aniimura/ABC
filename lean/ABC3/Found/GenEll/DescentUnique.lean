/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.BaseChangeRatTower

/-!
# [GenEll] Remark 1.5.1 —— **降下の一意性と持ち上げ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★同型の降下に要る 2 つ

`BaseChangeRatTower.lean` で**射**の降下は通った。★同型を降ろすには 2 つ足りない:

| 足りないもの | 本ファイル |
|---|---|
| 生成ファイバーで一致する 2 つの降下は有限段で一致する | ★`descent_unique_baseChangeRatTower` |
| 降下した射を `X′ ×_ℤ ℤ[1/i!]` へ持ち上げる | ★★`liftDescent` |

★★★どちらも既にある道具で出た——前者は mathlib の
`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`（単射側、Stacks 01ZC）、
後者は**引き戻しの普遍性**である。

## ★★★★持ち上げに整合性の仮定が要らなかった

`liftDescent` は `pullback.lift (pullback.fst …) g _` だが、
その `_`（`pullback.fst ≫ (Spec ℤ[1/i!] → Spec ℤ) = g ≫ f′`）は
★**`Spec ℤ` が終対象なので `specZIsTerminal.hom_ext` 1 行で出る**。
★★だから `g` が `Spec ℤ` 上の射であるという仮定すら要らない。

## ★★★★★残っている組み立て

同型 `α : X_ℚ ≅ X′_ℚ` から:

1. `α.hom` を降ろして `g : X ×_ℤ ℤ[1/i!] ⟶ X′` を得る（`exists_factor_baseChangeRatTower`）
2. `α.inv` を降ろして `h : X′ ×_ℤ ℤ[1/j!] ⟶ X` を得る
3. 両方を持ち上げる（`liftDescent`）
4. 共通の段 `k` で `G ≫ H` と `𝟙` を比べる——生成ファイバーで一致するので
   `descent_unique_baseChangeRatTower` が有限段での一致を出す
5. `IsIso` を結論する

★この組み立ては本ファイルには入っていない。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★降下の一意性 -/

/-- ★★★★★★**降下は有限段で一意** ——
生成ファイバーで一致する 2 つの降下は、有限段で一致する。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★mathlib の `Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`
（Stacks 01ZC の単射側）に `baseChangeRatTowerIsLimit` を流したもの。

★★仮定は `LocallyOfFiniteType f′`——**有限表示より弱い**（降下そのものは有限表示が要る）。 -/
theorem descent_unique_baseChangeRatTower {X X' : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    (f' : X' ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    [LocallyOfFiniteType f']
    {i j : ℕᵒᵖ} (a : (baseChangeRatTowerDiagram f).obj i ⟶ X')
    (ha : (baseChangeRatTowerToZ f).app i = a ≫ f')
    (b : (baseChangeRatTowerDiagram f).obj j ⟶ X')
    (hb : (baseChangeRatTowerToZ f).app j = b ≫ f')
    (hab : (baseChangeRatTowerCone f).π.app i ≫ a = (baseChangeRatTowerCone f).π.app j ≫ b) :
    ∃ (k : ℕᵒᵖ) (hik : k ⟶ i) (hjk : k ⟶ j),
      (baseChangeRatTowerDiagram f).map hik ≫ a = (baseChangeRatTowerDiagram f).map hjk ≫ b :=
  Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType
    (baseChangeRatTowerDiagram f) (baseChangeRatTowerToZ f) f'
    (baseChangeRatTowerCone f) (baseChangeRatTowerIsLimit f) a ha b hb hab

/-! ## ★★★★★★★★降下した射の持ち上げ -/

/-- ★★★★**降下した射を `X′ ×_ℤ ℤ[1/i!]` へ持ち上げる**（引き戻しの普遍性）。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★同型を降ろすには、`X ×_ℤ ℤ[1/i!] ⟶ X′` ではなく
`X ×_ℤ ℤ[1/i!] ⟶ X′ ×_ℤ ℤ[1/i!]` が要る——それを作る。

★★★整合性の仮定が要らない: `pullback.fst ≫ (Spec ℤ[1/i!] → Spec ℤ) = g ≫ f′` は
**`Spec ℤ` が終対象なので `specZIsTerminal.hom_ext` 1 行**で出る。 -/
noncomputable def liftDescent {X X' : Scheme.{0}}
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    (f' : X' ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    {i : ℕᵒᵖ} (g : (baseChangeRatTowerDiagram f).obj i ⟶ X') :
    (baseChangeRatTowerDiagram f).obj i ⟶ (baseChangeRatTowerDiagram f').obj i :=
  pullback.lift (pullback.fst (overRatTowerDiagram.obj i).hom f) g
    (specZIsTerminal.hom_ext _ _)

/-- ★持ち上げてから `X′` へ落とすと元の射。 -/
theorem liftDescent_snd {X X' : Scheme.{0}}
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    (f' : X' ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    {i : ℕᵒᵖ} (g : (baseChangeRatTowerDiagram f).obj i ⟶ X') :
    liftDescent f f' g ≫ pullback.snd (overRatTowerDiagram.obj i).hom f' = g :=
  pullback.lift_snd _ _ _

/-- ★★持ち上げは `Spec ℤ[1/i!]` への射影と両立する——**段を保つ**。 -/
theorem liftDescent_fst {X X' : Scheme.{0}}
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    (f' : X' ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    {i : ℕᵒᵖ} (g : (baseChangeRatTowerDiagram f).obj i ⟶ X') :
    liftDescent f f' g ≫ pullback.fst (overRatTowerDiagram.obj i).hom f'
      = pullback.fst (overRatTowerDiagram.obj i).hom f :=
  pullback.lift_fst _ _ _

/-! ## ★★★★図式の射と射影の関係 -/

/-- ★★**段を下げても `X` への射影は変わらない**。 -/
theorem baseChangeMap_snd {X : Scheme.{0}}
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    {i m : ℕᵒᵖ} (h : m ⟶ i) :
    (baseChangeRatTowerDiagram f).map h ≫ pullback.snd (overRatTowerDiagram.obj i).hom f
      = pullback.snd (overRatTowerDiagram.obj m).hom f := by
  simp only [baseChangeRatTowerDiagram, Functor.comp_map, Over.forget_map,
    Over.pullback_map_left]
  exact pullback.lift_snd _ _ _

/-- ★★**段を下げると `Spec ℤ[1/n!]` への射影は塔の射で送られる**。 -/
theorem baseChangeMap_fst {X : Scheme.{0}}
    (f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ)))
    {i m : ℕᵒᵖ} (h : m ⟶ i) :
    (baseChangeRatTowerDiagram f).map h ≫ pullback.fst (overRatTowerDiagram.obj i).hom f
      = pullback.fst (overRatTowerDiagram.obj m).hom f ≫ (overRatTowerDiagram.map h).left := by
  simp only [baseChangeRatTowerDiagram, Functor.comp_map, Over.forget_map,
    Over.pullback_map_left]
  exact pullback.lift_fst _ _ _

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 同型の降下の**組み立て**（上の 5 段）と、
因子 `D` の降下、`Σ` の外での conductor の一致が残っている。 -/

def descent_unique_baseChangeRatTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(降下は有限段で一意)",
    sectionId := "genell-rem-1-5-1" }

def liftDescent.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(降下した射を X′ ×_ℤ ℤ[1/i!] へ持ち上げる)",
    sectionId := "genell-rem-1-5-1" }

def liftDescent.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "baseChangeRatTowerIsLimit(X_ℚ = lim_n (X ×_ℤ ℤ[1/n!]))"
      (.inProject "ABC3" "ABC3.Found.GenEll.baseChangeRatTowerIsLimit") 9,
    .citation "[mathlib]" "Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType(単射側、Stacks 01ZC)"
      (.inMathlib "AlgebraicGeometry.Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType") 9,
    .implicitStep
      ("★★★持ち上げに整合性の仮定が要らなかった。pullback.fst ≫ (Spec ℤ[1/i!] → Spec ℤ)" ++
       " = g ≫ f′ は Spec ℤ が終対象なので specZIsTerminal.hom_ext 1 行で出る。" ++
       "だから g が Spec ℤ 上の射であるという仮定すら要らない") 9,
    .implicitStep
      ("★残っている組み立て: 同型 α : X_ℚ ≅ X′_ℚ から " ++
       "(1) α.hom を降ろす (2) α.inv を降ろす (3) 両方を持ち上げる " ++
       "(4) 共通の段で G ≫ H と 𝟙 を比べる(生成ファイバーで一致するので " ++
       "descent_unique_baseChangeRatTower が有限段での一致を出す) (5) IsIso を結論する") 9,
    .implicitStep
      ("★★そのあと因子 D の降下と Σ の外での conductor の一致を経て、" ++
       "LogCondSigma.lean の abs_logCond_sub_le_sum_log(既存)で Remark 1.5.1 が閉じる") 9 ]

def baseChangeMap_snd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(段を下げても X への射影は変わらない)",
    sectionId := "genell-rem-1-5-1" }

def baseChangeMap_fst.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(段を下げると Spec ℤ[1/n!] への射影は塔の射で送られる)",
    sectionId := "genell-rem-1-5-1" }

end ABC3.Found.GenEll
