/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.IsoDescent

/-!
# [GenEll] Remark 1.5.1 —— **切断の降下**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★何を作ったか

`IsoDescent.lean` で**モデルの同型**は降りた。★次は**因子 `D` が一緒に降りる**段である。

その土台になるのが**切断の降下**——`X_ℚ` 上で一致する 2 つの切断は、
**有限段 `X ×_ℤ ℤ[1/j!]` で既に一致している**:

    Γ(X_i, U) ∋ s, t     (X_ℚ で一致)  ⟹  ∃j ≤ i,  X_j で一致

★mathlib の `Scheme.exists_app_map_eq_map_of_isLimit`
（`AffineTransitionLimit.lean`、Stacks 01ZC 系）に `baseChangeRatTowerIsLimit` を流したもの。

## ★★★★★★有限個を一度に揃える

因子は**有限個の生成切断**で書けるので、1 本ずつではなく
**有限族を 1 つの段で一度に揃える**必要がある（`exists_bcMap_app_eq_forall`）。

★機構は余フィルター性（`IsCofiltered.inf_objs_exists`）だが、
★★`ℕᵒᵖ` が**前順序**であることが効いて、
「各 `a` ごとに得た道 `k a ≫ g a : S ⟶ i` が全部等しい」ことが
`Subsingleton.elim` で済む——**道を選び直す手間が消える**。

## ★★★★配管

`Γ(X, U)` は `X.presheaf.obj (op U)` であり、`X.presheaf` の型は
`TopCat.Presheaf CommRingCat X.toPresheafedSpace`（`(Opens X)ᵒᵖ ⥤ CommRingCat` の `def`）。
★そのため `ConcreteCategory.comp_apply` の `rw` が

    Note: The target expression is not type-correct under the `instances` transparency level

で落ちる。★★直し方は `rw` をやめて **`congrArg` で合成の外側を直接当てる**こと
（`app_map_comp_eq`）。★★★そのとき図式 `D` を**変数のまま**にしておくと
`D.obj l` が不透明になり、`baseChangeRatTowerDiagram` 特有の展開が起きない。

## ★残っている段

★因子 `D` の降下（本ファイルの切断の降下＋アフィン被覆で組む）、
★★`Σ` の外での conductor の一致。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★段を下げる射で切断の等式を運ぶ -/

/-- ★★★★**段を下げても切断の等式は保たれる**。

★図式 `D` を**変数のまま**にしてあるのが要点である——
`baseChangeRatTowerDiagram` を直接書くと `D.obj l` が展開されて
`instances transparency` のエラーになる（ファイル冒頭の配管を見よ）。 -/
theorem app_map_comp_eq {D : ℕᵒᵖ ⥤ Scheme.{0}} {i j l : ℕᵒᵖ} (g : j ⟶ i) (k : l ⟶ j)
    {U : (D.obj i).Opens} {s t : Γ(D.obj i, U)}
    (h : (D.map g).app U s = (D.map g).app U t) :
    (D.map (k ≫ g)).app U s = (D.map (k ≫ g)).app U t := by
  rw [Functor.map_comp, Scheme.Hom.comp_app]
  exact congrArg (ConcreteCategory.hom (Scheme.Hom.app (D.map k) (D.map g ⁻¹ᵁ U))) h

/-! ## ★★★★★★★★切断の降下 -/

/-- ★★★★★★★★**切断の降下** ——
`X_ℚ` 上で一致する 2 つの切断は、**有限段で既に一致している**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★mathlib の `Scheme.exists_app_map_eq_map_of_isLimit`（Stacks 01ZC 系）に
`baseChangeRatTowerIsLimit` を流したもの。

★★仮定 `hU` は `U` がコンパクトであること——原文の `X` は `ℤ`-固有なので、
その底変換の準コンパクト開集合はこれを満たす。 -/
theorem exists_bcMap_app_eq {X : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    {i : ℕᵒᵖ} {U : ((baseChangeRatTowerDiagram f).obj i).Opens}
    (hU : IsCompact (U : Set ((baseChangeRatTowerDiagram f).obj i)))
    (s t : Γ((baseChangeRatTowerDiagram f).obj i, U))
    (h : ((baseChangeRatTowerCone f).π.app i).app U s
       = ((baseChangeRatTowerCone f).π.app i).app U t) :
    ∃ (j : ℕᵒᵖ) (g : j ⟶ i),
      ((baseChangeRatTowerDiagram f).map g).app U s
        = ((baseChangeRatTowerDiagram f).map g).app U t :=
  exists_app_map_eq_map_of_isLimit _ _ (baseChangeRatTowerIsLimit f) hU s t h

/-- ★★★★★★★★**有限族を 1 つの段で一度に揃える**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★因子は有限個の生成切断で書けるので、1 本ずつでは足りない。

★★機構は余フィルター性（`IsCofiltered.inf_objs_exists`）だが、
`ℕᵒᵖ` が**前順序**であることが効いて、各 `a` ごとの道 `k a ≫ g a : S ⟶ i` が
**全部等しい**（`Subsingleton.elim`）——道を選び直す手間が消える。 -/
theorem exists_bcMap_app_eq_forall {X : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    {ι : Type} [Finite ι]
    {i : ℕᵒᵖ} {U : ((baseChangeRatTowerDiagram f).obj i).Opens}
    (hU : IsCompact (U : Set ((baseChangeRatTowerDiagram f).obj i)))
    (s t : ι → Γ((baseChangeRatTowerDiagram f).obj i, U))
    (h : ∀ a, ((baseChangeRatTowerCone f).π.app i).app U (s a)
            = ((baseChangeRatTowerCone f).π.app i).app U (t a)) :
    ∃ (j : ℕᵒᵖ) (g : j ⟶ i), ∀ a,
      ((baseChangeRatTowerDiagram f).map g).app U (s a)
        = ((baseChangeRatTowerDiagram f).map g).app U (t a) := by
  classical
  rcases isEmpty_or_nonempty ι with hι | hι
  · exact ⟨i, 𝟙 i, fun a => (hι.false a).elim⟩
  choose jj gg hgg using fun a => exists_bcMap_app_eq f hU (s a) (t a) (h a)
  cases nonempty_fintype ι
  obtain ⟨S, hS⟩ := IsCofiltered.inf_objs_exists (Finset.univ.image jj)
  have kk : ∀ a, S ⟶ jj a := fun a => (hS (Finset.mem_image_of_mem _ (Finset.mem_univ a))).some
  refine ⟨S, kk (Classical.arbitrary ι) ≫ gg _, fun a => ?_⟩
  have hg : kk (Classical.arbitrary ι) ≫ gg (Classical.arbitrary ι) = kk a ≫ gg a :=
    Subsingleton.elim _ _
  rw [hg]
  exact app_map_comp_eq (gg a) (kk a) (hgg a)

/-- ★★★★★★★★★**因子の降下が消費する形** ——
生成切断の族が `X_ℚ` で一致するなら、有限段で**生成するイデアルが一致する**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★これが原文の「因子 `D` も一緒に `ℤ[Σ^{-1}]` へ延びる」段の、
**アフィン局所での中身**である。
★★★残るのはアフィン被覆で貼る段（`IdealSheafData` の水準）である。 -/
theorem exists_bcMap_span_eq {X : Scheme.{0}}
    [CompactSpace X] [QuasiSeparatedSpace X]
    (f : X ⟶ Spec (CommRingCat.of ℤ))
    {ι : Type} [Finite ι]
    {i : ℕᵒᵖ} {U : ((baseChangeRatTowerDiagram f).obj i).Opens}
    (hU : IsCompact (U : Set ((baseChangeRatTowerDiagram f).obj i)))
    (s t : ι → Γ((baseChangeRatTowerDiagram f).obj i, U))
    (h : ∀ a, ((baseChangeRatTowerCone f).π.app i).app U (s a)
            = ((baseChangeRatTowerCone f).π.app i).app U (t a)) :
    ∃ (j : ℕᵒᵖ) (g : j ⟶ i),
      Ideal.span (Set.range fun a => ((baseChangeRatTowerDiagram f).map g).app U (s a))
        = Ideal.span (Set.range fun a => ((baseChangeRatTowerDiagram f).map g).app U (t a)) := by
  obtain ⟨j, g, hg⟩ := exists_bcMap_app_eq_forall f hU s t h
  exact ⟨j, g, by simp only [funext hg]⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` には
★因子 `D` の降下（アフィン被覆で貼る段）、★★`Σ` の外での conductor の一致、が残っている。 -/

def exists_bcMap_app_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(切断の降下——X_ℚ で一致する切断は有限段で一致する)",
    sectionId := "genell-rem-1-5-1" }

def exists_bcMap_span_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(因子の降下のアフィン局所での中身——貼る段は含まない)",
    sectionId := "genell-rem-1-5-1" }

def exists_bcMap_span_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.exists_app_map_eq_map_of_isLimit(切断の降下、Stacks 01ZC 系)"
      (.inMathlib "AlgebraicGeometry.exists_app_map_eq_map_of_isLimit") 9,
    .citation "[ABC3]" "baseChangeRatTowerIsLimit(X_ℚ = lim_n (X ×_ℤ ℤ[1/n!]))"
      (.inProject "ABC3" "ABC3.Found.GenEll.baseChangeRatTowerIsLimit") 9,
    .citation "[mathlib]" "IsCofiltered.inf_objs_exists(有限個の対象の共通の下界)"
      (.inMathlib "CategoryTheory.IsCofiltered.inf_objs_exists") 9,
    .citation "[mathlib]" "Scheme.IdealSheafData.comap(イデアル層の引き戻し)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap") 9,
    .implicitStep
      ("★★ℕᵒᵖ が前順序であること(hom が subsingleton)が効く——各 a ごとの道 " ++
       "k a ≫ g a : S ⟶ i が全部等しいので、道を選び直す手間が消える") 9,
    .implicitStep
      ("★★★残る段: アフィン被覆で貼って IdealSheafData の水準にする。" ++
       "★mathlib は comap / comap_comp / support_comap を持っているので、" ++
       "必要なのは「有限個のアフィン開で生成切断を取る」段である") 9,
    .implicitStep
      ("★★★★残る段: Σ の外での conductor の一致。" ++
       "★これが揃えば Skeleton/GenEll/Section1.lean の remark_1_5_1 が受けている" ++
       "仮定 hagree を証明で置き換えられる") 9 ]

end ABC3.Found.GenEll
