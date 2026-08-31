/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Remark151Sigma
import Mathlib.AlgebraicGeometry.ValuativeCriterion

/-!
# [GenEll] Remark 1.5.1 —— **固有性から点を延ばす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

## ★★★★★★★★原文が「where we recall that X is proper!」と叫ぶ箇所

原文は `x ∈ X(F)` から `Spec 𝓞_F ⟶ X` を作るのに**固有性**を使う。
★`Remark151Sigma.lean` の `remark_1_5_1_bdeq` が点の対応 `ePt` を
データとして受けているのは、この段をまだ持っていないからである。

## ★★★★★★機構は 2 段

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | **付値環**なら延びる（一意に） | ★**本ファイル**（mathlib の `ValuativeCriterion`） |
| 2 | Dedekind 環 `𝓞_F` へ貼り合わせる | ★★まだ |

★段 1 は mathlib の `IsProper.eq_valuativeCriterion`（Stacks 0BX5）に
`Spec ℤ` が終対象であることを合わせるだけで出る——**底での可換性が自動になる**。

★★段 2（各 `(𝓞_F)_v` で延ばして貼る）は Zariski 被覆ではないので、
`Spec 𝓞_F[1/m]` へ延ばしてから残りの有限個の素点を処理する必要がある。
★★★これは mathlib にも無く、ZMT／正規化の水準の作業である
（`ResearchPaper/mathlib-gap.json` の `dedekind-point-extension`）。
-/

namespace ABC3.Found.GenEll

open CategoryTheory AlgebraicGeometry

/-! ## ★★★★固有射は付値判定法を満たす -/

/-- ★★**固有射は付値判定法を満たす**（Stacks 0BX5 の片側）。

★mathlib の `IsProper.eq_valuativeCriterion` は 4 つの積で述べられているので、
成分を取り出しておく。 -/
theorem valuativeCriterion_of_isProper {X Y : Scheme.{0}} (f : X ⟶ Y) [h : IsProper f] :
    ValuativeCriterion f := by
  have : (ValuativeCriterion ⊓ @QuasiCompact ⊓ @QuasiSeparated ⊓ @LocallyOfFiniteType) f := by
    rw [← IsProper.eq_valuativeCriterion]; exact h
  exact this.1.1.1

/-! ## ★★★★★★★★付値環の点へ一意に延びる -/

/-- ★★★★★★★★**`ℤ` 上固有なら、分数体の点は付値環の点へ一意に延びる**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★★底が `Spec ℤ`（終対象）なので、付値可換四角形の**底での可換性が自動**になる
——`specZIsTerminal.hom_ext` 1 行。
★★★存在と一意性が同時に出るのは、`IsProper` が
`ValuativeCriterion`（＝存在 ⊓ 一意性）を含むからである。 -/
theorem exists_unique_extend_of_isProper {X' : Scheme.{0}}
    (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (R : Type) [CommRing R] [IsDomain R] [ValuationRing R]
    (K : Type) [Field K] [Algebra R K] [IsFractionRing R K]
    (xK : Spec (CommRingCat.of K) ⟶ X') :
    ∃! xR : Spec (CommRingCat.of R) ⟶ X',
      Spec.map (CommRingCat.ofHom (algebraMap R K)) ≫ xR = xK := by
  let S : ValuativeCommSq f' :=
    { R := R, K := K, i₁ := xK, i₂ := specZIsTerminal.from _
      commSq := ⟨specZIsTerminal.hom_ext _ _⟩ }
  obtain ⟨u⟩ := valuativeCriterion_of_isProper f' S
  refine ⟨u.default.l, u.default.fac_left, fun y hy => ?_⟩
  have huniq : (⟨y, hy, specZIsTerminal.hom_ext _ _⟩ : S.commSq.LiftStruct) = u.default :=
    u.uniq _
  exact congrArg CategoryTheory.CommSq.LiftStruct.l huniq

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 残っているのは**段 2**——
各 `(𝓞_F)_v` での延長を Dedekind 環 `𝓞_F` の上へ貼り合わせる段である。 -/

def exists_unique_extend_of_isProper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有性から点を延ばす——付値環の場合)",
    sectionId := "genell-def-1-5" }

def exists_unique_extend_of_isProper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsProper.eq_valuativeCriterion(Stacks 0BX5)"
      (.inMathlib "AlgebraicGeometry.IsProper.eq_valuativeCriterion") 8,
    .implicitStep
      ("★底が Spec ℤ(終対象)なので、付値可換四角形の底での可換性が自動になる" ++
       "——specZIsTerminal.hom_ext 1 行") 8,
    .implicitStep
      ("★★★残る段 2: 各 (𝓞_F)_v での延長を Dedekind 環 𝓞_F の上へ貼り合わせる。" ++
       "★Spec (𝓞_F)_v は Zariski 開ではないので素朴には貼れない。" ++
       "★★Spec 𝓞_F[1/m] へ延ばしてから残りの有限個の素点を処理する道になり、" ++
       "ZMT／正規化の水準の作業である(mathlib-gap.json の dedekind-point-extension)") 8 ]

end ABC3.Found.GenEll
