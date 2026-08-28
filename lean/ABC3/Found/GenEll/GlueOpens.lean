/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★開被覆に沿って射を貼る —— 段 E1d の道具（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★これは何か —— 段 E1d の貼り合わせ道具

段 E1d は「各チャートの射 `X_{s_i} ⟶ ℙᴺ_R`（`§9-816` の `chartToProj`）を貼って
`ψ : X ⟶ ℙᴺ_R` を作る」である。★その**貼り合わせの機構**が本ファイルである:

> 開集合の族 `U_i` が `X` を覆い、射 `f_i : U_i ⟶ Y` が重なりで一致するなら、
> **`X ⟶ Y` が一意に定まる**（`(U_i).ι ≫ g = f_i`）。

## ★★★機構 —— mathlib の 2 本を繋ぐだけ

| 道具 | 役割 |
|---|---|
| `Scheme.Cover.mkOfCovers` | 開集合の族から `OpenCover` を作る |
| `Scheme.Cover.glueMorphisms` | 被覆に沿って射を貼る（引き戻しでの一致が条件） |
| `isPullback_opens_inf` | ★**開集合 2 つの引き戻しは共通部分である** |

★3 本目が要点である。`glueMorphisms` の条件は `pullback (U_i).ι (U_j).ι` の言葉だが、
`isPullback_opens_inf` でそれを **`U_i ⊓ U_j`** に読み替えられる
——そうすると条件が「制限が一致する」という読める形になる。

## ★残っている段（明示）

★★**`chartToProj i` が重なりで一致すること**は本ファイルに無い。
それが段 E1d の**数学の側**であり、機構は `sectionRatio_agree`（段 D3）と同じ
「比の比は比」である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits

/-! ## ★開集合の族から `OpenCover` を作る -/

/-- ★**開集合の族が `X` を覆えば `OpenCover` になる**。 -/
noncomputable def opensCover {X : Scheme.{0}} {ι : Type} (U : ι → X.Opens)
    (hcov : (⨆ i, U i) = ⊤) : X.OpenCover :=
  Scheme.Cover.mkOfCovers ι (fun i => (U i : Scheme)) (fun i => (U i).ι) (by
    intro x
    have hx : x ∈ (⨆ i, U i) := by rw [hcov]; trivial
    obtain ⟨i, hi⟩ := TopologicalSpace.Opens.mem_iSup.1 hx
    exact ⟨i, ⟨x, hi⟩, rfl⟩)

/-! ## ★★★★★★★開被覆に沿って射を貼る -/

/-- ★★★★★★★**開被覆に沿って射を貼る** —— 段 E1d の機構。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

開集合の族 `U_i` が `X` を覆い、射 `f_i : U_i ⟶ Y` が**重なりで一致する**なら
`X ⟶ Y` が定まる。

★★条件を `U_i ⊓ U_j` の言葉で書けるのが要点である
——`glueMorphisms` の条件は `pullback (U_i).ι (U_j).ι` の言葉だが、
`isPullback_opens_inf`（mathlib）で共通部分に読み替えられる。 -/
noncomputable def glueOpens {X Y : Scheme.{0}} {ι : Type} (U : ι → X.Opens)
    (hcov : (⨆ i, U i) = ⊤) (f : ∀ i, ((U i : Scheme)) ⟶ Y)
    (hagree : ∀ i j, X.homOfLE (inf_le_left : U i ⊓ U j ≤ U i) ≫ f i
        = X.homOfLE (inf_le_right : U i ⊓ U j ≤ U j) ≫ f j) :
    X ⟶ Y := by
  refine Scheme.Cover.glueMorphisms (opensCover U hcov) f ?_
  intro i j
  show pullback.fst (U i).ι (U j).ι ≫ f i = pullback.snd (U i).ι (U j).ι ≫ f j
  have hpb := isPullback_opens_inf (U i) (U j)
  rw [← cancel_epi hpb.isoPullback.hom, ← Category.assoc, ← Category.assoc,
    hpb.isoPullback_hom_fst, hpb.isoPullback_hom_snd]
  exact hagree i j

/-- ★★★**貼った射はチャートの上で元の射に戻る**。 -/
theorem ι_glueOpens {X Y : Scheme.{0}} {ι : Type} (U : ι → X.Opens)
    (hcov : (⨆ i, U i) = ⊤) (f : ∀ i, ((U i : Scheme)) ⟶ Y)
    (hagree : ∀ i j, X.homOfLE (inf_le_left : U i ⊓ U j ≤ U i) ≫ f i
        = X.homOfLE (inf_le_right : U i ⊓ U j ≤ U j) ≫ f j) (i : ι) :
    (U i).ι ≫ glueOpens U hcov f hagree = f i :=
  Scheme.Cover.ι_glueMorphisms (opensCover U hcov) f _ i

/-! ## ★出典の紐付け(`.src`) -/

def opensCover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(開集合の族から OpenCover を作る)",
    sectionId := "genell-prop-1-4" }

def glueOpens.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(開被覆に沿って射を貼る——段 E1d の機構)",
    sectionId := "genell-prop-1-4" }

def ι_glueOpens.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(貼った射はチャートの上で元の射に戻る)",
    sectionId := "genell-prop-1-4" }

def glueOpens.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Cover.glueMorphisms / ι_glueMorphisms"
      (.inMathlib "AlgebraicGeometry.Scheme.Cover.glueMorphisms") 2,
    .citation "[mathlib]" "isPullback_opens_inf(開集合 2 つの引き戻しは共通部分)"
      (.inMathlib "AlgebraicGeometry.isPullback_opens_inf") 2,
    .citation "[mathlib]" "Scheme.Cover.mkOfCovers"
      (.inMathlib "AlgebraicGeometry.Scheme.Cover.mkOfCovers") 2,
    .implicitStep
      ("★★**chartToProj i が重なりで一致すること**は本ファイルに無い。" ++
       "それが段 E1d の**数学の側**であり、機構は sectionRatio_agree(段 D3)と同じ" ++
       "「比の比は比」である") 7 ]

end ABC3.Found.GenEll
