/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop16Proper
import ABC3.Found.GenEll.Prop16
import ABC3.Found.GenEll.HeightMetric

/-!
# [GenEll] Proposition 1.4, (ii)(iii) —— **固有性から但し書きが外れた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

## ★★原文が「思い出せ」と言っている当のもの

原文の (ii) の証明はこう書く ——

> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

★角括弧の中の「`X^arc` はコンパクトだった！」が、この段の**唯一の**根拠である。
`Found/Arakelov/UltraCompact.lean` の `compactSpace_arc` が固有性から直に
それを出すので、本ファイルは `Prop16.lean` の `ArcModel` 版から
**射影モデルの仮定を落とす**。

| 条 | 但し書きつきの版（既存） | ★本ファイル（但し書きなし） |
|---|---|---|
| (ii) | `prop_1_4_ii`（`ArcModel` を受ける） | ★`prop_1_4_ii_of_proper` |
| (iii) | `htArith_sub_abs_le`（`ArcModel` を受ける） | ★`htArith_sub_abs_le_of_proper` |

★★機構は `Prop16Proper.lean` と同一である（固有性 → `X^arc` コンパクト →
連続関数が有界 → 正規化で定数が `F` に依らなくなる）。

## ★★★(iv) は**この機構では外れない**——混同しない

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★(iv)（Northcott）の射影埋め込みは**固有性からではなく、主張の仮定 `L_ℚ ample` から**
出ている。したがって `NorthcottCoord.lean` の `northcott_of_projModel` の
「射影モデルを与えられたものとして」は、**原文の仮定そのもの**であって、
`compactSpace_arc` で消せる種類の但し書きではない。
★残るのは `ample ⟹ ある冪が very ample`（Serre）と、原文の吹き上げ `Y` の段である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

/-- ★★コンパクト空間上の連続関数は**両側**に有界。

★`exists_lower_bound_of_continuous`（`ArchBound.lean`）を `g` と `-g` に当てるだけ。
`Proposition 1.4, (iii)` は差の絶対値を測るので上下ともに要る。 -/
theorem exists_abs_bound_of_continuous {T : Type*} [TopologicalSpace T]
    [CompactSpace T] [Nonempty T] (g : T → ℝ) (hg : Continuous g) :
    ∃ C : ℝ, 0 ≤ C ∧ (∀ p, -C ≤ g p) ∧ ∀ p, g p ≤ C := by
  obtain ⟨C₁, hC₁, hlo⟩ := exists_lower_bound_of_continuous g hg
  obtain ⟨C₂, hC₂, hhi⟩ := exists_lower_bound_of_continuous (fun p => -g p) hg.neg
  refine ⟨max C₁ C₂, le_max_of_le_left hC₁, fun p => ?_, fun p => ?_⟩
  · have : -max C₁ C₂ ≤ -C₁ := by simp [neg_le_neg_iff]
    linarith [hlo p]
  · have h := hhi p
    have : C₂ ≤ max C₁ C₂ := le_max_right _ _
    simp only [neg_le_neg_iff] at h
    linarith

/-- ★★★★★★**[GenEll] Proposition 1.4, (ii)** —— `X` が固有なら `ht_D ≳ 0`。

原文 (GenEll p.6):
> by global sections [for instance, if the line bundle LQ is ample], then

★★★**原文の証明が名指しした根拠がそのまま補題になっている**:

> to a section of L ⊗M over X such that |t|L⊗M ≤1 on Xarc [where we recall that

1. `compactSpace_arc`（`UltraCompact.lean`）—— 固有性から `X^arc` がコンパクト
2. `exists_lower_bound_of_continuous` —— コンパクト空間上の連続関数は下に有界

★★`Prop16.lean` の `prop_1_4_ii` は同じ結論を `ArcModel`（＝射影モデル）から出していた。
本定理は**射影モデルを要らなくした**版である。

★**定数 `C` は `F` にも `x` にも依らない**ので、`X(ℚ̄)` 全体で `≳ 0` が成り立つ。 -/
theorem prop_1_4_ii_of_proper {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (complexPoints X)]
    (D : ArithCartier X)
    (hcont : @Continuous _ _ (arcTopology X) inferInstance D.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        -C ≤ htArith F D xF := by
  letI := arcTopology X
  haveI := compactSpace_arc hval
  obtain ⟨C, hC, hlo⟩ := exists_lower_bound_of_continuous D.green hcont
  refine ⟨C, hC, fun F _ _ xF => ?_⟩
  have h1 : (0:ℝ) ≤ degNormalized (idealADiv F (pullbackIdeal F D.divisor xF)) :=
    degNormalized_nonneg_of_isEffective F _ (idealADiv_isEffective F _)
  have h2 : -C ≤ (archADiv F D.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) :=
    archADiv_sum_div_finrank_ge F D.green xF C hC hlo
  have h3 : htArith F D xF
      = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
        + (archADiv F D.green xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) := by
    rw [htArith, degNormalized, degNormalized, deg_pullbackADiv]
    ring
  rw [h3]; linarith

/-- ★★★★★★**[GenEll] Proposition 1.4, (iii) のアルキメデス側** ——
`X` が固有なら、**因子が同じで計量が違う** 2 つの高さの差は一様に有界。

原文 (GenEll p.6):
> line bundle LQ on XQ. In particular, it makes sense to write [htL] or [htLQ] for the

★★`HeightMetric.lean` の `htArith_sub_abs_le` は同じ結論を `ArcModel` から出していた。
本定理は**射影モデルを要らなくした**版である。

★★★これが「`ht_L̄` の BD-類は `L_ℚ` の同型類にしか依らない」の**計量の側**である
（★因子の側＝`APic` の類にしか依らないことは `htArith_congr_of_pullbackAPic_eq`）。 -/
theorem htArith_sub_abs_le_of_proper {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (complexPoints X)]
    (D E : ArithCartier X) (hdiv : D.divisor = E.divisor)
    (hD : @Continuous _ _ (arcTopology X) inferInstance D.green)
    (hE : @Continuous _ _ (arcTopology X) inferInstance E.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        |htArith F D xF - htArith F E xF| ≤ C := by
  letI := arcTopology X
  haveI := compactSpace_arc hval
  obtain ⟨C₁, hC₁, hlo₁, hhi₁⟩ := exists_abs_bound_of_continuous D.green hD
  obtain ⟨C₂, hC₂, hlo₂, hhi₂⟩ := exists_abs_bound_of_continuous E.green hE
  refine ⟨C₁ + C₂, by linarith, ?_⟩
  intro F _ _ xF
  have h1 := htArith_eq_add F D xF
  have h2 := htArith_eq_add F E xF
  rw [hdiv] at h1
  have hd1lo := archADiv_sum_div_finrank_ge F D.green xF C₁ hC₁ hlo₁
  have hd1hi := archADiv_sum_div_finrank_le F D.green xF C₁ hhi₁
  have hd2lo := archADiv_sum_div_finrank_ge F E.green xF C₂ hC₂ hlo₂
  have hd2hi := archADiv_sum_div_finrank_le F E.green xF C₂ hhi₂
  rw [abs_le]
  refine ⟨by rw [h1, h2]; linarith, by rw [h1, h2]; linarith⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Proposition 1.4` は 4 条あり、
(iv) の射影埋め込み（`ample ⟹ very ample` と原文の吹き上げ `Y`）が残っているためである。
★本ファイルが外したのは (ii)(iii) の**射影モデルの但し書き**だけである。 -/

def exists_abs_bound_of_continuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(コンパクト空間上の連続関数は両側に有界)",
    sectionId := "genell-prop-1-4" }

def prop_1_4_ii_of_proper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (ii)(固有性から——射影モデルを要らなくした)",
    sectionId := "genell-prop-1-4" }

def prop_1_4_ii_of_proper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "compactSpace_arc(固有性から X^arc のコンパクト性——付値判定法)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 7,
    .citation "[ABC3]" "exists_lower_bound_of_continuous(コンパクト空間上の連続関数は下に有界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_lower_bound_of_continuous") 7,
    .implicitStep
      ("★★原文の証明が名指しした根拠(「X^arc はコンパクトだった！」)がそのまま補題になった。" ++
       "Prop16.lean の prop_1_4_ii は同じ結論を ArcModel(射影モデル)から出していたが、" ++
       "本定理は射影モデルを要らなくした") 7,
    .implicitStep
      ("★逸脱: L = O_X(D) との対応は含まない(htArith は ArithCartier から直に高さを作る)。" ++
       "また原文の「L_ℚ のある正冪が大域切断で生成される」という仮定は、" ++
       "因子表示では有限素点側が自動で非負になるため要らない") 7,
    .implicitStep
      "★Continuous D.green は仮定として受ける(GreenFn は任意の実数値関数のため)" 7 ]

def htArith_sub_abs_le_of_proper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iii)(計量の側——固有性から。射影モデルを要らなくした)",
    sectionId := "genell-prop-1-4" }

def htArith_sub_abs_le_of_proper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "compactSpace_arc(固有性から X^arc のコンパクト性)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 7,
    .citation "[ABC3]" "htArith_congr_of_pullbackAPic_eq(因子の側——APic の類にしか依らない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_congr_of_pullbackAPic_eq") 6,
    .implicitStep
      ("★本定理は (iii) の計量の側である。原文は (iii) を (i)(ii) から出すが" ++
       "(L ⊗ M⁻¹ に当てる)、因子表示では計量の差を直に測るほうが短い") 7 ]

end ABC3.Found.GenEll
