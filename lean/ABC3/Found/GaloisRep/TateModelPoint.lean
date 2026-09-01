/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateParamMap
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.PointTransport
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace

/-!
# 第 970 ブロック —— **★★★★★★★★★★★★★★★★★★★★Tate モデルの上の点と `j`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_torsion`（第 965）は Tate モデルの上の点 `P` と
`j` の一致 `hW′j` を受ける。★本ブロックは**その 2 つを大域のデータから作る**。

☆道は 5 段:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `(E_q) ⊗ Lv = C₀ • (E ⊗ Lv)` | `tateParamR_spec` ＋ `tateModel_baseChange`（第 944） |
| 2 | Vélu の商の楕円性 | `isElliptic_veluQuotient_vcPoint`（第 969） |
| 3 | `j` の一致（`C₀ • (E ⊗ Lv)` 側） | `j_map_velu_vcPoint`（第 967） |
| 4 | 点の位数は運べる | `addOrderOf_rhPoint` ＋ `addOrderOf_vcPoint`（在庫） |
| 5 | 曲線の等式に沿って点を運ぶ | `exists_point_image_eq`（第 968） |

★★配管の注意（実測）: 最後に `rw [himg, hbase]` を **`.j` の上で**やると
`IsElliptic` のせいで motive が壊れる。☆曲線の等式を先に作って
`j_congr_curve`（第 913）で渡すこと——第 944 と同じ穴である。

★★楕円性を `∃`/`∧` の中で運ぶと、後続の `.j` がそれをインスタンスとして
使えない。☆そこで結論を **`∀ _hell : …, j の等式`** の形にした。
呼ぶ側は `haveI` で入れてから使う。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll IsDedekindDomain NumberField WeierstrassCurve Finset
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★**大域の Vélu の商から、Tate モデルの上の
位数 `l` の点と `j` の一致を作る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 970）**——これで `minDeltaExp_eq_mul_of_torsion`（第 965）の
`P`・`hP`・`hP0`・`hW′j` がすべて揃う。 -/
theorem exists_point_j_tateModel {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    {l : ℕ} (hl : l.Prime) {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    ∃ P : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Point,
      l • P = 0 ∧ P ≠ 0 ∧
      ∀ _hell : (veluQuotientFull
          ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic,
        (E'.baseChange (p.adicCompletion L)).j
          = (veluQuotientFull
            ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j := by
  set Lv := p.adicCompletion L
  set R := p.adicCompletionIntegers L
  set φL : L →+* Lv := algebraMap L Lv with hφL
  set φR : R →+* Lv := algebraMap R Lv with hφR
  -- ★段 1: Tate モデルは `E ⊗ Lv` の変数変換
  obtain ⟨hq, C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase : (tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map φR
      = (C₀.map φR) • (E.map φL) :=
    tateModel_baseChange (E.baseChange Lv) h hCE
  -- ★段 2・3: 楕円性と `j` の一致（`C₀ • (E ⊗ Lv)` 側）
  haveI := isElliptic_veluQuotient_vcPoint φL (C₀.map φR) E E' hQ h2K hE'
  have hj := j_map_velu_vcPoint φL (C₀.map φR) E E' hQ h2K hE'
  -- ★段 4: 点の位数は運べる
  have hQ1 : addOrderOf (rhPoint φL E Q) = l := by rw [addOrderOf_rhPoint φL E Q, hQ]
  have hQ2 : addOrderOf (ABC3.Found.GenEll.vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q)) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q), hQ1]
  -- ★段 5: 曲線の等式に沿って点を運ぶ
  obtain ⟨P, hP, hP0, himg⟩ := exists_point_image_eq hbase hl _ hQ2
  refine ⟨P, hP, hP0, fun hell => ?_⟩
  haveI := hell
  have hcurve : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map φR)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull ((C₀.map φR) • (E.map φL))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q)))) := by
    rw [himg, hbase]
  rw [ABC3.Found.GenEll.j_congr_curve hcurve]
  exact hj

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_j_tateModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(大域の Vélu の商から Tate モデルの上の点と j の一致を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
