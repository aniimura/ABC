/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Meta.Claim

/-!
# 第 1241 ブロック —— **商の類の存在**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`quotLCyclicJ_spec` を使えるようにする

`quotLCyclicJ`（`EllModuliObjects`）は `∃ y, IsQuotClassJ x l y.1` のときだけ
意味を持つ（そうでなければ自分自身）。

☆位数 `l` の有理点 `Q` と、その Vélu の商が**楕円・半安定・乗法還元**であることから
その存在が出る。★個数の条件 `S.card + 1 = l` は第 1240。

## ★★逸脱の記録

☆商が**乗法還元を持つ**ことは仮説として受ける
（原文は「同種なので自動」と括弧で述べる。`VeluQuotOK` と同じ扱い）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep IsDedekindDomain NumberField ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★
**商の類は存在する**——★（第 1241）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆位数 `l` の有理点 `Q` と、その Vélu の商が楕円・半安定・乗法還元であることから
`IsQuotClassJ` の存在が出る。★個数の条件は第 1240。

★★逸脱: 商が乗法還元を持つことは仮説として受ける
（原文は「同種なので自動」と括弧で述べる）。 -/
theorem exists_isQuotClassJ (x : RealizedClass) {l : ℕ} (hl : l.Prime)
    (Q : x.rep.toSSCurve.W.toAffine.Point) (hQ : addOrderOf Q = l)
    (hell : (veluQuotientFull x.rep.toSSCurve.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
      SemistableAt p (veluQuotientFull x.rep.toSSCurve.W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))))
    (hmult : (quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))
      hell hss).HasMultRed) :
    ∃ y : RealizedClass, IsQuotClassJ x l y.1 := by
  refine ⟨⟨(quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).j,
    ⟨⟨quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss,
      hmult⟩, rfl⟩⟩, ?_⟩
  exact ⟨((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)), hell, hss,
    card_image_pointCoords_nsmul hl hQ, rfl⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def exists_isQuotClassJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の存在——商の乗法還元は仮説として受ける)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "card_image_pointCoords_nsmul(第 1240、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.card_image_pointCoords_nsmul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1241）**——`quotLCyclicJ_spec`（`EllModuliObjects`）を" ++
       "使えるようにする段である。☆★★逸脱: 商が**乗法還元を持つ**ことは" ++
       "仮説として受ける（原文は「同種なので自動」と括弧で述べる。" ++
       "`VeluQuotOK` と同じ扱い）。") 3 ]

end ABC3.Found.GenEll
