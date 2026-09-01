/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.Lemma35Concrete
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

/-! ## ★★★★★★★★★★★★商の類の `deg∞` -/

/-- ★★★★★★★★★★★★
**商の類の `deg∞` は `l` 倍**——★**無条件**（第 1242）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`degInfJ_eq`（類の値は曲線の `deg∞`）と
`degInfOf_eq_of_local`（局所の関係を足し上げる、在庫）を繋いだだけ。

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ` が要る形である
——局所の関係 `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` は `Lemma 3.2, (ii)`。 -/
theorem degInfJ_quot_eq (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ)
    (hloc : ∀ p : HeightOneSpectrum (𝓞 E.fld),
      minDeltaExp p (veluQuotientFull E.W S) = l * minDeltaExp p E.W) :
    degInfJ (quotSSCurve E S hell hss).j = (l : ℝ) * degInfJ E.j := by
  rw [degInfJ_eq, degInfJ_eq]
  exact ABC3.Found.GaloisRep.degInfOf_eq_of_local E.W (veluQuotientFull E.W S) l hloc

/-! ## ★★★★★★★★★★★★商の類の `ht^Falt` -/

/-- ★★★★★★★★★★★★
**商の類の `ht^Falt` の評価は曲線の評価そのもの**——★**無条件**（第 1243）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`faltingsHeightJ_eq`（類の値は曲線の `ht^Falt`、在庫）で置き換えるだけ。

★★★これが `EllModuliWitness` の `faltingsHeightJ_quotLCyclicJ` が要る形である
——同種写像の高さ評価 `hfalt` をそのまま類の言葉に移す。 -/
theorem faltingsHeightJ_quot_le (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ) (C₀ : ℝ)
    (hfalt : htFaltOf E.fld (veluQuotientFull E.W S)
      ≤ htFaltOf E.fld E.W + 2 * Real.log l + C₀) :
    faltingsHeightJ (quotSSCurve E S hell hss).j
      ≤ faltingsHeightJ E.j + 2 * Real.log l + C₀ := by
  rw [faltingsHeightJ_eq, faltingsHeightJ_eq]
  exact hfalt

/-! ## ★★★★★★★★★★★★★★★★`quotLCyclicJ` の水準へ -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `deg∞` は `l` 倍**——★（第 1244）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`quotLCyclicJ_spec`（在庫）で `S` を取り出し、第 1242 を当てる。
★局所の関係は `Lemma 3.2, (ii)` が与える（仮説として受ける）。

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ` そのものの形である。 -/
theorem degInfJ_quotLCyclicJ_of_local (x : RealizedClass) (l : ℕ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hloc : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
        minDeltaExp p (veluQuotientFull x.rep.toSSCurve.W S)
          = l * minDeltaExp p x.rep.toSSCurve.W) :
    degInfJ (quotLCyclicJ x l).cls = (l : ℝ) * degInfJ x.cls := by
  obtain ⟨S, hell, hss, hcard, hj⟩ := quotLCyclicJ_spec x l hex
  have hq := degInfJ_quot_eq x.rep.toSSCurve S hell hss l (hloc S hcard)
  rw [hj] at hq
  have hx : x.cls = x.rep.toSSCurve.j := (RealizedClass.rep_j x).symm
  rw [hx]
  exact hq

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def degInfJ_quotLCyclicJ_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の deg∞ は l 倍——局所の関係は Lemma 3.2, (ii) を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quot_le.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の ht^Falt の評価は曲線の評価そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quot_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の deg∞ は l 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

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
