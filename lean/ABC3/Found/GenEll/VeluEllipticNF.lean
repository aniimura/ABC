/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluEllipticField
import ABC3.Meta.Claim

/-!
# 第 1335 ブロック —— **数体の上で Vélu の商は楕円**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——埋め込みを 1 本取るだけ

第 1334 は「`ℂ` への埋め込みがあれば商は楕円」だった。
★数体は代数的で `ℂ` は代数閉体だから、埋め込みは `IsAlgClosed.lift` が与える。

☆これで `VeluQuotOK` の**楕円性の穴が閉じる**。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

open scoped Classical

/-- ★★★★★★**数体は `ℂ` に埋め込める**（`IsAlgClosed.lift`、在庫）。 -/
noncomputable def embedComplex (K : Type) [Field K] [NumberField K] : K →+* ℂ :=
  (IsAlgClosed.lift (R := ℚ) (S := K) (M := ℂ)).toRingHom

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**数体の上で `⟨Q⟩` による Vélu の商は楕円**——★**無条件**（第 1335）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `VeluQuotOK` の**楕円性の側そのもの**である。 -/
theorem isElliptic_veluQuotientFull_nsmul_nf (K : Type) [Field K] [NumberField K]
    [inst : DecidableEq K] (W : WeierstrassCurve K) [W.IsElliptic]
    {Q : W.toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    (veluQuotientFull W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic := by
  -- ★`Decidable` のインスタンスを古典的なものに揃える（配管）
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  exact isElliptic_veluQuotientFull_nsmul_of_embed W (embedComplex K) hl hQ

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**位数の仮定だけで Vélu の商は楕円**——★**無条件**（第 1335）。

☆`l = 0`（位数無限）のときは点集合が空で商は `W` 自身だから、
`0 < l` は要らない。★これが `VeluQuotOK` の楕円性の側にぴたりと嵌まる形である。 -/
theorem isElliptic_veluQuotientFull_nsmul_nf' (K : Type) [Field K] [NumberField K]
    [DecidableEq K] (W : WeierstrassCurve K) [W.IsElliptic]
    {Q : W.toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l) :
    (veluQuotientFull W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic := by
  rcases Nat.eq_zero_or_pos l with rfl | hl
  · rw [Finset.range_zero, Finset.erase_empty, Finset.image_empty, veluQuotientFull_empty]
    infer_instance
  · exact isElliptic_veluQuotientFull_nsmul_nf K W hl hQ

/-! ## ★出典の紐付け(`.src`) -/

def isElliptic_veluQuotientFull_nsmul_nf'.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数の仮定だけで Vélu の商は楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

def embedComplex.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数体は ℂ に埋め込める。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_nsmul_nf.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数体の上で ⟨Q⟩ による Vélu の商は楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_nsmul_nf.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isElliptic_veluQuotientFull_nsmul_of_embed(第 1334、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isElliptic_veluQuotientFull_nsmul_of_embed") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1335）**——`VeluQuotOK` の**楕円性の側そのもの**である。" ++
       "☆残るのは半安定性（良い素点）だけになった。") 3 ]

end ABC3.Found.GenEll
