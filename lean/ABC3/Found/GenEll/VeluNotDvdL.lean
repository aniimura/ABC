/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluBadPrimeAssembly
import ABC3.Meta.Claim

/-!
# 第 1407 ブロック —— **`p ∤ l` の側が閉じた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——良い素点と悪い素点を合わせる

* 良い素点（`0 ≤ jExp p E`）——第 1403 `semistableAt_veluQuot_goodPrime`
* 悪い素点（`jExp p E < 0`）——第 1406 `semistableAt_veluQuot_badPrime`

★★★これで **`semistableAt_veluQuotientFull` の `p ∤ l` の側が完全に閉じた**。

## ☆残る 2 つ

| # | 内容 | 状態 |
|---|---|---|
| 1 | `hcop`（`l ∤ jExp p E`）と `hc4L`（`c₄(E/C) ≠ 0`）を `VeluQuotOK` の側へ通す | ☆インタフェースの変更 |
| 2 | **`p ∣ l`** | ☆**形式群**が要る（`blocked-leaves.json` の `pDivLRevised2026_09_02`） |

☆`hcop` は消費側（`hdag_of_stableLine`）が `PrimeToLocalHeights l` として持っている。
`hc4L` は数学的には自動（同種な曲線は `j` の非整性を共有する）だが、
それを言うにはモジュラー多項式か同種不変性が要る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] `p ∤ l` では Vélu の商は半安定**——★（第 1407）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆良い素点（第 1403）と悪い素点（第 1406）を `jExp p E` の符号で場合分けするだけである。

★★★これで `semistableAt_veluQuotientFull` の **`p ∤ l` の側が完全に閉じた**。 -/
theorem semistableAt_veluQuot_of_not_dvd [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hc4L : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).c₄ ≠ 0) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  by_cases hj : 0 ≤ jExp p E
  · exact semistableAt_veluQuot_goodPrime p E hss hj hl hodd hlu Q hQ
  · exact semistableAt_veluQuot_badPrime p E hss (by omega) hl hodd hlu hcop Q hQ hc4L

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_of_not_dvd.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∤ l では Vélu の商は半安定。★l ∤ jExp p・c₄ ≠ 0)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_of_not_dvd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_goodPrime(第 1403、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_goodPrime") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_badPrime(第 1406、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_badPrime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1407）**——`semistableAt_veluQuotientFull` の " ++
       "**`p ∤ l` の側が完全に閉じた**。" ++
       "☆残るのは（1）`hcop`・`hc4L` を `VeluQuotOK` の側へ通すインタフェースの変更と" ++
       "（2）**`p ∣ l`**（形式群）である。") 17 ]

end ABC3.Found.GenEll
