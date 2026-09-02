/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSemistableAll
import ABC3.Found.GenEll.VeluGoodPrimeJ
import ABC3.Meta.Claim

/-!
# 第 1439 ブロック —— **★★★★★★★★★★★★★★★★★★★★★★★★
Vélu の商の半安定性が「`j` の整性」1 本になった**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1436 の仮定を**真のもの**に直す

第 1436 の `semistableAt_veluQuot_all_of_goodMem` が残していた仮定
「`p ∣ l` かつ良い素点のとき核の座標が `p` で整」は、
★**一般には成り立たない**——局所体の分岐が `e ≥ l−1` なら形式群の `l`-捩れ
`Ê(𝔪)[l]` が非自明になり、核が深い点を含みうる。

☆本ブロックは第 1438 を使って同じ場合を **`jExp p E′ ≥ 0`（`j(E′)` の整性）**に置き換える。
★★★こちらは**真**である——`E` が `p` で良還元なら同種な `E′` の `j` も整
（古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性から出る）。

| 場合 | 使うもの | 状態 |
|---|---|---|
| `p ∤ l` | 第 1417 | ★**無条件** |
| `p ∣ l`・悪い素点 | 第 1436 | ★**無条件** |
| `p ∣ l`・良い素点 | 第 1438 | ☆**`jExp p E′ ≥ 0`** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Vélu の商は半安定**——★**残る仮定は `j(E′)` の整性 1 本**（第 1439）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮定は

* `SemistableAt p E`
* `l` が素数で `5 ≤ l`（界面の側から取れる——第 1434 の測定）
* `addOrderOf Q = l`
* ★☆**`p ∣ l` かつ良い素点のときだけ** `0 ≤ jExp p E′`

の 4 つだけである。

★★★最後の 1 つは**真**の主張である——`E` が良還元なら同種な `E′` の `j` は整。
☆古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性、
あるいは形式群／Néron モデル（Néron–Ogg–Shafarevich）から出る。
★どちらも mathlib には無い（2026-09-02 確認）。 -/
theorem semistableAt_veluQuot_all_of_jExp [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E)
    {l : ℕ} (hl : l.Prime) (hl5 : 5 ≤ l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic]
    (hjE' : ¬ IsUnit ((l : primeSubring p)) → 0 ≤ jExp p E →
      0 ≤ jExp p (veluQuotientFull E
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hodd : l ≠ 2 := by omega
  by_cases hlu : IsUnit ((l : primeSubring p))
  · -- ★`p ∤ l` は第 1417 で無条件
    exact semistableAt_veluQuot_of_not_dvd_free p E hss hl hodd hlu Q hQ
  · -- ★`p ∣ l` なら `l ≥ 5` から `p ∤ 6`（第 1435・1438）
    have h2 := valAdd_2_eq_zero p hl hl5 hlu
    have h48 := valAdd_48_eq_zero p hl hl5 hlu
    have h864 := valAdd_864_eq_zero p hl hl5 hlu
    have h1728 := valAdd_1728_eq_zero p hl hl5 hlu
    by_cases hj : 0 ≤ jExp p E
    · exact semistableAt_veluQuot_goodPrime_of_jExp p E hss hj hl hodd h2 h48 h864 h1728
        Q hQ (hjE' hlu hj)
    · exact semistableAt_veluQuot_badPrime_all p E hss (by omega) hl hodd h48 h864 Q hQ

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_all_of_jExp.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は半安定。★残る仮定は良い素点 p ∣ l での j の整性 1 本)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_all_of_jExp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_of_not_dvd_free(第 1417、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_of_not_dvd_free") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_badPrime_all(第 1436、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_badPrime_all") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_goodPrime_of_jExp(第 1438-1439、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_goodPrime_of_jExp") 1,
    .citation "[Sil]" "The Arithmetic of Elliptic Curves, II.2 / Modular polynomials(j の整性)"
      (.absent
        ("`E` が `p` で良還元なら同種な `E′` の `j` も整である。" ++
         "★古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性から出るが、" ++
         "mathlib にモジュラー多項式は無い(2026-09-02 確認)。" ++
         "☆別の道は形式群／Néron モデル(Néron–Ogg–Shafarevich)で、これも無い")) 17,
    .implicitStep
      ("★★★★★**2026-09-02（第 1439）**——第 1436 が残していた仮定" ++
       "「核の座標が `p` で整」は**一般には偽**（分岐が `e ≥ l−1` なら形式群の " ++
       "`l`-捩れが非自明）だったので、第 1438 を使って **`jExp p E′ ≥ 0`** に置き換えた。" ++
       "☆こちらは真の主張である。" ++
       "★これで `Lemma 3.5` の半安定性に残るのは**この 1 本だけ**になった。") 17 ]

end ABC3.Found.GenEll
