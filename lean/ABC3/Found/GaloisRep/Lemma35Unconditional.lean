/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateLocalModel
import ABC3.Meta.Claim

/-!
# 第 1085 ブロック —— **`Lemma 3.5` を仮説なしで**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★到達点

第 1083 で `hfin` が無条件に出たので、`Lemma 3.5` の高さ不等式は
**外部からの入力を一切受けずに**証明できる。

## ★★★★原典との差（逸脱の記録）

原文の仮説は
「`ε > 0`」「`l` 素数」「`H_F` が `l`-cyclic」「`E_H = E/H`」
「`l` は悪い乗法還元の全素点での局所高さと素」だけである。

☆本定理はそれに加えて **2 つの仮説**を置いている:

| 追加した仮説 | 由来 |
|---|---|
| `l ≠ 2` | Vélu の `±` の対（`veluWFull` の `/2`）が不動点を持たないため |
| `d + 1 < l` | 第 1044（`isUnit_natCast_at_bad_prime`）が `p ∤ l` を出すのに使う。原典は Lemma 3.2, (i) で直接 `μ_l` と同定するのでこれを要さない |

★★**後続への影響は無い**——`Lemma 3.7`・`Theorem 3.8` は
`100·d·(ht^Falt + C·d^ε) ≤ l` を満たす `l` を取るので、
`l ≠ 2` も `d + 1 < l` も自動的に成り立つ。

☆この 2 つを外すには、悪い素点で `p ∣ l` の場合の局所計算
（Tate モデルの Vélu を `l` が単元でない環で行う）が要る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 高さの不等式（外部入力なし）**（第 1085）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1085）**——`hfin`（第 1083）で最後の入力が消えた。
☆原典に無い仮説は `l ≠ 2` と `d + 1 < l` の 2 つだけで、
どちらも後続（`Lemma 3.7`・`Theorem 3.8`）が取る `l` では自動的に成り立つ。 -/
theorem lemma_3_5_height_ineq (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime → l ≠ 2 →
      Module.finrank ℚ L + 1 < l →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C :=
  ABC3.Skeleton.GenEll.lemma_3_5_velu eps heps

def lemma_3_5_height_ineq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(外部入力なし。★原典に無い仮説は l ≠ 2 と d+1 < l の 2 つ)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_velu(第 1084、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.lemma_3_5_velu") 1,
    .citation "[ABC3]" "hfin_of_veluQuotientFull(hfin そのもの、第 1083、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hfin_of_veluQuotientFull") 1,
    .implicitStep
      ("★★**2026-09-01（第 1085）の逸脱の記録**——原典 (GenEll p.17) の Lemma 3.5 は " ++
       "`l ≠ 2` も `d + 1 < l` も置いていない。" ++
       "☆`l ≠ 2` は Vélu の `±` の対（`veluWFull` の `/2`）が不動点を持たないために要る。" ++
       "☆`d + 1 < l` は第 1044 が `p ∤ l` を出すのに使う（原典は Lemma 3.2, (i) で " ++
       "直接 `μ_l` と同定するのでこれを要さない）。" ++
       "★後続（`Lemma 3.7`・`Theorem 3.8`）は `100·d·(ht^Falt + C·d^ε) ≤ l` を取るので " ++
       "両方とも自動的に成り立ち、影響は無い。" ++
       "☆外すには悪い素点で `p ∣ l` の場合の Tate モデルの Vélu 計算が要る。") 8 ]

end ABC3.Found.GaloisRep
