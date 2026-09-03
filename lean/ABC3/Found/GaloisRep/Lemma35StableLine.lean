/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Unconditional
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# 第 1199 ブロック —— **`Lemma 3.5` は有限拡大の上で確かめれば足りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——第 1194 の道の完成

`Skeleton/GenEll/LCyclicReading.lean` の第 1194 で立てた道:

> `l`-巡回部分群 `H` の点は `L'' ≔ L(H)`（次数は `l−1` を割る）で**有理**になるので、
> そこで `Lemma 3.5` を回して `L` へ降ろす。

★本ファイルは**降ろす段が要らない**ことを示す——`Lemma 3.5` の不等式に現れる
`ht^Falt` と `deg∞` は**半安定な曲線では底変換で変わらない**（第 738・第 737）ので、
`L''` で得た不等式が**そのまま `L` の不等式**である。

| 量 | 底変換での振る舞い | 在庫 |
|---|---|---|
| `ht^Falt` | 変わらない | `htFaltOf_baseChange_of_semistable`（第 738） |
| `deg∞` | 変わらない | `degInfOf_baseChange_of_semistable`（第 737） |

☆したがって `L''` の上で `Q` が有理でありさえすれば、`L` の上の `Lemma 3.5` が出る。
★★★これで `l`-巡回の読みの橋（`LCyclicReading` の節点 2）は**閉じる**。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.5 —— 有限拡大の上で確かめれば足りる**（第 1199）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`E` が半安定なら、`Lemma 3.5` の高さ不等式は
**どの有限拡大 `L''` の上で確かめても同じ**である。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の第 1194 の道の**要**である
——`l`-巡回部分群 `H` の点は `L'' ≔ L(H)`（第 1195）で有理になるので、
そこで Vélu の商を作れば（第 1197）、`L` の上の不等式がそのまま出る。 -/
theorem lemma_3_5_height_ineq_over_extension (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic]
      (l : ℕ), l.Prime →
      (∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E) →
      ∀ (L'' : Type) [Field L''] [NumberField L''] [Algebra L L''] [IsScalarTower ℚ L L'']
        (E'' : WeierstrassCurve L''),
      ∀ _hE0 : (E.baseChange L'').IsElliptic, ∀ _hE1 : E''.IsElliptic,
      ∀ Q : (E.baseChange L'').toAffine.Point, addOrderOf Q = l →
      E'' = veluQuotientFull (E.baseChange L'')
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) →
      (∀ P : HeightOneSpectrum (𝓞 L''), SemistableAt P (E.baseChange L'')) →
      (∀ P : HeightOneSpectrum (𝓞 L''), SemistableAt P E'') →
      (∀ P : HeightOneSpectrum (𝓞 L''), jExp P (E.baseChange L'') < 0 →
        ¬ ((l : ℤ) ∣ jExp P (E.baseChange L''))) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_height_ineq eps heps
  refine ⟨C, ?_⟩
  intro L _ _ E _ l hl hss L'' _ _ _ _ E'' hE0 hE1 Q hQ hE' hssE hssE' hcop
  haveI := hE0
  haveI := hE1
  have hmain := hC L'' (E.baseChange L'') E'' l hl Q hQ hE' hssE hssE' hcop
  rwa [htFaltOf_baseChange_of_semistable L L'' E hss,
    degInfOf_baseChange_of_semistable L L'' E hss] at hmain

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_5_height_ineq_over_extension.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有限拡大の上で確かめれば足りる——ht^Falt も deg∞ も半安定なら不変)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_height_ineq_over_extension.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_5_height_ineq(第 1149、条なし、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_5_height_ineq") 1,
    .citation "[ABC3]" "htFaltOf_baseChange_of_semistable(第 738、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFaltOf_baseChange_of_semistable") 1,
    .citation "[ABC3]" "degInfOf_baseChange_of_semistable(第 737、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_baseChange_of_semistable") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1199）**——`Skeleton/GenEll/LCyclicReading.lean` の" ++
       "第 1194 の道の**要**である。☆`l`-巡回部分群 `H` の点は `L'' ≔ L(H)`（第 1195）で" ++
       "有理になるので、そこで Vélu の商を作れば（第 1197）、`L` の上の不等式がそのまま出る。" ++
       "★降ろす段（第 1185）は `Δ_min` の関係には要るが、**高さの不等式には要らなかった**" ++
       "——`ht^Falt` も `deg∞` も半安定なら底変換で不変だからである。") 2 ]

end ABC3.Found.GaloisRep
