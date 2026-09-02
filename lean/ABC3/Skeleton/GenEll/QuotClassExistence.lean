/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.VeluSemistable
import ABC3.Skeleton.GenEll.TateLocalModelK
import ABC3.Found.GenEll.ExtPoint
import ABC3.Found.GenEll.QuotClassExists
import ABC3.Found.GenEll.LCyclicPoint
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Meta.Claim

/-!
# 第 1348-1349 ブロック —— **商の類は存在する（`HasLCyclicJ` から）**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`IsQuotClassJ` の存在

`EllModuliWitness` の `degInfJ_quotLCyclicJ`（#2）に残っていた `sorry` は
**「`HasLCyclicJ` から `∃ y, IsQuotClassJ x l y.1` を出す」**一点だった。

☆材料は揃っている:

| # | 段 | 出どころ |
|---|---|---|
| 1 | `Gal`-安定な点 `Q`（`E.alg` の上） | `exists_stablePoint_of_hasLCyclicJ` |
| 2 | `L(H)` へ降ろし `ℂ` の中へ運ぶ | 第 1346 |
| 3 | `SSCurve` として持ち上げる | 第 1343 |
| 4 | 商は楕円 | 第 1335 |
| 5 | 商は半安定 | 第 1345（節点） |
| 6 | 商は乗法還元を持つ | 第 1348（下） |
-/

namespace ABC3.Skeleton.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★
**商も乗法還元を持つ**——★**無条件**（第 1348）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`E` の乗法還元の素点 `p` で Tate の関係
`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`（`isMuAtBadPrimes_of_veluQuotient_nodeg`、第 1141）
が効く。★`l ≠ 0` なので右辺は `0` でない。 -/
theorem hasMultRed_quotSSCurve (E : SSCurve) {l : ℕ} (hl : l.Prime)
    (hpr : E.PrimeToLocalHeights l) (hE : E.HasMultRed)
    (Q : E.W.toAffine.Point) (hQ : addOrderOf Q = l)
    (hell : (veluQuotientFull E.W
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld),
      SemistableAt p (veluQuotientFull E.W
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))) :
    (quotSSCurve E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).HasMultRed := by
  obtain ⟨p, hp⟩ := hE
  haveI := hell
  have hcop := E.not_dvd_jExp_of_primeToLocalHeights hl hpr
  have hmu := isMuAtBadPrimes_of_veluQuotient_nodeg E.W
    (veluQuotientFull E.W
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    hl Q hQ rfl E.ss hss hcop
  have hmd : minDeltaExp p E.W = max 0 (- jExp p E.W) :=
    minDeltaExp_eq_maxJ_of_semistable p E.W (E.ss p)
  have hlh : E.localHeightAt p = minDeltaExp p E.W := rfl
  rw [hlh] at hp
  have hneg : jExp p E.W < 0 := by
    by_cases hc : jExp p E.W < 0
    · exact hc
    · exact absurd (by rw [hmd]; omega) hp
  have hrel := hmu p hneg
  refine ⟨p, ?_⟩
  show minDeltaExp p (veluQuotientFull E.W
    (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ≠ 0
  rw [hrel]
  have hl0 : (l : ℤ) ≠ 0 := by exact_mod_cast hl.ne_zero
  exact mul_ne_zero hl0 hp

/-! ## ★★★★★★★★★★★★★★★★残っている段の測定（2026-09-02、第 1349）

☆**材料はすべて揃っている**が、最後の組み立てで `DecidableEq` のインスタンスが
食い違って通らなかった。★測ったことを書いておく。

| # | 段 | 状態 |
|---|---|---|
| 1 | `Gal`-安定な点 `Q`（`E.alg` の上） | ★`exists_stablePoint_of_hasLCyclicJ`（在庫） |
| 2 | `L(H)` へ降ろし `ℂ` の中へ運ぶ | ★第 1346（証明済み） |
| 3 | `SSCurve` として持ち上げる | ★第 1343（証明済み） |
| 4 | `j`・`PrimeToLocalHeights`・乗法還元の輸送 | ★第 1343・1347（証明済み） |
| 5 | 商は楕円 | ★第 1335（無条件） |
| 6 | 商は半安定 | ☆第 1345（節点） |
| 7 | 商は乗法還元 | ★第 1348（上、証明済み） |
| 8 | **組み立て** | ☆`DecidableEq` の合わせ込みが残る |

★★☆**(8) の中身**——`exists_ext_point_of_stable`（第 1346）は
`M₀` の `DecidableEq` を**古典的**に固定して返す（`exists_point_descent_of_stable` に合わせた）。
☆一方 `SSCurve` の語彙（`IsQuotClassJ`・`isElliptic_veluQuotientFull_nsmul_nf'`）は
`↥K` の**具体の** `Subtype.instDecidableEq` で合成される。
★橋は `addOrderOf_point_congr_dec`（第 1347、証明済み）にあるが、
`isElliptic_veluQuotientFull_nsmul_nf'` の側の暗黙引数の合わせ込みが要る。

☆**次にやること**——`exists_ext_point_of_stable` の結論から `letI` を外し、
証明の中で `addOrderOf_point_congr_dec` を 1 回使って具体のインスタンスに寄せる。
★そうすれば `SSCurve` の語彙とそのまま噛み合う。
-/

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def hasMultRed_quotSSCurve.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商も乗法還元を持つ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def hasMultRed_quotSSCurve.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1348）**——`IsQuotClassJ` の witness には " ++
       "`DegCurve`（＝`SSCurve` ＋ 乗法還元）が要るので、商の乗法還元を取る段である。") 2 ]

end ABC3.Skeleton.GenEll
