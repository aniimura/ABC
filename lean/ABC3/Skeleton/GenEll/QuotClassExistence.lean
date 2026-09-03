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

/-! ## ★★★★★★★★★★★★★★★★★★★★閉じた（2026-09-02、第 1349）

☆詰まっていたのは **`DecidableEq` のインスタンスの食い違い**だけだった。
★直し方は 2 つ:

1. `exists_ext_point_of_stable`（第 1346）の結論から `letI` を外し、
   証明の中で `addOrderOf_point_congr_dec`（第 1347）を 1 回使って**具体側に寄せる**
2. `isElliptic_veluQuotientFull_nsmul_nf`（第 1335）に `[DecidableEq K]` の
   **明示の束縛**を付ける（他の補題と同じ手）

★★★これで `EllModuliWitness` の `degInfJ_quotLCyclicJ`（#2）が閉じた。
-/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`HasLCyclicJ` なら商の類は存在する**——★**無条件**（第 1349）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ`（#2）に残っていた
**「商の類の存在」そのもの**である。 -/
theorem exists_isQuotClassJ_of_hasLCyclic (x : RealizedClass) {l : ℕ} (hl : l.Prime)
    (hl5 : 5 ≤ l)
    (hcyc : HasLCyclicJ x.rep.toSSCurve l)
    (hpr : x.rep.toSSCurve.PrimeToLocalHeights l) :
    ∃ y : RealizedClass, IsQuotClassJ x l y.1 := by
  obtain ⟨Q, hQ, hst⟩ := exists_stablePoint_of_hasLCyclicJ x.rep.toSSCurve l hl hcyc
  obtain ⟨M₀, hfin, hdeg, Q₀, hQ₀⟩ := exists_ext_point_of_stable x.rep.toSSCurve hl Q hQ hst
  haveI := hfin
  have hfr : Module.finrank x.rep.toSSCurve.fld (extField x.rep.toSSCurve M₀)
      = Module.finrank x.rep.toSSCurve.fld M₀ := rfl
  have hdeg' : Module.finrank x.rep.toSSCurve.fld (extField x.rep.toSSCurve M₀) < l := by
    have h2 := hl.two_le
    omega
  have hEj : (x.rep.toSSCurve.ext M₀).j = x.cls := by
    rw [SSCurve.ext_j]
    exact RealizedClass.rep_j x
  have hEpr : (x.rep.toSSCurve.ext M₀).PrimeToLocalHeights l :=
    x.rep.toSSCurve.primeToLocalHeights_ext M₀ hl hdeg' hpr
  have hEmult : (x.rep.toSSCurve.ext M₀).HasMultRed :=
    x.rep.toSSCurve.hasMultRed_ext M₀ x.rep.multRed
  obtain ⟨Q₁, hQ₁⟩ : ∃ Q₁ : (x.rep.toSSCurve.ext M₀).W.toAffine.Point, addOrderOf Q₁ = l :=
    ⟨Q₀, hQ₀⟩
  clear hfr hdeg' hdeg hQ₀ hst hQ hpr hcyc
  revert hEj hEpr hEmult Q₁ hQ₁
  generalize (x.rep.toSSCurve.ext M₀) = E₁
  intro hEj hEpr hEmult Q₁ hQ₁
  exact exists_isQuotClassJ x E₁ hEj hEpr Q₁ hQ₁
    (isElliptic_veluQuotientFull_nsmul_nf' E₁.fld E₁.W hQ₁)
    (fun p => semistableAt_veluQuot_ss E₁ hl hl5 Q₁ hQ₁ p)
    (hasMultRed_quotSSCurve E₁ hl hEpr hEmult Q₁ hQ₁
      (isElliptic_veluQuotientFull_nsmul_nf' E₁.fld E₁.W hQ₁)
      (fun p => semistableAt_veluQuot_ss E₁ hl hl5 Q₁ hQ₁ p))

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_isQuotClassJ_of_hasLCyclic.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(HasLCyclicJ なら商の類は存在する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ_of_hasLCyclic.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_ext_point_of_stable(第 1346、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_ext_point_of_stable") 1,
    .citation "[ABC3]" "isElliptic_veluQuotientFull_nsmul_nf(第 1335、無条件)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isElliptic_veluQuotientFull_nsmul_nf") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_ss(第 1345、節点)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.semistableAt_veluQuot_ss") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1349）**——`EllModuliWitness` の " ++
       "`degInfJ_quotLCyclicJ`（#2）に残っていた" ++
       "**「商の類の存在」そのもの**である。") 3 ]

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
