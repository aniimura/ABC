/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuHeadDenomFree
import ABC3.Meta.Claim

/-!
# 第 1103 ブロック —— **`Lemma 3.5` を条なしにする残り**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

`Found/GaloisRep/Lemma35Unconditional.lean` の `Lemma 3.5` は原典に無い仮説を
2 つ置いている（`l ≠ 2`・`d + 1 < l`）。本ファイルは **`d + 1 < l` を外す**残りを
節点に割り、進捗枠を置く。

☆`l ≠ 2` の方は別枠（`l = 2` の Vélu の分岐、20-40 ブロック）。

## ★★★★★★★★積んだ土台（`Found/`、`sorry` 0）

| 第 | 内容 |
|---|---|
| 1092 | `l^k·∑ f_i·inv(1−ζ^i)^k = ∑ f_i·∏_{j≠i}(1−ζ^j)^k` |
| 1097 | 6 種の頭項を分母なしの和に（橋） |
| 1098 | `DX`・`DY`・`D²X` の「対」を分母なしの和に（橋） |
| 1102 | **`E` 版 6 種の定義**（`S(η) = ∑ k·η^k` で無条件）と橋 |

★機構は在庫の `one_sub_mul_sum_nsmul`——`(1 − η)·S(η) = −l`、`IsUnit` 不要。

## ★★★★★★★★★★残り 5 節点（進捗枠 **0 / 5**）

| # | 節点 | 内容 | 重み |
|---|---|---|---|
| 1 | `sum_mu_dxpairE_zero` | `∑ DX` の `E` 版（`hu` なし） | 8 |
| 2 | `sum_mu_d2xpairE` | `∑ D²X` の `E` 版（`hu` なし） | 10 |
| 3 | `sum_veluV2E_eq_sum_tateDYpairE` | Vélu の `v` と `DY` の一致の `E` 版 | 6 |
| 4 | `c4_c6_velu_tateE` | `c₄`・`c₆` の恒等式の `E` 版 | 10 |
| 5 | `minDeltaExp_eq_mul_at_bad_prime_dvd` | `p ∣ l` の悪い素点で `Δ_min` が `l` 倍 | 12 |

☆総重み 46。★**本ファイルに型を置いたのは節点 1-2 だけである**（節点 3-5 は 1-2 が済んでから型を確定させる——`c₄`・`c₆` の `E` 版の形は節点 2 の結論に依存する）。★節点 1-2 が本体（μ_l の指標和を `ℤ[ζ_l]` の中で処理する段）。

### ★★★降下の道（第 1099-1100 で確定）

`E` 版の主張は `ℤ[ζ_l]` の中の恒等式（`q` の各次数ごと）に分かれる。
`ℤ[ζ_l]` は ℤ-捻れ自由なので `ℤ[ζ_l] → ℚ(ζ_l)` は単射で、
`ℚ(ζ_l)` では `1 − ζ^i ≠ 0` すなわち可逆だから**既存の `hu` つきの補題が使える**。
★単射性で降ろせば `hu` が消える。
☆「商体 `K` に移す」が通らない（第 1091）のは `q` の収束が壊れるからだが、
`q` の次数ごとに見れば収束は無関係である。
-/

namespace ABC3.Skeleton.GenEll

open Finset ABC3.Meta ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★★★**節点 1**——`∑ DX` の `E` 版（`hu` を受けない）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`sum_mu_dxpair_zero`（`TateODE.lean:144`）の分母なし版である。 -/
theorem sum_mu_dxpairE_zero {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) (h2 : (2 : R) ≠ 0) :
    ∑ i ∈ (range l).erase 0,
        (tateDXtermE l (ζ ^ i)
          + (l : R) ^ 3 * (tateDXtail (ζ ^ i) q hq
              - (tateDXterm (q * (ζ ^ i) ^ (l - 1))
                  + tateDXtail (q * (ζ ^ i) ^ (l - 1)) q hq))) = 0 := by
  sorry

def sum_mu_dxpairE_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 1——∑ DX の E 版。hu を受けない)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_dxpairE_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_tateDXterm(橋、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateDXterm") 1,
    .implicitStep
      ("★`ℤ[ζ_l]` は ℤ-捻れ自由なので `ℚ(ζ_l)` へ単射。" ++
       "☆`ℚ(ζ_l)` では `1 − ζ^i` が可逆なので既存の `sum_mu_dxpair_zero` が使え、" ++
       "単射性で降ろせば `hu` が消える。") 8 ]

/-- ★★★★★★★★★★**節点 2**——`∑ D²X` の `E` 版（`hu` を受けない）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`sum_mu_d2xpair`（`TateODE.lean:207`）の分母なし版。★`c4_velu_tate` の本体である。 -/
theorem sum_mu_d2xpairE {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    120 * ∑ i ∈ (range l).erase 0,
        (tateD2XtermE l (ζ ^ i)
          + (l : R) ^ 4 * (tateD2Xtail (ζ ^ i) q hq
              + (tateD2Xterm (q * (ζ ^ i) ^ (l - 1))
                  + tateD2Xtail (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = (l : R) ^ 4 * (((l : R) ^ 4 - 1)
          + 240 * ((l : R) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
              - evalAdic (sigmaSeries 3) q hq)) := by
  sorry

def sum_mu_d2xpairE.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(節点 2——∑ D²X の E 版。c4_velu_tate の本体)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_d2xpairE.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_tateD2Xterm(橋、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateD2Xterm") 1,
    .implicitStep "☆節点 1 と同じ降下（`ℤ[ζ_l] → ℚ(ζ_l)` の単射性）。" 10 ]

end ABC3.Skeleton.GenEll
