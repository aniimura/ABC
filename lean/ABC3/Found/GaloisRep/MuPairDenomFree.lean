/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuHeadDenomFree
import ABC3.Found.GaloisRep.TatePair
import ABC3.Meta.Claim

/-!
# 第 1118 ブロック —— **座標 `X`・`Y` の分母なし版**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

`tateXpair a w q hq` は「`a` での頭項＋尾」＋「`w` での頭項＋尾」−`2σ₁(q)` の形である。
☆`w = q·(ζ^i)^{l−1}` は `I` に入るので `1 − w` は単元で、そちら側に分母の問題は無い。
★問題になるのは **`a = ζ^i` での頭項だけ**である。

したがって `l²` を掛けた形（`Y` は `l³`）を、第 1102 の `E` 版で書き直せる:

    `tateXpairDF = tateXtermE + l²·(残り)`,  `tateYpairDF = tateYtermE + l³·(残り)`

☆定義に `Ring.inverse` が現れないので `p ∣ l` でも意味を持つ。
★これが節点 3（`∑ veluV2 = ∑ DY` の `E` 版）の土台である。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★**`l²·X(a,w,q)` の分母なし版**。 -/
noncomputable def tateXpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) : R :=
  tateXtermE l a
    + (l : R) ^ 2 * (tateXtail a q hq + (tateXterm w + tateXtail w q hq)
        - 2 * evalAdic (sigmaSeries 1) q hq)

/-- ★★★★**`l³·Y(a,w,q)` の分母なし版**。 -/
noncomputable def tateYpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) : R :=
  tateYtermE l a
    + (l : R) ^ 3 * (tateYtail a q hq - (tateXterm w + tateXtail w q hq)
        - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq)

section Bridge

variable [IsDomain R] {l : ℕ} {a : R}

/-- ★★★★★★**橋**——`hu` の下で `tateXpairDF = l²·tateXpair`。 -/
theorem natCast_pow_mul_tateXpair (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - a))
    (hpow : a ^ l = 1) (hsum : ∑ k ∈ range l, a ^ k = 0) (w q : R) (hq : q ∈ I) :
    (l : R) ^ 2 * tateXpair a w q hq = tateXpairDF l a w q hq := by
  rw [tateXpair, tateXpairDF, ← natCast_pow_mul_tateXterm hlu hu hpow hsum]
  ring

/-- ★★★★★★**橋**——`hu` の下で `tateYpairDF = l³·tateYpair`。 -/
theorem natCast_pow_mul_tateYpair (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - a))
    (hpow : a ^ l = 1) (hsum : ∑ k ∈ range l, a ^ k = 0) (w q : R) (hq : q ∈ I) :
    (l : R) ^ 3 * tateYpair a w q hq = tateYpairDF l a w q hq := by
  rw [tateYpair, tateYpairDF, ← natCast_pow_mul_tateYterm hlu hu hpow hsum]
  ring

end Bridge

/-! ## ★出典の紐付け(`.src`) -/

def tateXpairDF.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(座標 X の分母なし版——l² を掛けた形)",
    sectionId := "genell-def-3-3" }

def tateYpairDF.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(座標 Y の分母なし版——l³ を掛けた形)",
    sectionId := "genell-def-3-3" }

def tateXpairDF.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_tateXterm(頭項の橋、第 1102、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.natCast_pow_mul_tateXterm") 1,
    .implicitStep
      ("★★**2026-09-01（第 1118）**——節点 3（`∑ veluV2 = ∑ DY` の `E` 版）の土台。" ++
       "☆`veluV2 = 3x² + 2a₂x + a₄ − a₁y` は次数が混ざるので `l⁶` を掛ければ" ++
       "`3l²(l²x)² + 2a₂l⁴(l²x) + a₄l⁶ − a₁l³(l³y)` と書ける。") 8 ]

end ABC3.Found.GaloisRep
