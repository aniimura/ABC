/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuHeadDenomFree
import ABC3.Found.GaloisRep.TatePair
import ABC3.Found.GenEll.Velu
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

/-! ## ★★★★★★★★★★第 1119 —— Vélu の `v` の分母なし版

☆`veluV2 W x y = 3x² + 2a₂x + a₄ − a₁y` は **`x` と `y` で次数が混ざる**。
★`x` の頭の分母は 2 次、`y` は 3 次なので、`l⁶` を掛ければ両方払える:

    `l⁶·veluV2 = 3l²(l²x)² + 2a₂l⁴(l²x) + a₄l⁶ − a₁l³(l³y)`

☆右辺は `l²x`・`l³y` だけで書けているので、
第 1118 の `tateXpairDF`・`tateYpairDF` を入れれば**分母が消える**。 -/

/-- ★★★★**Vélu の `v` の分母なし版**——`X`・`Y` には `l²x`・`l³y` を入れる。 -/
noncomputable def veluV2DF (l : ℕ) (W : WeierstrassCurve R) (X Y : R) : R :=
  3 * (l : R) ^ 2 * X ^ 2 + 2 * W.a₂ * (l : R) ^ 4 * X + W.a₄ * (l : R) ^ 6
    - W.a₁ * (l : R) ^ 3 * Y

/-- ★★★★★★★★**無条件の代数恒等式**——`IsUnit` を一切使わない。 -/
theorem natCast_pow_mul_veluV2 (l : ℕ) (W : WeierstrassCurve R) (x y : R) :
    (l : R) ^ 6 * ABC3.Found.GenEll.veluV2 W x y
      = veluV2DF l W ((l : R) ^ 2 * x) ((l : R) ^ 3 * y) := by
  show (l : R) ^ 6 * ABC3.Found.GenEll.veluGx W x y = _
  rw [ABC3.Found.GenEll.veluGx, veluV2DF]
  ring

/-- ★★★★★★★★★★**Tate の座標での形**（`hu` の下で）。

★右辺は `Ring.inverse` を含まないので `p ∣ l` でも意味を持つ。 -/
theorem natCast_pow_mul_veluV2_tate [IsDomain R] {l : ℕ} {a : R}
    (hlu : IsUnit ((l : R))) (hu : IsUnit (1 - a))
    (hpow : a ^ l = 1) (hsum : ∑ k ∈ range l, a ^ k = 0)
    (W : WeierstrassCurve R) (w q : R) (hq : q ∈ I) :
    (l : R) ^ 6 * ABC3.Found.GenEll.veluV2 W (tateXpair a w q hq) (tateYpair a w q hq)
      = veluV2DF l W (tateXpairDF l a w q hq) (tateYpairDF l a w q hq) := by
  rw [natCast_pow_mul_veluV2, natCast_pow_mul_tateXpair hlu hu hpow hsum,
    natCast_pow_mul_tateYpair hlu hu hpow hsum]

def veluV2DF.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(VÃ©lu ã® v ã®åæ¯ãªãçââlâ¶ ãæãã¦ xã»y ã®æ¬¡æ°å·®ãå¸åãã)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_veluV2.needs : List ProofObligation :=
  [ .implicitStep
      ("★★**2026-09-01（第 1119）**——`l⁶·veluV2` を `l²x`・`l³y` だけで書いた。" ++
       "☆これは**無条件の代数恒等式**で、`ring` で閉じる。" ++
       "★第 1118 の `tateXpairDF`・`tateYpairDF` を入れれば分母が消える。") 5 ]

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
