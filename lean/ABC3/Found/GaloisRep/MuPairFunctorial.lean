/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuPairDenomFree
import ABC3.Found.GaloisRep.TateFunctorial
import ABC3.Meta.Claim

/-!
# 第 1125 ブロック —— **分母なし版は環準同型と可換する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

第 221（`TateFunctorial.lean`）は `map_tateXpair`・`map_tateYpair` を与えるが、
どちらも **`IsUnit (1 − a)` を要求する**——`Ring.inverse` が単元の上でしか
環準同型と可換しないからである。

★分母なし版（`tateXpairDF` など）は頭項が **`Ring.inverse` を含まない多項式**なので、
`a` 側の単元性を**要らない**。要るのは `w` 側だけで、そちらは `w ∈ I` から自動である。

| 部品 | 単元性 |
|---|---|
| `tateXtermE l a`（頭） | ★**不要**（第 1111 の `map_tateXtermE`） |
| `tateXtail a q hq`（尾） | ★**不要**（各項の引数が `I` に入る） |
| `tateXterm w`・尾（`w` 側） | ☆`w ∈ I` から `IsUnit (1 − w)` が自動 |

## ★★★★★★★★これで何が繋がるか

節点 3（`∑ veluV2 = ∑ DY` の `E` 版）の降下は

    `A₀ = PowerSeries ℤ[ζ_l]`  →（`PowerSeries.map`）→  `A₁ = PowerSeries ℚ(ζ_l)`

で行う。`A₁` では `l` が可逆だから `1 − ζ^i` も可逆で、**既存の `hu` つきの補題が使える**。
★本ブロックの `map_*DF` が「DF 形は `A₀ → A₁` を通り抜ける」を与えるので、
`PowerSeries.map` の単射性で `A₀` に降ろせる。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Meta

variable {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {J : Ideal R'}

/-! ## ★★★★`D` 項と `D` 尾の移送（第 221 に無かった分） -/

/-- ★★**`Df(t)` は単元の上で環準同型と可換**。 -/
theorem map_tateDXterm (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateDXterm t) = tateDXterm (φ t) := by
  rw [tateDXterm, tateDXterm, map_mul, map_mul, map_pow, map_ring_inverse φ ht,
    map_sub, map_one, map_add, map_one]

/-- ★★**`Dg(t)` は単元の上で環準同型と可換**。 -/
theorem map_tateDYterm (φ : R →+* R') {t : R} (ht : IsUnit (1 - t)) :
    φ (tateDYterm t) = tateDYterm (φ t) := by
  rw [tateDYterm, tateDYterm, map_mul, map_mul, map_pow, map_pow,
    map_ring_inverse φ ht, map_sub, map_one, map_add]
  simp [map_ofNat]

/-- ★★★★**`DX` の尾は環準同型で移る**。 -/
theorem map_tateDXtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (u q : R) (hq : q ∈ I) :
    φ (tateDXtail u q hq)
      = tateDXtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDXtail, map_adicSum φ hφ, tateDXtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateDXterm φ (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)),
    map_mul, map_pow]

/-- ★★★★**`DY` の尾は環準同型で移る**。 -/
theorem map_tateDYtail [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
    (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)
    (u q : R) (hq : q ∈ I) :
    φ (tateDYtail u q hq)
      = tateDYtail (φ u) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDYtail, map_adicSum φ hφ, tateDYtail]
  refine adicSum_congr _ _ fun n => ?_
  rw [map_tateDYterm φ (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)),
    map_mul, map_pow]

/-! ## ★★★★★★★★★★分母なし版は `a` 側の単元性を要らない -/

section DF

variable [IsAdicComplete I R] [IsAdicComplete J R'] (φ : R →+* R')
  (hφ : ∀ (n : ℕ) (x : R), x ∈ I ^ n → φ x ∈ J ^ n)

include hφ

/-- ★★★★★★★★**`tateXpairDF` は環準同型で移る**——★`a` 側の単元性は要らない。 -/
theorem map_tateXpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) (hw : IsUnit (1 - w)) :
    φ (tateXpairDF l a w q hq)
      = tateXpairDF l (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateXpairDF, tateXpairDF, map_add, map_tateXtermE, map_mul, map_pow, map_natCast,
    map_sub, map_add, map_add, map_tateXtail φ hφ a q hq, map_tateXterm φ hw,
    map_tateXtail φ hφ w q hq, map_mul, map_evalAdic φ hφ _ q hq, map_ofNat]

/-- ★★★★★★★★**`tateYpairDF` は環準同型で移る**——★`a` 側の単元性は要らない。 -/
theorem map_tateYpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) (hw : IsUnit (1 - w)) :
    φ (tateYpairDF l a w q hq)
      = tateYpairDF l (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateYpairDF, tateYpairDF, map_add, map_tateYtermE, map_mul, map_pow, map_natCast,
    map_add, map_sub, map_sub, map_add, map_add, map_tateYtail φ hφ a q hq,
    map_tateXterm φ hw, map_tateXtail φ hφ w q hq, map_tateYterm φ hw,
    map_tateYtail φ hφ w q hq, map_evalAdic φ hφ _ q hq]

/-- ★★★★★★★★**`tateDXpairDF` は環準同型で移る**——★`a` 側の単元性は要らない。 -/
theorem map_tateDXpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) (hw : IsUnit (1 - w)) :
    φ (tateDXpairDF l a w q hq)
      = tateDXpairDF l (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDXpairDF, tateDXpairDF, map_add, map_tateDXtermE, map_mul, map_pow, map_natCast,
    map_sub, map_tateDXtail φ hφ a q hq, map_add, map_tateDXterm φ hw,
    map_tateDXtail φ hφ w q hq]

/-- ★★★★★★★★**`tateDYpairDF` は環準同型で移る**——★`a` 側の単元性は要らない。 -/
theorem map_tateDYpairDF (l : ℕ) (a w q : R) (hq : q ∈ I) (hw : IsUnit (1 - w)) :
    φ (tateDYpairDF l a w q hq)
      = tateDYpairDF l (φ a) (φ w) (φ q) (by simpa using hφ 1 q (by simpa using hq)) := by
  rw [tateDYpairDF, tateDYpairDF, map_add, map_tateDYtermE, map_mul, map_pow, map_natCast,
    map_add, map_add, map_tateDYtail φ hφ a q hq, map_add, map_add,
    map_tateDXterm φ hw, map_tateDXtail φ hφ w q hq, map_tateDYterm φ hw,
    map_tateDYtail φ hφ w q hq]

end DF

/-- ★★★★★★**`veluV2DF` は環準同型で移る**——★無条件（多項式だけ）。 -/
theorem map_veluV2DF (φ : R →+* R') (l : ℕ) (W : WeierstrassCurve R) (X Y : R) :
    φ (veluV2DF l W X Y) = veluV2DF l (W.map φ) (φ X) (φ Y) := by
  rw [veluV2DF, veluV2DF, map_sub, map_add, map_add, map_mul, map_mul, map_mul,
    map_mul, map_mul, map_mul, map_pow, map_pow, map_pow, map_pow, map_natCast]
  simp [map_ofNat, WeierstrassCurve.map]

/-! ## ★出典の紐付け(`.src`) -/

def map_tateXpairDF.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(座標 X の分母なし版は環準同型で移る——a 側の単元性は要らない)",
    sectionId := "genell-def-3-3" }

def map_tateDYpairDF.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(DY の分母なし版は環準同型で移る——a 側の単元性は要らない)",
    sectionId := "genell-def-3-3" }

def map_tateXpairDF.needs : List ProofObligation :=
  [ .citation "[ABC3]" "map_tateXtermE(頭項の E 版は無条件に移る、第 1125、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_tateXtermE") 1,
    .citation "[ABC3]" "map_tateXtail(尾は移る、第 221、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.map_tateXtail") 1,
    .implicitStep
      ("★★**2026-09-01（第 1125）**——`map_tateXpair`（第 221）は " ++
       "`IsUnit (1 − a)` を要求するが、DF 形は要求しない。" ++
       "☆頭項が `Ring.inverse` を含まない多項式だからである。" ++
       "★これが `A₀ = PowerSeries ℤ[ζ_l] → A₁ = PowerSeries ℚ(ζ_l)` の" ++
       "降下を可能にする（`A₁` では `hu` が成り立つ）。") 4 ]

end ABC3.Found.GaloisRep
