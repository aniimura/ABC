/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateSeries
import ABC3.Meta.Claim

/-!
# 第 1110 ブロック —— **捻った σ 級数**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

第 1109 で決めたとおり、Tate の**尾**の万有な冪級数は

    `G(u) = ∑_N (∑_{m ∣ N} c_m·u^m)·X^N`

の形をしている（`c_m` は頭項の母関数の係数：`m`・`C(m,2)`・`m²`・`m·C(m,2)`・`m³`・`m⁴`）。

☆在庫の `sigmaSeries k = ∑_N σ_k(N)·X^N`（`TateSeries.lean:73`）は **`u = 1` の場合**である。
★本ファイルはそれを `u` で捻った版を置く。

## ★★★★これが要る理由

`Lemma 3.5` の分母払い（第 1092-1109）は、尾を
**係数が `ℤ[ζ_l]` に入る万有な冪級数**として書き直すことを要求する。
`u = ζ^i` を代入した `sigmaSeriesTwisted k (ζ^i)` がちょうどそれである。
☆そのうえで第 1107 の `evalAdicMap_eq_of_map_eq` を当てれば `hu` が消える。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {S : Type} [CommRing S]

/-- ★★★★**捻った σ 級数** `∑_N (∑_{m ∣ N} m^k·u^m)·X^N`。

☆`u = 1` なら在庫の `sigmaSeries k` に一致する（`coeff_sigmaSeriesTwisted_one`）。 -/
noncomputable def sigmaSeriesTwisted (k : ℕ) (u : S) : PowerSeries S :=
  PowerSeries.mk (fun N => if N = 0 then 0 else ∑ m ∈ N.divisors, (m : S) ^ k * u ^ m)

@[simp] theorem coeff_sigmaSeriesTwisted (k : ℕ) (u : S) (N : ℕ) :
    PowerSeries.coeff N (sigmaSeriesTwisted k u)
      = if N = 0 then 0 else ∑ m ∈ N.divisors, (m : S) ^ k * u ^ m := by
  rw [sigmaSeriesTwisted, PowerSeries.coeff_mk]

/-- ★★★**`u = 1` なら通常の σ である**。 -/
theorem coeff_sigmaSeriesTwisted_one (k N : ℕ) :
    PowerSeries.coeff N (sigmaSeriesTwisted k (1 : S))
      = if N = 0 then 0 else ((ArithmeticFunction.sigma k N : ℕ) : S) := by
  rw [coeff_sigmaSeriesTwisted]
  by_cases hN : N = 0
  · simp [hN]
  · simp only [if_neg hN, one_pow, mul_one]
    rw [ArithmeticFunction.sigma_apply]
    push_cast
    exact Finset.sum_congr rfl (fun m _ => rfl)

/-- ★★零番目の係数は `0`。 -/
@[simp] theorem coeff_zero_sigmaSeriesTwisted (k : ℕ) (u : S) :
    PowerSeries.coeff 0 (sigmaSeriesTwisted k u) = 0 := by
  simp

/-- ★★一番目の係数は `u`。 -/
theorem coeff_one_sigmaSeriesTwisted (k : ℕ) (u : S) :
    PowerSeries.coeff 1 (sigmaSeriesTwisted k u) = u := by
  rw [coeff_sigmaSeriesTwisted]
  simp

/-- ★★★**係数環の準同型で写る**——降下と特殊化に要る自然性。 -/
theorem map_sigmaSeriesTwisted {T : Type} [CommRing T] (ψ : S →+* T) (k : ℕ) (u : S) :
    PowerSeries.map ψ (sigmaSeriesTwisted k u) = sigmaSeriesTwisted k (ψ u) := by
  refine PowerSeries.ext (fun N => ?_)
  rw [PowerSeries.coeff_map, coeff_sigmaSeriesTwisted, coeff_sigmaSeriesTwisted]
  by_cases hN : N = 0
  · simp [hN]
  · simp only [if_neg hN, map_sum, map_mul, map_pow]
    exact Finset.sum_congr rfl (fun m _ => by rw [map_natCast])

/-! ## ★出典の紐付け(`.src`) -/

def sigmaSeriesTwisted.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(捻った σ 級数——Tate の尾の万有な形)",
    sectionId := "genell-def-3-3" }

def map_sigmaSeriesTwisted.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(捻った σ 級数は係数環の準同型で写る)",
    sectionId := "genell-def-3-3" }

def sigmaSeriesTwisted.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★**2026-09-01（第 1110）**——第 1109 で決めた万有な尾の形である。" ++
       "☆残るのは 6 種それぞれについて " ++
       "`evalAdicMap φ (sigmaSeriesTwisted k (ζ^i)) q hq = tail (ζ^i) q hq` を示す段。" ++
       "★`map_sigmaSeriesTwisted` が第 1107 の降下に要る自然性を与える。") 15 ]

end ABC3.Found.GaloisRep
