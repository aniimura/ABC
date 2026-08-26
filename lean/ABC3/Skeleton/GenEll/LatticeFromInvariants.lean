import ABC3.Found.GenEll.LatticeInvariance

/-!
# スケルトン —— **`(g₂,g₃)` は束を決める**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜこの節点を立てるのか——`htFalt` を固定する最後の一歩

第 348 で一意化(任意の `E/ℂ` が格子曲線と同型)が閉じ、
第 349-350 で **`archInv = ‖Δ_lat‖·covol⁶` は束だけで決まる**ところまで来た。

★★しかし `htFalt` を**曲線の関数**として定義するには、周期束の**選び方に依らない**ことが要る:

> `C • W = latticeCurve P` と `C' • W = latticeCurve P'` ⟹ `archInv P = archInv P'`

★★★`D := C'C⁻¹` の `u` について `g₂(P') = u⁻⁴g₂(P)`・`g₃(P') = u⁻⁶g₃(P)`、
すなわち `P'` と `scalePair P u` は**同じ `(g₂,g₃)`** を持つ。
★★★★したがってこの節点(`(g₂,g₃)` が束を決める)がそのまま
**`archInv` の一意性**を与え、界面の欠陥 #6(`htFalt` が固定されない)を塞ぐ鍵になる。

## ★★★★★★2026-08-26 の測定——**道具は mathlib にそろっている**

| 段 | 状態 |
|---|---|
| `℘` が同じ ⟹ 束が同じ | ✅ `Found/GenEll/LatticeInvariance.lean` `lattice_eq_of_weierstrassP_eq` |
| 束が同じ ⟹ `covol`・`archInv` が同じ | ✅ 同 `covol_congr`・`archInv_congr` |
| **`(g₂,g₃)` が同じ ⟹ `℘` が同じ** | ★**本節点** |

### ★★★★★冪級数の在庫(mathlib)

`hasFPowerSeriesAt_weierstrassPExcept` は `℘[L−0]` の `0` での冪級数係数が
**`(i+1)·sumInvPow 0 (i+2)`** であることを与え、`sumInvPow_zero : sumInvPow 0 = G` なので
係数は **`a_i = (i+1)·G(i+2)`** である。
★`G_eq_zero_of_odd` により奇数の `G` は消える。
★★`derivWeierstrassP_sq` は `℘'² = 4℘³ − g₂℘ − g₃` を与える。

### ★★★★★★段取り(係数の漸化式)

`f := ℘[L−0]` と置くと `℘ = f + z⁻²`(`weierstrassPExcept_def` を `l₀ = 0` で)。
微分方程式を 1 回微分して `2℘'` で割ると `℘'' = 6℘² − g₂/2` になり、`f` で書き直すと

    z²·f'' = 6·z²·f² + 12·f − (g₂/2)·z²

という**両辺 `0` で解析的な恒等式**になる。★係数を比べると

    (i−4)(i+3)·a_i = 6·Σ_{j+k=i−2} a_j a_k − (g₂/2)·[i = 2]

★★`i = 2` は `a₂ = 3G₄ = g₂/20`、`i = 4` は**係数が消える**——
`a₄ = 5G₆ = g₃/28` が 2 つ目の自由母数だからである。
★★★`i ≥ 5` では `(i−4)(i+3) ≠ 0` なので **`a_i` は `g₂, g₃` から順に決まる**。
★★★★したがって `(g₂,g₃)` が等しければ `f` の冪級数が一致し、`0` の近傍で `℘_L = ℘_{L'}`。
★★★★★あとは**一致の定理**で `(Λ ∪ Λ')ᶜ`(連結開集合)へ延ばし、極の位置で束が一致する。

★見積もり **4-8 ブロック**。

## ★★別の道(取らなかった)

古典的には `j` の `ℍ/Γ` 上での**単射性**——**valence 公式**(幅角の原理)——で出す。
★**2026-08-26 実測: mathlib に valence 公式は無い**。
★★上の冪級数の道は mathlib の在庫(`hasFPowerSeriesAt_weierstrassPExcept`)に直結しており、
留数計算・輪郭積分を一切使わない。★★★こちらを取る。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-- ★★★★★★★**`(g₂,g₃)` は束を決める**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これが `archInv` の一意性、ひいては `ht^Falt` の固定に効く。 -/
theorem lattice_eq_of_g₂_g₃_eq {L L' : PeriodPair}
    (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) : L.lattice = L'.lattice := by
  sorry

/-- ★★★★★★**アルキメデス不変量は `(g₂,g₃)` で決まる**——本節点の下流。 -/
theorem archInv_eq_of_g₂_g₃_eq {L L' : PeriodPair}
    (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) : archInv L = archInv L' :=
  archInv_congr (lattice_eq_of_g₂_g₃_eq h₂ h₃)

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def archInv_eq_of_g₂_g₃_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def archInv_eq_of_g₂_g₃_eq.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★本節点の直下流であり、追加の依存は無い——`Found/GenEll/LatticeInvariance.lean` の `archInv_congr`(束が同じなら `archInv` も同じ)に `lattice_eq_of_g₂_g₃_eq` を渡すだけである" 17,
    .otherPaper "GenEll" "Proposition 3.4(Faltings 高さと無限遠因子)" 17 ]

def lattice_eq_of_g₂_g₃_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def lattice_eq_of_g₂_g₃_eq.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★**2026-08-26 の測定(道具はそろっている)**: mathlib の `hasFPowerSeriesAt_weierstrassPExcept` は `℘[L−0]` の `0` での冪級数係数が `(i+1)·sumInvPow 0 (i+2)` であることを与え、`sumInvPow_zero : sumInvPow 0 = G` なので係数は `a_i = (i+1)·G(i+2)` である。★`G_eq_zero_of_odd` により奇数の `G` は消える。★★`derivWeierstrassP_sq` は `℘'² = 4℘³ − g₂℘ − g₃` を与える" 17,
    .implicitStep
      "★★★★★★**段取り(係数の漸化式)**: `f := ℘[L−0]` と置くと `℘ = f + z⁻²` で、微分方程式を 1 回微分して `2℘'` で割ると `℘'' = 6℘² − g₂/2`、`f` で書き直すと `z²f'' = 6z²f² + 12f − (g₂/2)z²` という両辺 `0` で解析的な恒等式になる。★係数を比べると `(i−4)(i+3)·a_i = 6·Σ_{j+k=i−2} a_j a_k − (g₂/2)·[i=2]`。★★`i=4` で係数が消えるのは `a₄ = 5G₆ = g₃/28` が 2 つ目の自由母数だからであり、`i ≥ 5` では `(i−4)(i+3) ≠ 0` なので `a_i` は `g₂,g₃` から順に決まる。★★★見積もり 4-8 ブロック" 17,
    .implicitStep
      "★★★★★**一致の定理で延ばす段**: 冪級数の一致は `0` の近傍でしか `℘_L = ℘_{L'}` を与えない。★`℘_L` は `Λᶜ` で、`℘_{L'}` は `Λ'ᶜ` で解析的なので、共通の定義域 `(Λ ∪ Λ')ᶜ` は開である。★★これが**連結**であること(`ℂ` から可算な離散集合を除いても連結)が要る。★★★そのうえで一致の定理を使い、最後に `Found/GenEll/LatticeInvariance.lean` の `lattice_eq_of_weierstrassP_eq`(束は `℘` の不連続点の集合)で束の一致に落とす" 17,
    .implicitStep
      "★★★★**なぜ要るか(下流)**: `C • W = latticeCurve P` と `C' • W = latticeCurve P'` のとき `D := C'C⁻¹` の `u` について `g₂(P') = u⁻⁴g₂(P)`・`g₃(P') = u⁻⁶g₃(P)`、すなわち `P'` と `scalePair P u` は同じ `(g₂,g₃)` を持つ。★本節点があれば `archInv P' = archInv (scalePair P u) = archInv P`(第 349 `archInv_scalePair`)となり、**`archInv` が曲線の関数になる**。★★これが界面の欠陥 #6(`htFalt` が固定されず `prop_3_4` が恒等的に成り立つ)を塞ぐ鍵である" 17,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, VI.3(格子は g₂, g₃ で決まる)"
      (.absent "mathlib の `Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` を `Injective|lattice_eq_of|valence` で検索して 0 件(2026-08-26 実測)。`Mathlib` 全体を `valence` で検索しても modular forms の valence 公式は無い") 17,
    .otherPaper "GenEll" "Proposition 3.4(Faltings 高さと無限遠因子)" 17 ]

end ABC3.Skeleton.GenEll
