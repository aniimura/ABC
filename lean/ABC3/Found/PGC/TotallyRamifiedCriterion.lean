import ABC3.Found.PGC.HarithPAdicSupply
import ABC3.Found.PGC.AdjoinIntegers

/-!
# [pGC] ★★★★`hvalK` は「値群の生成元 = 拡大次数」と同値で、
`adjoinIntegers` は `integerSubring` そのものだった

残っていた 2 本は `hnK`（次数）と `hvalK`（全分岐）である。
本ファイルはそのうち `hvalK` を ★**1 つの数の等式**に直し、
さらに★★**体の側の語彙への橋が実は目の前にあった**ことを測る。

## ★先に測ったこと（持ち場の証拠 1）

`TotallyRamifiedValueGroup.lean:430
 exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_totallyRamified` は
★**すでに目を通した**: これは跡びの列 `u` と `hu0`/`hult`/`hudvd`/`hbnd`/`hjump`/
`hbreak`/`hti` をまだ要求する**古い形**であり、★**`hnK` と `hvalK` はそこでも仮説**である。
⇒ ★本波の仕事は縮まなかった。木の**どの形でも** `hnK`/`hvalK` は仮説である。

## ★★出したもの 1 —— `hvalK` は数の等式

`c` を「`‖K^×‖` の生成元」（`TotallyRamified.exists_valSub_gen`、
`TotallyRamifiedValueGroup.lean:252` が与える）とすると、

    hvalK  ⟺  n ∣ c  ⟺  c = n　（`n = [M:K]`）

★後半は `c ≤ n`（`valIndex_le_finrank`）と `1 ≤ c` から。
`c ≤ n` の中身は `π^0, …, π^{c−1}` が `K` 上 1 次独立であること
（`TotallyRamified.card_le_finrank_of_ne_mod`）。
⇒ ★★**`hvalK` は「ノルムの言葉での完全分岐 `e = n`」そのもの**である。

## ★★★出したもの 2 —— #69 の境界はここでは越えられる（訂正）

`lean-idioms.md` #69 は「`adjoinField` / `adjoinIntegers` の境界は越えられない
（212 秒 timeout）」と言う。★**これは経路に固有である**。

    theorem adjoinIntegers_eq_integerSubring (K : PAdicLocalField p) (x : K.closure) :
        adjoinIntegers K x
          = IntegerNorm.integerSubring (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
      Subring.ext (fun _ => Iff.rfl)

★★**8.9 秒で通る。** 両方とも `{y | ‖y‖ ≤ 1}` を `Subring.mk` したものだからである
（`AdjoinIntegers.lean:70-91` と `IntegerSubringNorm.lean:100-113` を見比べると完全に同じ形）。
★#69 が止まったのは `Valued` / `isCompact_closedBall` を経由したときであり
（`AdjoinIntegers.lean:20-33` のモジュール docstring にそう書いてある）、
**部分環の同一性**ではない。
★本日 3 度目の「罠は層に固有」（#332 / #337 に続く）。

⇒ ★★★**`IntegerNorm.*` の道具一式（DVR / `maximalIdeal = (π)` /
`mem_lowerRamificationGroup_iff_norm` / 剰余体の標数 / Noether）は
`adjoinIntegers K x` にそのまま移せる。**

## ★★どこで止まったか（`file:line`）

残るのは ★**`ramificationIndex K x = c`**（体の側の `e` とノルムの側の `c` の一致）である。
* 定義は `UnramifiedExtension.lean:425 ramificationIndex` で
  `(IsLocalRing.maximalIdeal 𝒪[K.carrier]).ramificationIdx
   (IsLocalRing.maximalIdeal (adjoinIntegers K x))`。
* `IsTotallyRamifiedAdjoin K x` は `TotallyRamified.lean:54` で
  `inertiaDegree K x = 1`、`ramificationIndex_mul_inertiaDegree`（`UnramifiedExtension.lean:444`）
  で `e·f = [K(x):K]` なので ★`IsTotallyRamifiedAdjoin ⇒ e = [K(x):K]` は出る。
* ⇒ ★★**残る 1 ノードは `ramificationIdx = c`**、すなわち
  「`𝒪_K` の極大イデアルが `(π)^c` にちょうど入る」をノルムで言うことである。
  ★本波はそこに降りていない。

## ★前回の①③の現在（持ち場の依頼）

* ① **不分岐側（`f = p`）** —— ★**依然に生きている。本波も触っていない。**
  ちなみに `WildDepthFieldDescent.lean:141 axWildDescent_prime` の docstring は
  ★自分で「`∏_{k ∈ Icc 1 n} p = p^n` は非有界なので `AxLemma` は出ない。
  残っているのは**分岐による絞り込み**だけである」と書いており、
  ★これは本連鎖の出口（`axDecay`）が埋めようとしている穴と**一致する**。
* ③ **`IntermediateField` の層** —— ★★**本ファイルで 1 枚剥がれた**。
  部分環の同一性は越えられる。残るのはイデアルの側（`ramificationIdx`）である。

## 逸脱の記録

1. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace TotRamCriterion

open TotallyRamified

section Criterion

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★**値群の生成元 `c` は `[M:K]` 以下** —— ノルムだけの `e ≤ n`。

`π^0, …, π^{c−1}` は `‖K^×‖` の相異なる剰余類にあるので `K` 上 1 次独立
（`TotallyRamified.card_le_finrank_of_ne_mod`）。 -/
theorem valIndex_le_finrank [FiniteDimensional K M] {π : M} {c : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hc : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((c : ℤ) * m)) :
    c ≤ Module.finrank K M := by
  classical
  have hcard := card_le_finrank_of_ne_mod (K := K) (M := M) hπ0 (ne_of_lt hπ1) hc
    (ι := Fin c) (fun l => π ^ (l : ℕ)) (fun l => ((l : ℕ) : ℤ))
    (fun l => by rw [norm_pow, zpow_natCast])
    (fun i j hij => TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt
      (fun h => hij (Fin.ext h)))
  simpa using hcard

/-- ★★★**`hvalK` は `n ∣ c` と同値** —— `c` は `‖K^×‖` の生成元。 -/
theorem valK_iff_dvd {π : M} {c n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m)
    (hcgen : ∀ m : ℤ, (∃ d : K, d ≠ 0 ∧ ‖algebraMap K M d‖ = ‖π‖ ^ m) ↔ (c : ℤ) ∣ m) :
    (∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) ↔ (n : ℤ) ∣ (c : ℤ) := by
  constructor
  · intro hvalK
    obtain ⟨d, hd0, hd⟩ := (hcgen ((c : ℕ) : ℤ)).mpr dvd_rfl
    obtain ⟨m, hm⟩ := hvalK d hd0
    rw [hd] at hm
    have hexp : ((c : ℕ) : ℤ) = (n : ℤ) * m := (zpow_right_inj₀ hπ0 (ne_of_lt hπ1)).mp hm
    exact ⟨m, hexp⟩
  · rintro ⟨t, ht⟩ a ha
    have ha0 : algebraMap K M a ≠ 0 := (map_ne_zero (algebraMap K M)).mpr ha
    obtain ⟨m, hm⟩ := hval (algebraMap K M a) ha0
    obtain ⟨r, hr⟩ := (hcgen m).mp ⟨a, ha, hm⟩
    exact ⟨t * r, by rw [hm, hr, ht]; ring_nf⟩

/-- ★★★★**`hvalK` ⟺「値群の生成元が拡大次数に等しい」** ——
すなわち★**ノルムの言葉での「完全分岐」**。 -/
theorem valK_iff_valIndex_eq [FiniteDimensional K M] {π : M} {c n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n) (hc0 : 0 < c)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m)
    (hcgen : ∀ m : ℤ, (∃ d : K, d ≠ 0 ∧ ‖algebraMap K M d‖ = ‖π‖ ^ m) ↔ (c : ℤ) ∣ m) :
    (∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) ↔ c = n := by
  have hcle : c ≤ n := by
    rw [← hn]
    refine valIndex_le_finrank hπ0 hπ1 ?_
    intro a ha
    have ha0 : algebraMap K M a ≠ 0 := (map_ne_zero (algebraMap K M)).mpr ha
    obtain ⟨m, hm⟩ := hval (algebraMap K M a) ha0
    obtain ⟨r, hr⟩ := (hcgen m).mp ⟨a, ha, hm⟩
    exact ⟨r, by rw [hm, hr]⟩
  rw [valK_iff_dvd hπ0 hπ1 hval hcgen]
  constructor
  · intro hdvd
    have hdvdN : n ∣ c := by exact_mod_cast hdvd
    have := Nat.le_of_dvd hc0 hdvdN
    omega
  · intro hce
    rw [hce]

end Criterion

/-! ## §2 ★★測定 —— `adjoinIntegers` と `integerSubring` は同じものか -/

section Bridge

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**測定**: `adjoinIntegers K x`（`AdjoinIntegers.lean:70`）と
`IntegerNorm.integerSubring`（`IntegerSubringNorm.lean:100`）は
★**同じ部分環**である（どちらも `{y | ‖y‖ ≤ 1}`）。 -/
theorem adjoinIntegers_eq_integerSubring (K : PAdicLocalField p) (x : K.closure) :
    adjoinIntegers K x
      = IntegerNorm.integerSubring
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  Subring.ext (fun _ => Iff.rfl)

end Bridge


/-! ## `.src` と 公理 -/

def valK_iff_valIndex_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def adjoinIntegers_eq_integerSubring.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms valIndex_le_finrank
#print axioms valK_iff_dvd
#print axioms valK_iff_valIndex_eq
#print axioms adjoinIntegers_eq_integerSubring

end TotRamCriterion

end ABC3.Found.PGC
