import ABC3.Found.GenEll.BaseChange
import ABC3.Found.GenEll.DegMul
import Mathlib.NumberTheory.RamificationInertia.Ramification
import Mathlib.RingTheory.UniqueFactorizationDomain.Multiplicity

/-!
# [GenEll] Definition 1.2, (i) の足場 —— **イデアルの拡大は底変換である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★`Definition 1.2` の残り 4 段のうち、段 1

`PullbackBase.lean` で `x_K^* D = (x_F^* D)·𝓞_K` を取った。
★残るのはそれを **`ADiv` の言葉に翻訳する**ことである:

    `idealADiv K (J·𝓞_K) = baseChange F K (idealADiv F J)` の有限素点側

すなわち **`ord_w(J·𝓞_K) = e(w|v)·ord_v(J)`**。

## ★★★私の測り違いの訂正

一度「一般のイデアルの拡大の重複度公式は mathlib に無い」と書いた。
★★★**誤りだった。ある:**

    `NumberTheory/RamificationInertia/Ramification.lean`
    `emultiplicity_map_eq_ramificationIdx_mul`
      (h : I ≠ ⊥) … : emultiplicity w (I.map (algebraMap R S))
        = v.ramificationIdx w * emultiplicity v I

★`ramificationIdx` の語で検索し、**`emultiplicity` の語で検索しなかった**のが誤りである。

## ★★橋が 1 本要った

`idealADiv` は `Associates.count _ (Associates.mk J).factors` を使い、
mathlib の公式は `emultiplicity` で述べられている。
★**両者を直接繋ぐ補題は無い**(実測)。

★★**どちらも整除で特徴づけられる**ので橋は架かる:
`Associates.prime_pow_dvd_iff_le`(`p^n ∣ a ↔ n ≤ count`)と
`emultiplicity_eq_coe`(`p^n ∣ b ∧ ¬p^(n+1) ∣ b`)。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

/-! ## ★★橋 —— `Associates.count` と `emultiplicity` -/

/-- ★★**`Associates.count` は `emultiplicity` に等しい**。

★mathlib に無い(2026-08-17 実測)。
★★機構は「どちらも整除で特徴づけられる」ことだけ:
`Associates.prime_pow_dvd_iff_le` と `emultiplicity_eq_coe`。 -/
theorem emultiplicity_eq_associatesCount {α : Type*} [CommMonoidWithZero α]
    [UniqueFactorizationMonoid α] [DecidableEq (Associates α)]
    [DecidablePred (Irreducible : Associates α → Prop)] {p a : α} (hp : Prime p) (ha : a ≠ 0) :
    emultiplicity p a
      = Associates.count (Associates.mk p) (Associates.mk a).factors := by
  have hmk : Associates.mk a ≠ 0 := by simpa [Associates.mk_eq_zero] using ha
  have hirr : Irreducible (Associates.mk p) := by
    rw [Associates.irreducible_mk]; exact hp.irreducible
  rw [emultiplicity_eq_coe]
  refine ⟨?_, ?_⟩
  · rw [← Associates.mk_dvd_mk, Associates.mk_pow]
    exact (Associates.prime_pow_dvd_iff_le hmk hirr).2 le_rfl
  · rw [← Associates.mk_dvd_mk, Associates.mk_pow]
    intro hc
    have := (Associates.prime_pow_dvd_iff_le hmk hirr).1 hc
    omega

/-! ## ★★★イデアルの拡大の重複度 -/

section BaseChange

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

open scoped Classical in
/-- ★★★**イデアルを拡大すると重複度は `e(w|v)` 倍になる**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `ord_w(J·𝓞_K) = e(w|v)·ord_v(J)`

★★これが `Definition 1.2, (i)` の**段 1** である
——`idealADiv` を底変換と繋ぐ。

★機構は mathlib の `emultiplicity_map_eq_ramificationIdx_mul` に
本ファイルの橋(`emultiplicity_eq_associatesCount`)を当てるだけ。
★`LiesOver` インスタンスは `hosComap_liesOver`(`FinitePlaceRel.lean`)が供給する。 -/
theorem idealADiv_fin_map
    (J : Ideal (𝓞 F)) (hJ : J ≠ 0) (W : FinitePlace K)
    [W.asIdeal.LiesOver (hosComap F K W).asIdeal] :
    (idealADiv K (J.map (algebraMap (𝓞 F) (𝓞 K)))).fin W
      = (ramIdxOver F K W : ℤ) * (idealADiv F J).fin (hosComap F K W) := by
  set V := hosComap F K W with hV
  have hJK : J.map (algebraMap (𝓞 F) (𝓞 K)) ≠ 0 := by
    rw [Ideal.zero_eq_bot] at hJ ⊢
    exact Ideal.map_ne_bot_of_ne_bot hJ
  -- ★mathlib の公式(一般のイデアルについて)
  have hmul : emultiplicity W.asIdeal (J.map (algebraMap (𝓞 F) (𝓞 K)))
      = V.asIdeal.ramificationIdx W.asIdeal * emultiplicity V.asIdeal J := by
    refine Ideal.IsDedekindDomain.emultiplicity_map_eq_ramificationIdx_mul ?_
      V.irreducible W.irreducible W.ne_bot
    rwa [← Ideal.zero_eq_bot]
  -- ★橋で `Associates.count` に移し、ℕ の等式にする
  rw [emultiplicity_eq_associatesCount
      (W.asIdeal.prime_of_isPrime W.ne_bot W.isPrime) hJK,
    emultiplicity_eq_associatesCount
      (V.asIdeal.prime_of_isPrime V.ne_bot V.isPrime) hJ] at hmul
  have hcount : (Associates.mk W.asIdeal).count
        (Associates.mk (J.map (algebraMap (𝓞 F) (𝓞 K)))).factors
      = V.asIdeal.ramificationIdx W.asIdeal
        * (Associates.mk V.asIdeal).count (Associates.mk J).factors := by
    exact_mod_cast hmul
  rw [idealADiv_fin_apply K _ hJK, idealADiv_fin_apply F J hJ, hcount, ramIdxOver, ← hV]
  push_cast
  ring

end BaseChange

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
アルキメデス側の `ADiv` への翻訳と `X(ℚ̄)` の型の構成が残っている。 -/

def idealADiv_fin_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(イデアルの拡大が底変換になること——有限素点側のみ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
