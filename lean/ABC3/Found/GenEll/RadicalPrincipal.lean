import ABC3.Meta.Claim
import Mathlib.RingTheory.Radical.Basic
import Mathlib.RingTheory.Ideal.Operations

/-!
# [GenEll] Definition 1.5, (ii) の可換環論的核(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

## ★原文 (ii) の中身は何か

> 正規ネーター scheme `Z` の**正則部分**に含まれる有効 Cartier 因子 `E` について、
> `E_red` もまた有効 Cartier 因子である。

有効 Cartier = 局所的に**非零因子 1 つ**で生成、被約化 = **根基**である。
したがって局所的には

> **UFD において、主イデアルの根基はまた主イデアルである**

に帰着する(正則局所環は UFD——Auslander–Buchsbaum)。

## ★★mathlib が **TODO** と書いている場所である

`Mathlib/RingTheory/Radical/Basic.lean` の冒頭:

> ## TODO
> - Connect this notion with `Ideal.radical`. Particularly, for a principal ideal,
>   **`Ideal.radical (Ideal.span {a}) = Ideal.span {radical a}`**.

★本ファイルはそれを取る(`a ≠ 0` が要る——`a = 0` では左辺が `⊥`、右辺が `⊤` で破れる)。

## ★★これで (ii) はどこまで来たか

- **UFD の局所的な主張**: ★本ファイルで取れた
- **正則局所環は UFD**(Auslander–Buchsbaum): ★★**mathlib に無い**
  (`IsRegularLocalRing` は `RingTheory/RegularLocalRing/Defs.lean` 1 ファイルだけ、
  2026-08-17 実測)
- scheme への大域化: 未着手

★★**したがって `Definition 1.5, (ii)` は完了ではない。**
本ファイルが取ったのは**可換環論の側だけ**であり、`.src` も条つきである。
★**足りないものが 1 本の定理に絞れた**ことが本ファイルの成果である。

★なお `Spec(𝓞_F)`(原文が実際に使う場合)の局所環は DVR であり、
DVR は PID ゆえ UFD なので、**その場合には Auslander–Buchsbaum を経由せず本結果が使える**。
-/

namespace ABC3.Found.GenEll

open UniqueFactorizationMonoid

section RadicalPrincipal

variable {R : Type*} [CommRing R] [IsDomain R] [NormalizationMonoid R]
  [UniqueFactorizationMonoid R]

/-- ★`a` は `radical a` の**十分高いべき**を割る。

★指数は `normalizedFactors a` の**要素数**で足りる——
どの素因子の重複度もその数を超えないからである。
★★「最大重複度」を取り出す必要はなかった。 -/
theorem self_dvd_radical_pow {a : R} (ha : a ≠ 0) :
    a ∣ radical a ^ Multiset.card (normalizedFactors a) := by
  classical
  set n := Multiset.card (normalizedFactors a) with hn
  have hle : normalizedFactors a ≤ n • (normalizedFactors a).dedup := by
    rw [Multiset.le_iff_count]
    intro p
    rw [Multiset.count_nsmul]
    by_cases hp : p ∈ normalizedFactors a
    · rw [Multiset.count_dedup, if_pos hp, mul_one, hn]
      exact Multiset.count_le_card p _
    · rw [Multiset.count_eq_zero.2 hp]
      exact Nat.zero_le _
  have hprod : ((normalizedFactors a).prod) ∣ (n • (normalizedFactors a).dedup).prod :=
    Multiset.prod_dvd_prod_of_le hle
  rw [Multiset.prod_nsmul] at hprod
  have hrad : radical a = (normalizedFactors a).dedup.prod := by
    rw [radical, ← Finset.prod_val, primeFactors, Multiset.toFinset_val]
  rw [hrad]
  exact ((prod_normalizedFactors ha).dvd_iff_dvd_left).mp hprod

/-- ★★**mathlib の TODO** —— UFD において主イデアルの根基は主イデアルである。

> **`Ideal.radical (Ideal.span {a}) = Ideal.span {radical a}`**(`a ≠ 0`)

★`a ≠ 0` が要る: `a = 0` なら左辺は `⊥`(整域の冪零根基)、
右辺は `radical 0 = 1` より `⊤` で、**破れる**。
★★mathlib の TODO 文はこの条件を書いていない——**負の対照が 1 つ得られた**。 -/
theorem ideal_radical_span_singleton {a : R} (ha : a ≠ 0) :
    (Ideal.span {a} : Ideal R).radical = Ideal.span {radical a} := by
  apply le_antisymm
  · intro x hx
    obtain ⟨n, hn⟩ := hx
    rw [Ideal.mem_span_singleton] at hn ⊢
    rcases Nat.eq_zero_or_pos n with rfl | hnpos
    · -- `a ∣ 1` なら `a` は単元で `radical a = 1`
      rw [pow_zero] at hn
      rw [radical_of_isUnit (isUnit_of_dvd_one hn)]
      exact one_dvd x
    rcases eq_or_ne x 0 with rfl | hx0
    · exact dvd_zero _
    have hxn : x ^ n ≠ 0 := pow_ne_zero _ hx0
    have h1 : radical a ∣ radical (x ^ n) := radical_dvd_radical hn hxn
    rw [radical_pow x hnpos.ne'] at h1
    exact h1.trans radical_dvd_self
  · rw [Ideal.span_le, Set.singleton_subset_iff]
    exact ⟨Multiset.card (normalizedFactors a),
      Ideal.mem_span_singleton.2 (self_dvd_radical_pow ha)⟩

/-- ★★**`Definition 1.5, (ii)` の局所的な形** ——
UFD の主イデアルの根基は、**非零因子 1 つ**で生成される。

★有効 Cartier 因子の定義(局所的に非零因子 1 つで生成)がそのまま保たれる、
というのが原文 (ii) の主張であり、その局所版が本定理である。 -/
theorem exists_nonZeroDivisor_generator_radical {a : R} (ha : a ≠ 0) :
    ∃ b ∈ nonZeroDivisors R, (Ideal.span {a} : Ideal R).radical = Ideal.span {b} :=
  ⟨radical a, mem_nonZeroDivisors_of_ne_zero radical_ne_zero,
    ideal_radical_span_singleton ha⟩

/-- ★**負の対照** —— `a = 0` では `ideal_radical_span_singleton` が破れる。

`(span {0}).radical = ⊥` だが `span {radical 0} = span {1} = ⊤` である。
★★これが破れないなら `a ≠ 0` の仮定は効いていない。 -/
theorem ideal_radical_span_zero_ne :
    (Ideal.span {(0 : R)} : Ideal R).radical ≠ Ideal.span {radical (0 : R)} := by
  rw [radical_zero, Ideal.span_singleton_one, Ideal.span_singleton_eq_bot.2 rfl]
  intro hcon
  have h1 : (1 : R) ∈ (⊥ : Ideal R).radical := hcon ▸ Submodule.mem_top
  obtain ⟨n, hn⟩ := h1
  rw [one_pow] at hn
  exact one_ne_zero (Ideal.mem_bot.mp hn)

end RadicalPrincipal

/-! ## ★出典の紐付け(`.src`) -/

def ideal_radical_span_singleton.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(可換環論の核のみ——正則性は仮定していない)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
