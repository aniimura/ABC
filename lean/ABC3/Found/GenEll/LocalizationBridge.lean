/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ConductorDescent
import Mathlib.RingTheory.Localization.Ideal

/-!
# [GenEll] Remark 1.5.1 —— **局所化の橋**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★これが `Remark 1.5.1` の最後の 1 本

`ConductorDescent.lean` は「`X_n = X ×_ℤ ℤ[1/n!]` の点」について導手の一致を出した。
★原文は「`X` の点 `x` と `Σ` の外の素点 `v`」について述べている。
★★両者を繋ぐのが**本ファイル**である:

    I・J ⊆ 𝒪_F が 𝒪_F[1/N] で一致すれば、
    `N` を割らない `v` で **`ord_v(I) = ord_v(J)`**

★★★これは Dedekind 環の局所化の補題であり、**mathlib の欠落ではない**
——`FinitePlace F = IsDedekindDomain.HeightOneSpectrum 𝒪_F` なので、
mathlib の付値の道具（`intValuation_le_pow_iff_mem` ほか）がそのまま効く。

## ★★★★★★機構は 4 段

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | `m ∉ v` なら `v` の整付値は 1 | `intValuation_lt_one_iff_dvd` |
| 2 | `m·x ∈ v^n` かつ `m ∉ v` なら `x ∈ v^n` | `intValuation_le_pow_iff_mem`（付値の乗法性） |
| 3 | 局所化で `≤` が降りる: `J·A ≤ I·A` かつ `I ≤ v^n` なら `J ≤ v^n` | `mem_map_algebraMap_iff` ＋ 段 2 |
| 4 | ゆえに重複度が一致する | `Associates.prime_pow_dvd_iff_le` ＋ `Ideal.dvd_iff_le` |

★導手は**根基**を取るが、`IsLocalization.map_radical`（根基と局所化は可換）で吸収される。
-/

namespace ABC3.Found.GenEll

open IsDedekindDomain AlgebraicGeometry CategoryTheory

/-! ## ★★★★段 1-2 —— 付値で `v^n` への所属を測る -/

variable {R : Type*} [CommRing R] [IsDedekindDomain R]

/-- ★**`m ∉ v` なら `v` の整付値は 1**。

★`intValuation_le_one` と `intValuation_lt_one_iff_dvd` で挟む。 -/
theorem intValuation_eq_one_of_notMem (v : HeightOneSpectrum R) {m : R}
    (hm : m ∉ v.asIdeal) : v.intValuation m = 1 := by
  refine le_antisymm (v.intValuation_le_one m) (le_of_not_gt fun h => hm ?_)
  exact (Ideal.dvd_span_singleton).1 ((v.intValuation_lt_one_iff_dvd m).1 h)

/-- ★★**`v` の外の元を掛けても `v^n` への所属は変わらない**。

★機構は付値の乗法性——`v(m·x) = v(m)·v(x) = v(x)`。 -/
theorem mem_pow_of_mul_mem (v : HeightOneSpectrum R) {m x : R}
    (hm : m ∉ v.asIdeal) {n : ℕ} (h : m * x ∈ v.asIdeal ^ n) : x ∈ v.asIdeal ^ n := by
  rw [← v.intValuation_le_pow_iff_mem] at h ⊢
  rwa [map_mul, intValuation_eq_one_of_notMem v hm, one_mul] at h

/-! ## ★★★★★★段 3 —— 局所化で `≤` が降りる -/

/-- ★★★★★★**局所化で `v^n` への包含が降りる**。

`J·A ≤ I·A` かつ `I ≤ v^n` で、`M` が `v` を避けるなら `J ≤ v^n`。

★機構: `x ∈ J` なら `algebraMap x ∈ I·A` なので、
ある `m, c ∈ M` と `y ∈ I` で `c·m·x = c·y ∈ I ≤ v^n`。
★★`c·m ∈ M` は `v` の外にあるから、段 2 で `x ∈ v^n` が出る。 -/
theorem le_pow_of_map_le_map (M : Submonoid R) (A : Type*) [CommRing A] [Algebra R A]
    [IsLocalization M A] (v : HeightOneSpectrum R) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {I J : Ideal R} {n : ℕ} (hI : I ≤ v.asIdeal ^ n)
    (h : J.map (algebraMap R A) ≤ I.map (algebraMap R A)) : J ≤ v.asIdeal ^ n := by
  intro x hx
  have hxm : algebraMap R A x ∈ I.map (algebraMap R A) := h (Ideal.mem_map_of_mem _ hx)
  obtain ⟨⟨⟨y, hy⟩, ⟨m, hm⟩⟩, hEq⟩ := (IsLocalization.mem_map_algebraMap_iff M A).1 hxm
  simp only [← map_mul] at hEq
  obtain ⟨⟨c, hc⟩, hc2⟩ := (IsLocalization.eq_iff_exists M A).1 hEq
  simp only at hc2
  refine mem_pow_of_mul_mem v (m := c * m) (hM _ (M.mul_mem hc hm)) ?_
  have hxy : c * m * x = c * y := by rw [← hc2]; ring
  rw [hxy]
  exact hI (Ideal.mul_mem_left _ _ hy)

/-! ## ★★★★★★★★段 4 —— 重複度が一致する -/

open Classical in
/-- ★★★★★★★★**局所化で一致するイデアルは、避けた素点で同じ重複度を持つ**。

★機構: `v^n ∣ I ↔ I ≤ v^n`（`Ideal.dvd_iff_le`）と
`n ≤ count v I ↔ v^n ∣ I`（`Associates.prime_pow_dvd_iff_le`）で
重複度を包含に翻訳し、段 3 を両向きに使う。 -/
theorem count_eq_of_map_eq (M : Submonoid R) (A : Type*) [CommRing A] [Algebra R A]
    [IsLocalization M A] (v : HeightOneSpectrum R) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {I J : Ideal R} (hI : I ≠ 0) (hJ : J ≠ 0)
    (h : I.map (algebraMap R A) = J.map (algebraMap R A)) :
    (Associates.mk v.asIdeal).count (Associates.mk I).factors
      = (Associates.mk v.asIdeal).count (Associates.mk J).factors := by
  have hIa : Associates.mk I ≠ 0 := by simpa [Associates.mk_eq_zero] using hI
  have hJa : Associates.mk J ≠ 0 := by simpa [Associates.mk_eq_zero] using hJ
  have hirr : Irreducible (Associates.mk v.asIdeal) := Associates.irreducible_mk.2 v.irreducible
  have key : ∀ n : ℕ, (n ≤ (Associates.mk v.asIdeal).count (Associates.mk I).factors)
      ↔ (n ≤ (Associates.mk v.asIdeal).count (Associates.mk J).factors) := by
    intro n
    rw [← Associates.prime_pow_dvd_iff_le hIa hirr, ← Associates.prime_pow_dvd_iff_le hJa hirr,
      ← Associates.mk_pow, Associates.mk_le_mk_iff_dvd, Associates.mk_le_mk_iff_dvd,
      Ideal.dvd_iff_le, Ideal.dvd_iff_le]
    exact ⟨fun hle => le_pow_of_map_le_map M A v hM hle (le_of_eq h.symm),
      fun hle => le_pow_of_map_le_map M A v hM hle (le_of_eq h)⟩
  exact le_antisymm ((key _).1 le_rfl) ((key _).2 le_rfl)

/-! ## ★★★★★★★★★★到達点 —— `hagree` そのもの -/

variable (F : Type) [Field F] [NumberField F]

open Classical in
/-- ★★★★★★★★★**局所化で一致するイデアルは、避けた素点で同じ導手係数を与える**。

★根基を取る段は `IsLocalization.map_radical`（根基と局所化は可換）で吸収される。 -/
theorem idealADiv_radical_fin_eq_of_map_eq
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A]
    [IsLocalization M A] (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {I J : Ideal (NumberField.RingOfIntegers F)} (hI : I ≠ 0) (hJ : J ≠ 0)
    (h : I.map (algebraMap _ A) = J.map (algebraMap _ A)) :
    (idealADiv F I.radical).fin v = (idealADiv F J.radical).fin v := by
  have hIr : I.radical ≠ 0 := fun hc =>
    hI (le_bot_iff.1 (le_trans Ideal.le_radical (le_of_eq hc)))
  have hJr : J.radical ≠ 0 := fun hc =>
    hJ (le_bot_iff.1 (le_trans Ideal.le_radical (le_of_eq hc)))
  have hr : I.radical.map (algebraMap _ A) = J.radical.map (algebraMap _ A) := by
    rw [IsLocalization.map_radical (M := M), IsLocalization.map_radical (M := M), h]
  have hcount := count_eq_of_map_eq M A v hM hIr hJr hr
  unfold idealADiv ADiv.fin
  rw [dif_neg hIr, dif_neg hJr]
  simp only [Finsupp.ofSupportFinite_coe]
  exact_mod_cast hcount

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— `Σ` の外で導手の係数が一致する**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★**これが `Skeleton/GenEll/Section1.lean` の `remark_1_5_1` が仮定として受けている
`hagree` そのものである**——「`Σ` の外で」は「`M` が `v` を避ける」に対応する。

★★★★★橋の全体:
`ConductorDescent.lean` が `X_n` の点で因子の一致を出し、
本定理が「`𝒪_F[1/N]` で一致するなら `N` を割らない `v` で係数が一致する」を出す。 -/
theorem conductorADiv_fin_eq_of_map_eq
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A] [IsLocalization M A]
    (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {X X' : Scheme} (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hI : pullbackIdeal F D xF ≠ 0) (hJ : pullbackIdeal F D' xF' ≠ 0)
    (h : (pullbackIdeal F D xF).map (algebraMap _ A)
       = (pullbackIdeal F D' xF').map (algebraMap _ A)) :
    (conductorADiv F D xF).fin v = (conductorADiv F D' xF').fin v :=
  idealADiv_radical_fin_eq_of_map_eq F M A v hM hI hJ h

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
「`ConductorDescent.lean` の因子の一致を `𝒪_F[1/N]` 上のイデアルの一致に翻訳する」
配線だけである（`pullbackIdeal` の局所化との両立）。 -/

def count_eq_of_map_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(局所化の橋——避けた素点では重複度が一致する)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_map_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(Σ の外で導手の係数が一致する——hagree の中身)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_map_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "HeightOneSpectrum.intValuation_le_pow_iff_mem(v^n への所属は付値で測れる)"
      (.inMathlib "IsDedekindDomain.HeightOneSpectrum.intValuation_le_pow_iff_mem") 9,
    .citation "[mathlib]" "IsLocalization.mem_map_algebraMap_iff(局所化したイデアルの元)"
      (.inMathlib "IsLocalization.mem_map_algebraMap_iff") 9,
    .citation "[mathlib]" "IsLocalization.map_radical(根基と局所化は可換)"
      (.inMathlib "IsLocalization.map_radical") 9,
    .citation "[mathlib]" "Associates.prime_pow_dvd_iff_le(重複度を包含に翻訳)"
      (.inMathlib "Associates.prime_pow_dvd_iff_le") 9,
    .citation "[mathlib]" "Ideal.dvd_iff_le(Dedekind 環では割り切りは包含)"
      (.inMathlib "Ideal.dvd_iff_le") 9,
    .citation "[ABC3]" "conductorADiv_eq_of_square(X_n の点での因子の一致)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_eq_of_square") 9,
    .implicitStep
      ("★★原文の「Σ の外で」は、本定理では「M が v を避ける」に対応する。" ++
       "★M = powers (n! : 𝒪_F) と取れば、ch v ∉ Σ がちょうど hM になる") 9,
    .implicitStep
      ("★★★残る配線: ConductorDescent.lean の因子の一致(スキームの水準)を " ++
       "𝒪_F[1/N] 上のイデアルの一致(本定理の仮定 h)に翻訳する段。" ++
       "★pullbackIdeal が局所化と両立することを示せばよい") 9 ]

end ABC3.Found.GenEll
