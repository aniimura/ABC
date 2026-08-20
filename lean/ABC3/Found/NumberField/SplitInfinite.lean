/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitCount
import ABC3.Found.NumberField.PrimeDivisorsOfValues
import Mathlib.NumberTheory.NumberField.Ideal.KummerDedekind

/-!
# Dedekind–Kummer の 1 段 —— 根があれば剰余次数 1 の素がある(鎖 `cheb` の `cheb-spl-infinite`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★迂回路の最後の 1 段

鎖 `cheb` の `cheb-spl-infinite`(「Galois 拡大で完全分解する素数は無限個」)は
**全 Chebotarev を使わずに**出る:

| 段 | 内容 | 在庫 |
|---|---|---|
| 1 | `f` の値を割る素数は無限個(Schur) | `PrimeDivisorsOfValues.lean` |
| 2 | ★**`f` が mod `p` で根を持つ ⇒ 剰余次数 1 の素がある** | ★**本ファイル** |
| 3 | Galois で `e = f = 1` の素が 1 つあれば完全分解 | `SplitCount.lean` |

★2 は **Dedekind–Kummer**(mathlib の `NumberField/Ideal/KummerDedekind.lean`)そのもの。
`p ∤ exponent θ` のとき「`p` の上の素イデアル」と「`minpoly ℤ θ` mod `p` の
モニック既約因子」が 1 対 1 で、**剰余次数は因子の次数**である。
★根 `a` は 1 次因子 `X - a` を与えるから、対応する素の剰余次数は `1` である。

## ★★★底を `ℤ` にするために

`SplitCount.lean` の Galois の等式は底が数体である。★`ℚ` の整数環 `𝓞 ℚ` は
`ℤ` と**同型だが同一ではない**ので、`IsGaloisGroup Gal(K/ℚ) ℤ (𝓞 K)` を
自分で立てる必要がある。
★★中身は 1 行 —— **`ℤ` からの環準同型は一意**(`Subsingleton (ℤ →+* R)`)だから、
`𝓞 ℚ ≅ ℤ` を通した図式は自動的に可換である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `X_sub_C_mem_monicFactorsMod` | 根は 1 次のモニック既約因子を与える |
| `exists_inertiaDeg_eq_one_of_isRoot` | ★★根があれば剰余次数 1 の素がある |
| `isGaloisGroup_int` | `IsGaloisGroup Gal(K/ℚ) ℤ (𝓞 K)` |
| `splitsCompletely_of_isRoot` | ★★★★根が単純なら `p` は完全分解する |
-/

namespace ABC3.Found.NF

open _root_.NumberField Polynomial Ideal UniqueFactorizationMonoid
open scoped _root_.NumberField

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★1. 根はモニック既約因子を与える -/

omit [NumberField K] in
/-- ★**根 `a` は `X - a` という 1 次のモニック既約因子を与える**。

★`(ZMod p)[X]` は体上の 1 変数多項式環なので UFD で、`X - a` はモニックゆえ
normalize で不動、したがって `normalizedFactors` の元そのものになる。 -/
theorem X_sub_C_mem_monicFactorsMod (θ : 𝓞 K) {p : ℕ} [Fact p.Prime] {a : ZMod p}
    (ha : ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))).IsRoot a) :
    (X - C a) ∈ _root_.RingOfIntegers.monicFactorsMod θ p := by
  have hmon : ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))).Monic :=
    (minpoly.monic θ.isIntegral).map _
  have hg0 : ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))) ≠ 0 := hmon.ne_zero
  have hdvd : (X - C a) ∣ ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))) :=
    dvd_iff_isRoot.mpr ha
  obtain ⟨c, hc, hassoc⟩ :=
    exists_mem_normalizedFactors_of_dvd hg0 (irreducible_X_sub_C a) hdvd
  have hassoc' : Associated (X - C a) c := hassoc
  have hnorm : normalize c = c := normalize_normalized_factor c hc
  have hXnorm : normalize (X - C a) = X - C a := (monic_X_sub_C a).normalize_eq_self
  have : X - C a = c := by
    rw [← hXnorm, ← hnorm]
    exact normalize_eq_normalize hassoc'.dvd hassoc'.symm.dvd
  rw [this]
  exact Multiset.mem_toFinset.mpr hc

/-! ## ★2. 根があれば剰余次数 1 の素がある -/

/-- ★★★**Dedekind–Kummer の 1 段** —— `minpoly ℤ θ` が mod `p` で根 `a` を持てば、
`p` の上に**剰余次数 1** の素イデアルがある。分岐指数は `X - a` の重複度である。 -/
theorem exists_inertiaDeg_eq_one_of_isRoot (θ : 𝓞 K) {p : ℕ} [Fact p.Prime]
    (hp : ¬ p ∣ _root_.RingOfIntegers.exponent θ) {a : ZMod p}
    (ha : ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))).IsRoot a) :
    ∃ P : Ideal (𝓞 K), P.IsPrime ∧ P.LiesOver (span {(p : ℤ)}) ∧
      inertiaDeg (span {(p : ℤ)}) P = 1 ∧
      ramificationIdx (span {(p : ℤ)}) P
        = multiplicity (X - C a) ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))) := by
  have hmem := X_sub_C_mem_monicFactorsMod θ ha
  refine ⟨((_root_.NumberField.Ideal.primesOverSpanEquivMonicFactorsMod hp).symm
      ⟨X - C a, hmem⟩ : Ideal (𝓞 K)), ?_, ?_, ?_, ?_⟩
  · exact ((_root_.NumberField.Ideal.primesOverSpanEquivMonicFactorsMod hp).symm
      ⟨X - C a, hmem⟩).2.1
  · exact ((_root_.NumberField.Ideal.primesOverSpanEquivMonicFactorsMod hp).symm
      ⟨X - C a, hmem⟩).2.2
  · rw [_root_.NumberField.Ideal.inertiaDeg_primesOverSpanEquivMonicFactorsMod_symm_apply' hp hmem]
    simp
  · exact _root_.NumberField.Ideal.ramificationIdx_primesOverSpanEquivMonicFactorsMod_symm_apply'
      hp hmem

/-! ## ★3. 底を `ℤ` にする -/

/-- ★★**`ℤ` からの環準同型は一意** なので、`𝓞 ℚ ≅ ℤ` を通した図式は自動的に可換。 -/
theorem algebraMap_int_ringOfIntegersEquiv (a : 𝓞 ℚ) :
    algebraMap ℤ (𝓞 K) (Rat.ringOfIntegersEquiv a) = algebraMap (𝓞 ℚ) (𝓞 K) a := by
  have h : (algebraMap ℤ (𝓞 K))
      = (algebraMap (𝓞 ℚ) (𝓞 K)).comp (Rat.ringOfIntegersEquiv.symm : ℤ →+* 𝓞 ℚ) :=
    Subsingleton.elim _ _
  have h2 := congrFun (congrArg (fun f : ℤ →+* 𝓞 K => (f : ℤ → 𝓞 K)) h)
    (Rat.ringOfIntegersEquiv a)
  rw [h2]
  show algebraMap (𝓞 ℚ) (𝓞 K)
    (Rat.ringOfIntegersEquiv.symm (Rat.ringOfIntegersEquiv a)) = _
  rw [RingEquiv.symm_apply_apply]

/-- ★★**`Gal(K/ℚ)` は `ℤ ⊆ 𝓞 K` の Galois 群**。

★`IsGaloisGroup Gal(K/ℚ) (𝓞 ℚ) (𝓞 K)` は mathlib の instance で、
そこから `𝓞 ℚ ≅ ℤ` で移すだけである。 -/
instance isGaloisGroup_int [IsGalois ℚ K] : IsGaloisGroup (K ≃ₐ[ℚ] K) ℤ (𝓞 K) :=
  IsGaloisGroup.of_ringEquiv (K ≃ₐ[ℚ] K) (𝓞 ℚ) ℤ (𝓞 K)
    Rat.ringOfIntegersEquiv algebraMap_int_ringOfIntegersEquiv

/-! ## ★4. 完全分解 -/

/-- ★★★★★**Galois 拡大では、不分岐で剰余次数 1 の素が 1 つあれば `p` は完全分解する**
(底が `ℤ` の形)。

★Galois の基本等式 `g·e·f = |Gal(K/ℚ)| = [K:ℚ]` を当てるだけ。 -/
theorem ncard_primesOver_eq_finrank_of_int [IsGalois ℚ K] {p : ℕ} [Fact p.Prime]
    (P : Ideal (𝓞 K)) [P.IsPrime] [P.LiesOver (span {(p : ℤ)})]
    (he : ramificationIdx (span {(p : ℤ)}) P = 1)
    (hf : inertiaDeg (span {(p : ℤ)}) P = 1) :
    (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K := by
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := span {(p : ℤ)}) (by simp [NeZero.ne p]) (𝓞 K) (K ≃ₐ[ℚ] K)
  rw [Ideal.ramificationIdxIn_eq_ramificationIdx (span {(p : ℤ)}) P (K ≃ₐ[ℚ] K),
    Ideal.inertiaDegIn_eq_inertiaDeg (span {(p : ℤ)}) P (K ≃ₐ[ℚ] K), he, hf] at hG
  rw [IsGalois.card_aut_eq_finrank ℚ K] at hG
  simpa using hG

/-- ★★★★★★**根が単純なら `p` は完全分解する**(Galois の場合)。

★★これが迂回路 `cheb-spl-infinite` の骨である ——
Schur の補題(`infinite_primes_dvd_eval`)が「`f` の値を割る素数は無限個」を与え、
そのうち `p ∤ exponent θ` かつ根が単純なものについて、この定理が完全分解を与える。 -/
theorem splitsCompletely_of_isRoot [IsGalois ℚ K] (θ : 𝓞 K) {p : ℕ} [Fact p.Prime]
    (hp : ¬ p ∣ _root_.RingOfIntegers.exponent θ) {a : ZMod p}
    (ha : ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))).IsRoot a)
    (hsimple : multiplicity (X - C a) ((minpoly ℤ θ).map (Int.castRingHom (ZMod p))) = 1) :
    (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K := by
  obtain ⟨P, hP, hlies, hf, he⟩ := exists_inertiaDeg_eq_one_of_isRoot θ hp ha
  haveI := hP
  haveI := hlies
  exact ncard_primesOver_eq_finrank_of_int P (he.trans hsimple) hf

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Theorem 6.4, (iv)` が使う「完全分解する素数がある」の迂回路。 -/
def splitsCompletely_of_isRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Dedekind–Kummer による完全分解の判定(Chebotarev の迂回路)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
