import ABC3.Found.GaloisRep.TateUniversal

/-!
# Galois (G6) 第 224 ブロック —— **★★★★★★★★整除性を `A^n` と `W^n` に分ける**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★`(AW)^n ∣` は `A^n ∣` と `W^n ∣` に分かれる

第 223 で葉 (b) は「万有な環 `TateUniv` の中で `(AW)^n ∣ tateDefectTrunc n A W (AW)`」に
落ちた。本ブロックはこれを **`ℤ[A,W]` の中の二つの整除性**に分ける。

★★`A` と `W` は `ℤ[A,W]` の**素元**である。mathlib には `MvPolynomial.prime_X` が
無かったので、`MvPolynomial.finSuccEquiv` で `Polynomial` に移して `Polynomial.prime_X`・
`Polynomial.prime_C_iff` から作った。

★★★`A` と `W` は互いに素(`IsRelPrime`)——`A ∣ W` なら `A ↦ 0`、`W ↦ 1` の評価で
`0 ∣ 1` になってしまう。よって `IsRelPrime.mul_dvd` で

    A^n ∣ P  かつ  W^n ∣ P   ⟹   (AW)^n ∣ P

★`IsCoprime`(Bezout の意味)では**駄目**である——`(A, W)` は真のイデアルなので
`IsCoprime A W` は成り立たない。**`IsRelPrime`(共通の非単元因子が無い)を使う**。

## ★★★★局所化から分子へ

`TateUniv` の元 `x` は `IsLocalization.surj` で `x·(分母) = (分子)` と書ける。
分子が `(AW)^n` で割れれば、分母は単元なので `x` も `(AW)^n` で割れる(`dvd_pow_of_base`)。

## ★★★★★★★★到達点

> `ℤ[A,W]` の中で「分子が `A^n` で割れる」「分子が `W^n` で割れる」を示せば、
> **任意の完備環で Weierstrass 方程式が成り立つ**。

これが `tateDefect_eq_zero_of_base` である。

★★切り詰めた差も `(a,w)` の交換で不変(`tateDefectTrunc_symm`)なので、
`A` 側と `W` 側は**同じ議論の入れ替え**で済む。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `prime_univA` / `prime_univW` | ★★`A`・`W` は `ℤ[A,W]` の素元 |
| `isRelPrime_univA_univW` | ★★★`A` と `W` は互いに素 |
| `univQ_pow_dvd` | ★★★`A^n ∣ P ∧ W^n ∣ P → (AW)^n ∣ P` |
| `dvd_pow_of_base` | ★★★★局所化の整除性は分子の整除性から |
| `tateDefectTrunc_symm` | ★★★★切り詰めた差も交換で不変 |
| `tateDefect_eq_zero_of_base` | ★★★★★★★★**葉 (b) は二つの整除性に落ちた** |
-/

namespace ABC3.Found.GaloisRep

open MvPolynomial

/-! ## ★★`A`・`W` の素元性 -/

/-- ★★`A` は `ℤ[A,W]` の素元。

★mathlib に `MvPolynomial.prime_X` は無い。`finSuccEquiv` で `Polynomial` に移す。 -/
theorem prime_univA : Prime univA := by
  have hp : Prime (Polynomial.X : Polynomial (MvPolynomial (Fin 1) ℤ)) := Polynomial.prime_X
  rw [← MvPolynomial.finSuccEquiv_X_zero (R := ℤ) (n := 1)] at hp
  exact (MulEquiv.prime_iff (M := MvPolynomial (Fin 2) ℤ)
    (MvPolynomial.finSuccEquiv ℤ 1 : MvPolynomial (Fin 2) ℤ ≃ₐ[ℤ] _).toRingEquiv).1 hp

/-- ★★`W` は `ℤ[A,W]` の素元。 -/
theorem prime_univW : Prime univW := by
  have hp0 : Prime (Polynomial.X : Polynomial (MvPolynomial (Fin 0) ℤ)) := Polynomial.prime_X
  rw [← MvPolynomial.finSuccEquiv_X_zero (R := ℤ) (n := 0)] at hp0
  have h1 : Prime (X 0 : MvPolynomial (Fin 1) ℤ) :=
    (MulEquiv.prime_iff (M := MvPolynomial (Fin 1) ℤ)
      (MvPolynomial.finSuccEquiv ℤ 0 : MvPolynomial (Fin 1) ℤ ≃ₐ[ℤ] _).toRingEquiv).1 hp0
  have hp : Prime (Polynomial.C (X 0 : MvPolynomial (Fin 1) ℤ)) := Polynomial.prime_C_iff.2 h1
  rw [← MvPolynomial.finSuccEquiv_X_succ (R := ℤ) (n := 1) (j := 0)] at hp
  have hsucc : (0 : Fin 1).succ = (1 : Fin 2) := rfl
  rw [hsucc] at hp
  exact (MulEquiv.prime_iff (M := MvPolynomial (Fin 2) ℤ)
    (MvPolynomial.finSuccEquiv ℤ 1 : MvPolynomial (Fin 2) ℤ ≃ₐ[ℤ] _).toRingEquiv).1 hp

/-! ## ★★★互いに素であること -/

theorem not_univA_dvd_univW : ¬ (univA ∣ univW) := by
  rintro ⟨c, hc⟩
  have h := congrArg (tateEval (0 : ℤ) (1 : ℤ)) hc
  rw [tateEval_W, map_mul, tateEval_A] at h
  simp at h

/-- ★★★**`A` と `W` は互いに素**。

★`IsCoprime`(Bezout)では駄目——`(A,W)` は真のイデアル。`IsRelPrime` を使う。 -/
theorem isRelPrime_univA_univW : IsRelPrime univA univW := by
  intro z hzA hzW
  obtain ⟨c, hc⟩ := hzA
  rcases prime_univA.irreducible.isUnit_or_isUnit hc with h | h
  · exact h
  · exfalso
    obtain ⟨d, hd⟩ := isUnit_iff_exists_inv.1 h
    have hAz : univA ∣ z := ⟨d, by rw [hc, mul_assoc, hd, mul_one]⟩
    exact not_univA_dvd_univW (hAz.trans hzW)

/-- ★★★**`A^n ∣ P` かつ `W^n ∣ P` なら `(AW)^n ∣ P`**。 -/
theorem univQ_pow_dvd (n : ℕ) {P : TateBase} (hA : univA ^ n ∣ P) (hW : univW ^ n ∣ P) :
    univQ ^ n ∣ P := by
  rw [univQ, mul_pow]
  exact (isRelPrime_univA_univW.pow).mul_dvd hA hW

/-! ## ★★★★局所化から分子へ -/

/-- ★★★★**局所化での整除性は分子の整除性から出る**。 -/
theorem dvd_pow_of_base (n : ℕ) (x : TateUniv)
    (H : ∀ (P : TateBase) (d : tateDenoms),
      x * algebraMap TateBase TateUniv (d : TateBase) = algebraMap TateBase TateUniv P →
      univQ ^ n ∣ P) :
    (uA * uW) ^ n ∣ x := by
  obtain ⟨⟨P, d⟩, hPd⟩ := IsLocalization.surj tateDenoms x
  obtain ⟨c, hc⟩ := H P d hPd
  obtain ⟨v, hv⟩ := IsLocalization.map_units TateUniv d
  refine ⟨algebraMap TateBase TateUniv c * ↑v⁻¹, ?_⟩
  have hkey : x * (v : TateUniv) = (uA * uW) ^ n * algebraMap TateBase TateUniv c := by
    rw [hv, hPd, hc, map_mul, map_pow, uA_mul_uW]
  calc x = x * (v : TateUniv) * ↑v⁻¹ := by
        rw [mul_assoc, v.mul_inv, mul_one]
    _ = (uA * uW) ^ n * algebraMap TateBase TateUniv c * ↑v⁻¹ := by rw [hkey]
    _ = (uA * uW) ^ n * (algebraMap TateBase TateUniv c * ↑v⁻¹) := by ring

/-! ## ★★★★切り詰めの対称性 -/

section Symm

variable {R : Type} [CommRing R]

theorem tateXtrunc_symm (n : ℕ) (a w q : R) :
    tateXtrunc n w a q = tateXtrunc n a w q := by
  rw [tateXtrunc, tateXtrunc]
  ring

theorem tateYtrunc_swap (n : ℕ) (a w q : R) :
    tateYtrunc n w a q = -tateXtrunc n a w q - tateYtrunc n a w q := by
  rw [tateYtrunc, tateYtrunc, tateXtrunc]
  ring

/-- ★★★★**切り詰めた差も `(a,w)` の交換で不変**。 -/
theorem tateDefectTrunc_symm (n : ℕ) (a w q : R) :
    tateDefectTrunc n w a q = tateDefectTrunc n a w q := by
  rw [tateDefectTrunc, tateDefectTrunc, tateXtrunc_symm, tateYtrunc_swap]
  ring

end Symm

/-! ## ★★★★★★★★葉 (b) の最終的な帰着 -/

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**葉 (b) は万有な環の分子の二つの整除性に落ちた**。

`ℤ[A,W]` の中で「分子が `A^n` で割れる」「分子が `W^n` で割れる」を示せば、
**任意の完備環で Weierstrass 方程式が成り立つ**。

★★切り詰めた差は `(a,w)` の交換で不変(`tateDefectTrunc_symm`)なので、
二つは同じ議論の入れ替えで済む。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateDefect_eq_zero_of_base [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (HA : ∀ (n : ℕ) (P : TateBase) (d : tateDenoms),
      tateDefectTrunc n uA uW (uA * uW) * algebraMap TateBase TateUniv (d : TateBase)
        = algebraMap TateBase TateUniv P → univA ^ n ∣ P)
    (HW : ∀ (n : ℕ) (P : TateBase) (d : tateDenoms),
      tateDefectTrunc n uA uW (uA * uW) * algebraMap TateBase TateUniv (d : TateBase)
        = algebraMap TateBase TateUniv P → univW ^ n ∣ P) :
    tateDefect a w q hq = 0 := by
  refine tateDefect_eq_zero_of_univ a w q hq haw ha hw fun n => ?_
  refine dvd_pow_of_base n _ fun P d hPd => ?_
  exact univQ_pow_dvd n (HA n P d hPd) (HW n P d hPd)

/-! ## ★出典の紐付け(`.src`) -/

def prime_univA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——万有な環の素元)",
    sectionId := "genell-def-3-3" }

def univQ_pow_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——(AW)^n の整除性を A^n と W^n に分ける)",
    sectionId := "genell-def-3-3" }

def tateDefect_eq_zero_of_base.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (b) の分子の整除性への帰着)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
