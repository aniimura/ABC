import ABC3.Found.GaloisRep.CollUniversal

/-!
# Galois (G6) 第 252 ブロック —— **★★★★★★★★共線性を三つの整除性に落とす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★変数の素元性を一般の位置で取る

第 224 では `ℤ[A,W]` の `X 0`・`X 1` について、`finSuccEquiv` を 1 段・2 段と
使い分けて素元性を示した。★3 変数になるとこの手は `X 2` で 3 段になる。

★★★そこで **`Equiv.swap i 0` による `renameEquiv` で任意の変数を `X 0` に運ぶ**形に
一般化した(`prime_X_fin`)。`X 0` の素元性だけ `finSuccEquiv` で示せば、
あとは変数の名前替えで済む——変数の個数にも位置にもよらない。

## ★★★互いに素であること

`IsCoprime`(Bezout)では駄目で `IsRelPrime` を使う点は第 224 と同じ。
必要なのは `¬(U ∣ V)`、`¬(U ∣ W)`、`¬(V ∣ W)` の 3 つで、どれも
`collEval` に整数を代入して `1 = 0` を出すだけ。

## ★★★三つを合わせる

    U^n ∣ P,  V^n ∣ P,  W^n ∣ P   ⟹   (UVW)^n ∣ P          (`cQ_pow_dvd`)

`IsRelPrime.mul_left` で `(U^n V^n)` と `W^n` が互いに素になることを使う。

## ★★★★★★★★到達点

    collDefect u v w q = 0
      ⟸ 分子が `U^n`・`V^n`・`W^n` で割れる(万有な環の中で)

★葉 (b) では二つだったものが三つになっただけで、形は第 224 の
`tateDefect_eq_zero_of_base` と同じである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `prime_X_fin` | ★★★**任意の変数が素元**(名前替えで一般化) |
| `isRelPrime_of_prime_not_dvd` | ★★素元で割り切らないなら互いに素 |
| `cQ_pow_dvd` | ★★★三つの整除性から `(UVW)^n` |
| `coll_dvd_pow_of_base` | ★★★★局所化での整除性は分子から |
| `collDefect_eq_zero_of_base` | ★★★★★★★★**三つの整除性に落ちた** |
-/

namespace ABC3.Found.GaloisRep

open MvPolynomial

/-! ## ★★★変数の素元性 -/

theorem prime_X_zero (n : ℕ) : Prime (MvPolynomial.X 0 : MvPolynomial (Fin (n + 1)) ℤ) := by
  have hp : Prime (Polynomial.X : Polynomial (MvPolynomial (Fin n) ℤ)) := Polynomial.prime_X
  rw [← MvPolynomial.finSuccEquiv_X_zero (R := ℤ) (n := n)] at hp
  exact (MulEquiv.prime_iff (M := MvPolynomial (Fin (n + 1)) ℤ)
    (MvPolynomial.finSuccEquiv ℤ n : MvPolynomial (Fin (n + 1)) ℤ ≃ₐ[ℤ] _).toRingEquiv).1 hp

/-- ★★★**任意の変数が素元**——`Equiv.swap i 0` で `X 0` に運ぶ。

★第 224 のように変数ごとに `finSuccEquiv` の段数を変える必要がない。 -/
theorem prime_X_fin (n : ℕ) (i : Fin (n + 1)) :
    Prime (MvPolynomial.X i : MvPolynomial (Fin (n + 1)) ℤ) := by
  have hX : (MvPolynomial.renameEquiv ℤ (Equiv.swap i 0) :
      MvPolynomial (Fin (n + 1)) ℤ ≃ₐ[ℤ] MvPolynomial (Fin (n + 1)) ℤ)
      (MvPolynomial.X i) = MvPolynomial.X 0 := by
    show MvPolynomial.rename (Equiv.swap i 0) (MvPolynomial.X i) = _
    rw [MvPolynomial.rename_X, Equiv.swap_apply_left]
  refine (MulEquiv.prime_iff (M := MvPolynomial (Fin (n + 1)) ℤ)
    ((MvPolynomial.renameEquiv ℤ (Equiv.swap i 0) :
      MvPolynomial (Fin (n + 1)) ℤ ≃ₐ[ℤ] MvPolynomial (Fin (n + 1)) ℤ).toRingEquiv)).1 ?_
  rw [show ((MvPolynomial.renameEquiv ℤ (Equiv.swap i 0) :
      MvPolynomial (Fin (n + 1)) ℤ ≃ₐ[ℤ] MvPolynomial (Fin (n + 1)) ℤ).toRingEquiv)
      (MvPolynomial.X i) = MvPolynomial.X 0 from hX]
  exact prime_X_zero n

theorem prime_cU : Prime cU := prime_X_fin 2 0
theorem prime_cV : Prime cV := prime_X_fin 2 1
theorem prime_cW : Prime cW := prime_X_fin 2 2

/-! ## ★★★互いに素であること -/

/-- ★★素元 `a` が `b` を割らないなら `a` と `b` は互いに素。 -/
theorem isRelPrime_of_prime_not_dvd {a b : CollBase} (ha : Prime a) (hab : ¬ (a ∣ b)) :
    IsRelPrime a b := by
  intro z hza hzb
  obtain ⟨c, hc⟩ := hza
  rcases ha.irreducible.isUnit_or_isUnit hc with h | h
  · exact h
  · exfalso
    obtain ⟨d, hd⟩ := isUnit_iff_exists_inv.1 h
    have haz : a ∣ z := ⟨d, by rw [hc, mul_assoc, hd, mul_one]⟩
    exact hab (haz.trans hzb)

theorem not_cU_dvd_cV : ¬ (cU ∣ cV) := by
  rintro ⟨c, hc⟩
  have h := congrArg (collEval (0 : ℤ) 1 1) hc
  rw [collEval_V, map_mul, collEval_U] at h
  simp at h

theorem not_cU_dvd_cW : ¬ (cU ∣ cW) := by
  rintro ⟨c, hc⟩
  have h := congrArg (collEval (0 : ℤ) 1 1) hc
  rw [collEval_W, map_mul, collEval_U] at h
  simp at h

theorem not_cV_dvd_cW : ¬ (cV ∣ cW) := by
  rintro ⟨c, hc⟩
  have h := congrArg (collEval (1 : ℤ) 0 1) hc
  rw [collEval_W, map_mul, collEval_V] at h
  simp at h

theorem isRelPrime_cU_cV : IsRelPrime cU cV := isRelPrime_of_prime_not_dvd prime_cU not_cU_dvd_cV
theorem isRelPrime_cU_cW : IsRelPrime cU cW := isRelPrime_of_prime_not_dvd prime_cU not_cU_dvd_cW
theorem isRelPrime_cV_cW : IsRelPrime cV cW := isRelPrime_of_prime_not_dvd prime_cV not_cV_dvd_cW

/-! ## ★★★三つを合わせる -/

/-- ★★★**`U^n`・`V^n`・`W^n` で割れるなら `(UVW)^n` で割れる**。 -/
theorem cQ_pow_dvd (n : ℕ) {P : CollBase} (hU : cU ^ n ∣ P) (hV : cV ^ n ∣ P)
    (hW : cW ^ n ∣ P) : cQ ^ n ∣ P := by
  rw [cQ, mul_pow, mul_pow]
  have h1 : cU ^ n * cV ^ n ∣ P := (isRelPrime_cU_cV.pow).mul_dvd hU hV
  have h2 : IsRelPrime (cU ^ n * cV ^ n) (cW ^ n) :=
    (isRelPrime_cU_cW.pow).mul_left (isRelPrime_cV_cW.pow)
  exact h2.mul_dvd h1 hW

/-! ## ★★★★局所化から分子へ -/

/-- ★★★★**局所化での整除性は分子の整除性から出る**。 -/
theorem coll_dvd_pow_of_base (n : ℕ) (x : CollUniv)
    (H : ∀ (P : CollBase) (d : collDenoms),
      x * algebraMap CollBase CollUniv (d : CollBase) = algebraMap CollBase CollUniv P →
      cQ ^ n ∣ P) :
    (kU * kV * kW) ^ n ∣ x := by
  obtain ⟨⟨P, d⟩, hPd⟩ := IsLocalization.surj collDenoms x
  obtain ⟨c, hc⟩ := H P d hPd
  obtain ⟨v, hv⟩ := IsLocalization.map_units CollUniv d
  refine ⟨algebraMap CollBase CollUniv c * ↑v⁻¹, ?_⟩
  have hkey : x * (v : CollUniv) = (kU * kV * kW) ^ n * algebraMap CollBase CollUniv c := by
    rw [hv, hPd, hc, map_mul, map_pow, kQ_eq]
  calc x = x * (v : CollUniv) * ↑v⁻¹ := by rw [mul_assoc, v.mul_inv, mul_one]
    _ = (kU * kV * kW) ^ n * algebraMap CollBase CollUniv c * ↑v⁻¹ := by rw [hkey]
    _ = (kU * kV * kW) ^ n * (algebraMap CollBase CollUniv c * ↑v⁻¹) := by ring

/-! ## ★★★★★★★★三つの整除性への帰着 -/

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**共線性は万有な環の分子の三つの整除性に落ちた**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDefect_eq_zero_of_base [IsAdicComplete I R] (u v w q : R) (hq : q ∈ I)
    (huvw : u * v * w = q) (hcp : ∀ i, IsUnit (1 - collPts u v w i))
    (HU : ∀ (n : ℕ) (P : CollBase) (d : collDenoms),
      collDefectTrunc n kU kV kW (kU * kV * kW) * algebraMap CollBase CollUniv (d : CollBase)
        = algebraMap CollBase CollUniv P → cU ^ n ∣ P)
    (HV : ∀ (n : ℕ) (P : CollBase) (d : collDenoms),
      collDefectTrunc n kU kV kW (kU * kV * kW) * algebraMap CollBase CollUniv (d : CollBase)
        = algebraMap CollBase CollUniv P → cV ^ n ∣ P)
    (HW : ∀ (n : ℕ) (P : CollBase) (d : collDenoms),
      collDefectTrunc n kU kV kW (kU * kV * kW) * algebraMap CollBase CollUniv (d : CollBase)
        = algebraMap CollBase CollUniv P → cW ^ n ∣ P) :
    collDefect u v w q hq = 0 := by
  refine collDefect_eq_zero_of_univ u v w q hq huvw hcp fun n => ?_
  refine coll_dvd_pow_of_base n _ fun P d hPd => ?_
  exact cQ_pow_dvd n (HU n P d hPd) (HV n P d hPd) (HW n P d hPd)

/-! ## ★出典の紐付け(`.src`) -/

def prime_X_fin.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——変数の素元性、名前替えで一般化)",
    sectionId := "genell-def-3-3" }

def collDefect_eq_zero_of_base.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性を三つの整除性に落とす)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
