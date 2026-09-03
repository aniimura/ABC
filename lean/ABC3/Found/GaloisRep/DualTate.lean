/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DualAdic
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Found.GaloisRep.TateMultRed
import ABC3.Found.GaloisRep.TateSigma
import ABC3.Found.GaloisRep.TateCoeffTail

/-!
# Galois (G6) 第 851 ブロック —— **★★★★★★★★★★双対数で Tate 級数を微分する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★これは何か

第 850 で `R[ε]` が `I[ε]` 進完備であることを示した。本ブロックはその上で

    `tateXterm (mk t s) = mk (tateXterm t) (s(1+t)/(1−t)³)`
    `tateYterm (mk t s) = mk (tateYterm t) (s·t(2+t)/(1−t)⁴)`

を示す。★`s = t` を入れると `ε` 成分は `tateDXterm t`・`tateDYterm t`
（第 846）そのものであり、`s = −t` なら符号が反転する——
これが `X = f(a)+T(a)+f(w)+T(w)−2s₁` の `w` 側で符号が変わる理由である。

## ★★★★★計算の中身

`1 − mk t s = mk (1−t) (−s)` の逆元は `mk r (s r²)`（`r = (1−t)⁻¹`）である:

    `ε` 成分: `(1−t)·s r² + r·(−s) = s r((1−t)r − 1) = 0` ✓

これを使うと `(mk r (sr²))² = mk (r²) (2sr³)` なので

    `mk t s · mk (r²) (2sr³) = mk (t r²) (s(2t r³ + r²)) = mk (f(t)) (s(1+t)r³)`

（`r² = r³(1−t)` を使った）。★`s = t` で `t(1+t)r³ = Df(t)` ✓

## ★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ring_inverse_eq_of_mul_eq_one` | `x·y = 1` から `Ring.inverse x = y` |
| `isUnit_one_sub_dual` | `1 − mk t s` は単元 |
| `ringInverse_one_sub_dual` | ★★逆元の明示形 |
| `tateXterm_dual` / `tateYterm_dual` | ★★★★**項の微分** |
| `adicSum_dual` | ★★★**`adicSum` は成分ごと** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R]

/-! ## ★冪の成分 -/

namespace DualNum

@[simp] theorem re_pow (x : DualNum R) (n : ℕ) : (x ^ n).re = x.re ^ n := by
  induction n with
  | zero => simp
  | succ k ih => rw [pow_succ, re_mul, ih, pow_succ]

theorem eps_sq (x : DualNum R) : (x ^ 2).eps = 2 * x.re * x.eps := by
  rw [sq, eps_mul]
  ring

theorem eps_pow (x : DualNum R) (n : ℕ) :
    (x ^ (n + 1)).eps = ((n : R) + 1) * x.re ^ n * x.eps := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [pow_succ, eps_mul, ih, re_pow]
      push_cast
      ring

theorem eps_cube (x : DualNum R) : (x ^ 3).eps = 3 * x.re ^ 2 * x.eps := by
  have h3 : (3 : ℕ) = 2 + 1 := rfl
  rw [h3, pow_succ, eps_mul, eps_sq, re_pow]
  ring

end DualNum

/-! ## ★★★★★逆元 -/

theorem DualNum.re_one_sub_mk (t s : R) : (1 - DualNum.mk t s).re = 1 - t := by simp

theorem DualNum.eps_one_sub_mk (t s : R) : (1 - DualNum.mk t s).eps = -s := by simp

/-- ★★`1 − mk t s` は `1 − t` が単元なら単元。 -/
theorem isUnit_one_sub_dual {t s : R} (hu : IsUnit (1 - t)) :
    IsUnit (1 - DualNum.mk t s) := by
  rw [DualNum.isUnit_iff, DualNum.re_one_sub_mk]
  exact hu

/-- ★★★★★★**逆元の明示形** `(1 − mk t s)⁻¹ = mk r (s r²)`。 -/
theorem ringInverse_one_sub_dual {t s : R} (hu : IsUnit (1 - t)) :
    Ring.inverse (1 - DualNum.mk t s)
      = DualNum.mk (Ring.inverse (1 - t)) (s * Ring.inverse (1 - t) ^ 2) := by
  refine ring_inverse_eq_of_mul_eq_one (isUnit_one_sub_dual hu) ?_
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  ext
  · simp only [DualNum.re_mul, DualNum.re_one_sub_mk, DualNum.re_mk, DualNum.re_one]
    exact hr
  · simp only [DualNum.eps_mul, DualNum.re_one_sub_mk, DualNum.eps_one_sub_mk,
      DualNum.re_mk, DualNum.eps_mk, DualNum.eps_one]
    linear_combination (s * Ring.inverse (1 - t)) * hr

/-! ## ★★★★★★★★項の微分 -/

/-- ★★★★★★★★**`f(t + εs) = f(t) + ε·s(1+t)/(1−t)³`**。 -/
theorem tateXterm_dual {t s : R} (hu : IsUnit (1 - t)) :
    tateXterm (DualNum.mk t s)
      = DualNum.mk (tateXterm t) (s * (1 + t) * Ring.inverse (1 - t) ^ 3) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateXterm, ringInverse_one_sub_dual hu]
  ext
  · simp [tateXterm]
  · simp only [DualNum.eps_mul, DualNum.re_mk, DualNum.eps_mk, DualNum.re_pow,
      DualNum.eps_sq, DualNum.eps_mk]
    linear_combination (-(s * Ring.inverse (1 - t) ^ 2)) * hr

/-- ★★★★★★★★**`g(t + εs) = g(t) + ε·s·t(2+t)/(1−t)⁴`**。 -/
theorem tateYterm_dual {t s : R} (hu : IsUnit (1 - t)) :
    tateYterm (DualNum.mk t s)
      = DualNum.mk (tateYterm t) (s * t * (2 + t) * Ring.inverse (1 - t) ^ 4) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateYterm, ringInverse_one_sub_dual hu]
  ext
  · simp [tateYterm]
  · simp only [DualNum.eps_mul, DualNum.re_mk, DualNum.eps_mk, DualNum.re_pow,
      DualNum.eps_sq, DualNum.eps_cube]
    linear_combination (-(2 * s * t * Ring.inverse (1 - t) ^ 3)) * hr

/-! ## ★★★★★★`adicSum` は成分ごと -/

theorem DualNum.mk_pow (a : R) (n : ℕ) : (DualNum.mk a 0) ^ n = DualNum.mk (a ^ n) 0 := by
  induction n with
  | zero => ext <;> simp
  | succ k ih =>
      rw [pow_succ, ih]
      ext
      · simp [pow_succ]
      · simp

theorem DualNum.mk_zero_mul (a u v : R) :
    DualNum.mk a 0 * DualNum.mk u v = DualNum.mk (a * u) (a * v) := by
  ext <;> simp

theorem partialSum_dual (f g : ℕ → R) (n : ℕ) :
    partialSum (fun k => DualNum.mk (f k) (g k)) n
      = DualNum.mk (partialSum f n) (partialSum g n) := by
  induction n with
  | zero => ext <;> simp [partialSum]
  | succ m ih =>
      rw [partialSum, partialSum, partialSum, Finset.sum_range_succ,
        Finset.sum_range_succ, Finset.sum_range_succ]
      rw [show (∑ k ∈ Finset.range m, DualNum.mk (f k) (g k))
          = DualNum.mk (partialSum f m) (partialSum g m) from ih]
      ext <;> simp [partialSum]

theorem dual_mem_pow {I : Ideal R} {f g : ℕ → R} (hf : ∀ n, f n ∈ I ^ n)
    (hg : ∀ n, g n ∈ I ^ n) (n : ℕ) : DualNum.mk (f n) (g n) ∈ (dualIdeal I) ^ n := by
  rw [dualIdeal_pow]
  exact ⟨by simpa using hf n, by simpa using hg n⟩

/-- ★★★★★★**`adicSum` は成分ごとに取れる**。 -/
theorem adicSum_dual {I : Ideal R} [IsAdicComplete I R] (f g : ℕ → R)
    (hf : ∀ n, f n ∈ I ^ n) (hg : ∀ n, g n ∈ I ^ n) :
    adicSum (fun n => DualNum.mk (f n) (g n)) (dual_mem_pow hf hg)
      = DualNum.mk (adicSum f hf) (adicSum g hg) := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  rw [partialSum_dual, SModEq.sub_mem, mem_dual_smul_top]
  have h1 := adicSum_spec f hf n
  have h2 := adicSum_spec g hg n
  rw [SModEq.sub_mem, mem_smul_top_self] at h1 h2
  refine ⟨?_, ?_⟩
  · simpa using h1
  · simpa using h2

/-! ## ★★★★実部への埋め込み -/

/-- ★`R → R[ε]`、`x ↦ x + 0·ε`。 -/
def DualNum.inlHom : R →+* DualNum R where
  toFun x := DualNum.mk x 0
  map_one' := by ext <;> simp
  map_mul' a b := by ext <;> simp
  map_zero' := by ext <;> simp
  map_add' a b := by ext <;> simp

@[simp] theorem DualNum.inlHom_apply (x : R) : DualNum.inlHom x = DualNum.mk x 0 := rfl

theorem DualNum.intCast_eq (c : ℤ) : ((c : ℤ) : DualNum R) = DualNum.mk ((c : ℤ) : R) 0 := by
  rw [← DualNum.inlHom_apply, map_intCast]

theorem partialEval_dual (f : PowerSeries ℤ) (q : R) (N : ℕ) :
    partialEval f (DualNum.mk q 0) N = DualNum.mk (partialEval f q N) 0 := by
  simp only [partialEval]
  have hterm : ∀ n ∈ Finset.range N,
      ((PowerSeries.coeff n f : ℤ) : DualNum R) * (DualNum.mk q 0) ^ n
        = DualNum.inlHom (((PowerSeries.coeff n f : ℤ) : R) * q ^ n) := by
    intro n _
    rw [DualNum.intCast_eq, map_mul, map_pow, DualNum.inlHom_apply, DualNum.inlHom_apply]
  rw [Finset.sum_congr rfl hterm, ← map_sum, DualNum.inlHom_apply]

/-- ★★★★**`q` だけの級数は `ε` 成分を持たない**——`s₁` が微分で消える理由。 -/
theorem evalAdic_dual {I : Ideal R} [IsAdicComplete I R] (f : PowerSeries ℤ) (q : R)
    (hq : q ∈ I) :
    evalAdic f (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (evalAdic f q hq) 0 := by
  refine evalAdic_unique _ _ _ _ (fun n => ?_)
  rw [partialEval_dual, SModEq.sub_mem, mem_dual_smul_top]
  have h := evalAdic_spec f q hq n
  rw [SModEq.sub_mem, mem_smul_top_self] at h
  refine ⟨?_, ?_⟩
  · simpa using h
  · simp

/-! ## ★★★★★★★★尾の微分 -/

theorem tateXtail_dual_pos {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail (DualNum.mk u u) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateXtail u q hq) (tateDXtail u q hq) := by
  have hterm : ∀ n : ℕ,
      tateXterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u u)
        = DualNum.mk (tateXterm (q ^ (n + 1) * u)) (tateDXterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateXterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateDXterm]
  rw [tateXtail, tateXtail, tateDXtail,
    adicSum_congr (tateXtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateXtail_aux hq) (tateDXtail_aux hq)) hterm,
    adicSum_dual]

theorem tateXtail_dual_neg {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail (DualNum.mk u (-u)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateXtail u q hq) (-tateDXtail u q hq) := by
  have hmem : ∀ n : ℕ, (-1 : R) * tateDXterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateDXtail_aux hq n)
  have hterm : ∀ n : ℕ,
      tateXterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u (-u))
        = DualNum.mk (tateXterm (q ^ (n + 1) * u))
            ((-1 : R) * tateDXterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateXterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateDXterm]
    ext <;> simp <;> try ring
  rw [tateXtail, tateXtail, tateDXtail,
    adicSum_congr (tateXtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateXtail_aux hq) hmem) hterm,
    adicSum_dual _ _ (tateXtail_aux hq) hmem, adicSum_smul (-1 : R) _ (tateDXtail_aux hq)]
  ext <;> simp

theorem tateYtail_dual_pos {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail (DualNum.mk u u) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateYtail u q hq) (tateDYtail u q hq) := by
  have hterm : ∀ n : ℕ,
      tateYterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u u)
        = DualNum.mk (tateYterm (q ^ (n + 1) * u)) (tateDYterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateYterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateDYterm]
    ext <;> simp <;> try ring
  rw [tateYtail, tateYtail, tateDYtail,
    adicSum_congr (tateYtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateYtail_aux hq) (tateDYtail_aux hq)) hterm,
    adicSum_dual]

theorem tateYtail_dual_neg {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail (DualNum.mk u (-u)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateYtail u q hq) (-tateDYtail u q hq) := by
  have hmem : ∀ n : ℕ, (-1 : R) * tateDYterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateDYtail_aux hq n)
  have hterm : ∀ n : ℕ,
      tateYterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u (-u))
        = DualNum.mk (tateYterm (q ^ (n + 1) * u))
            ((-1 : R) * tateDYterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateYterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateDYterm]
    ext <;> simp <;> try ring
  rw [tateYtail, tateYtail, tateDYtail,
    adicSum_congr (tateYtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateYtail_aux hq) hmem) hterm,
    adicSum_dual _ _ (tateYtail_aux hq) hmem, adicSum_smul (-1 : R) _ (tateDYtail_aux hq)]
  ext <;> simp

/-! ## ★★★★★★★★★★`X`・`Y` の微分 -/

theorem DualNum.two_eq_mk : (2 : DualNum R) = DualNum.mk 2 0 := by
  have h : DualNum.inlHom (2 : R) = (2 : DualNum R) := map_ofNat DualNum.inlHom 2
  rw [← h, DualNum.inlHom_apply]

@[simp] theorem DualNum.re_ofNat (n : ℕ) [n.AtLeastTwo] :
    (OfNat.ofNat n : DualNum R).re = OfNat.ofNat n := by
  have h : DualNum.inlHom (OfNat.ofNat n : R) = (OfNat.ofNat n : DualNum R) := map_ofNat _ n
  rw [← h, DualNum.inlHom_apply, DualNum.re_mk]

@[simp] theorem DualNum.eps_ofNat (n : ℕ) [n.AtLeastTwo] :
    (OfNat.ofNat n : DualNum R).eps = 0 := by
  have h : DualNum.inlHom (OfNat.ofNat n : R) = (OfNat.ofNat n : DualNum R) := map_ofNat _ n
  rw [← h, DualNum.inlHom_apply, DualNum.eps_mk]

@[simp] theorem DualNum.re_two : (2 : DualNum R).re = 2 := by
  rw [DualNum.two_eq_mk, DualNum.re_mk]

@[simp] theorem DualNum.eps_two : (2 : DualNum R).eps = 0 := by
  rw [DualNum.two_eq_mk, DualNum.eps_mk]

/-- ★★★★★★★★★★**`X(a+εa, w−εw) = X + ε·DX`**。 -/
theorem tateXpair_dual {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateXpair (DualNum.mk a a) (DualNum.mk w (-w)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateXpair a w q hq) (tateDXpair a w q hq) := by
  rw [tateXpair, tateXterm_dual ha, tateXterm_dual hw, tateXtail_dual_pos a q hq,
    tateXtail_dual_neg w q hq, evalAdic_dual _ q hq]
  ext
  · simp [tateXpair]
  · simp [tateDXpair, tateDXterm]
    ring

/-- ★★★★★★★★★★**`Y(a+εa, w−εw) = Y + ε·DY`**。 -/
theorem tateYpair_dual {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateYpair (DualNum.mk a a) (DualNum.mk w (-w)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateYpair a w q hq) (tateDYpair a w q hq) := by
  rw [tateYpair, tateYterm_dual ha, tateXterm_dual hw, tateYterm_dual hw,
    tateYtail_dual_pos a q hq, tateXtail_dual_neg w q hq, tateYtail_dual_neg w q hq,
    evalAdic_dual _ q hq]
  ext
  · simp [tateYpair]
  · simp [tateDYpair, tateDXterm, tateDYterm]
    ring

theorem tateCurveAt_a4_dual {I : Ideal R} [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt (DualNum.mk q 0) (mk_mem_dualIdeal hq)).a₄
      = DualNum.mk ((tateCurveAt q hq).a₄) 0 := by
  rw [tateCurveAt_a₄, tateCurveAt_a₄, evalAdic_dual _ q hq]

theorem tateCurveAt_a6_dual {I : Ideal R} [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt (DualNum.mk q 0) (mk_mem_dualIdeal hq)).a₆
      = DualNum.mk ((tateCurveAt q hq).a₆) 0 := by
  rw [tateCurveAt_a₆, tateCurveAt_a₆, evalAdic_dual _ q hq]

/-! ## ★★★★★★★★★★★★★★★★★★★★ODE（掛けた形） -/

/-- ★★★★★★★★★★★★★★★★★★★★**Tate 曲線の ODE（`DX` を掛けた形）**:

    `DX · (DY − (3X² − Y + a₄)) = 0`

★`tate_equation`（`Y² + XY = X³ + a₄X + a₆`、証明済み）を `R[ε]` の中で
`(a+εa, w−εw, q)` に適用し、`ε` 成分を取ると

    `2Y·DY + X·DY + Y·DX = 3X²·DX + a₄·DX`

となり、`DX = 2Y + X`（第 846）を入れると `DX·(DY + Y − 3X² − a₄) = 0` である。
☆**2 の因子は一切出ない**——`DX` が非零因子であれば割れる。 -/
theorem tate_ode_mul {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair a w q hq
        * (tateDYpair a w q hq
            - (3 * tateXpair a w q hq ^ 2 - tateYpair a w q hq
                + (tateCurveAt q hq).a₄)) = 0 := by
  have haw' : DualNum.mk a a * DualNum.mk w (-w) = DualNum.mk q 0 := dualA_mul_dualW a w q haw
  have heq := tate_equation (I := dualIdeal I) (DualNum.mk a a) (DualNum.mk w (-w))
    (DualNum.mk q 0) (mk_mem_dualIdeal hq) haw'
    (isUnit_one_sub_dualA ha) (isUnit_one_sub_dualW hw)
  rw [tateXpair_dual a w q hq ha hw, tateYpair_dual a w q hq ha hw,
    tateCurveAt_a4_dual q hq, tateCurveAt_a6_dual q hq] at heq
  have heps := congrArg DualNum.eps heq
  simp only [DualNum.eps_add, DualNum.eps_mul, DualNum.eps_sq, DualNum.eps_cube,
    DualNum.re_pow, DualNum.re_mk, DualNum.eps_mk, DualNum.re_mul, DualNum.re_add,
    DualNum.eps_zero, mul_zero, add_zero, zero_add] at heps
  have hDX := tateDXpair_eq a w q hq ha hw
  linear_combination heps + tateDYpair a w q hq * hDX

/-! ## ★★★★★★★★★★微分でべき級数展開を上げる -/

/-- ★`re` を加法準同型として。 -/
def DualNum.reAddHom : DualNum R →+ R where
  toFun := DualNum.re
  map_zero' := DualNum.re_zero
  map_add' := DualNum.re_add

/-- ★`eps` を加法準同型として。 -/
def DualNum.epsAddHom : DualNum R →+ R where
  toFun := DualNum.eps
  map_zero' := DualNum.eps_zero
  map_add' := DualNum.eps_add

theorem DualNum.mk_sum {ι : Type*} (s : Finset ι) (f g : ι → R) :
    ∑ i ∈ s, DualNum.mk (f i) (g i)
      = DualNum.mk (∑ i ∈ s, f i) (∑ i ∈ s, g i) := by
  ext
  · show DualNum.reAddHom (∑ i ∈ s, DualNum.mk (f i) (g i)) = _
    rw [map_sum]
    simp [DualNum.reAddHom]
  · show DualNum.epsAddHom (∑ i ∈ s, DualNum.mk (f i) (g i)) = _
    rw [map_sum]
    simp [DualNum.epsAddHom]

theorem DualNum.mk_pow_self (t : R) (n : ℕ) :
    (DualNum.mk t t) ^ n = DualNum.mk (t ^ n) ((n : R) * t ^ n) := by
  induction n with
  | zero => ext <;> simp
  | succ k ih =>
      rw [pow_succ, ih]
      ext
      · simp [pow_succ]
      · simp [pow_succ]
        push_cast
        ring

theorem DualNum.natCast_eq (n : ℕ) : ((n : ℕ) : DualNum R) = DualNum.mk ((n : ℕ) : R) 0 := by
  rw [← DualNum.inlHom_apply, map_natCast]

/-- ★★★★★★**`Df(t+εt) = Df(t) + ε·D²f(t)`**。 -/
theorem tateDXterm_dual {t : R} (hu : IsUnit (1 - t)) :
    tateDXterm (DualNum.mk t t) = DualNum.mk (tateDXterm t) (tateD2Xterm t) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateDXterm_eq (isUnit_one_sub_dual hu), tateXterm_dual hu, tateYterm_dual hu]
  ext
  · simp [tateDXterm_eq hu]
  · simp only [DualNum.eps_add, DualNum.eps_mul, DualNum.eps_mk, DualNum.re_mk,
      DualNum.re_two, DualNum.eps_two, tateD2Xterm]
    linear_combination (-(t * (1 + t) * Ring.inverse (1 - t) ^ 3)) * hr

/-- ★★★★★★★★**`Df(t) = ∑_{n≥1} n² t^n`**——
`f(t) = ∑ n t^n` を双対数で微分して得る。 -/
theorem tateDXterm_eq_adicSum {I : Ideal R} [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateDXterm t
      = adicSum (fun n => (n : R) ^ 2 * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  have hut : IsUnit (1 - t) := isUnit_one_sub (I := I) ht
  have hmem : ∀ n : ℕ, (n : R) * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have hmem2 : ∀ n : ℕ, (n : R) ^ 2 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have htd' : DualNum.mk t t ∈ dualIdeal I := ⟨by simpa using ht, by simpa using ht⟩
  have h1 := tateXterm_eq_adicSum (I := dualIdeal I) htd'
  have hterm : ∀ n : ℕ, ((n : ℕ) : DualNum R) * (DualNum.mk t t) ^ n
      = DualNum.mk ((n : R) * t ^ n) ((n : R) ^ 2 * t ^ n) := by
    intro n
    rw [DualNum.natCast_eq, DualNum.mk_pow_self, DualNum.mk_zero_mul]
    ext
    · simp
    · simp
      ring
  rw [adicSum_congr _ (dual_mem_pow hmem hmem2) hterm,
    adicSum_dual _ _ hmem hmem2] at h1
  rw [tateXterm_dual hut] at h1
  simpa [tateDXterm] using congrArg DualNum.eps h1

/-- ★★★★★★★★**`D²f(t) = ∑_{n≥1} n³ t^n`**——もう一度微分する。 -/
theorem tateD2Xterm_eq_adicSum {I : Ideal R} [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateD2Xterm t
      = adicSum (fun n => (n : R) ^ 3 * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  have hut : IsUnit (1 - t) := isUnit_one_sub (I := I) ht
  have hmem2 : ∀ n : ℕ, (n : R) ^ 2 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have hmem3 : ∀ n : ℕ, (n : R) ^ 3 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have htd' : DualNum.mk t t ∈ dualIdeal I := ⟨by simpa using ht, by simpa using ht⟩
  have h1 := tateDXterm_eq_adicSum (I := dualIdeal I) htd'
  have hterm : ∀ n : ℕ, ((n : ℕ) : DualNum R) ^ 2 * (DualNum.mk t t) ^ n
      = DualNum.mk ((n : R) ^ 2 * t ^ n) ((n : R) ^ 3 * t ^ n) := by
    intro n
    rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
    ext
    · simp
    · simp
      ring
  rw [adicSum_congr _ (dual_mem_pow hmem2 hmem3) hterm,
    adicSum_dual _ _ hmem2 hmem3] at h1
  rw [tateDXterm_dual hut] at h1
  simpa using congrArg DualNum.eps h1

/-! ## ★★★★★★★★尾のべき級数展開も微分で上げる -/

theorem tateDXtail_dual_pos {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail (DualNum.mk u u) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateDXtail u q hq) (tateD2Xtail u q hq) := by
  have hterm : ∀ n : ℕ,
      tateDXterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u u)
        = DualNum.mk (tateDXterm (q ^ (n + 1) * u)) (tateD2Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateDXterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n))]
  rw [tateDXtail, tateDXtail, tateD2Xtail,
    adicSum_congr (tateDXtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateDXtail_aux hq) (tateD2Xtail_aux hq)) hterm,
    adicSum_dual]

/-- ★★★★★★**`T_{Df}(u) = ∑_n q^n ∑_{d∣n} d² u^d`**。 -/
theorem tateDXtail_eq_divisorSum {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 2 * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  have hmem1 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have hmem2 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 2 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have h1 := tateXtail_eq_divisorSum (I := dualIdeal I) (DualNum.mk u u)
    (DualNum.mk q 0) (mk_mem_dualIdeal hq)
  have hterm : ∀ n : ℕ,
      (DualNum.mk q 0) ^ n * ∑ d ∈ n.divisors, ((d : ℕ) : DualNum R) * (DualNum.mk u u) ^ d
        = DualNum.mk (q ^ n * ∑ d ∈ n.divisors, (d : R) * u ^ d)
            (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 2 * u ^ d) := by
    intro n
    have hin : ∀ d ∈ n.divisors, ((d : ℕ) : DualNum R) * (DualNum.mk u u) ^ d
        = DualNum.mk ((d : R) * u ^ d) ((d : R) ^ 2 * u ^ d) := by
      intro d _
      rw [DualNum.natCast_eq, DualNum.mk_pow_self, DualNum.mk_zero_mul]
      ext
      · simp
      · simp
        ring
    rw [Finset.sum_congr rfl hin, DualNum.mk_pow]
    rw [DualNum.mk_sum, DualNum.mk_zero_mul]
  rw [adicSum_congr _ (dual_mem_pow hmem1 hmem2) hterm, adicSum_dual _ _ hmem1 hmem2,
    tateXtail_dual_pos u q hq] at h1
  simpa using congrArg DualNum.eps h1

/-- ★★★★★★★★**`T_{D²f}(u) = ∑_n q^n ∑_{d∣n} d³ u^d`**——σ₃ が直に出る形。 -/
theorem tateD2Xtail_eq_divisorSum {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD2Xtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  have hmem2 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 2 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have hmem3 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have h1 := tateDXtail_eq_divisorSum (I := dualIdeal I) (DualNum.mk u u)
    (DualNum.mk q 0) (mk_mem_dualIdeal hq)
  have hterm : ∀ n : ℕ,
      (DualNum.mk q 0) ^ n * ∑ d ∈ n.divisors,
          ((d : ℕ) : DualNum R) ^ 2 * (DualNum.mk u u) ^ d
        = DualNum.mk (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 2 * u ^ d)
            (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * u ^ d) := by
    intro n
    have hin : ∀ d ∈ n.divisors, ((d : ℕ) : DualNum R) ^ 2 * (DualNum.mk u u) ^ d
        = DualNum.mk ((d : R) ^ 2 * u ^ d) ((d : R) ^ 3 * u ^ d) := by
      intro d _
      rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
      ext
      · simp
      · simp
        ring
    rw [Finset.sum_congr rfl hin, DualNum.mk_pow]
    rw [DualNum.mk_sum, DualNum.mk_zero_mul]
  rw [adicSum_congr _ (dual_mem_pow hmem2 hmem3) hterm, adicSum_dual _ _ hmem2 hmem3,
    tateDXtail_dual_pos u q hq] at h1
  simpa using congrArg DualNum.eps h1

/-! ## ★★★★★★★★`D²X` を微分して `D³X` -/

/-- ★★★★★★**`D²f(t + εs) = D²f(t) + ε·s(1+11t+11t²+t³)/(1−t)⁵`**。 -/
theorem tateD2Xterm_dual {t s : R} (hu : IsUnit (1 - t)) :
    tateD2Xterm (DualNum.mk t s)
      = DualNum.mk (tateD2Xterm t)
          (s * (1 + 11 * t + 11 * t ^ 2 + t ^ 3) * Ring.inverse (1 - t) ^ 5) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  have hn4 : (4 : DualNum R) = DualNum.mk 4 0 := by
    have h : DualNum.inlHom (4 : R) = (4 : DualNum R) := map_ofNat _ 4
    rw [← h, DualNum.inlHom_apply]
  rw [tateD2Xterm, ringInverse_one_sub_dual hu, hn4]
  ext
  · simp [tateD2Xterm]
  · simp only [DualNum.eps_mul, DualNum.re_mul, DualNum.re_mk, DualNum.eps_mk,
      DualNum.re_pow, DualNum.eps_sq, DualNum.re_add, DualNum.eps_add, DualNum.re_one,
      DualNum.eps_one, DualNum.re_ofNat, DualNum.eps_ofNat]
    have h4 : ((DualNum.mk (Ring.inverse (1 - t)) (s * Ring.inverse (1 - t) ^ 2)) ^ 4).eps
        = 4 * Ring.inverse (1 - t) ^ 3 * (s * Ring.inverse (1 - t) ^ 2) := by
      have e := DualNum.eps_pow
        (DualNum.mk (Ring.inverse (1 - t)) (s * Ring.inverse (1 - t) ^ 2)) 3
      simp only [DualNum.re_mk, DualNum.eps_mk] at e
      rw [e]
      push_cast
      ring
    rw [h4]
    linear_combination
      (-(s * (1 + 8 * t + 3 * t ^ 2) * Ring.inverse (1 - t) ^ 4)) * hr

theorem tateD2Xtail_dual_pos {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD2Xtail (DualNum.mk u u) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD2Xtail u q hq) (tateD3Xtail u q hq) := by
  have hterm : ∀ n : ℕ,
      tateD2Xterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u u)
        = DualNum.mk (tateD2Xterm (q ^ (n + 1) * u)) (tateD3Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateD2Xterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateD3Xterm]
  rw [tateD2Xtail, tateD2Xtail, tateD3Xtail,
    adicSum_congr (tateD2Xtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateD2Xtail_aux hq) (tateD3Xtail_aux hq)) hterm,
    adicSum_dual]

theorem tateD2Xtail_dual_neg {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD2Xtail (DualNum.mk u (-u)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD2Xtail u q hq) (-tateD3Xtail u q hq) := by
  have hmem : ∀ n : ℕ, (-1 : R) * tateD3Xterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateD3Xtail_aux hq n)
  have hterm : ∀ n : ℕ,
      tateD2Xterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u (-u))
        = DualNum.mk (tateD2Xterm (q ^ (n + 1) * u))
            ((-1 : R) * tateD3Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateD2Xterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateD3Xterm]
    ext <;> simp <;> try ring
  rw [tateD2Xtail, tateD2Xtail, tateD3Xtail,
    adicSum_congr (tateD2Xtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateD2Xtail_aux hq) hmem) hterm,
    adicSum_dual _ _ (tateD2Xtail_aux hq) hmem, adicSum_smul (-1 : R) _ (tateD3Xtail_aux hq)]
  ext <;> simp

/-- ★★★★★★★★**`D²X(a+εa, w−εw) = D²X + ε·D³X`**。 -/
theorem tateD2Xpair_dual {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateD2Xpair (DualNum.mk a a) (DualNum.mk w (-w)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD2Xpair a w q hq) (tateD3Xpair a w q hq) := by
  rw [tateD2Xpair, tateD2Xterm_dual ha, tateD2Xterm_dual hw,
    tateD2Xtail_dual_pos a q hq, tateD2Xtail_dual_neg w q hq]
  ext
  · simp [tateD2Xpair, tateD2Xterm]
  · simp [tateD3Xpair, tateD3Xterm]
    ring

/-! ## ★★★★★★★★`DX` を微分して `D²X` -/

/-- ★★★★**`Df(t + εs) = Df(t) + ε·s(1+4t+t²)/(1−t)⁴`**（一般の `s`）。 -/
theorem tateDXterm_dual' {t s : R} (hu : IsUnit (1 - t)) :
    tateDXterm (DualNum.mk t s)
      = DualNum.mk (tateDXterm t)
          (s * (1 + 4 * t + t ^ 2) * Ring.inverse (1 - t) ^ 4) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  rw [tateDXterm_eq (isUnit_one_sub_dual hu), tateXterm_dual hu, tateYterm_dual hu]
  ext
  · simp [tateDXterm_eq hu]
  · simp only [DualNum.eps_add, DualNum.eps_mul, DualNum.eps_mk, DualNum.re_mk,
      DualNum.re_two, DualNum.eps_two]
    linear_combination (-(s * (1 + t) * Ring.inverse (1 - t) ^ 3)) * hr

theorem tateDXtail_dual_neg {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail (DualNum.mk u (-u)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateDXtail u q hq) (-tateD2Xtail u q hq) := by
  have hmem : ∀ n : ℕ, (-1 : R) * tateD2Xterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateD2Xtail_aux hq n)
  have hterm : ∀ n : ℕ,
      tateDXterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u (-u))
        = DualNum.mk (tateDXterm (q ^ (n + 1) * u))
            ((-1 : R) * tateD2Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateDXterm_dual' (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateD2Xterm]
    ext <;> simp <;> try ring
  rw [tateDXtail, tateDXtail, tateD2Xtail,
    adicSum_congr (tateDXtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateDXtail_aux hq) hmem) hterm,
    adicSum_dual _ _ (tateDXtail_aux hq) hmem,
    adicSum_smul (-1 : R) _ (tateD2Xtail_aux hq)]
  ext <;> simp

/-- ★★★★★★★★**`DX(a+εa, w−εw) = DX + ε·D²X`**。 -/
theorem tateDXpair_dual {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair (DualNum.mk a a) (DualNum.mk w (-w)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateDXpair a w q hq) (tateD2Xpair a w q hq) := by
  rw [tateDXpair, tateDXterm_dual' ha, tateDXterm_dual' hw,
    tateDXtail_dual_pos a q hq, tateDXtail_dual_neg w q hq]
  ext
  · simp [tateDXpair, tateDXterm]
  · simp [tateD2Xpair, tateD2Xterm]
    ring

/-! ## ★★★★★★★★★★`D²X` と `D³X` の掛けた形 -/

/-- ★★★★★★**`DX·(D²X − (6X² + X + 2a₄)) = 0`**——`tate_ode_mul` の書き換え。 -/
theorem tate_d2x_mul {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair a w q hq
        * (tateD2Xpair a w q hq
            - (6 * tateXpair a w q hq ^ 2 + tateXpair a w q hq
                + 2 * (tateCurveAt q hq).a₄)) = 0 := by
  have h := tate_ode_mul a w q hq haw ha hw
  rw [tateDXpair_eq a w q hq ha hw] at h
  rw [tateD2Xpair_eq a w q hq ha hw, tateDXpair_eq a w q hq ha hw]
  linear_combination 2 * h

/-- ★★★★★★★★**`DX·(D³X − 12X·DX − DX) + D²X·(D²X − 6X² − X − 2a₄) = 0`**。

★`tate_d2x_mul` を双対数の点で評価し、`ε` 成分を取ったもの。 -/
theorem tate_d3x_mul {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair a w q hq
        * (tateD3Xpair a w q hq - 12 * tateXpair a w q hq * tateDXpair a w q hq
            - tateDXpair a w q hq)
      + tateD2Xpair a w q hq
        * (tateD2Xpair a w q hq - 6 * tateXpair a w q hq ^ 2 - tateXpair a w q hq
            - 2 * (tateCurveAt q hq).a₄) = 0 := by
  have haw' : DualNum.mk a a * DualNum.mk w (-w) = DualNum.mk q 0 := dualA_mul_dualW a w q haw
  have h := tate_d2x_mul (I := dualIdeal I) (DualNum.mk a a) (DualNum.mk w (-w))
    (DualNum.mk q 0) (mk_mem_dualIdeal hq) haw'
    (isUnit_one_sub_dualA ha) (isUnit_one_sub_dualW hw)
  have hn6 : (6 : DualNum R) = DualNum.mk 6 0 := by
    have hh : DualNum.inlHom (6 : R) = (6 : DualNum R) := map_ofNat _ 6
    rw [← hh, DualNum.inlHom_apply]
  rw [tateDXpair_dual a w q hq ha hw, tateD2Xpair_dual a w q hq ha hw,
    tateXpair_dual a w q hq ha hw, tateCurveAt_a4_dual q hq, hn6] at h
  have heps := congrArg DualNum.eps h
  simp only [DualNum.eps_mul, DualNum.eps_sub, DualNum.eps_add, DualNum.eps_mk,
    DualNum.re_mk, DualNum.re_mul, DualNum.re_sub, DualNum.re_add, DualNum.eps_sq,
    DualNum.re_pow, DualNum.eps_zero, DualNum.re_two, DualNum.eps_two] at heps
  linear_combination heps

/-! ## ★★★★★★★★`D³X` を微分して `D⁴X` -/

/-- ★★★★★★**`D³f(t+εs) = D³f(t) + ε·s(1+26t+66t²+26t³+t⁴)/(1−t)⁶`**。 -/
theorem tateD3Xterm_dual {t s : R} (hu : IsUnit (1 - t)) :
    tateD3Xterm (DualNum.mk t s)
      = DualNum.mk (tateD3Xterm t)
          (s * (1 + 26 * t + 66 * t ^ 2 + 26 * t ^ 3 + t ^ 4)
            * Ring.inverse (1 - t) ^ 6) := by
  have hr : (1 - t) * Ring.inverse (1 - t) = 1 := Ring.mul_inverse_cancel _ hu
  have hn11 : (11 : DualNum R) = DualNum.mk 11 0 := by
    have hh : DualNum.inlHom (11 : R) = (11 : DualNum R) := map_ofNat _ 11
    rw [← hh, DualNum.inlHom_apply]
  rw [tateD3Xterm, ringInverse_one_sub_dual hu, hn11]
  ext
  · simp [tateD3Xterm]
  · simp only [DualNum.eps_mul, DualNum.re_mul, DualNum.re_mk, DualNum.eps_mk,
      DualNum.re_pow, DualNum.eps_sq, DualNum.eps_cube, DualNum.re_add, DualNum.eps_add,
      DualNum.re_one, DualNum.eps_one]
    have h5 : ((DualNum.mk (Ring.inverse (1 - t)) (s * Ring.inverse (1 - t) ^ 2)) ^ 5).eps
        = 5 * Ring.inverse (1 - t) ^ 4 * (s * Ring.inverse (1 - t) ^ 2) := by
      have e := DualNum.eps_pow
        (DualNum.mk (Ring.inverse (1 - t)) (s * Ring.inverse (1 - t) ^ 2)) 4
      simp only [DualNum.re_mk, DualNum.eps_mk] at e
      rw [e]
      push_cast
      ring
    rw [h5]
    linear_combination
      (-(s * (1 + 22 * t + 33 * t ^ 2 + 4 * t ^ 3) * Ring.inverse (1 - t) ^ 5)) * hr

theorem tateD3Xtail_dual_pos {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD3Xtail (DualNum.mk u u) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD3Xtail u q hq) (tateD4Xtail u q hq) := by
  have hterm : ∀ n : ℕ,
      tateD3Xterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u u)
        = DualNum.mk (tateD3Xterm (q ^ (n + 1) * u)) (tateD4Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateD3Xterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateD4Xterm]
  rw [tateD3Xtail, tateD3Xtail, tateD4Xtail,
    adicSum_congr (tateD3Xtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateD3Xtail_aux hq) (tateD4Xtail_aux hq)) hterm,
    adicSum_dual]

theorem tateD3Xtail_dual_neg {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD3Xtail (DualNum.mk u (-u)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD3Xtail u q hq) (-tateD4Xtail u q hq) := by
  have hmem : ∀ n : ℕ, (-1 : R) * tateD4Xterm (q ^ (n + 1) * u) ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (tateD4Xtail_aux hq n)
  have hterm : ∀ n : ℕ,
      tateD3Xterm ((DualNum.mk q 0) ^ (n + 1) * DualNum.mk u (-u))
        = DualNum.mk (tateD3Xterm (q ^ (n + 1) * u))
            ((-1 : R) * tateD4Xterm (q ^ (n + 1) * u)) := by
    intro n
    rw [DualNum.mk_pow, DualNum.mk_zero_mul,
      tateD3Xterm_dual (isUnit_one_sub (I := I) (pow_succ_mul_mem_I hq n)), tateD4Xterm]
    ext <;> simp <;> try ring
  rw [tateD3Xtail, tateD3Xtail, tateD4Xtail,
    adicSum_congr (tateD3Xtail_aux (mk_mem_dualIdeal hq))
      (dual_mem_pow (tateD3Xtail_aux hq) hmem) hterm,
    adicSum_dual _ _ (tateD3Xtail_aux hq) hmem,
    adicSum_smul (-1 : R) _ (tateD4Xtail_aux hq)]
  ext <;> simp

/-- ★★★★★★★★**`D³X(a+εa, w−εw) = D³X + ε·D⁴X`**。 -/
theorem tateD3Xpair_dual {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateD3Xpair (DualNum.mk a a) (DualNum.mk w (-w)) (DualNum.mk q 0) (mk_mem_dualIdeal hq)
      = DualNum.mk (tateD3Xpair a w q hq) (tateD4Xpair a w q hq) := by
  rw [tateD3Xpair, tateD3Xterm_dual ha, tateD3Xterm_dual hw,
    tateD3Xtail_dual_pos a q hq, tateD3Xtail_dual_neg w q hq]
  ext
  · simp [tateD3Xpair, tateD3Xterm]
  · simp [tateD4Xpair, tateD4Xterm]
    ring

/-- ★★★★★★★★**`tate_d3x_mul` をもう一度微分した形**。

★`tate_d2x` と `tate_d3x` を入れると第 2・第 3 項が消え、
`DX ≠ 0` で割ると **`D⁴X = 12X·D²X + 12(DX)² + D²X`** が出る。 -/
theorem tate_d4x_mul {I : Ideal R} [IsAdicComplete I R] (a w q : R) (hq : q ∈ I)
    (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDXpair a w q hq
        * (tateD4Xpair a w q hq - 12 * tateXpair a w q hq * tateD2Xpair a w q hq
            - 12 * tateDXpair a w q hq ^ 2 - tateD2Xpair a w q hq)
      + 2 * tateD2Xpair a w q hq
        * (tateD3Xpair a w q hq - 12 * tateXpair a w q hq * tateDXpair a w q hq
            - tateDXpair a w q hq)
      + tateD3Xpair a w q hq
        * (tateD2Xpair a w q hq - 6 * tateXpair a w q hq ^ 2 - tateXpair a w q hq
            - 2 * (tateCurveAt q hq).a₄) = 0 := by
  have haw' : DualNum.mk a a * DualNum.mk w (-w) = DualNum.mk q 0 := dualA_mul_dualW a w q haw
  have h := tate_d3x_mul (I := dualIdeal I) (DualNum.mk a a) (DualNum.mk w (-w))
    (DualNum.mk q 0) (mk_mem_dualIdeal hq) haw'
    (isUnit_one_sub_dualA ha) (isUnit_one_sub_dualW hw)
  have hn6 : (6 : DualNum R) = DualNum.mk 6 0 := by
    have hh : DualNum.inlHom (6 : R) = (6 : DualNum R) := map_ofNat _ 6
    rw [← hh, DualNum.inlHom_apply]
  have hn12 : (12 : DualNum R) = DualNum.mk 12 0 := by
    have hh : DualNum.inlHom (12 : R) = (12 : DualNum R) := map_ofNat _ 12
    rw [← hh, DualNum.inlHom_apply]
  rw [tateDXpair_dual a w q hq ha hw, tateD2Xpair_dual a w q hq ha hw,
    tateD3Xpair_dual a w q hq ha hw, tateXpair_dual a w q hq ha hw,
    tateCurveAt_a4_dual q hq, hn6, hn12] at h
  have heps := congrArg DualNum.eps h
  simp only [DualNum.eps_mul, DualNum.eps_sub, DualNum.eps_add, DualNum.eps_mk,
    DualNum.re_mk, DualNum.re_mul, DualNum.re_sub, DualNum.re_add, DualNum.eps_sq,
    DualNum.re_pow, DualNum.eps_zero, DualNum.re_two, DualNum.eps_two] at heps
  linear_combination heps

/-! ## ★★★★★★`D³f`・`D⁴f` のべき級数展開 -/

/-- ★★★★**`D³f(t) = ∑_{n≥1} n⁴ t^n`**。 -/
theorem tateD3Xterm_eq_adicSum {I : Ideal R} [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateD3Xterm t
      = adicSum (fun n => (n : R) ^ 4 * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  have hut : IsUnit (1 - t) := isUnit_one_sub (I := I) ht
  have hmem3 : ∀ n : ℕ, (n : R) ^ 3 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have hmem4 : ∀ n : ℕ, (n : R) ^ 4 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have htd' : DualNum.mk t t ∈ dualIdeal I := ⟨by simpa using ht, by simpa using ht⟩
  have h1 := tateD2Xterm_eq_adicSum (I := dualIdeal I) htd'
  have hterm : ∀ n : ℕ, ((n : ℕ) : DualNum R) ^ 3 * (DualNum.mk t t) ^ n
      = DualNum.mk ((n : R) ^ 3 * t ^ n) ((n : R) ^ 4 * t ^ n) := by
    intro n
    rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
    ext
    · simp
    · simp
      ring
  rw [adicSum_congr _ (dual_mem_pow hmem3 hmem4) hterm,
    adicSum_dual _ _ hmem3 hmem4] at h1
  rw [tateD2Xterm_dual hut] at h1
  have h2 := congrArg DualNum.eps h1
  simp only [DualNum.eps_mk] at h2
  rw [← h2, tateD3Xterm]

/-- ★★★★**`D⁴f(t) = ∑_{n≥1} n⁵ t^n`**。 -/
theorem tateD4Xterm_eq_adicSum {I : Ideal R} [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    tateD4Xterm t
      = adicSum (fun n => (n : R) ^ 5 * t ^ n)
          (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)) := by
  have hut : IsUnit (1 - t) := isUnit_one_sub (I := I) ht
  have hmem4 : ∀ n : ℕ, (n : R) ^ 4 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have hmem5 : ∀ n : ℕ, (n : R) ^ 5 * t ^ n ∈ I ^ n :=
    fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow ht n)
  have htd' : DualNum.mk t t ∈ dualIdeal I := ⟨by simpa using ht, by simpa using ht⟩
  have h1 := tateD3Xterm_eq_adicSum (I := dualIdeal I) htd'
  have hterm : ∀ n : ℕ, ((n : ℕ) : DualNum R) ^ 4 * (DualNum.mk t t) ^ n
      = DualNum.mk ((n : R) ^ 4 * t ^ n) ((n : R) ^ 5 * t ^ n) := by
    intro n
    rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
    ext
    · simp
    · simp
      ring
  rw [adicSum_congr _ (dual_mem_pow hmem4 hmem5) hterm,
    adicSum_dual _ _ hmem4 hmem5] at h1
  rw [tateD3Xterm_dual hut] at h1
  have h2 := congrArg DualNum.eps h1
  simp only [DualNum.eps_mk] at h2
  rw [← h2, tateD4Xterm]

/-- ★★★★**`T_{D³f}(u) = ∑_n q^n ∑_{d∣n} d⁴ u^d`**。 -/
theorem tateD3Xtail_eq_divisorSum {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD3Xtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 4 * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  have hmem3 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have hmem4 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 4 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have h1 := tateD2Xtail_eq_divisorSum (I := dualIdeal I) (DualNum.mk u u)
    (DualNum.mk q 0) (mk_mem_dualIdeal hq)
  have hterm : ∀ n : ℕ,
      (DualNum.mk q 0) ^ n * ∑ d ∈ n.divisors,
          ((d : ℕ) : DualNum R) ^ 3 * (DualNum.mk u u) ^ d
        = DualNum.mk (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * u ^ d)
            (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 4 * u ^ d) := by
    intro n
    have hin : ∀ d ∈ n.divisors, ((d : ℕ) : DualNum R) ^ 3 * (DualNum.mk u u) ^ d
        = DualNum.mk ((d : R) ^ 3 * u ^ d) ((d : R) ^ 4 * u ^ d) := by
      intro d _
      rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
      ext
      · simp
      · simp
        ring
    rw [Finset.sum_congr rfl hin, DualNum.mk_pow, DualNum.mk_sum, DualNum.mk_zero_mul]
  rw [adicSum_congr _ (dual_mem_pow hmem3 hmem4) hterm, adicSum_dual _ _ hmem3 hmem4,
    tateD2Xtail_dual_pos u q hq] at h1
  simpa using congrArg DualNum.eps h1

/-- ★★★★★★**`T_{D⁴f}(u) = ∑_n q^n ∑_{d∣n} d⁵ u^d`**——σ₅ が直に出る形。 -/
theorem tateD4Xtail_eq_divisorSum {I : Ideal R} [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateD4Xtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 5 * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  have hmem4 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 4 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have hmem5 : ∀ n : ℕ, q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 5 * u ^ d ∈ I ^ n :=
    fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have h1 := tateD3Xtail_eq_divisorSum (I := dualIdeal I) (DualNum.mk u u)
    (DualNum.mk q 0) (mk_mem_dualIdeal hq)
  have hterm : ∀ n : ℕ,
      (DualNum.mk q 0) ^ n * ∑ d ∈ n.divisors,
          ((d : ℕ) : DualNum R) ^ 4 * (DualNum.mk u u) ^ d
        = DualNum.mk (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 4 * u ^ d)
            (q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 5 * u ^ d) := by
    intro n
    have hin : ∀ d ∈ n.divisors, ((d : ℕ) : DualNum R) ^ 4 * (DualNum.mk u u) ^ d
        = DualNum.mk ((d : R) ^ 4 * u ^ d) ((d : R) ^ 5 * u ^ d) := by
      intro d _
      rw [DualNum.natCast_eq, DualNum.mk_pow, DualNum.mk_pow_self, DualNum.mk_zero_mul]
      ext
      · simp
      · simp
        ring
    rw [Finset.sum_congr rfl hin, DualNum.mk_pow, DualNum.mk_sum, DualNum.mk_zero_mul]
  rw [adicSum_congr _ (dual_mem_pow hmem4 hmem5) hterm, adicSum_dual _ _ hmem4 hmem5,
    tateD3Xtail_dual_pos u q hq] at h1
  simpa using congrArg DualNum.eps h1

/-! ## ★★★★`c₄ = 1 + 240·s₃` -/

/-- ★★**べき級数の水準で `c₄ = 1 + 240·σ₃`**。

`b₂ = a₁² + 4a₂ = 1`、`b₄ = 2a₄ + a₁a₃ = 2a₄` なので
`c₄ = b₂² − 24b₄ = 1 − 48a₄ = 1 + 240 s₃`（`a₄ = −5s₃`）。 -/
theorem tateCurve_c4 : tateCurve.c₄ = 1 + 240 * sigmaSeries 3 := by
  rw [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  simp only [tateCurve, tateA4]
  have hC : (PowerSeries.C (-5 : ℤ)) = (-5 : PowerSeries ℤ) := by
    rw [show (-5 : ℤ) = -(5 : ℤ) from rfl, map_neg, map_ofNat]
  rw [hC]
  ring

/-- ★★★★**`c₄(E_q) = 1 + 240·s₃(q)`**。 -/
theorem tateCurveAt_c4_eq {I : Ideal R} [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).c₄ = 1 + 240 * evalAdic (sigmaSeries 3) q hq := by
  rw [tateCurveAt_c4, tateCurve_c4, map_add, map_mul, map_one,
    map_ofNat (evalAdicHom q hq) 240]
  rfl


theorem tateCurveAt_c6 {I : Ideal R} [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).c₆ = evalAdicHom q hq tateCurve.c₆ := by
  rw [tateCurveAt, WeierstrassCurve.map_c₆]

/-- ★★**べき級数の水準で `c₆ = −1 + 504·σ₅`**。

`b₂ = 1`、`b₄ = 2a₄`、`b₆ = 4a₆` なので `c₆ = −1 + 72a₄ − 864a₆` であり、
`864a₆ = 72·(12a₆) = −72(5s₃ + 7s₅)`（`twelve_mul_tateA6`）を入れると
`c₆ = −1 − 360s₃ + 360s₃ + 504s₅ = −1 + 504s₅`。★**割り算は一切要らない**。 -/
theorem tateCurve_c6 : tateCurve.c₆ = -1 + 504 * sigmaSeries 5 := by
  rw [WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆]
  simp only [tateCurve, tateA4]
  have hC : (PowerSeries.C (-5 : ℤ)) = (-5 : PowerSeries ℤ) := by
    rw [show (-5 : ℤ) = -(5 : ℤ) from rfl, map_neg, map_ofNat]
  have hC12 : (PowerSeries.C (12 : ℤ)) = (12 : PowerSeries ℤ) := map_ofNat _ 12
  have hC5 : (PowerSeries.C (5 : ℤ)) = (5 : PowerSeries ℤ) := map_ofNat _ 5
  have hC7 : (PowerSeries.C (7 : ℤ)) = (7 : PowerSeries ℤ) := map_ofNat _ 7
  have h12 := twelve_mul_tateA6
  rw [hC12, hC5, hC7] at h12
  rw [hC]
  linear_combination (-72 : PowerSeries ℤ) * h12

/-- ★★★★**`c₆(E_q) = −1 + 504·s₅(q)`**。 -/
theorem tateCurveAt_c6_eq {I : Ideal R} [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).c₆ = -1 + 504 * evalAdic (sigmaSeries 5) q hq := by
  rw [tateCurveAt_c6, tateCurve_c6, map_add, map_mul, map_neg, map_one,
    map_ofNat (evalAdicHom q hq) 504]
  rfl

end ABC3.Found.GaloisRep
