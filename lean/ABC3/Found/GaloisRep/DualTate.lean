/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DualAdic
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Found.GaloisRep.TateVelu

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

theorem eps_cube (x : DualNum R) : (x ^ 3).eps = 3 * x.re ^ 2 * x.eps := by
  have h3 : (3 : ℕ) = 2 + 1 := rfl
  rw [h3, pow_succ, eps_mul, eps_sq, re_pow]
  ring

end DualNum

/-! ## ★★★★★逆元 -/

theorem ring_inverse_eq_of_mul_eq_one {A : Type} [CommRing A] {x y : A} (hu : IsUnit x)
    (h : x * y = 1) : Ring.inverse x = y := by
  have h1 : Ring.inverse x * x = 1 := Ring.inverse_mul_cancel x hu
  calc Ring.inverse x = Ring.inverse x * (x * y) := by rw [h, mul_one]
    _ = Ring.inverse x * x * y := by ring
    _ = y := by rw [h1, one_mul]

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

end ABC3.Found.GaloisRep
