import ABC3.Found.GaloisRep.AdicEval
import Mathlib.RingTheory.PowerSeries.Trunc

/-!
# Galois (G6) 第 96 ブロック —— **★★★★完備環上の Tate 曲線**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★評価は環準同型である

第 95 ブロックで `ℤ⟦q⟧ → R`(完備 adic 環)の評価を作った。
★本ブロックでそれが**環準同型**であることを示す。

★★鍵は**打ち切り**(`PowerSeries.trunc`):

    部分和 = 打ち切り多項式の評価
    trunc n (f·g) と (trunc n f)(trunc n g) は **n 次未満で一致**
      ⟹ 差は `q^n` で割れる ⟹ `I^n` に入る

★★★mathlib の `trunc_trunc_mul_trunc` がその一致を与える。

## ★★★★★これで `E_q` が任意の完備 adic 環の上に立つ

    tateCurveAt q : WeierstrassCurve R    (q ∈ I, R は I-adic 完備)

★(G6) が要求する **Tate 一意化** `E(K) ≅ Kˣ/q^ℤ` は、この `E_q` の上で述べられる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `partialEval_eq_eval_trunc` | ★部分和 = 打ち切りの評価 |
| `eval_mem_pow_of_coeff_zero` | ★低次が消える多項式は `I^n` に値を取る |
| `partialEval_mul_smodEq` | ★★★部分和の乗法性(法 `I^n`) |
| `evalAdic_mul` / `evalAdic_C` / `evalAdic_one` | ★★環準同型の各条件 |
| `evalAdicHom` | ★★★**`ℤ⟦q⟧ →+* R`** |
| `tateCurveAt` | ★★★★**完備環上の Tate 曲線 `E_q`** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★打ち切り -/

theorem partialEval_eq_eval_trunc (f : PowerSeries ℤ) (q : R) (n : ℕ) :
    partialEval f q n = (PowerSeries.trunc n f).eval₂ (Int.castRingHom R) q := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    have h0 : PowerSeries.trunc 0 f = 0 := by
      ext m
      rw [PowerSeries.coeff_trunc]
      simp
    rw [h0, partialEval]
    simp
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    rw [Polynomial.eval₂_eq_sum_range' (Int.castRingHom R) (PowerSeries.natDegree_trunc_lt f m) q,
      partialEval]
    refine Finset.sum_congr rfl (fun k hk => ?_)
    rw [PowerSeries.coeff_trunc, if_pos (Finset.mem_range.1 hk)]
    rfl

/-- ★低次の係数が消える多項式は `q ∈ I` で `I^n` に値を取る。 -/
theorem eval_mem_pow_of_coeff_zero (p : Polynomial ℤ) (n : ℕ) (hp : ∀ k, k < n → p.coeff k = 0)
    {q : R} (hq : q ∈ I) : p.eval₂ (Int.castRingHom R) q ∈ I ^ n := by
  rw [Polynomial.eval₂_eq_sum_range' (Int.castRingHom R) (Nat.lt_succ_self p.natDegree) q]
  refine Submodule.sum_mem _ (fun k _ => ?_)
  rcases lt_or_ge k n with hk | hk
  · rw [hp k hk]
    simp
  · have hpow : q ^ k = q ^ n * q ^ (k - n) := by
      rw [← pow_add]
      congr 1
      omega
    rw [hpow, ← mul_assoc]
    exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n))

/-! ## ★★★乗法性 -/

/-- ★★★部分和は法 `I^n` で乗法的である。 -/
theorem partialEval_mul_smodEq (f g : PowerSeries ℤ) {q : R} (hq : q ∈ I) (n : ℕ) :
    partialEval (f * g) q n ≡ partialEval f q n * partialEval g q n
      [SMOD (I ^ n • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem, partialEval_eq_eval_trunc, partialEval_eq_eval_trunc,
    partialEval_eq_eval_trunc, ← Polynomial.eval₂_mul, ← Polynomial.eval₂_sub]
  have hcoeff : ∀ k, k < n →
      (PowerSeries.trunc n (f * g)
        - (PowerSeries.trunc n f) * (PowerSeries.trunc n g)).coeff k = 0 := by
    intro k hk
    rw [Polynomial.coeff_sub, PowerSeries.coeff_trunc, if_pos hk]
    have h1 := PowerSeries.trunc_trunc_mul_trunc (n := n) f g
    have h2 := congrArg (fun p : Polynomial ℤ => p.coeff k) h1
    simp only [PowerSeries.coeff_trunc, if_pos hk] at h2
    have h3 : Polynomial.coeff ((PowerSeries.trunc n f) * (PowerSeries.trunc n g)) k
        = PowerSeries.coeff k
            ((((PowerSeries.trunc n f) : Polynomial ℤ) : PowerSeries ℤ)
              * (((PowerSeries.trunc n g) : Polynomial ℤ) : PowerSeries ℤ)) := by
      rw [← Polynomial.coe_mul, Polynomial.coeff_coe]
    rw [h3, h2]
    ring
  have hmem := eval_mem_pow_of_coeff_zero
    (PowerSeries.trunc n (f * g) - (PowerSeries.trunc n f) * (PowerSeries.trunc n g)) n hcoeff hq
  simpa using hmem

/-- ★★★**評価は乗法的である**。 -/
theorem evalAdic_mul [IsAdicComplete I R] (f g : PowerSeries ℤ) (q : R) (hq : q ∈ I) :
    evalAdic (f * g) q hq = evalAdic f q hq * evalAdic g q hq := by
  refine evalAdic_unique _ _ _ _ (fun n => ?_)
  refine (partialEval_mul_smodEq f g hq n).trans ?_
  have ha := evalAdic_spec f q hq n
  have hb := evalAdic_spec g q hq n
  rw [SModEq.sub_mem] at ha hb ⊢
  have hexp : partialEval f q n * partialEval g q n - evalAdic f q hq * evalAdic g q hq
      = (partialEval f q n - evalAdic f q hq) * partialEval g q n
        + evalAdic f q hq * (partialEval g q n - evalAdic g q hq) := by ring
  rw [hexp]
  refine Submodule.add_mem _ ?_ ?_
  · simpa using Ideal.mul_mem_right _ _ (by simpa using ha)
  · simpa using Ideal.mul_mem_left _ _ (by simpa using hb)

/-- ★定数の値。 -/
theorem evalAdic_C [IsAdicComplete I R] (c : ℤ) (q : R) (hq : q ∈ I) :
    evalAdic (PowerSeries.C c) q hq = (c : R) := by
  refine evalAdic_unique _ _ _ _ (fun n => ?_)
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    have htop : (I ^ 0 • (⊤ : Submodule R R)) = ⊤ := by simp
    rw [SModEq.sub_mem, htop]
    exact Submodule.mem_top
  · have hp : partialEval (PowerSeries.C c) q n = (c : R) := by
      obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
      rw [partialEval, Finset.sum_range_succ']
      have hz : ∀ x ∈ Finset.range m,
          ((PowerSeries.coeff (x + 1) (PowerSeries.C c) : ℤ) : R) * q ^ (x + 1) = 0 := by
        intro x _
        rw [PowerSeries.coeff_C, if_neg (Nat.succ_ne_zero x)]
        simp
      rw [Finset.sum_eq_zero hz]
      simp
    rw [hp]

theorem evalAdic_one [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    evalAdic 1 q hq = 1 := by
  have h1 : (1 : PowerSeries ℤ) = PowerSeries.C 1 := by simp
  rw [h1, evalAdic_C]
  simp

/-- ★★★**評価は環準同型である** `ℤ⟦q⟧ →+* R`。 -/
noncomputable def evalAdicHom [IsAdicComplete I R] (q : R) (hq : q ∈ I) : PowerSeries ℤ →+* R where
  toFun := fun f => evalAdic f q hq
  map_one' := evalAdic_one q hq
  map_mul' := fun f g => evalAdic_mul f g q hq
  map_zero' := by
    have h0 : (0 : PowerSeries ℤ) = PowerSeries.C 0 := by simp
    have hc := evalAdic_C (I := I) 0 q hq
    rw [← h0] at hc
    simpa using hc
  map_add' := fun f g => evalAdic_add f g q hq

/-! ## ★★★★完備環上の Tate 曲線 -/

/-- ★★★★**完備環上の Tate 曲線** `E_q`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★(G6) の **Tate 一意化** `E(K) ≅ Kˣ/q^ℤ` はこの曲線の上で述べられる。 -/
noncomputable def tateCurveAt [IsAdicComplete I R] (q : R) (hq : q ∈ I) : WeierstrassCurve R :=
  tateCurve.map (evalAdicHom q hq)

/-! ## ★出典の紐付け(`.src`) -/

def tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(完備環上の Tate 曲線 E_q)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
