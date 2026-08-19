import ABC3.Found.GaloisRep.TateJ

/-!
# Galois (G6) 第 100 ブロック —— **★★★★★Tate 母数の存在(j 反転)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★`1/j` が反転できた

第 99 ブロックで `Δ/c₄³ = X·(単元)` を得た。★本ブロックはそれを**反転**する:

    任意の `t ∈ I` に対し、`q ∈ I` が在って `(X·u)(q) = t`

★★機構は**不動点反復**である。`u` の逆 `v`(冪級数として)を取ると
`q = t·v(q)` の不動点を探せばよく、

    q₀ = 0,   q_{n+1} = t · v(q_n)

★★★連続する 2 項の差は `I^{n+1}` に入る(`t ∈ I` と評価の連続性から)ので
`IsPrecomplete` で極限が取れ、`IsHausdorff` で不動点方程式が成り立つ。

★★★★**位相を使っていない**——`IsAdicComplete` だけで完結する。

## ★★これで Tate 母数が作れる

`j` が与えられたとき、`t := 1/j` に対する `q` が **Tate 母数**である。
★(G6) に残るのは **一意化定理** `E_q(K) ≅ Kˣ/q^ℤ` だけになる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `partialEval_sub_mem` / `evalAdic_sub_mem` | ★★★評価の連続性 |
| `tateIter` / `tateIter_succ_sub` / `tateIter_cauchy` | ★★不動点反復と Cauchy 性 |
| `exists_fixedPoint` | ★★★★**不動点の存在** |
| `exists_evalAdic_eq` | ★★★★★**`X·u` は `I` の上で全射** |
| `exists_tateParam` | ★★★★★**Tate 母数の存在** |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★評価の連続性 -/

/-- ★★部分和は Lipschitz——`a ≡ b (mod I^k)` なら値も `I^k` で近い。 -/
theorem partialEval_sub_mem (f : PowerSeries ℤ) {a b : R} (k : ℕ) (hab : a - b ∈ I ^ k) (N : ℕ) :
    partialEval f a N - partialEval f b N ∈ I ^ k := by
  rw [partialEval, partialEval, ← Finset.sum_sub_distrib]
  refine Submodule.sum_mem _ (fun n _ => ?_)
  have hdiff : ((PowerSeries.coeff n f : ℤ) : R) * a ^ n - ((PowerSeries.coeff n f : ℤ) : R) * b ^ n
      = ((PowerSeries.coeff n f : ℤ) : R) * (a ^ n - b ^ n) := by ring
  rw [hdiff]
  exact Ideal.mul_mem_left _ _ (Ideal.mem_of_dvd _ (sub_dvd_pow_sub_pow a b n) hab)

/-- ★★★**評価は連続である**。 -/
theorem evalAdic_sub_mem [IsAdicComplete I R] (f : PowerSeries ℤ) {a b : R} (ha : a ∈ I)
    (hb : b ∈ I) (k : ℕ) (hab : a - b ∈ I ^ k) :
    evalAdic f a ha - evalAdic f b hb ∈ I ^ k := by
  have hsa := evalAdic_spec f a ha k
  have hsb := evalAdic_spec f b hb k
  rw [SModEq.sub_mem] at hsa hsb
  have hsa' : partialEval f a k - evalAdic f a ha ∈ I ^ k := by simpa using hsa
  have hsb' : partialEval f b k - evalAdic f b hb ∈ I ^ k := by simpa using hsb
  have hpart := partialEval_sub_mem (I := I) f k hab k
  have hexp : evalAdic f a ha - evalAdic f b hb
      = (partialEval f a k - partialEval f b k)
        - (partialEval f a k - evalAdic f a ha) + (partialEval f b k - evalAdic f b hb) := by
    ring
  rw [hexp]
  exact Submodule.add_mem _ (Submodule.sub_mem _ hpart hsa') hsb'

/-! ## ★★不動点反復 -/

/-- ★不動点反復 `x ↦ t·v(x)`。 -/
noncomputable def tateIter [IsAdicComplete I R] (v : PowerSeries ℤ) {t : R} (ht : t ∈ I) :
    {x : R // x ∈ I} → {x : R // x ∈ I} :=
  fun x => ⟨t * evalAdic v x.1 x.2, Ideal.mul_mem_right _ _ ht⟩

/-- ★★連続する 2 項の差は `I^{n+1}` に入る。 -/
theorem tateIter_succ_sub [IsAdicComplete I R] (v : PowerSeries ℤ) {t : R} (ht : t ∈ I)
    (x0 : {x : R // x ∈ I}) : ∀ n : ℕ,
      (((tateIter v ht)^[n + 1] x0).1 - ((tateIter v ht)^[n] x0).1) ∈ I ^ (n + 1) := by
  intro n
  induction n with
  | zero =>
    have h1 : I ^ (0 + 1) = I := by simp
    rw [h1]
    exact Submodule.sub_mem _ ((tateIter v ht)^[0 + 1] x0).2 ((tateIter v ht)^[0] x0).2
  | succ m ih =>
    have hstep : (((tateIter v ht)^[m + 2] x0).1 - ((tateIter v ht)^[m + 1] x0).1)
        = t * (evalAdic v ((tateIter v ht)^[m + 1] x0).1 ((tateIter v ht)^[m + 1] x0).2
              - evalAdic v ((tateIter v ht)^[m] x0).1 ((tateIter v ht)^[m] x0).2) := by
      rw [Function.iterate_succ_apply' (f := tateIter v ht) (n := m + 1),
        Function.iterate_succ_apply' (f := tateIter v ht) (n := m)]
      simp only [tateIter]
      ring
    rw [hstep]
    have hdiff := evalAdic_sub_mem (I := I) v
      ((tateIter v ht)^[m + 1] x0).2 ((tateIter v ht)^[m] x0).2 (m + 1) ih
    have hmul : t * (evalAdic v ((tateIter v ht)^[m + 1] x0).1 ((tateIter v ht)^[m + 1] x0).2
        - evalAdic v ((tateIter v ht)^[m] x0).1 ((tateIter v ht)^[m] x0).2) ∈ I * I ^ (m + 1) :=
      Ideal.mul_mem_mul ht hdiff
    simpa [pow_succ'] using hmul

/-- ★★反復列は `I` 進 Cauchy である。 -/
theorem tateIter_cauchy [IsAdicComplete I R] (v : PowerSeries ℤ) {t : R} (ht : t ∈ I)
    (x0 : {x : R // x ∈ I}) {m n : ℕ} (hmn : m ≤ n) :
    ((tateIter v ht)^[m] x0).1 ≡ ((tateIter v ht)^[n] x0).1
      [SMOD (I ^ m • ⊤ : Submodule R R)] := by
  induction n, hmn using Nat.le_induction with
  | base => rfl
  | succ n hmn ih =>
    refine ih.trans ?_
    rw [SModEq.sub_mem]
    have h := tateIter_succ_sub v ht x0 n
    have hle : I ^ (n + 1) ≤ I ^ m := Ideal.pow_le_pow_right (by omega)
    simpa using neg_mem (hle h)

/-- ★★★★**不動点の存在**——`q = t·v(q)` を満たす `q ∈ I` がある。 -/
theorem exists_fixedPoint [IsAdicComplete I R] (v : PowerSeries ℤ) {t : R} (ht : t ∈ I) :
    ∃ (q : R) (hq : q ∈ I), q = t * evalAdic v q hq := by
  set x0 : {x : R // x ∈ I} := ⟨0, I.zero_mem⟩ with hx0
  set f : ℕ → R := fun n => ((tateIter v ht)^[n] x0).1 with hf
  obtain ⟨L, hL⟩ := IsPrecomplete.prec' f (fun {m n} hmn => tateIter_cauchy v ht x0 hmn)
  have hLI : L ∈ I := by
    have h1 := hL 1
    rw [SModEq.sub_mem] at h1
    have h2 : f 1 - L ∈ I := by simpa using h1
    have h3 : f 1 ∈ I := ((tateIter v ht)^[1] x0).2
    have h4 : L = f 1 - (f 1 - L) := by ring
    rw [h4]
    exact Submodule.sub_mem _ h3 h2
  refine ⟨L, hLI, ?_⟩
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => ?_)
  rw [SModEq.sub_mem]
  have hn1 := hL (n + 1)
  rw [SModEq.sub_mem] at hn1
  have hfn1 : f (n + 1) - L ∈ I ^ (n + 1) := by simpa using hn1
  have hfnL : f n - L ∈ I ^ n := by
    have hx := hL n
    rw [SModEq.sub_mem] at hx
    simpa using hx
  have hstep : ((tateIter v ht)^[n + 1] x0).1
      = t * evalAdic v ((tateIter v ht)^[n] x0).1 ((tateIter v ht)^[n] x0).2 := by
    rw [Function.iterate_succ_apply' (f := tateIter v ht) (n := n)]
    rfl
  have hvdiff : evalAdic v ((tateIter v ht)^[n] x0).1 ((tateIter v ht)^[n] x0).2
      - evalAdic v L hLI ∈ I ^ n :=
    evalAdic_sub_mem (I := I) v ((tateIter v ht)^[n] x0).2 hLI n hfnL
  have hmul : t * (evalAdic v ((tateIter v ht)^[n] x0).1 ((tateIter v ht)^[n] x0).2
      - evalAdic v L hLI) ∈ I ^ (n + 1) := by
    have hm := Ideal.mul_mem_mul ht hvdiff
    simpa [pow_succ'] using hm
  have hexp : L - t * evalAdic v L hLI
      = -(f (n + 1) - L) + t * (evalAdic v ((tateIter v ht)^[n] x0).1
          ((tateIter v ht)^[n] x0).2 - evalAdic v L hLI) := by
    show L - t * evalAdic v L hLI
      = -(((tateIter v ht)^[n + 1] x0).1 - L) + _
    rw [hstep]
    ring
  have hfinal : L - t * evalAdic v L hLI ∈ I ^ (n + 1) := by
    rw [hexp]
    exact Submodule.add_mem _ (neg_mem hfn1) hmul
  have hle : I ^ (n + 1) ≤ I ^ n := Ideal.pow_le_pow_right (Nat.le_succ n)
  simpa using hle hfinal

/-! ## ★★★★★反転 -/

/-- ★★★★★**`X·u` は `I` の上で全射である**。 -/
theorem exists_evalAdic_eq [IsAdicComplete I R] (u v : PowerSeries ℤ) (huv : v * u = 1)
    {t : R} (ht : t ∈ I) :
    ∃ (q : R) (hq : q ∈ I), evalAdic (PowerSeries.X * u) q hq = t := by
  obtain ⟨q, hq, hfix⟩ := exists_fixedPoint v ht
  refine ⟨q, hq, ?_⟩
  rw [evalAdic_mul, evalAdic_X]
  calc q * evalAdic u q hq = (t * evalAdic v q hq) * evalAdic u q hq := by rw [← hfix]
    _ = t * (evalAdic v q hq * evalAdic u q hq) := by ring
    _ = t * evalAdic (v * u) q hq := by rw [evalAdic_mul]
    _ = t * evalAdic 1 q hq := by rw [huv]
    _ = t := by rw [evalAdic_one]; ring

/-- ★★★★★**Tate 母数の存在**——`1/j` の値を指定すると `q` が定まる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが古典的な「`j(q) = 1/q + 744 + ⋯` の反転」である。 -/
theorem exists_tateParam [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    ∃ (q : R) (hq : q ∈ I),
      evalAdic tateCurve.Δ q hq = t * evalAdic (tateCurve.c₄ ^ 3) q hq := by
  obtain ⟨w, hw, hΔ⟩ := tateInvJ_eq_X_mul_unit
  obtain ⟨wu, hwu⟩ := hw
  have hinv : (↑wu⁻¹ : PowerSeries ℤ) * w = 1 := by
    rw [← hwu]
    simp
  obtain ⟨q, hq, hqe⟩ := exists_evalAdic_eq w (↑wu⁻¹) hinv ht
  refine ⟨q, hq, ?_⟩
  rw [hΔ, evalAdic_mul, hqe]

/-! ## ★出典の紐付け(`.src`) -/

def exists_tateParam.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数の存在——j の反転)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
