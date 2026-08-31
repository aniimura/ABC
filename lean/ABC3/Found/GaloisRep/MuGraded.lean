/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicFinsetSum
import ABC3.Meta.Claim

/-!
# μ-等級付き `I` 進級数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`tateModel_of_quot_mu`（葉 1）の**手順 3**では
`v = ∑_ζ (3X(ζ)² + a₄ − Y(ζ))` のように **`X(ζ)` の積**が現れる。

★`X(ζ^i)`・`Y(ζ^i)` はどちらも

    `∑_n (∑_{a<l} A n a · (ζ^i)^a)`   （`I` 進和）

の形に書ける（`ζ` の指数は `ζ^l = 1` により **`mod l` で揃う**）。
★★本ファイルはこの形を `muEval` として名前を付け、

* `μ_l∖{1}` 上の和が `ζ`-free になること（`sum_mu_muEval`）
* 積がふたたび同じ形になること（`muEval_mul`）

を与える。☆これで手順 3 が**有限個の係数計算**に落ちる。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**μ-等級付き `I` 進級数の値**——`z = ζ^i` を代入した値。 -/
noncomputable def muEval [IsAdicComplete I R] (l : ℕ) (A : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (z : R) : R :=
  adicSum (I := I) (fun n => ∑ a ∈ range l, A n a * z ^ a)
    (fun n => Submodule.sum_mem _ (fun a _ => Ideal.mul_mem_right _ _ (hA n a)))

/-- ★★★★★★★★★★★★★★★★★★★★
**μ-等級付き級数を `μ_l∖{1}` 上で足すと `ζ` が消える**。

    `∑_{i≠0} muEval l A (ζ^i) = adicSum (n ↦ ∑_{a<l} A n a · ([l ∣ a]·l − 1))`

★★`a` は `0 ≤ a < l` なので `l ∣ a ⇔ a = 0`——★右辺は
`l·A n 0 − ∑_{a<l} A n a` である。 -/
theorem sum_mu_muEval [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (A : ℕ → ℕ → R) (hA : ∀ n a, A n a ∈ I ^ n) :
    ∑ i ∈ (range l).erase 0, muEval l A hA (ζ ^ i)
      = adicSum (I := I)
          (fun n => ∑ a ∈ range l, A n a * ((if l ∣ a then (l : R) else 0) - 1))
          (fun n => Submodule.sum_mem _
            (fun a _ => Ideal.mul_mem_right _ _ (hA n a))) := by
  classical
  simp only [muEval]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => ∑ a ∈ range l, A n a * (ζ ^ i) ^ a)
      (fun i n => Submodule.sum_mem _
        (fun a _ => Ideal.mul_mem_right _ _ (hA n a)))]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [← Finset.mul_sum, ← sum_mu_pow_erase_zero hl hζ a]
  exact congrArg _ (Finset.sum_congr rfl (fun i _ => by rw [← pow_mul]))

/-- ★★★★★★`a < l` のとき `l ∣ a ⇔ a = 0` なので、和は簡単になる。 -/
theorem sum_mu_muEval' [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (A : ℕ → ℕ → R) (hA : ∀ n a, A n a ∈ I ^ n) :
    ∑ i ∈ (range l).erase 0, muEval l A hA (ζ ^ i)
      = adicSum (I := I)
          (fun n => (l : R) * A n 0 - ∑ a ∈ range l, A n a)
          (fun n => Submodule.sub_mem _
            (Ideal.mul_mem_left _ _ (hA n 0))
            (Submodule.sum_mem _ (fun a _ => hA n a))) := by
  classical
  rw [sum_mu_muEval hl hζ A hA]
  refine adicSum_congr _ _ (fun n => ?_)
  have hsplit : ∑ a ∈ range l, A n a * ((if l ∣ a then (l : R) else 0) - 1)
      = ∑ a ∈ range l, (A n a * (if l ∣ a then (l : R) else 0) - A n a) := by
    exact Finset.sum_congr rfl (fun a _ => by ring)
  rw [hsplit, Finset.sum_sub_distrib]
  congr 1
  refine Finset.sum_eq_single 0 (fun a ha hne => ?_) (fun h => absurd (Finset.mem_range.2 hl.pos) h)
    |>.trans ?_
  · rw [if_neg (fun hdvd => hne (Nat.eq_zero_of_dvd_of_lt hdvd (Finset.mem_range.1 ha))), mul_zero]
  · rw [if_pos (dvd_zero l), mul_comm]

/-! ## ★★★★★★★★★★★★★★★★積 -/

/-- ★`z^l = 1` なら指数は `mod l` で済む。 -/
theorem pow_mod_eq {z : R} {l : ℕ} (hl : 0 < l) (hz : z ^ l = 1) (m : ℕ) :
    z ^ m = z ^ (m % l) := by
  conv_lhs => rw [← Nat.div_add_mod m l]
  rw [pow_add, pow_mul, hz, one_pow, one_mul]

/-- ★★★★★★**畳み込み**——`q`-次数は Cauchy 積、`ζ`-次数は `mod l` で足す。 -/
noncomputable def muConv (l : ℕ) (A B : ℕ → ℕ → R) : ℕ → ℕ → R :=
  fun n c => ∑ k ∈ range (n + 1), ∑ a ∈ range l, ∑ b ∈ range l,
    if (a + b) % l = c then A k a * B (n - k) b else 0

theorem muConv_mem {l : ℕ} {A B : ℕ → ℕ → R} (hA : ∀ n a, A n a ∈ I ^ n)
    (hB : ∀ n a, B n a ∈ I ^ n) (n c : ℕ) : muConv l A B n c ∈ I ^ n := by
  refine Submodule.sum_mem _ (fun k hk => Submodule.sum_mem _ (fun a _ =>
    Submodule.sum_mem _ (fun b _ => ?_)))
  have hle : k ≤ n := Nat.lt_succ_iff.1 (Finset.mem_range.1 hk)
  by_cases h : (a + b) % l = c
  · rw [if_pos h]
    have hmem : A k a * B (n - k) b ∈ I ^ k * I ^ (n - k) :=
      Ideal.mul_mem_mul (hA k a) (hB (n - k) b)
    have he : k + (n - k) = n := by omega
    rwa [← pow_add, he] at hmem
  · rw [if_neg h]
    exact Submodule.zero_mem _

/-- ★★★★★★★★★★★★★★★★★★★★★★
**μ-等級付き級数の積はふたたび μ-等級付き**。

★★これで `X(ζ)²`・`X(ζ)Y(ζ)` などがすべて同じ枠に乗る。 -/
theorem muEval_mul [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) (A B : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (hB : ∀ n a, B n a ∈ I ^ n)
    {z : R} (hz : z ^ l = 1) :
    muEval l A hA z * muEval l B hB z
      = muEval l (muConv l A B) (muConv_mem hA hB) z := by
  classical
  simp only [muEval]
  rw [adicSum_mul]
  refine adicSum_congr _ _ (fun n => ?_)
  have hinner : ∀ k : ℕ,
      (∑ a ∈ range l, A k a * z ^ a) * (∑ b ∈ range l, B (n - k) b * z ^ b)
        = ∑ c ∈ range l, (∑ a ∈ range l, ∑ b ∈ range l,
            (if (a + b) % l = c then A k a * B (n - k) b else 0)) * z ^ c := by
    intro k
    calc (∑ a ∈ range l, A k a * z ^ a) * (∑ b ∈ range l, B (n - k) b * z ^ b)
        = ∑ a ∈ range l, ∑ b ∈ range l, (A k a * B (n - k) b) * z ^ ((a + b) % l) := by
          rw [Finset.sum_mul_sum]
          refine Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun b _ => ?_))
          rw [← pow_mod_eq hl hz (a + b), pow_add]
          ring
      _ = ∑ a ∈ range l, ∑ b ∈ range l, ∑ c ∈ range l,
            (if (a + b) % l = c then (A k a * B (n - k) b) * z ^ c else 0) := by
          refine Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun b _ => ?_))
          rw [Finset.sum_eq_single ((a + b) % l)]
          · rw [if_pos rfl]
          · intro c _ hne
            rw [if_neg (Ne.symm hne)]
          · intro h
            exact absurd (Finset.mem_range.2 (Nat.mod_lt _ hl)) h
      _ = ∑ c ∈ range l, ∑ a ∈ range l, ∑ b ∈ range l,
            (if (a + b) % l = c then (A k a * B (n - k) b) * z ^ c else 0) := by
          have hswap : ∀ a : ℕ,
              (∑ b ∈ range l, ∑ c ∈ range l,
                (if (a + b) % l = c then (A k a * B (n - k) b) * z ^ c else 0))
              = ∑ c ∈ range l, ∑ b ∈ range l,
                (if (a + b) % l = c then (A k a * B (n - k) b) * z ^ c else 0) :=
            fun a => Finset.sum_comm
          rw [Finset.sum_congr rfl (fun a _ => hswap a), Finset.sum_comm]
      _ = ∑ c ∈ range l, (∑ a ∈ range l, ∑ b ∈ range l,
            (if (a + b) % l = c then A k a * B (n - k) b else 0)) * z ^ c := by
          refine Finset.sum_congr rfl (fun c _ => ?_)
          rw [Finset.sum_mul]
          refine Finset.sum_congr rfl (fun a _ => ?_)
          rw [Finset.sum_mul]
          refine Finset.sum_congr rfl (fun b _ => ?_)
          by_cases h : (a + b) % l = c <;> simp [h]
  rw [Finset.sum_congr rfl (fun k _ => hinner k), Finset.sum_comm]
  refine Finset.sum_congr rfl (fun c _ => ?_)
  rw [muConv, ← Finset.sum_mul]

/-! ## ★★★★★★★★和・定数倍・`ζ`-free な項 -/

theorem muEval_add [IsAdicComplete I R] {l : ℕ} (A B : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (hB : ∀ n a, B n a ∈ I ^ n) (z : R) :
    muEval l A hA z + muEval l B hB z
      = muEval l (fun n a => A n a + B n a)
          (fun n a => Submodule.add_mem _ (hA n a) (hB n a)) z := by
  classical
  simp only [muEval]
  rw [← adicSum_add]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl (fun a _ => by ring)

theorem muEval_smul [IsAdicComplete I R] {l : ℕ} (c : R) (A : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (z : R) :
    c * muEval l A hA z
      = muEval l (fun n a => c * A n a) (fun n a => Ideal.mul_mem_left _ _ (hA n a)) z := by
  classical
  simp only [muEval]
  rw [← adicSum_smul c (fun n => ∑ a ∈ range l, A n a * z ^ a)
    (fun n => Submodule.sum_mem _ (fun a _ => Ideal.mul_mem_right _ _ (hA n a)))]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl (fun a _ => by ring)

/-- ★★★★★★★★**`ζ`-free な `I` 進和も μ-等級付きと見られる**。 -/
theorem adicSum_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) (f : ℕ → R)
    (hf : ∀ n, f n ∈ I ^ n) (z : R) :
    adicSum f hf
      = muEval l (fun n a => if a = 0 then f n else 0)
          (fun n a => by
            by_cases h : a = 0
            · simpa [h] using hf n
            · simpa [h] using Submodule.zero_mem (I ^ n)) z := by
  classical
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.sum_eq_single 0]
  · simp
  · intro a _ hne
    simp [hne]
  · intro h
    exact absurd (Finset.mem_range.2 hl) h

/-! ## ★★★★★★★★★★★★★★★★★★Tate の尾を μ-等級付きに直す -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★
**`tateXtail(z, q)` は μ-等級付き級数である**（`z^l = 1` のとき）。

★係数は `A n a = q^n · ∑_{d ∣ n, d ≡ a (mod l)} d`。 -/
theorem tateXtail_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) :
    tateXtail z q hq
      = muEval l (fun n a => q ^ n * ∑ d ∈ n.divisors.filter (fun d => d % l = a), (d : R))
          (fun n a => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) z := by
  classical
  rw [tateXtail_eq_divisorSum z q hq]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  have hfib : ∑ a ∈ range l, ∑ d ∈ n.divisors.filter (fun d => d % l = a),
        (d : R) * z ^ d
      = ∑ d ∈ n.divisors, (d : R) * z ^ d :=
    Finset.sum_fiberwise_of_maps_to (fun d _ => Finset.mem_range.2 (Nat.mod_lt _ hl)) _
  rw [← hfib, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [mul_assoc, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl (fun d hd => ?_)
  have hda : d % l = a := (Finset.mem_filter.1 hd).2
  rw [pow_mod_eq hl hz d, hda]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★**`Y` 側の尾も同じ**。 -/
theorem tateYtail_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) :
    tateYtail z q hq
      = muEval l (fun n a => q ^ n * ∑ d ∈ n.divisors.filter (fun d => d % l = a),
            ((d.choose 2 : ℕ) : R))
          (fun n a => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) z := by
  classical
  rw [tateYtail_eq_divisorSum z q hq]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  have hfib : ∑ a ∈ range l, ∑ d ∈ n.divisors.filter (fun d => d % l = a),
        ((d.choose 2 : ℕ) : R) * z ^ d
      = ∑ d ∈ n.divisors, ((d.choose 2 : ℕ) : R) * z ^ d :=
    Finset.sum_fiberwise_of_maps_to (fun d _ => Finset.mem_range.2 (Nat.mod_lt _ hl)) _
  rw [← hfib, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [mul_assoc, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl (fun d hd => ?_)
  have hda : d % l = a := (Finset.mem_filter.1 hd).2
  rw [pow_mod_eq hl hz d, hda]

open Finset in
/-- ★★★★★★★★★★★★★★★★**`tateXterm(t)`（`t ∈ I`）も μ-等級付き**。

`t = q·z^{m}` の形のときに使う。 -/
theorem tateXterm_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) (m : ℕ) :
    tateXterm (q * z ^ m)
      = muEval l (fun n a => if (n * m) % l = a then (n : R) * q ^ n else 0)
          (fun n a => by
            by_cases h : (n * m) % l = a
            · simpa [h] using Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)
            · simpa [h] using Submodule.zero_mem (I ^ n)) z := by
  classical
  rw [tateXterm_eq_adicSum (Ideal.mul_mem_right _ _ hq)]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.sum_eq_single ((n * m) % l)]
  · rw [if_pos rfl, mul_pow, ← pow_mul, ← pow_mod_eq hl hz (n * m)]
    ring
  · intro a _ hne
    rw [if_neg (fun h => hne h.symm), zero_mul]
  · intro h
    exact absurd (Finset.mem_range.2 (Nat.mod_lt _ hl)) h

/-! ## ★★★★★★★★★★★★定数（`q`-次数 0）の項 -/

/-- ★★★★★★**0 次に集中した `I` 進和はその項そのもの**。 -/
theorem adicSum_concentrated [IsAdicComplete I R] (c : R) :
    adicSum (I := I) (fun n => if n = 0 then c else 0)
        (fun n => by
          by_cases h : n = 0
          · simpa [h] using Submodule.mem_top (x := c)
          · simpa [h] using Submodule.zero_mem (I ^ n)) = c := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  cases n with
  | zero =>
      simp only [partialSum, Finset.range_zero, Finset.sum_empty]
      rw [SModEq.sub_mem]
      simpa using Submodule.mem_top (x := (0 : R) - c)
  | succ m =>
      have hps : partialSum (fun n => if n = 0 then c else 0) (m + 1) = c := by
        simp only [partialSum]
        rw [Finset.sum_eq_single 0]
        · simp
        · intro b _ hb
          simp [hb]
        · intro h
          exact absurd (Finset.mem_range.2 (Nat.succ_pos m)) h
      rw [hps]

/-- ★★★★★★★★★★★★★★**`ζ` の多項式は μ-等級付き（`q`-次数 0）**。 -/
theorem muEval_const [IsAdicComplete I R] {l : ℕ} (c : ℕ → R) (z : R) :
    ∑ a ∈ range l, c a * z ^ a
      = muEval (I := I) l (fun n a => if n = 0 then c a else 0)
          (fun n a => by
            by_cases h : n = 0
            · simpa [h] using Submodule.mem_top (x := c a)
            · simpa [h] using Submodule.zero_mem (I ^ n)) z := by
  classical
  simp only [muEval]
  rw [← adicSum_concentrated (I := I) (∑ a ∈ range l, c a * z ^ a)]
  refine adicSum_congr _ _ (fun n => ?_)
  by_cases h : n = 0
  · simp [h]
  · simp [h]

/-! ## ★★★★★★★★★★★★★★★★`1/(1−ζ)` の環版と `tateXterm(ζ)` -/

/-- ★★★★★★★★★★★★**`Ring.inverse (1−η) = −(1/l)·∑_{k<l} kη^k`**（環版）。 -/
theorem ringInverse_one_sub_eq [IsDomain R] {l : ℕ} (hlu : IsUnit ((l : R))) {η : R}
    (hu : IsUnit (1 - η)) (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    Ring.inverse (1 - η) = -(Ring.inverse ((l : R))) * ∑ k ∈ range l, (k : R) * η ^ k := by
  have hcore : (1 - η) * ∑ k ∈ range l, (k : R) * η ^ k = -(l : R) :=
    one_sub_mul_sum_nsmul (R := R) hpow hsum
  have hmul : (1 - η) * (-(Ring.inverse ((l : R))) * ∑ k ∈ range l, (k : R) * η ^ k) = 1 := by
    rw [show (1 - η) * (-(Ring.inverse ((l : R))) * ∑ k ∈ range l, (k : R) * η ^ k)
        = -(Ring.inverse ((l : R))) * ((1 - η) * ∑ k ∈ range l, (k : R) * η ^ k) from by ring,
      hcore, neg_mul_neg, Ring.inverse_mul_cancel _ hlu]
  obtain ⟨u, hu'⟩ := hu
  rw [← hu', Ring.inverse_unit]
  refine Units.inv_eq_of_mul_eq_one_right ?_
  rw [hu']
  exact hmul

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★
**`tateXterm(ζ) = ζ/(1−ζ)²` は `ζ` の多項式**（`l` が単元のとき）。

★係数は `(1/l)²·∑_{k,m<l, k+m+1 ≡ a (mod l)} k·m`。
★★これで**定数項も μ-等級付きの枠に乗る**。 -/
theorem tateXterm_zeta_eq_poly [IsDomain R] {l : ℕ} (hl : 0 < l) (hlu : IsUnit ((l : R)))
    {z : R} (hu : IsUnit (1 - z)) (hpow : z ^ l = 1) (hsum : ∑ k ∈ range l, z ^ k = 0) :
    tateXterm z
      = ∑ a ∈ range l, ((Ring.inverse ((l : R))) ^ 2 *
          ∑ k ∈ range l, ∑ m ∈ range l,
            (if (k + m + 1) % l = a then (k : R) * (m : R) else 0)) * z ^ a := by
  classical
  have hcore : ∑ a ∈ range l, (∑ k ∈ range l, ∑ m ∈ range l,
          (if (k + m + 1) % l = a then (k : R) * (m : R) else 0)) * z ^ a
      = ∑ k ∈ range l, ∑ m ∈ range l, ((k : R) * (m : R)) * z ^ ((k + m + 1) % l) := by
    have h1 : ∀ a ∈ range l, (∑ k ∈ range l, ∑ m ∈ range l,
            (if (k + m + 1) % l = a then (k : R) * (m : R) else 0)) * z ^ a
        = ∑ k ∈ range l, ∑ m ∈ range l,
            (if (k + m + 1) % l = a then ((k : R) * (m : R)) * z ^ a else 0) := by
      intro a _
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl (fun k _ => ?_)
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl (fun m _ => ?_)
      by_cases h : (k + m + 1) % l = a <;> simp [h]
    rw [Finset.sum_congr rfl h1, Finset.sum_comm]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [Finset.sum_eq_single ((k + m + 1) % l)]
    · rw [if_pos rfl]
    · intro c _ hne
      rw [if_neg (Ne.symm hne)]
    · intro h
      exact absurd (Finset.mem_range.2 (Nat.mod_lt _ hl)) h
  rw [tateXterm, ringInverse_one_sub_eq hlu hu hpow hsum]
  have hLHS : z * (-(Ring.inverse ((l : R))) * ∑ k ∈ range l, (k : R) * z ^ k) ^ 2
      = (Ring.inverse ((l : R))) ^ 2 *
          ∑ k ∈ range l, ∑ m ∈ range l, ((k : R) * (m : R)) * z ^ ((k + m + 1) % l) := by
    have h2 : z * (-(Ring.inverse ((l : R))) * ∑ k ∈ range l, (k : R) * z ^ k) ^ 2
        = (Ring.inverse ((l : R))) ^ 2 *
          (((∑ k ∈ range l, (k : R) * z ^ k) * ∑ m ∈ range l, (m : R) * z ^ m) * z) := by
      ring
    rw [h2, Finset.sum_mul_sum, Finset.sum_mul]
    congr 1
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [← pow_mod_eq hl hpow (k + m + 1)]
    rw [pow_add, pow_add, pow_one]
    ring
  rw [hLHS, ← hcore, Finset.mul_sum]
  exact Finset.sum_congr rfl (fun a _ => by ring)

/-! ## ★★★★★★★★★★★★`w = q·z^m` 側の尾 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★**`tateXtail(q·z^m, q)` も μ-等級付き**。 -/
theorem tateXtail_qz_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) (m : ℕ) :
    tateXtail (q * z ^ m) q hq
      = muEval l (fun n a => q ^ n * ∑ d ∈ n.divisors.filter (fun d => (m * d) % l = a),
            (d : R) * q ^ d)
          (fun n a => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) z := by
  classical
  rw [tateXtail_eq_divisorSum (q * z ^ m) q hq]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  have hfib : ∑ a ∈ range l, ∑ d ∈ n.divisors.filter (fun d => (m * d) % l = a),
        (d : R) * (q * z ^ m) ^ d
      = ∑ d ∈ n.divisors, (d : R) * (q * z ^ m) ^ d :=
    Finset.sum_fiberwise_of_maps_to (fun d _ => Finset.mem_range.2 (Nat.mod_lt _ hl)) _
  rw [← hfib, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [mul_assoc, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl (fun d hd => ?_)
  have hda : (m * d) % l = a := (Finset.mem_filter.1 hd).2
  rw [mul_pow, ← pow_mul, pow_mod_eq hl hz (m * d), hda]
  ring

open Finset in
/-- ★★★★★★★★★★★★★★★★**`tateYtail(q·z^m, q)` も μ-等級付き**。 -/
theorem tateYtail_qz_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) (m : ℕ) :
    tateYtail (q * z ^ m) q hq
      = muEval l (fun n a => q ^ n * ∑ d ∈ n.divisors.filter (fun d => (m * d) % l = a),
            ((d.choose 2 : ℕ) : R) * q ^ d)
          (fun n a => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) z := by
  classical
  rw [tateYtail_eq_divisorSum (q * z ^ m) q hq]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  have hfib : ∑ a ∈ range l, ∑ d ∈ n.divisors.filter (fun d => (m * d) % l = a),
        ((d.choose 2 : ℕ) : R) * (q * z ^ m) ^ d
      = ∑ d ∈ n.divisors, ((d.choose 2 : ℕ) : R) * (q * z ^ m) ^ d :=
    Finset.sum_fiberwise_of_maps_to (fun d _ => Finset.mem_range.2 (Nat.mod_lt _ hl)) _
  rw [← hfib, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [mul_assoc, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl (fun d hd => ?_)
  have hda : (m * d) % l = a := (Finset.mem_filter.1 hd).2
  rw [mul_pow, ← pow_mul, pow_mod_eq hl hz (m * d), hda]
  ring

open Finset in
/-- ★★★★★★★★★★★★★★★★**`tateYterm(q·z^m)` も μ-等級付き**。 -/
theorem tateYterm_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l)
    {z : R} (hz : z ^ l = 1) (q : R) (hq : q ∈ I) (m : ℕ) :
    tateYterm (q * z ^ m)
      = muEval l (fun n a => if (n * m) % l = a then ((n.choose 2 : ℕ) : R) * q ^ n else 0)
          (fun n a => by
            by_cases h : (n * m) % l = a
            · simpa [h] using Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)
            · simpa [h] using Submodule.zero_mem (I ^ n)) z := by
  classical
  rw [tateYterm_eq_adicSum (Ideal.mul_mem_right _ _ hq)]
  simp only [muEval]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.sum_eq_single ((n * m) % l)]
  · rw [if_pos rfl, mul_pow, ← pow_mul, ← pow_mod_eq hl hz (n * m)]
    ring
  · intro a _ hne
    rw [if_neg (fun h => hne h.symm), zero_mul]
  · intro h
    exact absurd (Finset.mem_range.2 (Nat.mod_lt _ hl)) h

/-! ## ★★★★★★★★係数の合同と差 -/

theorem muEval_congr [IsAdicComplete I R] {l : ℕ} (A B : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (hB : ∀ n a, B n a ∈ I ^ n)
    (h : ∀ n a, A n a = B n a) (z : R) :
    muEval l A hA z = muEval l B hB z := by
  simp only [muEval]
  exact adicSum_congr _ _ (fun n => Finset.sum_congr rfl (fun a _ => by rw [h n a]))

theorem muEval_sub [IsAdicComplete I R] {l : ℕ} (A B : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (hB : ∀ n a, B n a ∈ I ^ n) (z : R) :
    muEval l A hA z - muEval l B hB z
      = muEval l (fun n a => A n a - B n a)
          (fun n a => Submodule.sub_mem _ (hA n a) (hB n a)) z := by
  classical
  simp only [muEval]
  rw [← adicSum_sub]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl (fun a _ => by ring)

/-- ★★★★★★★★★★★★**`s₁(q)` は `ζ`-free な μ-等級付き**。 -/
theorem sigmaOne_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) (q : R) (hq : q ∈ I) :
    evalAdic (sigmaSeries 1) q hq
      = muEval (I := I) l
          (fun n a => if a = 0 then q ^ n * ∑ d ∈ n.divisors, (d : R) else 0)
          (fun n a => by
            by_cases h : a = 0
            · simpa [h] using Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
            · simpa [h] using Submodule.zero_mem (I ^ n)) q := by
  classical
  rw [← tateXtail_one q hq, tateXtail_eq_divisorSum 1 q hq]
  rw [adicSum_eq_muEval (l := l) hl (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) * (1 : R) ^ d)
    (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) q]
  refine muEval_congr _ _ _ _ (fun n a => ?_) q
  by_cases h : a = 0
  · simp [h]
  · simp [h]

/-! ## ★出典の紐付け(`.src`) -/

def pow_mod_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(z^l = 1 なら指数は mod l で済む。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muConv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の畳み込み)",
    sectionId := "genell-lemma-3-2" }

def muConv_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(畳み込みの係数も I^n に入る。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の積はふたたび μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の定数倍。★無条件)",
    sectionId := "genell-lemma-3-2" }

def adicSum_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ-free な I 進和も μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateXtail_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateXtail は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateYtail_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateYtail は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateXterm_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateXterm(q·z^m) は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def adicSum_concentrated.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(0 次に集中した I 進和はその項。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_const.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ の多項式は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def ringInverse_one_sub_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Ring.inverse(1−η) を η の多項式で書く。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateXterm_zeta_eq_poly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateXterm(ζ) は ζ の多項式。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateXtail_qz_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateXtail(q·z^m,q) は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateYtail_qz_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateYtail(q·z^m,q) は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateYterm_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(tateYterm(q·z^m) は μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(係数が等しければ値も等しい。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の差。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sigmaOne_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(s₁(q) は ζ-free な μ-等級付き。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き I 進級数の値)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の μ_l 和は ζ-free。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_muEval'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l 和は l·A n 0 − ∑_a A n a。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
