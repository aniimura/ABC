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
