import ABC3.Found.GaloisRep.TateTailBound

/-!
# Galois (G6) 第 228 ブロック —— **★★★★★ℂ の上の二重和**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★残る穴——`s₁` の q 展開を ℂ の上で

第 227 で尾の評価が取れた。残るのは `X` の定義に入っている `−2 s₁(q)` の突き合わせである:

    ∑_{m≥1} f(qᵐ) = ∑_{N≥1} σ₁(N) qᴺ

★形式側では第 111 ブロック(`tateXtail_one`)で取れているが、**ℂ の上では別に要る**。
`evalAdic` の議論は `I` 進完備環でしか使えないからである。

## ★★★★二重和に開く

`‖t‖ < 1` なら `f(t) = ∑_{d≥0}(d+1)t^{d+1}`(幾何級数の微分)。よって

    ∑_{m≥1} f(qᵐ) = ∑_{(m,d) ∈ ℕ×ℕ} (d+1)·q^{(m+1)(d+1)}

★★**二重族の可和性**が要る。`(m+1)(d+1) ≥ m + (d+1)` なので

    ‖(d+1)q^{(m+1)(d+1)}‖ ≤ ‖q‖^m · (d+1)‖q‖^{d+1}

と**積の形**に分かれ、`Summable.mul_of_nonneg` で押さえられる(`summable_double_family`)。

## ★★次の段

あとは `N = (m+1)(d+1)` で括り直せば `∑_{d|N} d = σ₁(N)` が出る。
`Equiv.sigmaFiberEquiv` と `Summable.tsum_sigma`、そして
`Nat.divisorsAntidiagonal`(`(a,b) ↦ ab = N` の有限集合)が道具である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hasSum_nat_mul_geometric` | ★`f(t) = ∑_{d≥0} d·tᵈ`(ℂ の上) |
| `hasSum_nat_succ_mul_geometric` | ★`f(t) = ∑_{d≥0}(d+1)t^{d+1}` |
| `add_succ_le_mul_succ` | ★`a + (b+1) ≤ (a+1)(b+1)` |
| `summable_double_family` | ★★★★**二重族は可和** |
| `tsum_double_family_eq` | ★★★★★**二重和の内側は `f(q^{m+1})`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★幾何級数の微分 -/

/-- ★`f(t) = ∑_{d≥0} d·tᵈ`(ℂ の上)。 -/
theorem hasSum_nat_mul_geometric {t : ℂ} (ht : ‖t‖ < 1) :
    HasSum (fun d : ℕ => (d : ℂ) * t ^ d) (tateXterm t) := by
  rw [tateXterm_field]
  exact hasSum_coe_mul_geometric_of_norm_lt_one ht

/-- ★`f(t) = ∑_{d≥0} (d+1)·t^{d+1}`——`d = 0` の項は 0 なのでずらせる。 -/
theorem hasSum_nat_succ_mul_geometric {t : ℂ} (ht : ‖t‖ < 1) :
    HasSum (fun d : ℕ => ((d + 1 : ℕ) : ℂ) * t ^ (d + 1)) (tateXterm t) := by
  have h := (hasSum_nat_add_iff' (f := fun d : ℕ => (d : ℂ) * t ^ d) 1).2
    (hasSum_nat_mul_geometric ht)
  simpa using h

theorem summable_nat_succ_mul_geometric_real {r : ℝ} (hr0 : 0 ≤ r) (hr : r < 1) :
    Summable fun d : ℕ => ((d + 1 : ℕ) : ℝ) * r ^ (d + 1) := by
  have h := (summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) 1 (by rwa [Real.norm_eq_abs,
    abs_of_nonneg hr0]))
  have h2 := (summable_nat_add_iff 1).2 h
  refine h2.congr fun d => ?_
  push_cast
  ring

/-! ## ★★★★二重族の可和性 -/

theorem add_succ_le_mul_succ (a b : ℕ) : a + (b + 1) ≤ (a + 1) * (b + 1) := by
  calc a + (b + 1) = 0 + (a + b + 1) := by ring
    _ ≤ a * b + (a + b + 1) := Nat.add_le_add_right (Nat.zero_le _) _
    _ = (a + 1) * (b + 1) := by ring

set_option maxHeartbeats 1000000 in
/-- ★★★★**二重族 `(m,d) ↦ (d+1)q^{(m+1)(d+1)}` は可和**。

★`(m+1)(d+1) ≥ m + (d+1)` なので**積の形**に分かれ、`Summable.mul_of_nonneg` で押さえられる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem summable_double_family (q : ℂ) (hq : ‖q‖ < 1) :
    Summable fun p : ℕ × ℕ => ((p.2 + 1 : ℕ) : ℂ) * q ^ ((p.1 + 1) * (p.2 + 1)) := by
  have hr0 : (0 : ℝ) ≤ ‖q‖ := norm_nonneg q
  have hg : Summable fun p : ℕ × ℕ =>
      (‖q‖ ^ p.1) * (((p.2 + 1 : ℕ) : ℝ) * ‖q‖ ^ (p.2 + 1)) :=
    Summable.mul_of_nonneg (summable_geometric_of_lt_one hr0 hq)
      (summable_nat_succ_mul_geometric_real hr0 hq)
      (fun m => by positivity) (fun d => by positivity)
  refine Summable.of_norm_bounded hg fun p => ?_
  rw [norm_mul, norm_pow, Complex.norm_natCast]
  have hpow : ‖q‖ ^ ((p.1 + 1) * (p.2 + 1)) ≤ ‖q‖ ^ (p.1 + (p.2 + 1)) :=
    pow_le_pow_of_le_one hr0 hq.le (add_succ_le_mul_succ p.1 p.2)
  have hd : (0 : ℝ) ≤ ((p.2 + 1 : ℕ) : ℝ) := by positivity
  calc ((p.2 + 1 : ℕ) : ℝ) * ‖q‖ ^ ((p.1 + 1) * (p.2 + 1))
      ≤ ((p.2 + 1 : ℕ) : ℝ) * ‖q‖ ^ (p.1 + (p.2 + 1)) := mul_le_mul_of_nonneg_left hpow hd
    _ = ‖q‖ ^ p.1 * (((p.2 + 1 : ℕ) : ℝ) * ‖q‖ ^ (p.2 + 1)) := by
        rw [pow_add]
        ring

/-! ## ★★★★★二重和の内側 -/

/-- ★★★★★**二重和の内側は `f(q^{m+1})`**——`∑_{m≥1}f(qᵐ)` が二重和に開く。 -/
theorem tsum_double_family_eq (q : ℂ) (hq : ‖q‖ < 1) :
    (∑' p : ℕ × ℕ, ((p.2 + 1 : ℕ) : ℂ) * q ^ ((p.1 + 1) * (p.2 + 1)))
      = ∑' m : ℕ, tateXterm (q ^ (m + 1)) := by
  rw [(summable_double_family q hq).tsum_prod]
  refine tsum_congr fun m => ?_
  have hqm : ‖q ^ (m + 1)‖ < 1 := norm_pow_succ_lt_one hq m
  have h := (hasSum_nat_succ_mul_geometric hqm).tsum_eq
  rw [← h]
  refine tsum_congr fun d => ?_
  rw [← pow_mul]

/-! ## ★出典の紐付け(`.src`) -/

def summable_double_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——二重族の可和性)",
    sectionId := "genell-def-3-3" }

def tsum_double_family_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——s1 の二重和への展開)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
