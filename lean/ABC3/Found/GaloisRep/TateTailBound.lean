import ABC3.Found.GaloisRep.TateSplitNat

/-!
# Galois (G6) 第 227 ブロック —— **★★★★★尾の評価**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★切り詰めとの差を測る

第 226 で解析側の `X`・`Y` が「頭 + `ℕ` の尾」の形になった。形式側の切り詰め
(第 222 の `tateXtrunc`)は**尾を `n` 項で止めたもの**である。両者の差は
**残りの尾**であり、これが `O(‖q‖^{n+1})` であることを示す。

★★★**`‖q‖ ≤ 1/2` を仮定に置く**。第 218 の評価 `‖f(t)‖ ≤ 4‖t‖` は `‖t‖ ≤ 1/2` で
しか使えないが、`‖q‖ ≤ 1/2` かつ `‖c‖ ≤ 1` なら `‖q^{m+1}c‖ ≤ 1/2` が出る。
★`q → 0` を見るのだから `‖q‖ ≤ 1/2` は制限にならない。

## ★★★★★一般の尾の評価

各項が `K·r^{m+1}` で押さえられ `r ≤ 1/2` なら:

    ‖(∑_{m≥0} f(m)) − (∑_{m<n} f(m))‖ = ‖∑_{m≥0} f(m+n)‖ ≤ K·r^{n+1}/(1−r) ≤ 2K·r^{n+1}

★`1/(1−r) ≤ 2` は `r ≤ 1/2` から出る。これが `norm_tsum_sub_partialSum_le` である。

| 項 | `‖t‖ ≤ 1/2` での評価 | 尾の評価 |
|---|---|---|
| `f(t) = t/(1−t)²` | `≤ 4‖t‖`(第 218) | `8‖q‖^{n+1}` |
| `g(t) = t²/(1−t)³` | `≤ 8‖t‖`(`2g = h − f` から) | `16‖q‖^{n+1}` |

## ★★★`s₁` の側

`X` の定義には `−2 s₁(q)` が入っている。解析側では `∑_{n∈ℤ} f(q^{−n})` の形で現れるが、
`f(q^{−m}) = f(qᵐ)`(反転不変)かつ `f(q⁰) = f(1) = 0` なので

    ∑_{n∈ℤ} f(q^{−n}) = 2·∑_{m≥1} f(qᵐ)

★★`f(1) = 0` は Lean の `Ring.inverse 0 = 0` の規約から出る——`n = 0` の段が
**自動的に落ちる**のはここでも効いている(第 218 の定数項 `1/12` と同じ仕掛け)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `norm_tsum_sub_partialSum_le` | ★★★★★**一般の尾の評価** |
| `norm_tateYterm_le_of_small` | ★`‖g(t)‖ ≤ 8‖t‖` |
| `norm_pow_succ_mul_le_half` | ★`‖q^{m+1}c‖ ≤ 1/2` |
| `norm_tateXtail_sub_partialSum_le` | ★★★★★`X` の尾は `8‖q‖^{n+1}` 以下 |
| `norm_tateYtail_sub_partialSum_le` | ★★★★★`Y` の尾は `16‖q‖^{n+1}` 以下 |
| `tsum_int_tateXterm_one` | ★★★★★`s₁` の側の ℤ 和は `ℕ` 側の 2 倍 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★★★★一般の尾の評価 -/

/-- ★★★★★**尾の一般評価**——各項が `K·r^{m+1}` で押さえられ `r ≤ 1/2` なら、
`n` 次で切り詰めた尾は `2K·r^{n+1}` 以下。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_tsum_sub_partialSum_le (f : ℕ → ℂ) (r K : ℝ) (hr0 : 0 ≤ r) (hr : r ≤ 1 / 2)
    (hK : 0 ≤ K) (hf : Summable f) (hb : ∀ m, ‖f m‖ ≤ K * r ^ (m + 1)) (n : ℕ) :
    ‖(∑' m : ℕ, f m) - partialSum f n‖ ≤ 2 * K * r ^ (n + 1) := by
  have hr1 : r < 1 := by linarith
  have hshift := hf.sum_add_tsum_nat_add n
  have heq : (∑' m : ℕ, f m) - partialSum f n = ∑' m : ℕ, f (m + n) := by
    rw [partialSum, ← hshift]
    ring
  rw [heq]
  have hgsum : Summable fun m : ℕ => K * r ^ (n + 1) * r ^ m :=
    ((summable_geometric_of_lt_one hr0 hr1).mul_left _)
  have hbb : ∀ m : ℕ, ‖f (m + n)‖ ≤ K * r ^ (n + 1) * r ^ m := by
    intro m
    refine (hb (m + n)).trans (le_of_eq ?_)
    rw [show m + n + 1 = (n + 1) + m by omega, pow_add]
    ring
  have hns : Summable fun m : ℕ => ‖f (m + n)‖ :=
    Summable.of_norm_bounded hgsum (fun m => by rw [norm_norm]; exact hbb m)
  refine (norm_tsum_le_tsum_norm hns).trans ?_
  refine (hns.tsum_le_tsum hbb hgsum).trans ?_
  rw [(summable_geometric_of_lt_one hr0 hr1).tsum_mul_left,
    tsum_geometric_of_lt_one hr0 hr1]
  have hinv : (1 - r)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]
    linarith
  have hpos : 0 ≤ K * r ^ (n + 1) := by positivity
  nlinarith

/-! ## ★項ごとの評価 -/

/-- ★`‖t‖ ≤ 1/2` なら `‖g(t)‖ ≤ 8‖t‖`(`2g = h − f` から)。 -/
theorem norm_tateYterm_le_of_small {t : ℂ} (ht : ‖t‖ ≤ 1 / 2) :
    ‖tateYterm t‖ ≤ 8 * ‖t‖ := by
  have hd := norm_tateDterm_le_of_small ht
  have hx := norm_tateXterm_le_of_small ht
  have hy : tateYterm t = (tateDterm t - tateXterm t) / 2 := tateYterm_eq_half t
  rw [hy, norm_div]
  have hsub : ‖tateDterm t - tateXterm t‖ ≤ ‖tateDterm t‖ + ‖tateXterm t‖ := norm_sub_le _ _
  have h2 : ‖(2 : ℂ)‖ = 2 := by simp
  rw [h2]
  linarith

theorem norm_pow_succ_mul_le_half {q c : ℂ} (hq : ‖q‖ ≤ 1 / 2) (hc : ‖c‖ ≤ 1) (m : ℕ) :
    ‖q ^ (m + 1) * c‖ ≤ 1 / 2 := by
  rw [norm_mul, norm_pow]
  have h1 : ‖q‖ ^ (m + 1) ≤ (1 / 2 : ℝ) ^ (m + 1) :=
    pow_le_pow_left₀ (norm_nonneg q) hq (m + 1)
  have h2 : ((1 : ℝ) / 2) ^ (m + 1) ≤ (1 / 2 : ℝ) ^ 1 :=
    pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega)
  have h3 : ‖q‖ ^ (m + 1) ≤ 1 / 2 := by
    rw [pow_one] at h2
    linarith
  nlinarith [norm_nonneg c, pow_nonneg (norm_nonneg q) (m + 1)]

theorem norm_pow_succ_lt_one {q : ℂ} (hq : ‖q‖ < 1) (m : ℕ) : ‖q ^ (m + 1)‖ < 1 := by
  rw [norm_pow]
  exact pow_lt_one₀ (norm_nonneg q) hq (Nat.succ_ne_zero m)

/-! ## ★★★★★尾の評価 -/

/-- ★★★★★**`X` の尾の評価**——`‖q‖ ≤ 1/2` なら切り詰めとの差は `8‖q‖^{n+1}` 以下。 -/
theorem norm_tateXtail_sub_partialSum_le (q c : ℂ) (hq : ‖q‖ ≤ 1 / 2) (hc : ‖c‖ ≤ 1) (n : ℕ) :
    ‖(∑' m : ℕ, tateXterm (q ^ (m + 1) * c))
        - partialSum (fun m => tateXterm (q ^ (m + 1) * c)) n‖ ≤ 8 * ‖q‖ ^ (n + 1) := by
  have hq1 : ‖q‖ < 1 := by linarith
  have hb : ∀ m : ℕ, ‖tateXterm (q ^ (m + 1) * c)‖ ≤ 4 * ‖q‖ ^ (m + 1) := by
    intro m
    refine (norm_tateXterm_le_of_small (norm_pow_succ_mul_le_half hq hc m)).trans ?_
    rw [norm_mul, norm_pow]
    nlinarith [pow_nonneg (norm_nonneg q) (m + 1), norm_nonneg c]
  have h := norm_tsum_sub_partialSum_le (fun m => tateXterm (q ^ (m + 1) * c)) ‖q‖ 4
    (norm_nonneg q) hq (by norm_num) (summable_nat_tateXterm_succ q c hq1) hb n
  linarith [h]

/-- ★★★★★**`Y` の尾の評価**。 -/
theorem norm_tateYtail_sub_partialSum_le (q c : ℂ) (hq : ‖q‖ ≤ 1 / 2) (hc : ‖c‖ ≤ 1) (n : ℕ) :
    ‖(∑' m : ℕ, tateYterm (q ^ (m + 1) * c))
        - partialSum (fun m => tateYterm (q ^ (m + 1) * c)) n‖ ≤ 16 * ‖q‖ ^ (n + 1) := by
  have hq1 : ‖q‖ < 1 := by linarith
  have hb : ∀ m : ℕ, ‖tateYterm (q ^ (m + 1) * c)‖ ≤ 8 * ‖q‖ ^ (m + 1) := by
    intro m
    refine (norm_tateYterm_le_of_small (norm_pow_succ_mul_le_half hq hc m)).trans ?_
    rw [norm_mul, norm_pow]
    nlinarith [pow_nonneg (norm_nonneg q) (m + 1), norm_nonneg c]
  have h := norm_tsum_sub_partialSum_le (fun m => tateYterm (q ^ (m + 1) * c)) ‖q‖ 8
    (norm_nonneg q) hq (by norm_num) (summable_nat_tateYterm_succ q c hq1) hb n
  linarith [h]

/-! ## ★★★`s₁` の側 -/

/-- ★★★★★**`s₁` の側の ℤ 和は `ℕ` 側の 2 倍**——`f(q^{−m}) = f(qᵐ)`、`f(1) = 0`。 -/
theorem tsum_int_tateXterm_one (q : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) :
    (∑' n : ℤ, tateXterm (q ^ (-n))) = 2 * ∑' m : ℕ, tateXterm (q ^ (m + 1)) := by
  have hpos : ∀ m : ℕ, tateXterm (q ^ (-((m : ℤ) + 1))) = tateXterm (q ^ (m + 1)) := by
    intro m
    have hne0 : q ^ (m + 1) ≠ 0 := pow_ne_zero _ hq0
    have hne1 : q ^ (m + 1) ≠ 1 := ne_one_of_norm_lt_one (norm_pow_succ_lt_one hq m)
    rw [show (-((m : ℤ) + 1)) = -((m + 1 : ℕ) : ℤ) by push_cast; ring, zpow_neg, zpow_natCast,
      tateXterm_inv hne0 hne1]
  have hneg : ∀ m : ℕ, tateXterm (q ^ (-(-((m : ℤ) + 1)))) = tateXterm (q ^ (m + 1)) := by
    intro m
    rw [neg_neg, show ((m : ℤ) + 1) = ((m + 1 : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
  have hs : Summable fun m : ℕ => tateXterm (q ^ (m + 1)) := by
    refine (summable_nat_tateXterm_succ q 1 hq).congr fun m => ?_
    rw [mul_one]
  have h1 : Summable fun m : ℕ => tateXterm (q ^ (-((m : ℤ) + 1))) := by
    refine hs.congr fun m => ?_
    rw [hpos m]
  have h2 : Summable fun m : ℕ => tateXterm (q ^ (-(-((m : ℤ) + 1)))) := by
    refine hs.congr fun m => ?_
    rw [hneg m]
  have key := tsum_of_add_one_of_neg_add_one (f := fun n : ℤ => tateXterm (q ^ (-n))) h1 h2
  rw [key]
  simp only [hpos, hneg, neg_zero, zpow_zero, tateXterm_one_complex]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def norm_tsum_sub_partialSum_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——尾の一般評価)",
    sectionId := "genell-def-3-3" }

def norm_tateXtail_sub_partialSum_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の尾の評価)",
    sectionId := "genell-def-3-3" }

def tsum_int_tateXterm_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——s1 の側の ℤ 和)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
