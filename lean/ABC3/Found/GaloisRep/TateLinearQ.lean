import ABC3.Found.GaloisRep.TateAnnulus

/-!
# Galois (G6) 第 270 ブロック —— **★★★★★★`a₆` の 1 次の項が `q` を決める**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★制約 `a·w = q` を回復する道具

第 269 の環帯の解 `(a,w)` は制約 `a·w = q` を満たすとは限らない。回復の筋は

> `(x,y)` が `E_q` 上で `X(a,w,q) = x`、`Y(a,w,q) = y` なら `defect(a,w,q) = 0`。
> 一方 `defect(a,w,a·w) = 0`(葉 (b))。
> `defect` の `q` についての**主要項が `q − aw`** なので `q = aw`。

★その「主要項」を作るのが本ブロックである。

## ★★★★★★評価の 1 次の項

冪級数 `f`(定数項 0)について

    f(a) − f(b) − c₁·(a − b) ∈ I^{k+1}      (`a − b ∈ I^k`、`a, b ∈ I`)

(`evalAdic_sub_linear_mem`)。★2 次以上の項 `aⁿ − bⁿ` は
`(a−b)a^{n−1} + (a^{n−1}−b^{n−1})b` と分けるだけで `I^{k+1}` に入る
(`pow_sub_pow_mem_succ`)。**係数を触らずに済む**。

★配線:`partialEval` を `N = k+2` まで取り、`Finset.sum_erase_add` で `n = 1` の項だけ
抜き出す。残差 `evalAdic − partialEval ∈ I^{k+2} ⊆ I^{k+1}` は自動である。

## ★★`a₆` の 1 次の係数は `−1`

`coeff 1 tateA6 = −(5σ₃(1) + 7σ₅(1))/12 = −12/12 = −1`(在庫、第 `TateDelta` ブロック)。
したがって

    a₆(q) − a₆(q') + (q − q') ∈ I^{k+1}

★**`a₆` が `q` を 1 次で拾う**——これが制約の回復の芯である。
`a₄(q) = −5s₃(q)` の側は `∈ I^k` で十分(`defect` の中で `X ∈ I` が掛かるから)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pow_sub_pow_mem_succ` | ★★★`aⁿ − bⁿ ∈ I^{k+1}`(`n ≥ 2`) |
| `evalAdic_sub_linear_mem` | ★★★★★★**評価の 1 次の項** |
| `tateCurveAt_a6_sub_linear` | ★★★★★★**`a₆` の 1 次の項は `−q`** |
| `tateCurveAt_a4_sub_mem` | ★★`a₄` の差 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★2 次以上の差 -/

/-- ★★★**`aⁿ − bⁿ` は `n ≥ 2` なら 1 つ位が上がる**。 -/
theorem pow_sub_pow_mem_succ {k n : ℕ} {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (hab : a - b ∈ I ^ k) (hn : 2 ≤ n) : a ^ n - b ^ n ∈ I ^ (k + 1) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have hm : 1 ≤ m := by omega
  have h1 : a ^ (m + 1) - b ^ (m + 1) = (a - b) * a ^ m + (a ^ m - b ^ m) * b := by ring
  have ham : a ^ m ∈ I := by
    obtain ⟨j, rfl⟩ : ∃ j, m = j + 1 := ⟨m - 1, by omega⟩
    rw [pow_succ]
    exact Ideal.mul_mem_left _ _ ha
  have hdiff : a ^ m - b ^ m ∈ I ^ k :=
    Ideal.mem_of_dvd _ (sub_dvd_pow_sub_pow a b m) hab
  rw [h1, pow_succ]
  exact Ideal.add_mem _ (Ideal.mul_mem_mul hab ham) (Ideal.mul_mem_mul hdiff hb)

/-! ## ★★★★★★評価の 1 次の項 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★**評価の 1 次の項**——`f(a) − f(b) ≡ c₁·(a − b)`(位が 1 つ上がる)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem evalAdic_sub_linear_mem [IsAdicComplete I R] (f : PowerSeries ℤ)
    (hf0 : PowerSeries.coeff 0 f = 0) {a b : R} (ha : a ∈ I) (hb : b ∈ I) (k : ℕ)
    (hab : a - b ∈ I ^ k) :
    evalAdic f a ha - evalAdic f b hb - ((PowerSeries.coeff 1 f : ℤ) : R) * (a - b)
      ∈ I ^ (k + 1) := by
  set N := k + 2 with hN
  have hsa := evalAdic_spec (I := I) f a ha N
  have hsb := evalAdic_spec (I := I) f b hb N
  rw [SModEq.sub_mem] at hsa hsb
  have hsa' : partialEval f a N - evalAdic f a ha ∈ I ^ N := by simpa using hsa
  have hsb' : partialEval f b N - evalAdic f b hb ∈ I ^ N := by simpa using hsb
  have hup : ∀ x : R, x ∈ I ^ N → x ∈ I ^ (k + 1) := fun x hx =>
    Ideal.pow_le_pow_right (by omega) hx
  have hpart : partialEval f a N - partialEval f b N
      - ((PowerSeries.coeff 1 f : ℤ) : R) * (a - b) ∈ I ^ (k + 1) := by
    rw [partialEval, partialEval, ← Finset.sum_sub_distrib]
    have h1mem : (1 : ℕ) ∈ Finset.range N := by simp [hN]
    have hsplit : ∑ n ∈ Finset.range N,
        (((PowerSeries.coeff n f : ℤ) : R) * a ^ n - ((PowerSeries.coeff n f : ℤ) : R) * b ^ n)
        - ((PowerSeries.coeff 1 f : ℤ) : R) * (a - b)
        = ∑ n ∈ (Finset.range N).erase 1,
          (((PowerSeries.coeff n f : ℤ) : R) * a ^ n
            - ((PowerSeries.coeff n f : ℤ) : R) * b ^ n) := by
      rw [← Finset.sum_erase_add _ _ h1mem]
      simp only [pow_one]
      ring
    rw [hsplit]
    refine Ideal.sum_mem _ fun n hn => ?_
    have hn1 : n ≠ 1 := (Finset.mem_erase.1 hn).1
    rcases Nat.eq_zero_or_pos n with rfl | hpos
    · simp [hf0]
    · have hn2 : 2 ≤ n := by omega
      have hd : ((PowerSeries.coeff n f : ℤ) : R) * a ^ n
          - ((PowerSeries.coeff n f : ℤ) : R) * b ^ n
          = ((PowerSeries.coeff n f : ℤ) : R) * (a ^ n - b ^ n) := by ring
      rw [hd]
      exact Ideal.mul_mem_left _ _ (pow_sub_pow_mem_succ ha hb hab hn2)
  have hexp : evalAdic f a ha - evalAdic f b hb - ((PowerSeries.coeff 1 f : ℤ) : R) * (a - b)
      = -(partialEval f a N - evalAdic f a ha) + (partialEval f b N - evalAdic f b hb)
        + (partialEval f a N - partialEval f b N
          - ((PowerSeries.coeff 1 f : ℤ) : R) * (a - b)) := by ring
  rw [hexp]
  exact Ideal.add_mem _ (Ideal.add_mem _ (neg_mem (hup _ hsa')) (hup _ hsb')) hpart

/-! ## ★★★★★★`a₆` の 1 次の項 -/

theorem tateCurveAt_a4_eq [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₄ = evalAdic tateA4 q hq := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]

theorem tateCurveAt_a6_eq [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).a₆ = evalAdic tateA6 q hq := by
  simp [tateCurveAt, tateCurve, WeierstrassCurve.map, evalAdicHom]

/-- ★★★★★★**`a₆` の 1 次の項は `−q`**——これが制約 `a·w = q` の回復の芯である。 -/
theorem tateCurveAt_a6_sub_linear [IsAdicComplete I R] {q q' : R} (hq : q ∈ I) (hq' : q' ∈ I)
    (k : ℕ) (hqq : q - q' ∈ I ^ k) :
    (tateCurveAt q hq).a₆ - (tateCurveAt q' hq').a₆ + (q - q') ∈ I ^ (k + 1) := by
  have h := evalAdic_sub_linear_mem (I := I) tateA6 (by simp [tateA6]) hq hq' k hqq
  rw [coeff_one_tateA6] at h
  rw [tateCurveAt_a6_eq, tateCurveAt_a6_eq]
  have he : evalAdic tateA6 q hq - evalAdic tateA6 q' hq' + (q - q')
      = evalAdic tateA6 q hq - evalAdic tateA6 q' hq' - ((-1 : ℤ) : R) * (q - q') := by
    push_cast
    ring
  rw [he]
  exact h

theorem tateCurveAt_a4_sub_mem [IsAdicComplete I R] {q q' : R} (hq : q ∈ I) (hq' : q' ∈ I)
    (k : ℕ) (hqq : q - q' ∈ I ^ k) :
    (tateCurveAt q hq).a₄ - (tateCurveAt q' hq').a₄ ∈ I ^ k := by
  rw [tateCurveAt_a4_eq, tateCurveAt_a4_eq]
  exact evalAdic_sub_mem tateA4 hq hq' k hqq

/-! ## ★★★★★★★★形式逆関数定理（`I` 進版） -/

/-- ★★★★★★★★**`f(q) = q + O(q²)` なら `q ↦ f(q)` は `I` の上で単射**。

★道は `evalAdic_sub_linear_mem` だけである: `a − b ∈ I^k` なら

    `f(a) − f(b) − 1·(a − b) ∈ I^{k+1}`

なので `f(a) = f(b)` なら `a − b ∈ I^{k+1}`。帰納してすべての `k` で成り立ち、
Hausdorff で `a = b` 。

☆これが `Lemma 3.5` の残る義務のひとつ
「`q ↦ j(E_q)` の単射性」の中身である——`1/j = Δ/c₄³ = q + O(q²)`。 -/
theorem evalAdic_injective_of_coeff_one [IsAdicComplete I R] (f : PowerSeries ℤ)
    (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = 1)
    {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (h : evalAdic f a ha = evalAdic f b hb) : a = b := by
  have key : ∀ k : ℕ, a - b ∈ I ^ k := by
    intro k
    induction k with
    | zero => simp
    | succ n ih =>
        have hm := evalAdic_sub_linear_mem f hf0 ha hb n ih
        rw [h, sub_self, hf1] at hm
        simp only [Int.cast_one, one_mul, zero_sub] at hm
        exact neg_mem_iff.1 hm
  have hz : a - b = 0 := by
    refine IsHausdorff.haus' (I := I) (a - b) (fun n => ?_)
    rw [SModEq.sub_mem, sub_zero]
    simpa using key n
  exact sub_eq_zero.1 hz

/-! ## ★★★★★★★★★★形式逆関数定理——存在の側 -/

/-- ★★★★★★★★★★**`f(q) = t` は解ける**（`f = X + 高次`、`t ∈ I`）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★第 875 の `evalAdic_injective_of_coeff_one` の**対**である——
形式逆関数定理の「全射」の側。

☆道は逐次近似である。`q₀ = 0`、`q_{n+1} = q_n + (t − f(q_n))` とおくと

    `f(q_n) − t ∈ Iⁿ`

が帰納的に成り立つ（`evalAdic_sub_linear_mem` を 1 回使うだけ）。
`IsPrecomplete` で極限 `L` を取り、`IsHausdorff` で `f(L) = t` を結ぶ。

★これがあれば、`v_p(j) < 0` なる任意の `j` に対して
**Tate 母数 `q` が作れる**——分裂性を問わずに。 -/
theorem evalAdic_surjective_of_coeff_one [IsAdicComplete I R] (f : PowerSeries ℤ)
    (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = 1)
    {t : R} (ht : t ∈ I) :
    ∃ (q : R) (hq : q ∈ I), evalAdic f q hq = t := by
  classical
  -- ★逐次近似の列（`I` の元であることを抱えたまま）
  set step : {x : R // x ∈ I} → {x : R // x ∈ I} := fun x =>
    ⟨x.1 + (t - evalAdic f x.1 x.2),
      add_mem x.2 (sub_mem ht (evalAdic_mem f hf0 x.1 x.2))⟩ with hstep
  set g : ℕ → {x : R // x ∈ I} := fun n => step^[n] ⟨0, I.zero_mem⟩ with hg
  have hgsucc : ∀ n, g (n + 1) = step (g n) := by
    intro n
    rw [hg]
    exact Function.iterate_succ_apply' step n _
  -- ★不変量: `f(q_n) − t ∈ Iⁿ`
  have hinv : ∀ n, evalAdic f (g n).1 (g n).2 - t ∈ I ^ n := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
        have hdiff : (g (n + 1)).1 - (g n).1 ∈ I ^ n := by
          rw [hgsucc n, hstep]
          simpa using neg_mem ih
        have hlin := evalAdic_sub_linear_mem f hf0 (g (n + 1)).2 (g n).2 n hdiff
        rw [hf1] at hlin
        simp only [Int.cast_one, one_mul] at hlin
        have hval : (g (n + 1)).1 - (g n).1 = t - evalAdic f (g n).1 (g n).2 := by
          rw [hgsucc n, hstep]
          ring
        rw [hval] at hlin
        have : evalAdic f (g (n + 1)).1 (g (n + 1)).2 - t
            = evalAdic f (g (n + 1)).1 (g (n + 1)).2 - evalAdic f (g n).1 (g n).2
              - (t - evalAdic f (g n).1 (g n).2) := by ring
        rw [this]
        exact hlin
  -- ★Cauchy
  have hstepdiff : ∀ n, (g (n + 1)).1 - (g n).1 ∈ I ^ n := by
    intro n
    rw [hgsucc n, hstep]
    simpa using neg_mem (hinv n)
  have hcauchy : ∀ {m n : ℕ}, m ≤ n → (g n).1 - (g m).1 ∈ I ^ m := by
    intro m n hmn
    induction n with
    | zero =>
        have : m = 0 := Nat.le_zero.1 hmn
        subst this
        simp
    | succ k ih =>
        rcases Nat.lt_or_ge m (k + 1) with h | h
        · have hmk : m ≤ k := Nat.lt_succ_iff.1 h
          have h1 := ih hmk
          have h2 : (g (k + 1)).1 - (g k).1 ∈ I ^ m :=
            Ideal.pow_le_pow_right hmk (hstepdiff k)
          have : (g (k + 1)).1 - (g m).1
              = ((g (k + 1)).1 - (g k).1) + ((g k).1 - (g m).1) := by ring
          rw [this]
          exact add_mem h2 h1
        · have : m = k + 1 := le_antisymm hmn h
          subst this
          simp
  -- ★極限
  obtain ⟨L, hL⟩ := IsPrecomplete.prec' (I := I) (M := R) (fun n => (g n).1) (by
    intro m n hmn
    rw [SModEq.sub_mem]
    simpa using neg_mem (hcauchy hmn))
  have hLmem : L ∈ I := by
    have h1 := hL 1
    rw [SModEq.sub_mem] at h1
    have h2 : (g 1).1 - L ∈ I := by simpa using h1
    have : L = (g 1).1 - ((g 1).1 - L) := by ring
    rw [this]
    exact sub_mem (g 1).2 h2
  refine ⟨L, hLmem, ?_⟩
  -- ★`f(L) − t ∈ Iⁿ` すべての `n` で
  have hall : ∀ n : ℕ, evalAdic f L hLmem - t ∈ I ^ n := by
    intro n
    have hLn : L - (g n).1 ∈ I ^ n := by
      have := hL n
      rw [SModEq.sub_mem] at this
      simpa using neg_mem (by simpa using this : (g n).1 - L ∈ I ^ n)
    have hlin := evalAdic_sub_linear_mem f hf0 hLmem (g n).2 n hLn
    rw [hf1] at hlin
    simp only [Int.cast_one, one_mul] at hlin
    have hup : evalAdic f L hLmem - evalAdic f (g n).1 (g n).2 - (L - (g n).1) ∈ I ^ n :=
      Ideal.pow_le_pow_right (by omega) hlin
    have : evalAdic f L hLmem - t
        = (evalAdic f L hLmem - evalAdic f (g n).1 (g n).2 - (L - (g n).1))
          + (L - (g n).1) + (evalAdic f (g n).1 (g n).2 - t) := by ring
    rw [this]
    exact add_mem (add_mem hup hLn) (hinv n)
  have hz : evalAdic f L hLmem - t = 0 := by
    refine IsHausdorff.haus' (I := I) _ (fun n => ?_)
    rw [SModEq.sub_mem, sub_zero]
    simpa using hall n
  exact sub_eq_zero.1 hz

/-! ## ★出典の紐付け(`.src`) -/

def evalAdic_surjective_of_coeff_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—形式逆関数定理の存在の側。★無条件)",
    sectionId := "genell-def-3-3" }

def evalAdic_sub_linear_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——評価の 1 次の項)",
    sectionId := "genell-def-3-3" }

def tateCurveAt_a6_sub_linear.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——a6 の 1 次の項は -q)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
