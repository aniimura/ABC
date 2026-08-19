import ABC3.Found.GaloisRep.TateExpand

/-!
# Galois (G6) 第 109 ブロック —— **★★★★★`I` 進級数の並べ替え**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★葉 (b) には**二重和の入れ替え**が要る

第 108 ブロックで `f(t) = ∑_{n≥0} n tⁿ` を得た。★これを `t = qᵐu` に入れると

    ∑_{m≥1} f(qᵐu) = ∑_{m≥1} ∑_{d≥1} d·q^{md}·u^d
                   = ∑_{N≥1} (∑_{d|N} d·u^d) q^N      ← ★★入れ替え

が古典的な q 展開である。★★★この入れ替えを支える道具を積む。

## ★★★★道具は 3 つ

| 道具 | 内容 |
|---|---|
| `adicSum_sub_partialSum_of_tail` | ★★★**尾の細かい評価**——添字 `N` 以降が `I^k` なら尾も `I^k` |
| `adicSum_reindex_mul` | ★★★★**添字の付け替え** `d ↦ m·d` |
| `adicSum_fubini` | ★★★★★**二重和の入れ替え** |

★`adicSum` の既定の評価は「和 − 部分和 ∈ `I^N`」しか与えないので、
`m·d` のような**粗い添字**を扱うには足りない。
★★`adicSum_sub_partialSum_of_tail` はそこを埋める——
「添字 `N` 以降の項がすべて `I^k` に入る」という**項ごとの情報**から尾を評価する。

## ★★★★★Fubini の形

対角より下が消える(`n ≤ m ⟹ c m n = 0`)二重級数について

    ∑_m (∑_n c m n)  =  ∑_n (∑_{m<n} c m n)

★これは `∑_{m≥1} f(qᵐu)` の各項が `q^{m·d}`(`m·d ≥ m`)しか出さないことに対応する。
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★尾の細かい評価 -/

theorem partialSum_sub_partialSum {a : ℕ → R} {m n : ℕ} (hmn : m ≤ n) :
    partialSum a n - partialSum a m = ∑ k ∈ Finset.Ico m n, a k := by
  rw [partialSum, partialSum, Finset.range_eq_Ico, Finset.range_eq_Ico,
    ← Finset.sum_Ico_consecutive a (Nat.zero_le m) hmn]
  ring

/-- ★★★**尾の細かい評価**——添字 `N` 以降の項がすべて `I^k` に入るなら、
和と部分和の差も `I^k` に入る。

★既定の評価(`I^N`)では足りない場面のための道具である。 -/
theorem adicSum_sub_partialSum_of_tail [IsAdicComplete I R] (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n)
    (N k : ℕ) (hk : ∀ n, N ≤ n → a n ∈ I ^ k) :
    adicSum a ha - partialSum a N ∈ I ^ k := by
  set M := max N k with hM
  have hNM : N ≤ M := le_max_left _ _
  have hkM : k ≤ M := le_max_right _ _
  have hspec := adicSum_spec a ha M
  rw [SModEq.sub_mem] at hspec
  have h1 : partialSum a M - adicSum a ha ∈ I ^ M := by simpa using hspec
  have h1' : adicSum a ha - partialSum a M ∈ I ^ k :=
    Ideal.pow_le_pow_right hkM (by simpa using neg_mem h1)
  have h2 : partialSum a M - partialSum a N ∈ I ^ k := by
    rw [partialSum_sub_partialSum hNM]
    exact Submodule.sum_mem _ (fun j hj => hk j (Finset.mem_Ico.1 hj).1)
  have h3 : adicSum a ha - partialSum a N
      = (adicSum a ha - partialSum a M) + (partialSum a M - partialSum a N) := by ring
  rw [h3]
  exact Submodule.add_mem _ h1' h2

/-! ## ★★★★並べ替えの判定法 -/

/-- ★★★★**並べ替えの判定法(合同版)**——各 `k` について、
十分先の項が `I^k` に入り、そこまでの部分和が `I^k` を法として一致するなら、和は等しい。 -/
theorem adicSum_eq_of_partial' [IsAdicComplete I R] (a b : ℕ → R)
    (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n)
    (h : ∀ k : ℕ, ∃ Na Nb : ℕ, (∀ n, Na ≤ n → a n ∈ I ^ k) ∧ (∀ n, Nb ≤ n → b n ∈ I ^ k)
      ∧ partialSum a Na - partialSum b Nb ∈ I ^ k) :
    adicSum a ha = adicSum b hb := by
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun k => ?_)
  rw [SModEq.sub_mem]
  obtain ⟨Na, Nb, hA, hB, hAB⟩ := h k
  have h1 : adicSum a ha - partialSum a Na ∈ I ^ k :=
    adicSum_sub_partialSum_of_tail a ha Na k hA
  have h2 : adicSum b hb - partialSum b Nb ∈ I ^ k :=
    adicSum_sub_partialSum_of_tail b hb Nb k hB
  have hsub : adicSum a ha - adicSum b hb
      = (adicSum a ha - partialSum a Na) - (adicSum b hb - partialSum b Nb)
        + (partialSum a Na - partialSum b Nb) := by ring
  have h4 : adicSum a ha - adicSum b hb ∈ I ^ k := by
    rw [hsub]
    exact Submodule.add_mem _ (Submodule.sub_mem _ h1 h2) hAB
  simpa using h4

/-- ★★★★**並べ替えの判定法(一致版)**。 -/
theorem adicSum_eq_of_partial [IsAdicComplete I R] (a b : ℕ → R)
    (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n)
    (h : ∀ k : ℕ, ∃ Na Nb : ℕ, (∀ n, Na ≤ n → a n ∈ I ^ k) ∧ (∀ n, Nb ≤ n → b n ∈ I ^ k)
      ∧ partialSum a Na = partialSum b Nb) :
    adicSum a ha = adicSum b hb := by
  refine adicSum_eq_of_partial' a b ha hb (fun k => ?_)
  obtain ⟨Na, Nb, hA, hB, hAB⟩ := h k
  refine ⟨Na, Nb, hA, hB, ?_⟩
  rw [hAB, sub_self]
  exact Submodule.zero_mem _

/-! ## ★★★★添字の付け替え -/

theorem sum_range_mul_dvd (m : ℕ) (hm : 1 ≤ m) (g : ℕ → R) (k : ℕ) :
    ∑ n ∈ Finset.range (m * k), (if m ∣ n then g (n / m) else 0)
      = ∑ d ∈ Finset.range k, g d := by
  induction k with
  | zero => simp
  | succ j ih =>
    have hsplit : m * (j + 1) = m * j + m := by ring
    rw [hsplit, Finset.sum_range_add, ih, Finset.sum_range_succ]
    congr 1
    refine (Finset.sum_eq_single 0 ?_ ?_).trans ?_
    · intro i hi hi0
      have him : i < m := Finset.mem_range.1 hi
      have hnd : ¬ (m ∣ (m * j + i)) := by
        intro h
        have hdi : m ∣ i := (Nat.dvd_add_right ⟨j, by ring⟩).mp h
        have := Nat.le_of_dvd (Nat.pos_of_ne_zero hi0) hdi
        omega
      simp [hnd]
    · intro h0
      exact absurd (Finset.mem_range.2 hm) h0
    · have hd : m ∣ (m * j + 0) := ⟨j, by ring⟩
      rw [if_pos hd]
      congr 1
      rw [add_zero]
      exact Nat.mul_div_cancel_left j (by omega)

theorem reindex_mem (m : ℕ) {g : ℕ → R} (hg : ∀ d, g d ∈ I ^ (m * d)) (n : ℕ) :
    (if m ∣ n then g (n / m) else 0) ∈ I ^ n := by
  by_cases h : m ∣ n
  · rw [if_pos h]
    have hmem : g (n / m) ∈ I ^ (m * (n / m)) := hg (n / m)
    rwa [Nat.mul_div_cancel' h] at hmem
  · rw [if_neg h]
    exact Submodule.zero_mem _

/-- ★★★★**添字の付け替え**——`g d ∈ I^{m·d}` なら
`∑_d g d` は `k = m·d` で並べ直したものに等しい。 -/
theorem adicSum_reindex_mul [IsAdicComplete I R] (m : ℕ) (hm : 1 ≤ m) (g : ℕ → R)
    (hg : ∀ d, g d ∈ I ^ (m * d)) :
    adicSum g (fun d => Ideal.pow_le_pow_right (Nat.le_mul_of_pos_left d (by omega)) (hg d))
      = adicSum (fun n => if m ∣ n then g (n / m) else 0) (reindex_mem m hg) := by
  refine adicSum_eq_of_partial _ _ _ _ (fun k => ⟨k, m * k, ?_, ?_, ?_⟩)
  · intro n hn
    exact Ideal.pow_le_pow_right hn
      (Ideal.pow_le_pow_right (Nat.le_mul_of_pos_left n (by omega)) (hg n))
  · intro n hn
    have hk : k ≤ n := le_trans (Nat.le_mul_of_pos_left k (by omega)) hn
    exact Ideal.pow_le_pow_right hk (reindex_mem m hg n)
  · rw [partialSum, partialSum]
    exact (sum_range_mul_dvd m hm g k).symm

/-! ## ★★★★★Fubini -/

theorem adicSum_mem_pow_of_vanish [IsAdicComplete I R] (c : ℕ → ℕ → R)
    (hc : ∀ m n, c m n ∈ I ^ n) (hz : ∀ m n, n ≤ m → c m n = 0) (m : ℕ) :
    adicSum (c m) (hc m) ∈ I ^ m := by
  have hp : partialSum (c m) (m + 1) = 0 := by
    rw [partialSum]
    exact Finset.sum_eq_zero (fun n hn => hz m n (by
      have := Finset.mem_range.1 hn
      omega))
  have hspec := adicSum_spec (c m) (hc m) (m + 1)
  rw [SModEq.sub_mem] at hspec
  have h1 : partialSum (c m) (m + 1) - adicSum (c m) (hc m) ∈ I ^ (m + 1) := by simpa using hspec
  rw [hp] at h1
  have h2 : adicSum (c m) (hc m) ∈ I ^ (m + 1) := by simpa using neg_mem h1
  exact Ideal.pow_le_pow_right (Nat.le_succ m) h2

theorem diag_mem (c : ℕ → ℕ → R) (hc : ∀ m n, c m n ∈ I ^ n) (n : ℕ) :
    (∑ m ∈ Finset.range n, c m n) ∈ I ^ n :=
  Submodule.sum_mem _ (fun m _ => hc m n)

/-- ★★★★★**`I` 進 Fubini**——対角より下が消える二重級数は順序を入れ替えられる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem adicSum_fubini [IsAdicComplete I R] (c : ℕ → ℕ → R)
    (hc : ∀ m n, c m n ∈ I ^ n) (hz : ∀ m n, n ≤ m → c m n = 0) :
    adicSum (fun m => adicSum (c m) (hc m)) (adicSum_mem_pow_of_vanish c hc hz)
      = adicSum (fun n => ∑ m ∈ Finset.range n, c m n) (diag_mem c hc) := by
  refine adicSum_eq_of_partial' _ _ _ _ (fun k => ⟨k, k, ?_, ?_, ?_⟩)
  · intro m hm
    exact Ideal.pow_le_pow_right hm (adicSum_mem_pow_of_vanish c hc hz m)
  · intro n hn
    exact Ideal.pow_le_pow_right hn (diag_mem c hc n)
  · have hext : ∀ n ∈ Finset.range k,
        ∑ m ∈ Finset.range n, c m n = ∑ m ∈ Finset.range k, c m n := by
      intro n hn
      have hnk : n ≤ k := le_of_lt (Finset.mem_range.1 hn)
      refine Finset.sum_subset
        (fun x hx => Finset.mem_range.2 (lt_of_lt_of_le (Finset.mem_range.1 hx) hnk)) ?_
      intro m _ hmn
      exact hz m n (by
        have := Finset.mem_range.not.1 hmn
        omega)
    have hkey : partialSum (fun n => ∑ m ∈ Finset.range n, c m n) k
        = ∑ m ∈ Finset.range k, partialSum (c m) k := by
      rw [partialSum, Finset.sum_congr rfl hext, Finset.sum_comm]
      exact Finset.sum_congr rfl (fun m _ => rfl)
    have hkey2 : partialSum (fun m => adicSum (c m) (hc m)) k
        = ∑ m ∈ Finset.range k, adicSum (c m) (hc m) := rfl
    rw [hkey, hkey2, ← Finset.sum_sub_distrib]
    refine Submodule.sum_mem _ (fun m _ => ?_)
    have hspec := adicSum_spec (c m) (hc m) k
    rw [SModEq.sub_mem] at hspec
    have h1 : partialSum (c m) k - adicSum (c m) (hc m) ∈ I ^ k := by simpa using hspec
    simpa using neg_mem h1

/-! ## ★出典の紐付け(`.src`) -/

def adicSum_fubini.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 級数の二重和の入れ替え——葉 (b) の道具)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
