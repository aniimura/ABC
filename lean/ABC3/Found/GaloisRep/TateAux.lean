import ABC3.Found.GaloisRep.TateAddEquiv

/-!
# Galois (G6) 第 299 ブロック —— **★★★★★★★★補助母数は取れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> 有限個(6 個の点と 3 個の平方)を避ける類 `e` が存在する(`exists_aux`)

★★★倍化(第 298)に要る補助母数がこれである。

## ★★★★★★具体的な無限族 `1 + q^{n+1}`

抽象的な「`Kˣ/q^ℤ` は無限」を使わず、**明示的な族**を取る:

    e_n := [1 + q^{n+1}]        (`n : ℕ`)

★どれも `R` の単元なので付値は 0(`vAdd_auxUnit`)。

## ★★★★★★族は単射、平方は高々 2 対 1

| 主張 | 根拠 |
|---|---|
| `e_n = e_m ⟹ n = m` | 付値 0 なので `q^k` の因子は `k = 0`、あとは `q^n` の単射性 |
| `e_n² = e_m² ⟹ 1 + q^{n+1} = ±(1 + q^{m+1})` | 体だから `(A−B)(A+B) = 0` |
| 平方の逆像は高々 2 点 | 3 点あれば `−` の枝が 2 回 ⟹ `q^{m+1} = q^{l+1}` で矛盾 |

★★★★**「符号の自由度が 2 つある」ことが、そのまま「高々 2 対 1」になっている。**

## ★★★★数え上げ

`n ∈ {0, …, 12}` のうち
* 6 個の点に当たる `n` は高々 6 個(単射だから)
* 3 個の平方に当たる `n` は高々 `2 × 3 = 6` 個(2 対 1 だから)

★合わせて高々 12 個。**13 個あるので 1 つは残る**。
★★抽象的な濃度の議論を避け、`Finset.card` の不等式だけで済ませた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pow_inj_of_mem` | ★★`q^n` は `n` について単射 |
| `auxUnit`・`auxCls` | ★★★補助母数の族 |
| `auxCls_inj` | ★★★★★★族は単射 |
| `auxUnit_sq_eq`・`auxCls_sq_not_three` | ★★★★★★平方は高々 2 対 1 |
| `exists_aux` | ★★★★★★★★**補助母数は取れる** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]

/-! ## ★★`q^n` は `n` について単射 -/

/-- ★★`q^n` は `n` について単射(`q ≠ 0`、`q` は非単元)。 -/
theorem pow_inj_of_mem (q : R) (hq : q ∈ I) (hq0 : q ≠ 0) {n m : ℕ}
    (h : q ^ n = q ^ m) : n = m := by
  by_contra hne
  rcases Nat.lt_or_ge n m with hlt | hge
  · obtain ⟨k, rfl⟩ : ∃ k, m = n + (k + 1) := ⟨m - n - 1, by omega⟩
    rw [pow_add] at h
    have h2 : q ^ n * (1 - q ^ (k + 1)) = 0 := by linear_combination h
    have hu : IsUnit (1 - q ^ (k + 1)) :=
      isUnit_one_sub (Ideal.pow_mem_of_mem I hq (k + 1) (Nat.succ_pos k))
    rcases mul_eq_zero.1 h2 with h3 | h3
    · exfalso
      rcases Nat.eq_zero_or_pos n with rfl | hpos
      · simp at h3
      · exact hq0 ((pow_eq_zero_iff hpos.ne').1 h3)
    · rw [h3] at hu
      exact not_isUnit_zero hu
  · obtain ⟨k, rfl⟩ : ∃ k, n = m + (k + 1) := ⟨n - m - 1, by omega⟩
    rw [pow_add] at h
    have h2 : q ^ m * (1 - q ^ (k + 1)) = 0 := by linear_combination -h
    have hu : IsUnit (1 - q ^ (k + 1)) :=
      isUnit_one_sub (Ideal.pow_mem_of_mem I hq (k + 1) (Nat.succ_pos k))
    rcases mul_eq_zero.1 h2 with h3 | h3
    · exfalso
      rcases Nat.eq_zero_or_pos m with rfl | hpos
      · simp at h3
      · exact hq0 ((pow_eq_zero_iff hpos.ne').1 h3)
    · rw [h3] at hu
      exact not_isUnit_zero hu

/-! ## ★★★補助母数の族 -/

section Aux

variable {K : Type} [Field K] [Algebra R K]

theorem isUnit_one_add_pow (S : TateSetup R I K) (n : ℕ) : IsUnit (1 + S.q ^ (n + 1)) :=
  isUnit_add_mem (I := I) isUnit_one (Ideal.pow_mem_of_mem I S.hq (n + 1) (Nat.succ_pos n))

/-- ★★★補助母数の族 `1 + q^{n+1}`。 -/
noncomputable def auxUnit (S : TateSetup R I K) (n : ℕ) : Kˣ :=
  Units.mk0 (algebraMap R K (1 + S.q ^ (n + 1)))
    (((isUnit_one_add_pow S n).map (algebraMap R K)).ne_zero)

noncomputable def auxCls (S : TateSetup R I K) (n : ℕ) : Kˣ ⧸ Subgroup.zpowers S.Q :=
  QuotientGroup.mk (auxUnit S n)

theorem vAdd_auxUnit (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t) (n : ℕ) :
    vAdd S.v (auxUnit S n) = 0 :=
  vAdd_eq_zero_of_isUnit S hvR (1 + S.q ^ (n + 1)) (isUnit_one_add_pow S n) _

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**補助母数の族は単射**。 -/
theorem auxCls_inj (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hq0 : S.q ≠ 0) {n m : ℕ} (h : auxCls S n = auxCls S m) : n = m := by
  rw [auxCls, auxCls, QuotientGroup.eq] at h
  obtain ⟨k, hk⟩ := h
  have hk' : S.Q ^ k = (auxUnit S n)⁻¹ * auxUnit S m := hk
  have hv : vAdd S.v ((auxUnit S n)⁻¹ * auxUnit S m) = 0 := by
    rw [vAdd_mul, vAdd_inv, vAdd_auxUnit S hvR n, vAdd_auxUnit S hvR m]
    omega
  have hvk : (k : ℤ) * vAdd S.v S.Q = 0 := by
    rw [← vAdd_zpow, hk']
    exact hv
  have hQ := S.hQ
  have hk0 : k = 0 := by
    rcases mul_eq_zero.1 hvk with h1 | h1
    · exact h1
    · omega
  rw [hk0, zpow_zero] at hk'
  have heq : auxUnit S n = auxUnit S m := inv_mul_eq_one.1 hk'.symm
  have h3 : algebraMap R K (1 + S.q ^ (n + 1)) = algebraMap R K (1 + S.q ^ (m + 1)) :=
    congrArg (fun u : Kˣ => (u : K)) heq
  have h4 : (1 : R) + S.q ^ (n + 1) = 1 + S.q ^ (m + 1) := S.hinj h3
  have h5 : S.q ^ (n + 1) = S.q ^ (m + 1) := by linear_combination h4
  have h6 := pow_inj_of_mem (I := I) S.q S.hq hq0 h5
  omega

/-! ## ★★★★★★平方は高々 2 対 1 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**平方が一致するのは `±` の 2 通りだけ**。 -/
theorem auxUnit_sq_eq (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    {n m : ℕ} (h : auxCls S n * auxCls S n = auxCls S m * auxCls S m) :
    (1 : R) + S.q ^ (n + 1) = 1 + S.q ^ (m + 1)
      ∨ (1 : R) + S.q ^ (n + 1) = -(1 + S.q ^ (m + 1)) := by
  have h2 : (QuotientGroup.mk (auxUnit S n * auxUnit S n) : Kˣ ⧸ Subgroup.zpowers S.Q)
      = QuotientGroup.mk (auxUnit S m * auxUnit S m) := by
    rw [QuotientGroup.mk_mul, QuotientGroup.mk_mul]
    exact h
  rw [QuotientGroup.eq] at h2
  obtain ⟨k, hk⟩ := h2
  have hk' : S.Q ^ k = (auxUnit S n * auxUnit S n)⁻¹ * (auxUnit S m * auxUnit S m) := hk
  have hv : vAdd S.v ((auxUnit S n * auxUnit S n)⁻¹ * (auxUnit S m * auxUnit S m)) = 0 := by
    rw [vAdd_mul, vAdd_inv, vAdd_mul, vAdd_mul, vAdd_auxUnit S hvR n, vAdd_auxUnit S hvR m]
    omega
  have hvk : (k : ℤ) * vAdd S.v S.Q = 0 := by
    rw [← vAdd_zpow, hk']
    exact hv
  have hQ := S.hQ
  have hk0 : k = 0 := by
    rcases mul_eq_zero.1 hvk with h1 | h1
    · exact h1
    · omega
  rw [hk0, zpow_zero] at hk'
  have heq : auxUnit S n * auxUnit S n = auxUnit S m * auxUnit S m := inv_mul_eq_one.1 hk'.symm
  have h3 : algebraMap R K ((1 + S.q ^ (n + 1)) * (1 + S.q ^ (n + 1)))
      = algebraMap R K ((1 + S.q ^ (m + 1)) * (1 + S.q ^ (m + 1))) := by
    rw [map_mul, map_mul]
    exact congrArg (fun u : Kˣ => (u : K)) heq
  have h4 : (1 + S.q ^ (n + 1)) * (1 + S.q ^ (n + 1))
      = (1 + S.q ^ (m + 1)) * (1 + S.q ^ (m + 1)) := S.hinj h3
  have h5 : ((1 + S.q ^ (n + 1)) - (1 + S.q ^ (m + 1)))
      * ((1 + S.q ^ (n + 1)) + (1 + S.q ^ (m + 1))) = 0 := by linear_combination h4
  rcases mul_eq_zero.1 h5 with h6 | h6
  · exact Or.inl (sub_eq_zero.1 h6)
  · exact Or.inr (by linear_combination h6)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**平方の逆像は高々 2 点**。 -/
theorem auxCls_sq_not_three (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hq0 : S.q ≠ 0) {n m l : ℕ} (hnm : n ≠ m) (hnl : n ≠ l) (hml : m ≠ l)
    (h1 : auxCls S n * auxCls S n = auxCls S m * auxCls S m)
    (h2 : auxCls S n * auxCls S n = auxCls S l * auxCls S l) : False := by
  have hA : ∀ {a b : ℕ}, (1 : R) + S.q ^ (a + 1) = 1 + S.q ^ (b + 1) → a = b := by
    intro a b h
    have h5 : S.q ^ (a + 1) = S.q ^ (b + 1) := by linear_combination h
    have h6 := pow_inj_of_mem (I := I) S.q S.hq hq0 h5
    omega
  rcases auxUnit_sq_eq S hvR h1 with e1 | e1
  · exact hnm (hA e1)
  rcases auxUnit_sq_eq S hvR h2 with e2 | e2
  · exact hnl (hA e2)
  · refine hml (hA ?_)
    linear_combination e1 - e2

/-! ## ★★★★★★★★補助母数は取れる -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**補助母数は取れる**——有限個の条件を避ける `e`。

★13 個の候補のうち、悪いものは高々 12 個。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_aux (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hq0 : S.q ≠ 0) (A B : Finset (Kˣ ⧸ Subgroup.zpowers S.Q))
    (hA : A.card ≤ 6) (hB : B.card ≤ 3) :
    ∃ n : ℕ, auxCls S n ∉ A ∧ auxCls S n * auxCls S n ∉ B := by
  classical
  by_contra hcon
  push_neg at hcon
  set T := Finset.range 13 with hT
  set T1 := T.filter (fun n => auxCls S n ∈ A) with hT1
  set T2 := T.filter (fun n => auxCls S n * auxCls S n ∈ B) with hT2
  have hcover : T ⊆ T1 ∪ T2 := by
    intro n hn
    by_cases h : auxCls S n ∈ A
    · exact Finset.mem_union_left _ (Finset.mem_filter.2 ⟨hn, h⟩)
    · exact Finset.mem_union_right _ (Finset.mem_filter.2 ⟨hn, hcon n h⟩)
  have hc1 : T1.card ≤ 6 := by
    refine le_trans (Finset.card_le_card_of_injOn (auxCls S) ?_ ?_) hA
    · intro n hn
      exact (Finset.mem_filter.1 hn).2
    · intro a _ b _ hab
      exact auxCls_inj S hvR hq0 hab
  have hc2 : T2.card ≤ 6 := by
    have hfib : ∀ a ∈ T2.image (fun n => auxCls S n * auxCls S n),
        (T2.filter (fun n => auxCls S n * auxCls S n = a)).card ≤ 2 := by
      intro a _
      by_contra hgt
      push_neg at hgt
      obtain ⟨x, y, z, hx, hy, hz, hxy, hxz, hyz⟩ := Finset.two_lt_card_iff.1 hgt
      have ex := (Finset.mem_filter.1 hx).2
      have ey := (Finset.mem_filter.1 hy).2
      have ez := (Finset.mem_filter.1 hz).2
      exact auxCls_sq_not_three S hvR hq0 hxy hxz hyz (ex.trans ey.symm) (ex.trans ez.symm)
    have h1 := Finset.card_le_mul_card_image (f := fun n => auxCls S n * auxCls S n) T2 2 hfib
    have h2 : (T2.image (fun n => auxCls S n * auxCls S n)).card ≤ B.card := by
      refine Finset.card_le_card ?_
      intro a ha
      obtain ⟨n, hn, rfl⟩ := Finset.mem_image.1 ha
      exact (Finset.mem_filter.1 hn).2
    omega
  have h3 : T.card ≤ T1.card + T2.card :=
    (Finset.card_le_card hcover).trans (Finset.card_union_le _ _)
  rw [hT, Finset.card_range] at h3
  omega

end Aux

/-! ## ★出典の紐付け(`.src`) -/

def auxCls.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——補助母数の族)",
    sectionId := "genell-def-3-3" }

def exists_aux.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——補助母数は取れる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
