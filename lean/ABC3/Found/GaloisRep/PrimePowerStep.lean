import ABC3.Found.GaloisRep.TorsionGroup

/-!
# Galois (G1) 第 68 ブロック —— **★★★★★素数冪の段の部品**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★数え上げを使わない道が見つかった

§9-400 では「位数 `p^k` の元の個数」と「悪い `b` の個数」を比べる予定だったが、
★**第一同型定理だけで済む**ことが分かった:

    #A = #range(p^{k−1}·) · #ker(p^{k−1}·) = #range · p^{2k−2} = p^{2k}
      ⟹ #range(p^{k−1}·) = p²
    range(p^{k−1}·) ⊆ ker(p·)  かつ  #ker(p·) = p²
      ⟹ **range(p^{k−1}·) = ker(p·) = A[p]**

★★したがって `A[p]` の元はすべて `p^{k−1}b` の形に書ける
——`⟨a⟩` の外の `c ∈ A[p]` を取れば、その持ち上げ `b` が求めるものである。

## ★★★機構

| 段 | 内容 |
|---|---|
| `card_range_mul_card_ker` | ★第一同型定理 + Lagrange |
| `range_eq_ker` | ★★★**`p^{k−1}A = A[p]`** |
| `exists_high_order` | ★`p^{k−1}•a ≠ 0` となる `a` の存在 |
| `addOrderOf_eq_prime_pow` | ★位数が `p^k` であること |

## ★★本ブロックで取れるもの

上の 4 つ。★★★残るのは `c ∈ A[p] \ ⟨a⟩` の存在と独立性の検証である。
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★**`m` 倍写像**。 -/
def nsmulHom (A : Type u) [AddCommGroup A] (m : ℕ) : A →+ A :=
  AddMonoidHom.mk' (fun x => m • x) (by intro a b; simp [smul_add])

theorem nsmulHom_apply {A : Type u} [AddCommGroup A] (m : ℕ) (x : A) :
    nsmulHom A m x = m • x := rfl

theorem mem_ker_nsmulHom {A : Type u} [AddCommGroup A] {m : ℕ} {x : A} :
    x ∈ (nsmulHom A m).ker ↔ m • x = 0 := Iff.rfl

/-- ★**第一同型定理 + Lagrange**。 -/
theorem card_range_mul_card_ker {A : Type u} [AddCommGroup A] [Finite A] (m : ℕ) :
    Nat.card (nsmulHom A m).range * Nat.card (nsmulHom A m).ker = Nat.card A := by
  rw [AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup (nsmulHom A m).ker,
    Nat.card_congr (QuotientAddGroup.quotientKerEquivRange (nsmulHom A m)).toEquiv]

/-- ★★★**`p^{k−1}A = A[p]`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★個数が両方 `p²` で片方が他方に含まれるから。 -/
theorem range_eq_ker {A : Type u} [AddCommGroup A] [Finite A] (p k : ℕ) (hp : 1 < p) (hk : 1 ≤ k)
    (hexp : ∀ x : A, (p ^ k) • x = 0)
    (hcardA : Nat.card A = (p ^ k) ^ 2)
    (hcardk1 : Nat.card (nsmulHom A (p ^ (k - 1))).ker = (p ^ (k - 1)) ^ 2)
    (hcard1 : Nat.card (nsmulHom A p).ker = p ^ 2) :
    (nsmulHom A (p ^ (k - 1))).range = (nsmulHom A p).ker := by
  have hle : (nsmulHom A (p ^ (k - 1))).range ≤ (nsmulHom A p).ker := by
    rintro _ ⟨y, rfl⟩
    rw [AddMonoidHom.mem_ker, nsmulHom_apply, nsmulHom_apply, smul_smul,
      show p * p ^ (k - 1) = p ^ k by rw [← pow_succ']; congr 1; omega]
    exact hexp y
  have hrng := card_range_mul_card_ker (A := A) (p ^ (k - 1))
  rw [hcardk1, hcardA] at hrng
  have hcr : Nat.card (nsmulHom A (p ^ (k - 1))).range = p ^ 2 := by
    have hsplit : ((p ^ k) ^ 2 : ℕ) = p ^ 2 * (p ^ (k - 1)) ^ 2 := by
      rw [← mul_pow, ← pow_succ']
      congr 2
      omega
    rw [hsplit] at hrng
    exact Nat.eq_of_mul_eq_mul_right (by positivity) hrng
  exact AddSubgroup.eq_of_le_of_card_ge hle (by rw [hcr, hcard1])

/-- ★**`p^{k−1}•a ≠ 0` となる `a` がある**。 -/
theorem exists_high_order {A : Type u} [AddCommGroup A] [Finite A] (p k : ℕ) (hk : 1 ≤ k)
    (hcardA : Nat.card A = (p ^ k) ^ 2)
    (hcardk1 : Nat.card (nsmulHom A (p ^ (k - 1))).ker = (p ^ (k - 1)) ^ 2)
    (hp : 1 < p) :
    ∃ a : A, (p ^ (k - 1)) • a ≠ 0 := by
  by_contra hcon
  push_neg at hcon
  have htop : (nsmulHom A (p ^ (k - 1))).ker = ⊤ := by
    ext x; simp only [AddSubgroup.mem_top, iff_true, AddMonoidHom.mem_ker, nsmulHom_apply]
    exact hcon x
  rw [htop] at hcardk1
  rw [Nat.card_congr (AddSubgroup.topEquiv (G := A)).toEquiv, hcardA] at hcardk1
  have hlt : ((p ^ (k - 1)) ^ 2 : ℕ) < (p ^ k) ^ 2 := by
    refine Nat.pow_lt_pow_left ?_ (by norm_num)
    exact Nat.pow_lt_pow_right hp (by omega)
  omega

/-- ★**位数が `p^k` であること**。 -/
theorem addOrderOf_eq_prime_pow {A : Type u} [AddCommGroup A] (p k : ℕ) (hp : p.Prime)
    (hk : 1 ≤ k) (a : A) (h1 : (p ^ k) • a = 0) (h2 : (p ^ (k - 1)) • a ≠ 0) :
    addOrderOf a = p ^ k := by
  have hdvd : addOrderOf a ∣ p ^ k := addOrderOf_dvd_of_nsmul_eq_zero h1
  obtain ⟨j, hj, hord⟩ := (Nat.dvd_prime_pow hp).1 hdvd
  rcases Nat.lt_or_ge j k with hjk | hjk
  · exfalso
    refine h2 ?_
    have hdd : addOrderOf a ∣ p ^ (k - 1) := by rw [hord]; exact pow_dvd_pow p (by omega)
    exact addOrderOf_dvd_iff_nsmul_eq_zero.1 hdd
  · rw [hord]; congr 1; omega

/-! ## ★出典の紐付け(`.src`) -/

def range_eq_ker.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(素数冪捩れの構造——p^{k−1}A = A[p])",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
