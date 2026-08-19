import ABC3.Found.GaloisRep.PrimePowerStep

/-!
# Galois (G1) 第 69 ブロック —— **★★★★★★★素数冪の構造定理 `A ≃+ (ℤ/p^k)²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★数え上げなしで独立な 2 元が取れた

| 段 | 内容 |
|---|---|
| ★ | `p^{k−1}•a ≠ 0` となる `a`(第 68)——位数は `p^k` |
| ★★ | `⟨a⟩ ∩ A[p] ⊆ ⟨p^{k−1}a⟩` で個数 `p`、`#A[p] = p²` ⟹ **外に `c` がある** |
| ★★★ | `A[p] = p^{k−1}A`(第 68)⟹ `c = p^{k−1}b` と書ける |
| ★★★★ | `j•b ∈ ⟨a⟩ ⟹ p^k ∣ j`——Bezout で `gcd(j,p^k)•b ∈ ⟨a⟩`、`c ∈ ⟨a⟩` に反する |

★★★★★これで `i•a + j•b = 0 ⟹ p^k ∣ i かつ p^k ∣ j` が出て、
第 67 ブロックの `addEquiv_of_indep` により **`(ℤ/p^k)² ≃+ A`** ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `zsmul_mem_of_p_smul` | ★`⟨a⟩ ∩ A[p] ⊆ ⟨p^{k−1}a⟩` |
| `addOrderOf_high` | ★`p^{k−1}a` の位数は `p` |
| `exists_outside` | ★★`⟨a⟩` の外の `c ∈ A[p]` |
| `dvd_of_zsmul_mem` | ★★★★独立性の核 |
| `exists_indep_prime_pow` | ★★★★★★独立な 2 元 |
| `addEquiv_prime_pow` | ★★★★★★★**`A ≃+ (ℤ/p^k)²`** |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★`⟨a⟩` の `p` 捩れは `⟨p^{k−1}a⟩` に入る。 -/
theorem zsmul_mem_of_p_smul {A : Type u} [AddCommGroup A] (p k : ℕ) (hp : 0 < p) (hk : 1 ≤ k)
    (a : A) (hord : addOrderOf a = p ^ k) (i : ℤ) (hx : p • (i • a) = 0) :
    ∃ t : ℤ, t • ((p ^ (k - 1)) • a) = i • a := by
  have hz : ((p : ℤ) * i) • a = 0 := by
    rw [mul_zsmul]
    simpa using hx
  rw [← addOrderOf_dvd_iff_zsmul_eq_zero, hord] at hz
  have hpk : (((p : ℕ) ^ k : ℕ) : ℤ) = (p : ℤ) * (p : ℤ) ^ (k - 1) := by
    push_cast
    rw [← pow_succ']
    congr 1
    omega
  rw [hpk] at hz
  have hpne : (p : ℤ) ≠ 0 := by exact_mod_cast hp.ne'
  have hdvd : ((p : ℤ) ^ (k - 1)) ∣ i := (mul_dvd_mul_iff_left hpne).1 hz
  obtain ⟨t, ht⟩ := hdvd
  refine ⟨t, ?_⟩
  rw [ht, mul_comm, mul_zsmul, ← Int.natCast_pow, natCast_zsmul]

/-- ★`p^{k−1}a` の位数は `p`。 -/
theorem addOrderOf_high {A : Type u} [AddCommGroup A] (p k : ℕ) (hp : p.Prime) (hk : 1 ≤ k)
    (a : A) (hord : addOrderOf a = p ^ k) (hne : (p ^ (k - 1)) • a ≠ 0) :
    addOrderOf ((p ^ (k - 1)) • a) = p := by
  have hz : p • ((p ^ (k - 1)) • a) = 0 := by
    rw [smul_smul, show p * p ^ (k - 1) = p ^ k by rw [← pow_succ']; congr 1; omega, ← hord]
    exact addOrderOf_nsmul_eq_zero a
  have hdvd : addOrderOf ((p ^ (k - 1)) • a) ∣ p := addOrderOf_dvd_of_nsmul_eq_zero hz
  rcases (Nat.Prime.eq_one_or_self_of_dvd hp _ hdvd) with h1 | h1
  · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h1) hne
  · exact h1

/-- ★★**`⟨a⟩` の外の `c ∈ A[p]` がある**——`#A[p] = p² > p`。 -/
theorem exists_outside {A : Type u} [AddCommGroup A] [Finite A] (p k : ℕ) (hp : p.Prime)
    (hk : 1 ≤ k) (a : A) (hord : addOrderOf a = p ^ k) (hne : (p ^ (k - 1)) • a ≠ 0)
    (hcard1 : Nat.card (nsmulHom A p).ker = p ^ 2) :
    ∃ c : A, p • c = 0 ∧ c ∉ AddSubgroup.zmultiples a := by
  by_contra hcon
  push_neg at hcon
  have hle : (nsmulHom A p).ker ≤ AddSubgroup.zmultiples ((p ^ (k - 1)) • a) := by
    intro x hx
    have hx' : p • x = 0 := hx
    obtain ⟨i, hi⟩ := hcon x hx'
    have hi' : i • a = x := hi
    obtain ⟨t, ht⟩ := zsmul_mem_of_p_smul p k hp.pos hk a hord i (by rw [hi']; exact hx')
    refine ⟨t, ?_⟩
    show t • ((p ^ (k - 1)) • a) = x
    rw [ht, hi']
  have hcle : Nat.card (nsmulHom A p).ker
      ≤ Nat.card (AddSubgroup.zmultiples ((p ^ (k - 1)) • a)) :=
    Nat.card_le_card_of_injective (AddSubgroup.inclusion hle) (AddSubgroup.inclusion_injective hle)
  rw [hcard1, Nat.card_zmultiples, addOrderOf_high p k hp hk a hord hne] at hcle
  have hlt : p < p ^ 2 := by nlinarith [hp.two_le]
  omega

/-- ★★★★**独立性の核**——`j•b ∈ ⟨a⟩` なら `p^k ∣ j`。

★Bezout で `gcd(j, p^k)•b ∈ ⟨a⟩`、`gcd = p^s` で `s < k` なら `p^{k−1}b ∈ ⟨a⟩` に至る。 -/
theorem dvd_of_zsmul_mem {A : Type u} [AddCommGroup A] (p k : ℕ) (hp : p.Prime) (hk : 1 ≤ k)
    (a b : A) (hb : (p ^ k) • b = 0)
    (hcout : ((p ^ (k - 1)) • b) ∉ AddSubgroup.zmultiples a)
    (j : ℤ) (hmem : j • b ∈ AddSubgroup.zmultiples a) :
    ((p ^ k : ℕ) : ℤ) ∣ j := by
  by_contra hnd
  set d : ℕ := Int.gcd j ((p ^ k : ℕ) : ℤ) with hd
  have hdr : (d : ℤ) ∣ ((p ^ k : ℕ) : ℤ) := Int.gcd_dvd_right j _
  have hdvd : d ∣ p ^ k := Int.natCast_dvd_natCast.mp hdr
  have hdl : (d : ℤ) ∣ j := Int.gcd_dvd_left j _
  obtain ⟨s, hs, hds⟩ := (Nat.dvd_prime_pow hp).1 hdvd
  have hsk : s < k := by
    rcases Nat.lt_or_ge s k with h | h
    · exact h
    · exfalso
      apply hnd
      have hde : d = p ^ k := by rw [hds]; congr 1; omega
      rw [hde] at hdl
      exact_mod_cast hdl
  obtain ⟨u, v, huv⟩ : ∃ u v : ℤ, (d : ℤ) = j * u + ((p ^ k : ℕ) : ℤ) * v :=
    ⟨Int.gcdA j _, Int.gcdB j _, Int.gcd_eq_gcd_ab j _⟩
  have hpkb : ((p ^ k : ℕ) : ℤ) • b = 0 := by rw [natCast_zsmul]; exact hb
  have hdb : ((d : ℕ) : ℤ) • b ∈ AddSubgroup.zmultiples a := by
    rw [huv, add_zsmul, mul_comm j u, mul_zsmul, mul_comm ((p ^ k : ℕ) : ℤ) v, mul_zsmul,
      hpkb, zsmul_zero, add_zero]
    exact AddSubgroup.zsmul_mem _ hmem u
  apply hcout
  have hcomp : (p ^ (k - 1) : ℕ) • b = ((p ^ (k - 1 - s) : ℕ) : ℤ) • (((d : ℕ) : ℤ) • b) := by
    rw [hds, smul_smul, ← Nat.cast_mul, ← pow_add,
      show k - 1 - s + s = k - 1 by omega, natCast_zsmul]
  rw [hcomp]
  exact AddSubgroup.zsmul_mem _ hdb _

/-- ★★★★★★**素数冪の場合の独立な 2 元**。 -/
theorem exists_indep_prime_pow {A : Type u} [AddCommGroup A] [Finite A] (p k : ℕ) (hp : p.Prime)
    (hk : 1 ≤ k)
    (hexp : ∀ x : A, (p ^ k) • x = 0)
    (hcardA : Nat.card A = (p ^ k) ^ 2)
    (hcardk1 : Nat.card (nsmulHom A (p ^ (k - 1))).ker = (p ^ (k - 1)) ^ 2)
    (hcard1 : Nat.card (nsmulHom A p).ker = p ^ 2) :
    ∃ a b : A, (p ^ k) • a = 0 ∧ (p ^ k) • b = 0 ∧
      ∀ i j : ℤ, i • a + j • b = 0 → ((p ^ k : ℕ) : ℤ) ∣ i ∧ ((p ^ k : ℕ) : ℤ) ∣ j := by
  obtain ⟨a, ha⟩ := exists_high_order p k hk hcardA hcardk1 hp.one_lt
  have hord : addOrderOf a = p ^ k := addOrderOf_eq_prime_pow p k hp hk a (hexp a) ha
  obtain ⟨c, hc0, hcout⟩ := exists_outside p k hp hk a hord ha hcard1
  have hcr : c ∈ (nsmulHom A (p ^ (k - 1))).range := by
    rw [range_eq_ker p k hp.one_lt hk hexp hcardA hcardk1 hcard1]
    exact hc0
  obtain ⟨b, hb⟩ := hcr
  have hb' : (p ^ (k - 1)) • b = c := hb
  refine ⟨a, b, hexp a, hexp b, ?_⟩
  intro i j hij
  have hjb : j • b ∈ AddSubgroup.zmultiples a := by
    refine ⟨-i, ?_⟩
    show (-i) • a = j • b
    rw [neg_zsmul, neg_eq_iff_add_eq_zero]
    exact hij
  have hjd : ((p ^ k : ℕ) : ℤ) ∣ j :=
    dvd_of_zsmul_mem p k hp hk a b (hexp b) (by rw [hb']; exact hcout) j hjb
  refine ⟨?_, hjd⟩
  have hjb0 : j • b = 0 := by
    obtain ⟨t, ht⟩ := hjd
    rw [ht, mul_comm, mul_zsmul, natCast_zsmul, hexp b, zsmul_zero]
  have hia : i • a = 0 := by rw [hjb0, add_zero] at hij; exact hij
  rw [← hord]
  exact addOrderOf_dvd_iff_zsmul_eq_zero.2 hia

/-- ★★★★★★★**素数冪の場合の構造定理** `A ≃+ (ℤ/p^k)²`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem addEquiv_prime_pow {A : Type u} [AddCommGroup A] [Finite A] (p k : ℕ) (hp : p.Prime)
    (hk : 1 ≤ k)
    (hexp : ∀ x : A, (p ^ k) • x = 0)
    (hcardA : Nat.card A = (p ^ k) ^ 2)
    (hcardk1 : Nat.card (nsmulHom A (p ^ (k - 1))).ker = (p ^ (k - 1)) ^ 2)
    (hcard1 : Nat.card (nsmulHom A p).ker = p ^ 2) :
    Nonempty ((ZMod (p ^ k) × ZMod (p ^ k)) ≃+ A) := by
  obtain ⟨a, b, ha, hb, hind⟩ := exists_indep_prime_pow p k hp hk hexp hcardA hcardk1 hcard1
  have h1 : 1 ≤ p ^ k := Nat.one_le_pow _ _ hp.pos
  exact addEquiv_of_indep (p ^ k) h1 a b ha hb hind hcardA

/-! ## ★出典の紐付け(`.src`) -/

def addEquiv_prime_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(素数冪捩れの構造定理)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
