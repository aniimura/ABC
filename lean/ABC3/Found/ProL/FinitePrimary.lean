/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.Torsion
import Mathlib.Data.Nat.Factorization.Basic
import Mathlib.RingTheory.Coprime.Lemmas
import Mathlib.FieldTheory.Finite.Basic

/-!
# 有限アーベル群の準素射影 —— チェーン `prol` の葉 `finite-primary`

★`ResearchPaper/frdi-decomposition.json` の `prol` チェーンの葉。
最終目標は [FrdI] `Definition 2.8, (ii)`(副有限アーベル群の pro-`l` 分解)である。

## ★何を作るか

有限アーベル群 `A` と素数 `p` について、**標準的な `p`-準素射影**

  `π_p : A →* A`,  `π_p x = x ^ (ordCompl[p] |A|) ^ φ(ordProj[p] |A|)`

を作り、次の 4 つを示す:

| 性質 | 宣言 |
|---|---|
| `π_p x` は `p`-準素成分に入る | `primProj_mem_primaryComponent` |
| `π_p` は `p`-準素成分の上で恒等 | `primProj_eq_self_of_mem` |
| `π_p` は `r ≠ p` の準素成分を潰す | `primProj_eq_one_of_mem_ne` |
| ★**`x = ∏_{p ∣ |A|} π_p x`** | `prod_primProj` |
| ★★**自然性 `f (π_p x) = π_p (f x)`** | `primProj_naturality` |

## ★★測定 —— mathlib にあったもの・無かったもの(2026-08-18)

★`CommGroup.primaryComponent`(`Mathlib/GroupTheory/Torsion.lean:376`)は**ある**。
`mem_primaryComponent_iff_orderOf` / `primaryComponent.disjoint` / `primaryComponent.isPGroup` も。
★★**無いのは分解そのもの**(`⨆_p primaryComponent = ⊤`)と**射影**である。
`Sylow.directProductOfNormal` は同型を与えるが射影の自然性は出さない。

★★★**設計の要**: 指数に **Euler の定理**を使ったこと。
`m := ordCompl[p] n` は `q := ordProj[p] n` と互いに素なので `m^{φ(q)} ≡ 1 [MOD q]`。
これで **Bezout を使わずに** `π_p` が書ける。
★積公式だけは中国剰余定理が要る(`sum_primExp_modEq_one`、ℤ の中で
`Finset.prod_dvd_of_coprime` を使う)。

★★**自然性は指数の等式では出ない**(`primExp` は群ごとに違う)。
**分解の一意性**(`p`-成分で恒等・他を潰す)から出る——それが `primProj_naturality` の証明である。
-/

namespace ABC3.Found.ProL

open CommGroup

/-! ## ★1. 指数 -/

/-- `p`-準素射影に使う指数。`m := ordCompl[p] n` の `φ(ordProj[p] n)` 乗。 -/
def primExp (n p : ℕ) : ℕ :=
  (n / p ^ n.factorization p) ^ (p ^ n.factorization p).totient

/-- `φ(q) ≥ 1` なので `m ∣ primExp`。 -/
theorem ordCompl_dvd_primExp {n p : ℕ} (hp : p.Prime) :
    (n / p ^ n.factorization p) ∣ primExp n p := by
  have hq : 0 < p ^ n.factorization p := pow_pos hp.pos _
  have hN : 0 < (p ^ n.factorization p).totient := Nat.totient_pos.mpr hq
  exact dvd_pow_self _ hN.ne'

/-- ★**Euler の定理** —— `m ^ φ(q) ≡ 1 [MOD q]`。 -/
theorem primExp_modEq_one {n p : ℕ} (hn : n ≠ 0) (hp : p.Prime) :
    primExp n p ≡ 1 [MOD p ^ n.factorization p] :=
  Nat.ModEq.pow_totient ((Nat.coprime_ordCompl hp hn).symm.pow_right _)

/-- ★`r ≠ p` が素数なら `ordProj[r] n ∣ ordCompl[p] n`。 -/
theorem ordProj_dvd_ordCompl {n p r : ℕ} (hp : p.Prime) (hr : r.Prime) (hne : r ≠ p) :
    r ^ n.factorization r ∣ n / p ^ n.factorization p := by
  have hcop : Nat.Coprime (r ^ n.factorization r) (p ^ n.factorization p) :=
    Nat.Coprime.pow _ _ ((Nat.coprime_primes hr hp).mpr hne)
  refine hcop.dvd_of_dvd_mul_left ?_
  rw [Nat.ordProj_mul_ordCompl_eq_self]
  exact Nat.ordProj_dvd n r

/-- 素因数分解の積。 -/
theorem prod_ordProj {n : ℕ} (hn : n ≠ 0) :
    ∏ p ∈ n.primeFactors, p ^ n.factorization p = n := by
  conv_rhs => rw [← Nat.prod_factorization_pow_eq_self hn]
  rw [Finsupp.prod, Nat.support_factorization]

/-- ★★**中国剰余定理の段** —— `∑_{p ∣ n} primExp n p ≡ 1 [MOD n]`。

★各素数冪 `q_r` を法にすると、`r` の項だけが `1` で他は `0` になる
(`q_r ∣ m_p ∣ primExp n p`)。素数冪は互いに素で積が `n` なので、ℤ の中で
`Finset.prod_dvd_of_coprime` が使える。 -/
theorem sum_primExp_modEq_one {n : ℕ} (hn : n ≠ 0) :
    (∑ p ∈ n.primeFactors, primExp n p) ≡ 1 [MOD n] := by
  have key : ((n : ℤ)) ∣ (∑ p ∈ n.primeFactors, (primExp n p : ℤ)) - 1 := by
    have hprod : ∏ r ∈ n.primeFactors, ((r ^ n.factorization r : ℕ) : ℤ) = (n : ℤ) := by
      rw [← Nat.cast_prod, prod_ordProj hn]
    rw [← hprod]
    refine Finset.prod_dvd_of_coprime ?_ ?_
    · intro r hr s hs hne
      have hrp : r.Prime := Nat.prime_of_mem_primeFactors hr
      have hsp : s.Prime := Nat.prime_of_mem_primeFactors hs
      show IsCoprime ((r ^ n.factorization r : ℕ) : ℤ) ((s ^ n.factorization s : ℕ) : ℤ)
      rw [Nat.isCoprime_iff_coprime]
      exact Nat.Coprime.pow _ _ ((Nat.coprime_primes hrp hsp).mpr hne)
    · intro r hr
      have hrp : r.Prime := Nat.prime_of_mem_primeFactors hr
      have hsplit : (∑ p ∈ n.primeFactors, (primExp n p : ℤ))
          = (primExp n r : ℤ) + ∑ p ∈ n.primeFactors.erase r, (primExp n p : ℤ) :=
        (Finset.add_sum_erase _ _ hr).symm
      have hre : (primExp n r : ℤ) + (∑ p ∈ n.primeFactors.erase r, (primExp n p : ℤ)) - 1
          = ((primExp n r : ℤ) - 1) + ∑ p ∈ n.primeFactors.erase r, (primExp n p : ℤ) := by ring
      rw [hsplit, hre]
      refine dvd_add (dvd_sub_comm.mp (primExp_modEq_one hn hrp).dvd) (Finset.dvd_sum ?_)
      intro p hp
      have hpp : p.Prime := Nat.prime_of_mem_primeFactors (Finset.mem_of_mem_erase hp)
      have hne : r ≠ p := fun h => (Finset.ne_of_mem_erase hp) h.symm
      have hd : r ^ n.factorization r ∣ primExp n p :=
        (ordProj_dvd_ordCompl hpp hrp hne).trans (ordCompl_dvd_primExp hpp)
      exact_mod_cast Int.natCast_dvd_natCast.mpr hd
  rw [Nat.modEq_iff_dvd]
  push_cast
  exact dvd_sub_comm.mp key

/-! ## ★2. 射影 -/

/-- ★★**`p`-準素射影** `π_p : A →* A`。 -/
noncomputable def primProj (A : Type*) [CommGroup A] [Finite A] (p : ℕ) : A →* A :=
  powMonoidHom (primExp (Nat.card A) p)

@[simp] theorem primProj_apply {A : Type*} [CommGroup A] [Finite A] (p : ℕ) (x : A) :
    primProj A p x = x ^ primExp (Nat.card A) p := rfl

section
variable {A : Type*} [CommGroup A] [Finite A]

theorem primProj_eq_self_of_dvd_ordProj {p : ℕ} (hp : p.Prime) {y : A}
    (hy : orderOf y ∣ p ^ (Nat.card A).factorization p) : primProj A p y = y := by
  have hn : Nat.card A ≠ 0 := Nat.card_pos.ne'
  have h := (primExp_modEq_one hn hp).of_dvd hy
  calc primProj A p y = y ^ primExp (Nat.card A) p := rfl
    _ = y ^ 1 := pow_eq_pow_iff_modEq.mpr h
    _ = y := pow_one y

theorem primProj_eq_one_of_dvd_ordCompl {p : ℕ} (hp : p.Prime) {y : A}
    (hy : orderOf y ∣ Nat.card A / p ^ (Nat.card A).factorization p) : primProj A p y = 1 :=
  orderOf_dvd_iff_pow_eq_one.mp (hy.trans (ordCompl_dvd_primExp hp))

theorem orderOf_dvd_ordProj_of_mem {p : ℕ} (hp : p.Prime) {y : A}
    (hy : y ∈ primaryComponent A p) : orderOf y ∣ p ^ (Nat.card A).factorization p := by
  haveI : Fact p.Prime := ⟨hp⟩
  obtain ⟨j, hj⟩ := mem_primaryComponent_iff_orderOf.mp hy
  have hdvd : p ^ j ∣ Nat.card A := hj ▸ orderOf_dvd_natCard y
  have hle : j ≤ (Nat.card A).factorization p :=
    (hp.pow_dvd_iff_le_factorization Nat.card_pos.ne').mp hdvd
  rw [hj]
  exact pow_dvd_pow p hle

theorem orderOf_dvd_ordCompl_of_mem {p r : ℕ} (hp : p.Prime) (hr : r.Prime) (hne : r ≠ p) {y : A}
    (hy : y ∈ primaryComponent A r) :
    orderOf y ∣ Nat.card A / p ^ (Nat.card A).factorization p := by
  haveI : Fact r.Prime := ⟨hr⟩
  obtain ⟨j, hj⟩ := mem_primaryComponent_iff_orderOf.mp hy
  have hcop : Nat.Coprime (orderOf y) (p ^ (Nat.card A).factorization p) := by
    rw [hj]; exact Nat.Coprime.pow _ _ ((Nat.coprime_primes hr hp).mpr hne)
  refine hcop.dvd_of_dvd_mul_left ?_
  rw [Nat.ordProj_mul_ordCompl_eq_self]
  exact orderOf_dvd_natCard y

/-- ★`π_p x` は `p`-準素成分に入る。 -/
theorem primProj_mem_primaryComponent {p : ℕ} (hp : p.Prime) (x : A) :
    primProj A p x ∈ primaryComponent A p := by
  refine ⟨(Nat.card A).factorization p, ?_⟩
  show (x ^ primExp (Nat.card A) p) ^ p ^ (Nat.card A).factorization p = 1
  rw [← pow_mul]
  refine orderOf_dvd_iff_pow_eq_one.mp ((orderOf_dvd_natCard x).trans ?_)
  conv_lhs => rw [← Nat.ordProj_mul_ordCompl_eq_self (Nat.card A) p]
  rw [mul_comm (primExp (Nat.card A) p)]
  exact Nat.mul_dvd_mul_left _ (ordCompl_dvd_primExp hp)

/-- ★`π_p` は `p`-準素成分の上で恒等。 -/
theorem primProj_eq_self_of_mem {p : ℕ} (hp : p.Prime) {y : A}
    (hy : y ∈ primaryComponent A p) : primProj A p y = y :=
  primProj_eq_self_of_dvd_ordProj hp (orderOf_dvd_ordProj_of_mem hp hy)

/-- ★**負の対照** —— `π_p` は恒等写像ではない: `r ≠ p` の準素成分を潰す。 -/
theorem primProj_eq_one_of_mem_ne {p r : ℕ} (hp : p.Prime) (hr : r.Prime) (hne : r ≠ p) {y : A}
    (hy : y ∈ primaryComponent A r) : primProj A p y = 1 :=
  primProj_eq_one_of_dvd_ordCompl hp (orderOf_dvd_ordCompl_of_mem hp hr hne hy)

/-- ★★★**準素分解** —— `x = ∏_{p ∣ |A|} π_p x`。 -/
theorem prod_primProj (x : A) : ∏ p ∈ (Nat.card A).primeFactors, primProj A p x = x := by
  have hn : Nat.card A ≠ 0 := Nat.card_pos.ne'
  have h1 : ∏ p ∈ (Nat.card A).primeFactors, primProj A p x
      = x ^ (∑ p ∈ (Nat.card A).primeFactors, primExp (Nat.card A) p) := by
    rw [← Finset.prod_pow_eq_pow_sum]; rfl
  rw [h1]
  calc x ^ (∑ p ∈ (Nat.card A).primeFactors, primExp (Nat.card A) p)
      = x ^ 1 :=
        pow_eq_pow_iff_modEq.mpr ((sum_primExp_modEq_one hn).of_dvd (orderOf_dvd_natCard x))
    _ = x := pow_one x

/-- `|A|` の素因数でない `p` については `π_p = 1`。 -/
theorem primProj_eq_one_of_not_mem {p : ℕ} (h : p ∉ (Nat.card A).primeFactors) (x : A) :
    primProj A p x = 1 := by
  have hf : (Nat.card A).factorization p = 0 := by
    by_contra hc
    exact h (Nat.support_factorization (Nat.card A) ▸ Finsupp.mem_support_iff.mpr hc)
  show x ^ primExp (Nat.card A) p = 1
  rw [show primExp (Nat.card A) p = Nat.card A by simp [primExp, hf]]
  exact pow_card_eq_one'

end

/-! ## ★3. 自然性 -/

section
variable {A B : Type*} [CommGroup A] [CommGroup B]

theorem map_mem_primaryComponent (f : A →* B) {q : ℕ} {y : A} (hy : y ∈ primaryComponent A q) :
    f y ∈ primaryComponent B q := by
  obtain ⟨k, hk⟩ := hy
  exact ⟨k, by rw [← map_pow, hk, map_one]⟩

end

section
variable {A B : Type*} [CommGroup A] [Finite A] [CommGroup B] [Finite B]

/-- ★★★**準素射影の自然性** —— `f (π_p x) = π_p (f x)`。

★★指数 `primExp` は群ごとに違うのに等式が成り立つ。理由は
`π_p` が**準素分解に沿った標準的な射影**だからで、証明は分解の一意性
(`p`-成分で恒等・他の成分を潰す)から出る。**指数の等式からは出ない。** -/
theorem primProj_naturality (f : A →* B) {p : ℕ} (hp : p.Prime) (x : A) :
    f (primProj A p x) = primProj B p (f x) := by
  by_cases hmem : p ∈ (Nat.card A).primeFactors
  · have hx : f x = ∏ q ∈ (Nat.card A).primeFactors, f (primProj A q x) := by
      rw [← map_prod, prod_primProj]
    rw [hx, map_prod, Finset.prod_eq_single p]
    · exact (primProj_eq_self_of_mem hp
        (map_mem_primaryComponent f (primProj_mem_primaryComponent hp x))).symm
    · intro q hq hne
      exact primProj_eq_one_of_mem_ne hp (Nat.prime_of_mem_primeFactors hq) hne
        (map_mem_primaryComponent f (primProj_mem_primaryComponent
          (Nat.prime_of_mem_primeFactors hq) x))
    · intro h; exact absurd hmem h
  · rw [primProj_eq_one_of_not_mem hmem x, map_one]
    have hord : orderOf (f x) ∣ Nat.card A := (orderOf_map_dvd f x).trans (orderOf_dvd_natCard x)
    have hnp : ¬ (p ∣ Nat.card A) := fun hd =>
      hmem (Nat.mem_primeFactors.mpr ⟨hp, hd, Nat.card_pos.ne'⟩)
    have hcop : Nat.Coprime (orderOf (f x)) (p ^ (Nat.card B).factorization p) :=
      Nat.Coprime.pow_right _ ((Nat.Prime.coprime_iff_not_dvd hp).mpr
        (fun hd => hnp (hd.trans hord))).symm
    refine (primProj_eq_one_of_dvd_ordCompl hp ?_).symm
    refine hcop.dvd_of_dvd_mul_left ?_
    rw [Nat.ordProj_mul_ordCompl_eq_self]
    exact orderOf_dvd_natCard (f x)

end

end ABC3.Found.ProL
