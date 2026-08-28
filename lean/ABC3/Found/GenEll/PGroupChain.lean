/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.TameRamification
import ABC3.Meta.Claim

/-!
# ★★★★★★段 EC4 —— p-群を次数 `p` の段に割る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★これは何か —— "elementary claim" の群論の側

`Proposition 1.7` の証明にある **elementary claim**

> 素数 `p` と正整数 `d` を固定すると、正整数 `n` が存在して、
> `[L:K] ≤ d` なる `Q_p` の有限拡大の有限 Galois 拡大 `L/K` について `𝔡_{L/K} ⊇ p^n·O_L`

の証明で、原文は

> since Gal(L/K) is a [necessarily solvable!] p-group of order ≤ d,
> it suffices to consider the case where [L:K] = p

と書く。★本ファイルはその**群論の側**を取る。

## ★★★機構 —— Sylow を入れ子に使う

★★★★**測定（2026-08-28）**: mathlib に `IsPGroup.exists_le_subgroup_card_eq` も
`IsPGroup.exists_normal_subgroup_card_eq` も**無い**。
あるのは `Sylow.exists_subgroup_card_pow_prime`
（`p^n ∣ |G|` なら位数 `p^n` の部分群が存在する）だけである。

★これを**部分群の中で使い直す**と、指数 `p` の段が 1 つ取れる:

    `|H| = p^{k+1}` ⟹ `∃ K ≤ H, |K| = p^k`

★★正規性は要らない——欲しいのは**中間体の塔**なので、部分群で足りる。

## ★残っている段（明示）

★★★Galois 対応で**中間体**に翻訳する段（`IsGalois.intermediateFieldEquivSubgroup`）と、
different の塔公式（`§9-DifferentKummer` の `pow_mem_differentIdeal_tower`）で継ぐ段。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★指数 `p` の部分群 -/

/-- ★★★★**位数 `p^{k+1}` の群は指数 `p` の部分群を持つ**。

★`Sylow.exists_subgroup_card_pow_prime` に `p^k ∣ p^{k+1}` を渡すだけである。 -/
theorem exists_subgroup_card_pow_of_pow_succ (p : ℕ) [Fact p.Prime] {G : Type*} [Group G]
    [Finite G] {k : ℕ} (hG : Nat.card G = p ^ (k + 1)) :
    ∃ H : Subgroup G, Nat.card H = p ^ k ∧ H.index = p := by
  obtain ⟨H, hH⟩ := Sylow.exists_subgroup_card_pow_prime (G := G) p
    (n := k) (by rw [hG]; exact pow_dvd_pow p (Nat.le_succ k))
  refine ⟨H, hH, ?_⟩
  have hmul := Subgroup.card_mul_index H
  rw [hH, hG] at hmul
  have hp : 0 < p := (Fact.out : p.Prime).pos
  have h2 : p ^ k * H.index = p ^ k * p := by rw [hmul, pow_succ]
  exact Nat.eq_of_mul_eq_mul_left (pow_pos hp k) h2

/-- ★★★★★★**部分群の中でも 1 段降りられる**。

★これが「`p`-群を次数 `p` の塔に割る」帰納の 1 歩である
——`Sylow` を**部分群の中で使い直す**のが要点。 -/
theorem exists_le_subgroup_card (p : ℕ) [Fact p.Prime] {G : Type*} [Group G] [Finite G]
    (H : Subgroup G) {k : ℕ} (hH : Nat.card H = p ^ (k + 1)) :
    ∃ K : Subgroup G, K ≤ H ∧ Nat.card K = p ^ k := by
  obtain ⟨K', hK', -⟩ := exists_subgroup_card_pow_of_pow_succ p (G := H) hH
  exact ⟨K'.map H.subtype, Subgroup.map_subtype_le K',
    (Subgroup.card_subtype H K').trans hK'⟩

/-! ## ★★★★★★★Galois 対応 —— 中間体へ翻訳する -/

open IntermediateField in
/-- ★★★★★★★**次数 `p^{k+1}` の Galois 拡大は次数 `p` の中間体を持つ**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★`exists_le_subgroup_card` を Galois 対応（`fixedField`）で中間体に翻訳したものである。
★★これで different の塔公式（`pow_mem_differentIdeal_tower`）に渡せる形になる。 -/
theorem exists_intermediateField_finrank (p : ℕ) [Fact p.Prime]
    (K L : Type) [Field K] [Field L] [Algebra K L] [FiniteDimensional K L] [IsGalois K L]
    {k : ℕ} (h : Module.finrank K L = p ^ (k + 1)) :
    ∃ M : IntermediateField K L, Module.finrank M L = p ^ k ∧ Module.finrank K M = p := by
  have hcard : Nat.card (L ≃ₐ[K] L) = p ^ (k + 1) := by
    rw [IsGalois.card_aut_eq_finrank, h]
  obtain ⟨H, hH⟩ : ∃ H : Subgroup (L ≃ₐ[K] L), Nat.card H = p ^ k :=
    ⟨_, (Sylow.exists_subgroup_card_pow_prime (G := (L ≃ₐ[K] L)) p
      (n := k) (by rw [hcard]; exact pow_dvd_pow p (Nat.le_succ k))).choose_spec⟩
  refine ⟨fixedField H, ?_, ?_⟩
  · rw [finrank_fixedField_eq_card H, hH]
  · have hmul := Module.finrank_mul_finrank K (fixedField H) L
    rw [finrank_fixedField_eq_card H, hH, h, pow_succ, mul_comm (p ^ k) p] at hmul
    have hp : 0 < p := (Fact.out : p.Prime).pos
    exact Nat.eq_of_mul_eq_mul_right (pow_pos hp k) hmul

/-! ## ★★段数は `d` だけで決まる -/

/-- ★★**次数が `d` 以下なら段数は `log_p d` 以下**。

★これが「`n` が `p` と `d` だけで決まる」ことの中身である
——原文の `Fix a prime number p and a positive integer d. Then there exists a positive
integer n` の `n` の一様性はここから来る。 -/
theorem exponent_le_log (p d k : ℕ) (hp : 2 ≤ p) (hd : d ≠ 0) (h : p ^ k ≤ d) :
    k ≤ Nat.log p d :=
  (Nat.le_log_iff_pow_le hp hd).mpr h

/-! ## ★★★★★正規な位数 `p` の部分群 -/

/-- ★**中心に含まれる部分群は正規である**。 -/
theorem normal_of_le_center {G : Type*} [Group G] (H : Subgroup G)
    (h : H ≤ Subgroup.center G) : H.Normal := by
  constructor
  intro n hn g
  have hc : ∀ x : G, x * n = n * x := (Subgroup.mem_center_iff).1 (h hn)
  rw [hc g, mul_assoc, mul_inv_cancel, mul_one]
  exact hn

/-- ★★**自明でない有限 `p`-群は位数 `p` の元を持つ**。 -/
theorem exists_orderOf_eq_prime (p : ℕ) [Fact p.Prime] {G : Type} [Group G]
    [Fintype G] [Nontrivial G] (hp : IsPGroup p G) : ∃ x : G, orderOf x = p := by
  obtain ⟨n, hn⟩ := (IsPGroup.iff_card).1 hp
  rw [Nat.card_eq_fintype_card] at hn
  have h1 : 1 < Fintype.card G := Fintype.one_lt_card
  have hn0 : n ≠ 0 := by rintro rfl; simp at hn; omega
  refine exists_prime_orderOf_dvd_card p ?_
  rw [hn]
  exact dvd_pow_self p hn0

/-- ★★★★★**自明でない有限 `p`-群は位数 `p` の**正規**部分群を持つ**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★中心が自明でない（`IsPGroup.center_nontrivial`）ので、
中心の中に位数 `p` の元を取り、その生成する部分群を取ればよい。
★★**正規性が要るのは、塔の各段を Galois にしたいから**である
——`§9-894` の `exists_intermediateField_finrank` は部分群だけでよかったが、
different の塔を帰納で回すには中間体も Galois でなければならない。 -/
theorem exists_normal_subgroup_card_prime (p : ℕ) [Fact p.Prime] {G : Type} [Group G]
    [Fintype G] [Nontrivial G] (hp : IsPGroup p G) :
    ∃ H : Subgroup G, H.Normal ∧ Nat.card H = p := by
  haveI : Fintype (Subgroup.center G) := Fintype.ofFinite _
  haveI : Nontrivial (Subgroup.center G) := hp.center_nontrivial
  obtain ⟨z, hz⟩ := exists_orderOf_eq_prime p (hp.to_subgroup (Subgroup.center G))
  refine ⟨(Subgroup.zpowers ((Subgroup.center G).subtype z)), ?_, ?_⟩
  · refine normal_of_le_center _ ?_
    rw [Subgroup.zpowers_le]
    exact z.2
  · rw [Nat.card_zpowers,
      orderOf_injective (Subgroup.center G).subtype (Subgroup.subtype_injective _) z]
    exact hz

def exists_normal_subgroup_card_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——p-群は位数 p の正規部分群を持つ)",
    sectionId := "genell-prop-1-7" }

/-! ## ★★★★★★★★★★Galois な中間体の存在 -/

open IntermediateField in
/-- ★★★★★★★★★★**次数 `p^{k+1}` の Galois 拡大は、`[L:M] = p` かつ `M/K` も Galois な中間体を持つ**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★位数 `p` の**正規**部分群 `N`（`exists_normal_subgroup_card_prime`）を取り、
`M ≔ L^N` とすると `[L:M] = |N| = p`、`[K:M] = p^k`で、
`N` が正規なので `M/K` も Galois である。
★★**これで different の塔を帰納で回せる形になった**
——`EC2`（次数 `p` の段）を `L/M` に、帰納法の仮定を `M/K` に当て、
`EC3`（塔）で継ぐ。 -/
theorem exists_intermediateField_galois (p : ℕ) [Fact p.Prime]
    (K L : Type) [Field K] [Field L] [Algebra K L] [FiniteDimensional K L] [IsGalois K L]
    {k : ℕ} (h : Module.finrank K L = p ^ (k + 1)) :
    ∃ M : IntermediateField K L, Module.finrank M L = p ∧ Module.finrank K M = p ^ k
      ∧ IsGalois K M := by
  have hcard : Nat.card (L ≃ₐ[K] L) = p ^ (k + 1) := by
    rw [IsGalois.card_aut_eq_finrank, h]
  haveI : Fintype (L ≃ₐ[K] L) := Fintype.ofFinite _
  haveI : Nontrivial (L ≃ₐ[K] L) := by
    refine Finite.one_lt_card_iff_nontrivial.mp ?_
    rw [hcard]
    exact Nat.one_lt_pow (Nat.succ_ne_zero k) (Fact.out : p.Prime).one_lt
  obtain ⟨N, hN, hcardN⟩ := exists_normal_subgroup_card_prime p (IsPGroup.of_card hcard)
  haveI := hN
  refine ⟨fixedField N, ?_, ?_, IsGalois.of_fixedField_normal_subgroup N⟩
  · rw [finrank_fixedField_eq_card N, hcardN]
  · have hmul := Module.finrank_mul_finrank K (fixedField N) L
    rw [finrank_fixedField_eq_card N, hcardN, h, pow_succ] at hmul
    have hp : 0 < p := (Fact.out : p.Prime).pos
    exact Nat.eq_of_mul_eq_mul_right hp hmul

def exists_intermediateField_galois.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——Galois な中間体の存在)",
    sectionId := "genell-prop-1-7" }

/-! ## ★★★★★★★★指数の一様化 -/

/-- ★**指数を上げても入ったまま** —— `x^m ∈ I` かつ `m ≤ n` なら `x^n ∈ I`。 -/
theorem pow_mem_of_le {R : Type*} [CommRing R] (I : Ideal R) (x : R) {m n : ℕ}
    (h : m ≤ n) (hm : x ^ m ∈ I) : x ^ n ∈ I := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le h
  rw [pow_add]
  exact Ideal.mul_mem_right _ _ hm

/-- ★★★★★★★★**`n` は `p` と `d` だけで決まる**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★塔の段数 `s` は `p^s ≤ d` を満たすので `s ≤ log_p d`であり（`exponent_le_log`）、
各段で `p^p` が出るので塔全体で `p^{s·p}` となる。
★★指数は上げても入ったままなので（`pow_mem_of_le`）、
**`n ≔ p·log_p d` という `L` にも `K` にも依らない値**で揃えられる。 -/
theorem pow_mem_uniform {R : Type*} [CommRing R] (I : Ideal R) (p d s N : ℕ)
    (hp : 2 ≤ p) (hd : d ≠ 0) (hs : p ^ s ≤ d) (hN : Nat.log p d * p ≤ N)
    (h : (p : R) ^ (s * p) ∈ I) : (p : R) ^ N ∈ I := by
  have hle : s ≤ Nat.log p d := (Nat.le_log_iff_pow_le hp hd).mpr hs
  exact pow_mem_of_le I _ (le_trans (Nat.mul_le_mul_right p hle) hN) h

def pow_mem_uniform.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——n は p と d だけで決まる)",
    sectionId := "genell-prop-1-7" }

/-! ## ★出典の紐付け(`.src`) -/

def exists_intermediateField_finrank.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——次数 p^{k+1} の Galois 拡大は次数 p の中間体を持つ)",
    sectionId := "genell-prop-1-7" }

def exponent_le_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——段数は log_p d 以下)",
    sectionId := "genell-prop-1-7" }

def exists_subgroup_card_pow_of_pow_succ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——位数 p^{k+1} の群は指数 p の部分群を持つ)",
    sectionId := "genell-prop-1-7" }

def exists_le_subgroup_card.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——部分群の中でも 1 段降りられる)",
    sectionId := "genell-prop-1-7" }

def exists_le_subgroup_card.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Sylow.exists_subgroup_card_pow_prime(p^n ∣ |G| なら位数 p^n の部分群)"
      (.inMathlib "Sylow.exists_subgroup_card_pow_prime") 1,
    .implicitStep
      ("★★★★測定(2026-08-28): mathlib に IsPGroup.exists_le_subgroup_card_eq も " ++
       "IsPGroup.exists_normal_subgroup_card_eq も**無い**。" ++
       "あるのは Sylow.exists_subgroup_card_pow_prime だけで、" ++
       "これを**部分群の中で使い直す**と指数 p の段が 1 つ取れる") 2,
    .implicitStep
      ("★★正規性は要らない——欲しいのは**中間体の塔**なので部分群で足りる") 1,
    .implicitStep
      ("★★★残るのは Galois 対応で**中間体**に翻訳する段" ++
       "(IsGalois.intermediateFieldEquivSubgroup)と、" ++
       "different の塔公式(§9-DifferentKummer の pow_mem_differentIdeal_tower)で継ぐ段") 4 ]

end ABC3.Found.GenEll
