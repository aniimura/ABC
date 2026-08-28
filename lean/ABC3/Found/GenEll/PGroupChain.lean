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
