/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Lemma65ii

/-!
# [FrdI] Theorem 6.4, (iii) —— `deg(Ψ^rlf) ∈ ℚ>0`(`Found`、`sorry` 無し)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> Next, we observe that assertion (iii) follows formally from assertion (ii), by

## ★★原文の論証(p.116)

原文はこう書く。

> If deg(Ψrlf) ∉ Q>0, then one verifies immediately that there exist three
> nonarchimedean valuations w1, w3, w5 ∈ V(L1) lying over primes p1, p3, p5,
> … such that p1, …, p6 are distinct. But this implies that
> log(p1)/log(p2), log(p3)/log(p4), log(p5)/log(p6) ∈ (deg(Ψrlf))⁻¹ · Q>0

★すなわち **`Lemma 6.5, (ii)`(six exponentials)への当て込み**である。

## ★★★★測定の訂正(2026-08-25)—— スケルトンの主張が原文とずれていた

`Skeleton/FrdI/Thm64Deg.lean` の旧 `deg_rat_of_six_exp` は仮定を

```
∀ p 素数, ∃ q, c · log p = log q
```

と置いていたが、これは原文の形ではない。原文が使うのは

```
∀ p 素数, ∃ p' 素数, ∃ λ ∈ ℚ>0, log p / log p' = c⁻¹ · λ
```

——**素数から素数への対応**(`Ψ` が誘導する `V(L₁) ≃ V(L₂)`)であり、
その対応が `c⁻¹ · ℚ>0` の比を与える、という形である。★本ファイルは
**原文の形**で書き直したうえで閉じた。

## ★★★★★「one verifies immediately」の中身は 3 段だった

| 段 | 中身 |
|---|---|
| (a) | `c ∉ ℚ` なら `f p ≠ p`(`f p = p` なら `1 = c⁻¹·λ` で `c = λ ∈ ℚ`) |
| (b) | `f` は**単射**(`f p = f q` なら `log p / log q = λ_p/λ_q ∈ ℚ` で `Lemma 6.5, (i)`) |
| (c) | 単射なので、除外すべき素数は各段で有限個。素数は無限個なので 3 組取れる |

★★(b) が要点である —— **単射性は仮定ではなく、`c` の無理性から出る**。
これが無いと「6 つが相異なる」を選べない。★原文の「immediately」の中身はここ。
-/

namespace ABC3.Found.FrdI

open ABC3.Meta

/-! ## ★1. 素数の対数の比が有理数なら同じ素数 -/

/-- ★★**`log p / log q ∈ ℚ` なら `p = q`** —— `Lemma 6.5, (i)` の言い換え。 -/
theorem prime_eq_of_log_div_rat {p q : Nat.Primes} {r : ℚ}
    (h : Real.log ((p : ℕ) : ℝ) / Real.log ((q : ℕ) : ℝ) = (r : ℝ)) : p = q := by
  by_contra hne
  have hq : Real.log ((q : ℕ) : ℝ) ≠ 0 :=
    (Real.log_pos (by exact_mod_cast q.2.one_lt)).ne'
  have hmul : Real.log ((p : ℕ) : ℝ) = (r : ℝ) * Real.log ((q : ℕ) : ℝ) := by
    field_simp at h; linarith [h]
  exact log_prime_ne_rat_mul hne r hmul

/-! ## ★2. 3 組の比から矛盾 —— `Lemma 6.5, (ii)` への当て込み -/

/-- ★★★★**3 組の比が `c⁻¹·ℚ>0` に入れば矛盾**(6 素数が相異なるとき)。

★原文の
`log(p1)/log(p2), log(p3)/log(p4), log(p5)/log(p6) ∈ (deg(Ψrlf))⁻¹ · Q>0`
がそのまま `Lemma 6.5, (ii)` の仮定になる(`λ₁ = l1/l3`、`λ₂ = l1/l5`)。 -/
theorem not_six_primes_ratio (c : ℝ) (hc : 0 < c)
    (P : Fin 6 → Nat.Primes) (hP : Function.Injective P)
    (l1 l3 l5 : ℚ) (hl1 : 0 < l1) (hl3 : 0 < l3) (hl5 : 0 < l5)
    (h1 : Real.log ((P 0 : ℕ) : ℝ) / Real.log ((P 1 : ℕ) : ℝ) = c⁻¹ * (l1 : ℝ))
    (h3 : Real.log ((P 2 : ℕ) : ℝ) / Real.log ((P 3 : ℕ) : ℝ) = c⁻¹ * (l3 : ℝ))
    (h5 : Real.log ((P 4 : ℕ) : ℝ) / Real.log ((P 5 : ℕ) : ℝ) = c⁻¹ * (l5 : ℝ)) :
    False := by
  refine lemma_6_5_ii P hP (l1 / l3) (l1 / l5) (by positivity) (by positivity) ?_ ?_
  · rw [h1, h3]; push_cast; field_simp
  · rw [h1, h5]; push_cast; field_simp

/-! ## ★3. 選び方 —— 単射性から「除外は有限個」 -/

/-- ★★**単射なら、除外すべき素数は有限個**なので外に取れる。

★`T` と `f ⁻¹ T` の合併が有限、素数は無限。 -/
theorem exists_prime_avoiding (f : Nat.Primes → Nat.Primes)
    (hf : Function.Injective f) (T : Finset Nat.Primes) :
    ∃ p : Nat.Primes, p ∉ T ∧ f p ∉ T := by
  have hfin : ((T : Set Nat.Primes) ∪ f ⁻¹' (T : Set Nat.Primes)).Finite :=
    (T.finite_toSet).union (Set.Finite.preimage hf.injOn T.finite_toSet)
  obtain ⟨p, hp⟩ := hfin.infinite_compl.nonempty
  simp only [Set.mem_compl_iff, Set.mem_union, Set.mem_preimage, Finset.mem_coe, not_or] at hp
  exact ⟨p, hp.1, hp.2⟩

/-- ★**6 つ組の単射性**(15 本の相異なりから)。 -/
theorem six_injective (f : Nat.Primes → Nat.Primes) (a b d : Nat.Primes)
    (e1 : f a ≠ a) (e2 : b ≠ a) (e3 : b ≠ f a) (e4 : f b ≠ a) (e5 : f b ≠ f a)
    (e6 : f b ≠ b) (e7 : d ≠ a) (e8 : d ≠ f a) (e9 : d ≠ b) (e10 : d ≠ f b)
    (e11 : f d ≠ a) (e12 : f d ≠ f a) (e13 : f d ≠ b) (e14 : f d ≠ f b) (e15 : f d ≠ d) :
    Function.Injective ![a, f a, b, f b, d, f d] := by
  intro i j hij
  fin_cases i <;> fin_cases j <;>
    first
      | rfl
      | exact absurd hij (by first | assumption | exact Ne.symm (by assumption))

/-! ## ★4. 本体 -/

/-- ★★★★★★**[FrdI] Theorem 6.4, (iii)** —— `deg(Ψ^rlf) ∈ ℚ>0`。

★仮定は原文の形そのもの —— `Ψ` が誘導する素数の対応 `f` と、
その比が `c⁻¹ · ℚ>0` に入ること(これが `(ii)` の `deg` の定義から出る段)。

★★中身は「`c` が無理数だと仮定すると `f` は不動点なしかつ**単射**になり、
素数が無限個あるので 6 つ相異なる素数が取れて `Lemma 6.5, (ii)` に矛盾する」。 -/
theorem deg_rat_of_pairing
    (c : ℝ) (hc : 0 < c)
    (f : Nat.Primes → Nat.Primes) (l : Nat.Primes → ℚ)
    (hl : ∀ p, 0 < l p)
    (hval : ∀ p : Nat.Primes,
      Real.log ((p : ℕ) : ℝ) / Real.log ((f p : ℕ) : ℝ) = c⁻¹ * (l p : ℝ)) :
    ∃ q : ℚ, 0 < q ∧ c = (q : ℝ) := by
  classical
  by_contra hirr
  have hcirr : ∀ r : ℚ, c ≠ (r : ℝ) := by
    intro r hr
    have hr0 : 0 < r := by rw [hr] at hc; exact_mod_cast hc
    exact hirr ⟨r, hr0, hr⟩
  have hlogpos : ∀ p : Nat.Primes, 0 < Real.log ((p : ℕ) : ℝ) := fun p =>
    Real.log_pos (by exact_mod_cast p.2.one_lt)
  -- ★(a) 不動点は無い
  have hne : ∀ p, f p ≠ p := by
    intro p hp
    have h := hval p
    rw [hp, div_self (hlogpos p).ne'] at h
    have hcl : c = (l p : ℝ) := by field_simp at h; linarith [h]
    exact hcirr (l p) hcl
  -- ★(b) `f` は単射 —— これが原文の「immediately」の中身
  have hinj : Function.Injective f := by
    intro p q hpq
    have hp := hval p
    have hq := hval q
    rw [hpq] at hp
    have hlfq : Real.log ((f q : ℕ) : ℝ) ≠ 0 := (hlogpos (f q)).ne'
    rw [div_eq_iff hlfq] at hp hq
    have hlqne : (l q : ℝ) ≠ 0 := by exact_mod_cast (hl q).ne'
    refine prime_eq_of_log_div_rat (r := l p / l q) ?_
    rw [hp, hq]
    push_cast
    field_simp
  -- ★(c) 3 組を選ぶ
  obtain ⟨a⟩ : Nonempty Nat.Primes := ⟨⟨2, Nat.prime_two⟩⟩
  obtain ⟨b, hb1, hb2⟩ := exists_prime_avoiding f hinj ({a, f a} : Finset Nat.Primes)
  obtain ⟨d, hd1, hd2⟩ := exists_prime_avoiding f hinj ({a, f a, b, f b} : Finset Nat.Primes)
  simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hb1 hb2 hd1 hd2
  exact not_six_primes_ratio c hc ![a, f a, b, f b, d, f d]
    (six_injective f a b d (hne a) hb1.1 hb1.2 hb2.1 hb2.2 (hne b)
      hd1.1 hd1.2.1 hd1.2.2.1 hd1.2.2.2 hd2.1 hd2.2.1 hd2.2.2.1 hd2.2.2.2 (hne d))
    (l a) (l b) (l d) (hl a) (hl b) (hl d) (hval a) (hval b) (hval d)

/-! ## ★5. `(iii)` の「In particular」—— `V(L₁) ≃ V(L₂)` は ℚ の素点を保つ

原文 (FrdI p.115):
> denote by Spec(L2)], A2 = (Ψpf)un-tr(A1), then the bijection
-/

/-- ★★★★**`(iii)` の「In particular」** —— `deg` が有理数だと分かった後は、
素数の対応 `f` は**恒等**になる(＝ `v₁` の上の `ℚ` の素点と `v₂` の上のそれが一致)。

★`log p / log (f p) = c⁻¹ · λ` で `c ∈ ℚ` なら右辺は有理数なので、
`Lemma 6.5, (i)` から `p = f p`。 -/
theorem pairing_eq_of_rat
    (c : ℝ) (q : ℚ) (hcq : c = (q : ℝ)) (_hq : q ≠ 0)
    (f : Nat.Primes → Nat.Primes) (l : Nat.Primes → ℚ)
    (hval : ∀ p : Nat.Primes,
      Real.log ((p : ℕ) : ℝ) / Real.log ((f p : ℕ) : ℝ) = c⁻¹ * (l p : ℝ)) :
    ∀ p, f p = p := by
  intro p
  refine (prime_eq_of_log_div_rat (p := p) (q := f p) (r := q⁻¹ * l p) ?_).symm
  rw [hval p, hcq]
  push_cast
  ring

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def deg_rat_of_pairing.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — deg(Ψ^rlf) ∈ ℚ>0",
    sectionId := "frdi-thm-6-4" }

def deg_rat_of_pairing.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Lemma 6.5, (ii)(six exponentials、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.FrdI.lemma_6_5_ii") 116,
    .citation "[ABC3]" "Lemma 6.5, (i)(素数の対数の ℚ 上一次独立性)"
      (.inProject "ABC3" "ABC3.Found.FrdI.log_primes_linearIndependent") 116,
    .derivation
      "c が無理数なら素数の対応 f は不動点なしかつ単射。素数は無限個なので 6 つ相異なる素数が取れる" 116,
    .implicitStep
      "★原文は「one verifies immediately that there exist three nonarchimedean valuations」で畳んでいる。中身は f の単射性である" 116 ]

def not_six_primes_ratio.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — 3 組の比が (deg)⁻¹·ℚ>0 に入れば矛盾",
    sectionId := "frdi-thm-6-4" }

def not_six_primes_ratio.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_6_5_ii"
      (.inProject "ABC3" "ABC3.Found.FrdI.lemma_6_5_ii") 116 ]

def prime_eq_of_log_div_rat.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — log p / log q ∈ ℚ なら p = q",
    sectionId := "frdi-thm-6-4" }

def prime_eq_of_log_div_rat.needs : List ProofObligation :=
  [ .citation "[ABC3]" "log_prime_ne_rat_mul(Lemma 6.5, (i) の系)"
      (.inProject "ABC3" "ABC3.Found.FrdI.log_prime_ne_rat_mul") 116 ]

def exists_prime_avoiding.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — 除外は有限個なので素数を外に取れる",
    sectionId := "frdi-thm-6-4" }

def exists_prime_avoiding.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Set.Finite.preimage(単射なら逆像は有限)"
      (.inMathlib "Set.Finite.preimage") 116,
    .derivation "素数は無限個なので、有限集合の補集合は空でない" 116 ]

def six_injective.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iii) — 6 つ組の単射性",
    sectionId := "frdi-thm-6-4" }

def six_injective.needs : List ProofObligation :=
  [ .derivation "15 本の相異なりから Fin 6 の 36 通りを潰す" 116 ]

def pairing_eq_of_rat.src : Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (iii) — V(L₁) ≃ V(L₂) は ℚ の素点の上下を保つ",
    sectionId := "frdi-thm-6-4" }

def pairing_eq_of_rat.needs : List ProofObligation :=
  [ .citation "[ABC3]" "prime_eq_of_log_div_rat"
      (.inProject "ABC3" "ABC3.Found.FrdI.prime_eq_of_log_div_rat") 115,
    .derivation "deg が有理数なら比も有理数なので、Lemma 6.5, (i) で p = f p" 115 ]

end ABC3.Found.FrdI
