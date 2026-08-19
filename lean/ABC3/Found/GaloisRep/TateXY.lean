import ABC3.Found.GaloisRep.AdicSeries

/-!
# Galois (G6) 第 105 ブロック —— **★★★★★Tate 級数 `X(u,q)`・`Y(u,q)` の項と尾**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★葉 (a)——両側和を片側 2 本に畳む

Tate 一意化の座標は**両側**和である:

    X(u,q) = ∑_{n∈ℤ} qⁿu/(1−qⁿu)² − 2s₁(q)

★ここで `f(t) = t/(1−t)²` と置くと、**`f(1/t) = f(t)`**(`tateXterm_inv`)なので

    ∑_{n∈ℤ} f(qⁿu) = f(u) + ∑_{m≥1} f(qᵐu) + ∑_{m≥1} f(qᵐu⁻¹)

——★★**片側 2 本に畳める**。★★★`qᵐu` も `qᵐu⁻¹` も `I^m` に入るので、
第 104 ブロックの `adicSum` がそのまま使える。

## ★★★★完備 adic 環では `1 − x` が単元

分母 `(1−qⁿu)` が可逆であることが要る。★これは

    x ∈ I,  IsAdicComplete I R  ⟹  IsUnit (1 − x)

であり、逆元は等比級数 `∑ xⁿ` である(`mul_adicSum_geom`)。
★★mathlib には無かった(2026-08-20 実測)ので自分で積んだ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mul_adicSum_geom` | ★★★`(1−x)·∑xⁿ = 1` |
| `isUnit_one_sub` | ★★★★**完備 adic 環では `1−x` は単元** |
| `tateXterm` / `tateYterm` | ★`f(t) = t/(1−t)²`、`g(t) = t²/(1−t)³` |
| `mul_tateXterm` / `mul_tateYterm` | ★★分母を払った形 |
| `tateXtail` / `tateYtail` | ★★★**片側の尾 `∑_{m≥1} f(qᵐu)`** |
| `tateXtail_rec` / `tateYtail_rec` | ★★★漸化式 `T(u) = f(qu) + T(qu)` |
| `tateXterm_inv` | ★★★★**`f(1/t) = f(t)`** |
| `tateYterm_add_inv` | ★★★★**`g(t) + g(1/t) = −f(t)`** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★完備 adic 環では `1 − x` が単元 -/

/-- ★★★**`(1−x)·∑ xⁿ = 1`**。 -/
theorem mul_adicSum_geom [IsAdicComplete I R] (x : R) (hx : x ∈ I) :
    (1 - x) * adicSum (fun n => x ^ n) (fun n => Ideal.pow_mem_pow hx n) = 1 := by
  set S := adicSum (fun n => x ^ n) (fun n => Ideal.pow_mem_pow hx n) with hSdef
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun n => ?_)
  rw [SModEq.sub_mem]
  set P := partialSum (fun n => x ^ n) n with hPdef
  have hspec := adicSum_spec (fun n => x ^ n) (fun n => Ideal.pow_mem_pow hx n) n
  rw [SModEq.sub_mem] at hspec
  have h1 : P - S ∈ I ^ n := by simpa using hspec
  have h1' : S - P ∈ I ^ n := by simpa using neg_mem h1
  have h2 : (1 - x) * P = 1 - x ^ n := by
    rw [hPdef, partialSum]
    exact mul_neg_geom_sum x n
  have h3 : (1 - x) * S - 1 = (1 - x) * (S - P) - x ^ n := by
    rw [mul_sub, h2]
    ring
  have h4 : (1 - x) * S - 1 ∈ I ^ n := by
    rw [h3]
    exact Submodule.sub_mem _ (Ideal.mul_mem_left _ _ h1') (Ideal.pow_mem_pow hx n)
  simpa using h4

/-- ★★★★**完備 adic 環では `1 − x`(`x ∈ I`)は単元である**。

★mathlib には無い(2026-08-20 実測)。★★分母 `(1−qⁿu)` の可逆性がこれで出る。 -/
theorem isUnit_one_sub [IsAdicComplete I R] {x : R} (hx : x ∈ I) : IsUnit (1 - x) :=
  IsUnit.of_mul_eq_one _ (mul_adicSum_geom x hx)

/-! ## ★級数の項 -/

/-- ★`f(t) = t/(1−t)²` —— `X` の項。 -/
noncomputable def tateXterm (t : R) : R := t * Ring.inverse (1 - t) ^ 2

/-- ★`g(t) = t²/(1−t)³` —— `Y` の項。 -/
noncomputable def tateYterm (t : R) : R := t ^ 2 * Ring.inverse (1 - t) ^ 3

theorem tateXterm_mem_pow {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateXterm t ∈ I ^ k :=
  Ideal.mul_mem_right _ _ ht

theorem tateYterm_mem_pow {k : ℕ} {t : R} (ht : t ∈ I ^ k) : tateYterm t ∈ I ^ k := by
  rw [tateYterm, sq]
  exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_right _ _ ht)

/-- ★★`(1−t)²·f(t) = t` —— `f` が本当に `t/(1−t)²` であること。 -/
theorem mul_tateXterm [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    (1 - t) ^ 2 * tateXterm t = t := by
  have hu : IsUnit (1 - t) := isUnit_one_sub ht
  rw [tateXterm,
    show (1 - t) ^ 2 * (t * Ring.inverse (1 - t) ^ 2)
      = t * ((1 - t) * Ring.inverse (1 - t)) ^ 2 by ring,
    Ring.mul_inverse_cancel _ hu, one_pow, mul_one]

/-- ★★`(1−t)³·g(t) = t²`。 -/
theorem mul_tateYterm [IsAdicComplete I R] {t : R} (ht : t ∈ I) :
    (1 - t) ^ 3 * tateYterm t = t ^ 2 := by
  have hu : IsUnit (1 - t) := isUnit_one_sub ht
  rw [tateYterm,
    show (1 - t) ^ 3 * (t ^ 2 * Ring.inverse (1 - t) ^ 3)
      = t ^ 2 * ((1 - t) * Ring.inverse (1 - t)) ^ 3 by ring,
    Ring.mul_inverse_cancel _ hu, one_pow, mul_one]

/-! ## ★★★片側の尾 -/

theorem tateXtail_aux {u q : R} (hq : q ∈ I) (n : ℕ) :
    tateXterm (q ^ (n + 1) * u) ∈ I ^ n :=
  Ideal.pow_le_pow_right (Nat.le_succ n)
    (tateXterm_mem_pow (Ideal.mul_mem_right u _ (Ideal.pow_mem_pow hq (n + 1))))

theorem tateYtail_aux {u q : R} (hq : q ∈ I) (n : ℕ) :
    tateYterm (q ^ (n + 1) * u) ∈ I ^ n :=
  Ideal.pow_le_pow_right (Nat.le_succ n)
    (tateYterm_mem_pow (Ideal.mul_mem_right u _ (Ideal.pow_mem_pow hq (n + 1))))

/-- ★★★**`∑_{m≥1} f(qᵐu)`** —— `X(u,q)` の片側の尾。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tateXtail [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : R :=
  adicSum (fun n => tateXterm (q ^ (n + 1) * u)) (tateXtail_aux hq)

/-- ★★★**`∑_{m≥1} g(qᵐu)`** —— `Y(u,q)` の片側の尾。 -/
noncomputable def tateYtail [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : R :=
  adicSum (fun n => tateYterm (q ^ (n + 1) * u)) (tateYtail_aux hq)

theorem tateXtail_mem [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : tateXtail u q hq ∈ I := by
  refine adicSum_mem _ _ ?_
  simpa using tateXterm_mem_pow (I := I) (k := 1)
    (by simpa using Ideal.mul_mem_right u I hq)

theorem tateYtail_mem [IsAdicComplete I R] (u q : R) (hq : q ∈ I) : tateYtail u q hq ∈ I := by
  refine adicSum_mem _ _ ?_
  simpa using tateYterm_mem_pow (I := I) (k := 1)
    (by simpa using Ideal.mul_mem_right u I hq)

/-- ★★★**尾の漸化式** `T(u) = f(qu) + T(qu)`。 -/
theorem tateXtail_rec [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail u q hq = tateXterm (q * u) + tateXtail (q * u) q hq := by
  rw [tateXtail, adicSum_shift]
  congr 1
  · norm_num
  · exact adicSum_congr _ _
      (fun n => by rw [show q ^ (n + 1 + 1) * u = q ^ (n + 1) * (q * u) by ring])

theorem tateYtail_rec [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail u q hq = tateYterm (q * u) + tateYtail (q * u) q hq := by
  rw [tateYtail, adicSum_shift]
  congr 1
  · norm_num
  · exact adicSum_congr _ _
      (fun n => by rw [show q ^ (n + 1 + 1) * u = q ^ (n + 1) * (q * u) by ring])

/-! ## ★★★★両側和が畳める理由 -/

section Field

variable {K : Type} [Field K]

theorem tateXterm_eq_div {t : K} (ht1 : t ≠ 1) : tateXterm t = t / (1 - t) ^ 2 := by
  have h1 : (1 : K) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  rw [tateXterm, Ring.inverse_eq_inv', div_eq_mul_inv, inv_pow]

theorem tateYterm_eq_div {t : K} (ht1 : t ≠ 1) : tateYterm t = t ^ 2 / (1 - t) ^ 3 := by
  have h1 : (1 : K) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  rw [tateYterm, Ring.inverse_eq_inv', div_eq_mul_inv, inv_pow]

theorem one_sub_inv_ne_zero {t : K} (ht1 : t ≠ 1) : (1 : K) - t⁻¹ ≠ 0 := by
  rw [sub_ne_zero]
  intro h
  exact ht1 (by rw [← inv_inv t, ← h, inv_one])

/-- ★★★★**`f(1/t) = f(t)`** —— 両側和が片側 2 本に畳める理由。 -/
theorem tateXterm_inv {t : K} (ht : t ≠ 0) (ht1 : t ≠ 1) :
    tateXterm t⁻¹ = tateXterm t := by
  have h1 : (1 : K) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  have h2 : (1 : K) - t⁻¹ ≠ 0 := one_sub_inv_ne_zero ht1
  have hne : t⁻¹ ≠ 1 := fun h => h2 (by rw [h]; ring)
  rw [tateXterm_eq_div hne, tateXterm_eq_div ht1,
    div_eq_div_iff (pow_ne_zero 2 h2) (pow_ne_zero 2 h1)]
  field_simp
  ring

theorem tateYterm_inv {t : K} (ht : t ≠ 0) (ht1 : t ≠ 1) :
    tateYterm t⁻¹ = - (t / (1 - t) ^ 3) := by
  have h1 : (1 : K) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  have h2 : (1 : K) - t⁻¹ ≠ 0 := one_sub_inv_ne_zero ht1
  have hne : t⁻¹ ≠ 1 := fun h => h2 (by rw [h]; ring)
  rw [tateYterm_eq_div hne, neg_div', div_eq_div_iff (pow_ne_zero 3 h2) (pow_ne_zero 3 h1)]
  field_simp
  ring

/-- ★★★★**`g(t) + g(1/t) = −f(t)`** —— 原典の `Y(u) + Y(1/u) = −X(u)`。 -/
theorem tateYterm_add_inv {t : K} (ht : t ≠ 0) (ht1 : t ≠ 1) :
    tateYterm t + tateYterm t⁻¹ = - tateXterm t := by
  have h1 : (1 : K) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  rw [tateYterm_inv ht ht1, tateYterm_eq_div ht1, tateXterm_eq_div ht1]
  field_simp
  ring

end Field

/-! ## ★出典の紐付け(`.src`) -/

def tateXtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化の級数 X(u,q)・Y(u,q)——葉 (a))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
