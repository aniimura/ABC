import ABC3.Found.GaloisRep.TateNormalForm

/-!
# Galois (G6) 第 308 ブロック —— **★★★★★★★★標準形に残る 1 径数で `a₄` を合わせる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> Tate 標準形の曲線は、**残る 1 径数**を動かして `a₄` を**任意の `𝔪` の値に合わせられる**
> (`exists_residual_match_a4`)

★★★これが G6 第 3 段の芯である。

## ★★★★★★★標準形はまだ 1 径数残っている

`a₁ = 1`・`a₂ = a₃ = 0` を保つ変数変換 `(u, r, s, t)` は

    u = 1 + 2s,    3r = s + s²,    r = −2t

を満たすものだけである。★標数を割らずに済ませるには **`t` で母数化**するのがよい:
`t ∈ 𝔪` を与えると `r := −2t`、`s` は **`s + s² = −6t`**(Artin–Schreier、第 307)で取れる。

## ★★★★★★★★★`u² = 1 − 24t`——`s` が消える

`u = 1 + 2s` と `s + s² = −6t` から

    u² = 1 + 4(s + s²) = 1 − 24t

★★★★**`s` が消えて `t` だけの式になる**。これで `a₄` の変換則が

    (1 − 24t)² · a₄' = a₄ − t + 12t²

という **`t` の 2 次式**になり、Hensel で持ち上げた `s` の連続性を論じずに済む。
★★★★★ここが本ブロックの要点である——**`s` を消せると分かった時点で 1 変数になった。**

## ★★★★★★縮小写像は 1 変数

`a₄' = A` を解くと

    t = a₄ + 12t² − (1 − 24t)²A =: F(t)

★`F : 𝔪 → 𝔪` で、差は `F(x) − F(y) = (x−y)·(12(x+y) + 24(2−24x−24y)A)` と
**因数分解できる**(`ring`)。★★★括弧が `𝔪` に入る(`x, y, A ∈ 𝔪`)ので縮小写像であり、
第 100 の `exists_fixedPoint_of_contraction` がそのまま効く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `residual_a1`・`residual_a2`・`residual_a3` | ★★★★標準形が保たれる |
| `residual_a4_eq` | ★★★★★★★`(1−24t)²a₄' = a₄ − t + 12t²` |
| `exists_residual_param` | ★★★★★★縮小写像で `t` を取る |
| `exists_residual_match_a4` | ★★★★★★★★**`a₄` を合わせられる** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

section Residual

variable {R : Type} [CommRing R]

/-- ★★★★残る 1 径数の変換でも `a₁ = 1`。 -/
theorem residual_a1 (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (u : Rˣ) (s t : R)
    (hu : ((u : R)) = 1 + 2 * s) :
    ((⟨u, -2 * t, s, t⟩ : VariableChange R) • E).a₁ = 1 := by
  rw [full_a1, h1, ← hu]
  rw [← Units.val_mul, inv_mul_cancel]
  rfl

/-- ★★★★残る 1 径数の変換でも `a₂ = 0`(`s + s² = −6t` のとき)。 -/
theorem residual_a2 (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (u : Rˣ) (s t : R)
    (hs : s + s ^ 2 = -6 * t) :
    ((⟨u, -2 * t, s, t⟩ : VariableChange R) • E).a₂ = 0 := by
  rw [full_a2, h1, h2]
  have hz : (0 : R) - s * 1 + 3 * (-2 * t) - s ^ 2 = 0 := by linear_combination -hs
  rw [hz, mul_zero]

/-- ★★★★残る 1 径数の変換でも `a₃ = 0`。 -/
theorem residual_a3 (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h3 : E.a₃ = 0) (u : Rˣ) (s t : R) :
    ((⟨u, -2 * t, s, t⟩ : VariableChange R) • E).a₃ = 0 := by
  rw [full_a3, h1, h3]
  have hz : (0 : R) + (-2 * t) * 1 + 2 * t = 0 := by ring
  rw [hz, mul_zero]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**`a₄` の変換則**——`s` が消えて `t` の 2 次式になる。 -/
theorem residual_a4_eq (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (h3 : E.a₃ = 0)
    (u : Rˣ) (s t : R) (hu : ((u : R)) = 1 + 2 * s) (hs : s + s ^ 2 = -6 * t) :
    (1 - 24 * t) ^ 2 * ((⟨u, -2 * t, s, t⟩ : VariableChange R) • E).a₄
      = E.a₄ - t + 12 * t ^ 2 := by
  rw [full_a4, h1, h2, h3]
  have hu2 : (1 - 24 * t) = ((u : R)) ^ 2 := by
    rw [hu]
    linear_combination -4 * hs
  have huv : ((u⁻¹ : Rˣ) : R) * ((u : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, inv_mul_cancel]
    rfl
  rw [hu2]
  have hexp : (((u : R)) ^ 2) ^ 2 * (((u⁻¹ : Rˣ) : R) ^ 4
      * (E.a₄ - s * 0 + 2 * (-2 * t) * 0 - (t + (-2 * t) * s) * 1 + 3 * (-2 * t) ^ 2 - 2 * s * t))
      = ((((u⁻¹ : Rˣ) : R)) * ((u : R))) ^ 4 * (E.a₄ - t + 12 * t ^ 2) := by
    ring
  rw [hexp, huv, one_pow, one_mul]

end Residual

section Local

variable {R : Type} [CommRing R] [IsLocalRing R]

/-- ★★★★残る 1 径数の変換でも `a₆ ∈ 𝔪`。 -/
theorem residual_a6_mem (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (h3 : E.a₃ = 0)
    (h4 : E.a₄ ∈ IsLocalRing.maximalIdeal R) (h6 : E.a₆ ∈ IsLocalRing.maximalIdeal R)
    (u : Rˣ) (s t : R) (ht : t ∈ IsLocalRing.maximalIdeal R) :
    ((⟨u, -2 * t, s, t⟩ : VariableChange R) • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  rw [full_a6, h1, h2, h3]
  refine Ideal.mul_mem_left _ _ ?_
  have hz : E.a₆ + (-2 * t) * E.a₄ + (-2 * t) ^ 2 * 0 + (-2 * t) ^ 3 - t * 0 - t ^ 2
      - (-2 * t) * t * 1
      = E.a₆ + (-2 * E.a₄) * t + (-8 * t * t) * t + t * t := by ring
  rw [hz]
  exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.add_mem _ h6 (Ideal.mul_mem_left _ _ ht))
    (Ideal.mul_mem_left _ _ ht)) (Ideal.mul_mem_left _ _ ht)

end Local

section Complete

variable {R : Type} [CommRing R] [IsLocalRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]

/-- ★★★★★★**縮小写像で残る径数 `t` を取る**。 -/
theorem exists_residual_param (a A : R) (ha : a ∈ IsLocalRing.maximalIdeal R)
    (hA : A ∈ IsLocalRing.maximalIdeal R) :
    ∃ t ∈ IsLocalRing.maximalIdeal R, a + 12 * t ^ 2 - (1 - 24 * t) ^ 2 * A = t := by
  refine exists_fixedPoint_of_contraction (I := IsLocalRing.maximalIdeal R)
    (fun t => a + 12 * t ^ 2 - (1 - 24 * t) ^ 2 * A) ?_ ?_ ?_
  · intro x hx
    exact Ideal.sub_mem _ (Ideal.add_mem _ ha (Ideal.mul_mem_left _ _
      (by rw [pow_two]; exact Ideal.mul_mem_left _ _ hx))) (Ideal.mul_mem_left _ _ hA)
  · have hz : a + 12 * (0 : R) ^ 2 - (1 - 24 * 0) ^ 2 * A = a - A := by ring
    show a + 12 * (0 : R) ^ 2 - (1 - 24 * 0) ^ 2 * A ∈ _
    rw [hz]
    exact Ideal.sub_mem _ ha hA
  · intro x hx y hy k hxy
    have hfac : (a + 12 * x ^ 2 - (1 - 24 * x) ^ 2 * A)
        - (a + 12 * y ^ 2 - (1 - 24 * y) ^ 2 * A)
        = (x - y) * (12 * (x + y) + 24 * (2 - 24 * x - 24 * y) * A) := by ring
    show (a + 12 * x ^ 2 - (1 - 24 * x) ^ 2 * A) - (a + 12 * y ^ 2 - (1 - 24 * y) ^ 2 * A) ∈ _
    rw [hfac, pow_succ]
    refine Ideal.mul_mem_mul hxy ?_
    exact Ideal.add_mem _ (Ideal.mul_mem_left _ _ (Ideal.add_mem _ hx hy))
      (Ideal.mul_mem_left _ _ hA)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**標準形のまま `a₄` を任意の `𝔪` の値に合わせられる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_residual_match_a4 (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0)
    (h3 : E.a₃ = 0) (h4 : E.a₄ ∈ IsLocalRing.maximalIdeal R)
    (h6 : E.a₆ ∈ IsLocalRing.maximalIdeal R) (A : R) (hA : A ∈ IsLocalRing.maximalIdeal R) :
    ∃ C : VariableChange R, (C • E).a₁ = 1 ∧ (C • E).a₂ = 0 ∧ (C • E).a₃ = 0
      ∧ (C • E).a₄ = A ∧ (C • E).a₆ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨t, ht, hfix⟩ := exists_residual_param E.a₄ A h4 hA
  obtain ⟨s, hs, hseq0⟩ := exists_artin_schreier (-6 * t) (Ideal.mul_mem_left _ _ ht)
  have hseq : s + s ^ 2 = -6 * t := by linear_combination hseq0
  have hu : IsUnit (1 + 2 * s) := isUnit_one_add_mem (Ideal.mul_mem_left _ _ hs)
  have huval : ((hu.unit : Rˣ) : R) = 1 + 2 * s := hu.unit_spec
  refine ⟨⟨hu.unit, -2 * t, s, t⟩, residual_a1 E h1 _ s t huval,
    residual_a2 E h1 h2 _ s t hseq, residual_a3 E h1 h3 _ s t, ?_,
    residual_a6_mem E h1 h2 h3 h4 h6 _ s t ht⟩
  have hkey := residual_a4_eq E h1 h2 h3 hu.unit s t huval hseq
  have hunit24 : IsUnit ((1 - 24 * t) ^ 2) := by
    refine IsUnit.pow 2 ?_
    have hz : (1 : R) - 24 * t = 1 + (-24 * t) := by ring
    rw [hz]
    exact isUnit_one_add_mem (Ideal.mul_mem_left _ _ ht)
  have heq : (1 - 24 * t) ^ 2 * ((⟨hu.unit, -2 * t, s, t⟩ : VariableChange R) • E).a₄
      = (1 - 24 * t) ^ 2 * A := by
    rw [hkey]
    linear_combination hfix
  exact hunit24.mul_left_cancel heq

end Complete

/-! ## ★出典の紐付け(`.src`) -/

def residual_a4_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(標準形に残る 1 径数の a₄ 変換則)",
    sectionId := "genell-def-3-3" }

def exists_residual_match_a4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(残る 1 径数で a₄ を合わせる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
