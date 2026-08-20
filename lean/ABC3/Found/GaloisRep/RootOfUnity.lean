import ABC3.Found.GaloisRep.D2Bridge

/-!
# Galois (G5) 第 176 ブロック —— **★★★★★`τ(g)/g` は定数である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★Weil 対の値が `F` に落ちること

D2 が閉じたので、`P ∈ E[n]` に対して `g_P ∈ F(W)` で `g_P^n = μ f_P` なるものが取れる。
★Weil 対は

    e_n(P, Q) := τ_Q(g_P) / g_P        (Q ∈ E[n]、τ_Q は平行移動)

で定める。★★これが **`F` の元**であることが定義の前提であり、本ブロックで示す。

### ★★★機構は 2 段

1. **`(τ_Q(g_P)/g_P)^n = 1`** —— 第 168 の `τ ∘ μ = μ` から
   `τ_Q(g_P)^n = τ_Q(g_P^n) = τ_Q(μ f_P) = μ f_P = g_P^n`(`pow_aut_div_eq_one`)。
2. **1 の `n` 乗根は定数** —— `h^n = 1` なら `h` は `X^n − 1` の根なので
   `F[W]` 上整。★第 137 の**整閉性**から `h ∈ F[W]`、そこで `h^n = 1` は単元、
   第 128 の**単元は定数**から `h ∈ F`(`const_of_pow_eq_one`)。

★★★**第 137(整閉性)と第 128(単元は定数)がここで効いた。**
どちらも D2 のために積んだものである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `const_of_pow_eq_one` | ★★★★★**1 の `n` 乗根は定数** |
| `pow_aut_div_eq_one` | ★★★★`(τ g / g)^n = 1` |
| `exists_const_aut_div` | ★★★★★**`τ(g)/g` は定数** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★**1 の `n` 乗根は定数である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`h^n = 1` なら `h` は `X^n − 1` の根なので `F[W]` 上整。第 137 の整閉性で
`h ∈ F[W]`、そこで単元なので第 128 から定数。 -/
theorem const_of_pow_eq_one (h2 : IsUnit (2 : F)) {h : W.FunctionField} {n : ℕ} (hn : 1 ≤ n)
    (hpow : h ^ n = 1) : ∃ c : F, c ≠ 0 ∧ h = algebraMap F W.FunctionField c := by
  haveI := isIntegrallyClosed_coordinateRing h2 W
  have hint : IsIntegral W.CoordinateRing h := by
    refine ⟨Polynomial.X ^ n - 1, ?_, ?_⟩
    · exact (Polynomial.monic_X_pow n).sub_of_left (by
        refine lt_of_le_of_lt (Polynomial.degree_one_le) ?_
        rw [Polynomial.degree_X_pow]
        exact_mod_cast (by omega : (0 : ℕ) < n))
    · simp [Polynomial.eval₂_sub, Polynomial.eval₂_pow, hpow]
  obtain ⟨u, hu⟩ := IsIntegrallyClosed.isIntegral_iff.1 hint
  have hun : u ^ n = 1 := by
    refine IsFractionRing.injective W.CoordinateRing W.FunctionField ?_
    rw [map_pow, hu, hpow, map_one]
  have hmul : u * u ^ (n - 1) = 1 := by
    rw [← pow_succ', show n - 1 + 1 = n by omega]
    exact hun
  have huu : IsUnit u := ⟨⟨u, u ^ (n - 1), hmul, by rw [mul_comm]; exact hmul⟩, rfl⟩
  obtain ⟨c, hc0, hcu⟩ := isUnit_coordinateRing huu
  refine ⟨c, hc0, ?_⟩
  rw [← hu, hcu, ← IsScalarTower.algebraMap_apply]

/-- ★★★★`τ` が `z = g^n` を固定すれば `(τ g / g)^n = 1`。 -/
theorem pow_aut_div_eq_one {τ : W.FunctionField ≃ₐ[F] W.FunctionField}
    {g z : W.FunctionField} {n : ℕ}
    (hg : g ^ n = z) (hz : z ≠ 0) (hτz : τ z = z) :
    (τ g / g) ^ n = 1 := by
  rw [div_pow, ← map_pow, hg, hτz, div_self hz]

/-- ★★★★★**`τ(g)/g` は定数である**——Weil 対の値が `F` に落ちる根拠。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_const_aut_div (h2 : IsUnit (2 : F))
    {τ : W.FunctionField ≃ₐ[F] W.FunctionField}
    {g z : W.FunctionField} {n : ℕ} (hn : 1 ≤ n)
    (hg : g ^ n = z) (hz : z ≠ 0) (hτz : τ z = z) :
    ∃ c : F, c ≠ 0 ∧ τ g / g = algebraMap F W.FunctionField c :=
  const_of_pow_eq_one W h2 hn (pow_aut_div_eq_one W hg hz hτz)

/-! ## ★出典の紐付け(`.src`) -/

def const_of_pow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——1 の n 乗根が定数であること)",
    sectionId := "genell-thm-3-8" }

def exists_const_aut_div.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——τ(g)/g が定数であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
