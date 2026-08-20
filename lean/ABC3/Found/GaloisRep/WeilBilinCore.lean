import ABC3.Found.GaloisRep.AddRelation

/-!
# Galois (G5) 第 183 ブロック —— **★★★★★★★双線型性の中核**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★3 つのブロックが合流する

| 材料 | 出どころ |
|---|---|
| `f₁ f₂ = u · h^n · f₃` | 第 182 `elem_relation_of_add` |
| `μ̃`(`μ` の関数体への延長)と `τ ∘ μ̃ = μ̃` | 第 181 |
| `(g₁g₂/(μ̃(h)·g₃))^n = 定数 ⟹ 定数` | 第 177 `const_of_pow_eq_const` |
| `g₁g₂ = ζ·z·g₃ ⟹ (τg₁/g₁)(τg₂/g₂) = τg₃/g₃` | 第 181 `aut_div_mul` |

★これらを繋ぐと **`e_n(P₁,Q)·e_n(P₂,Q) = e_n(P₁+P₂,Q)`** の中核が出る。

### ★★★★★機構

1. `μ̃` を `f₁f₂ = u·h^n·f₃` に当てて `μf₁·μf₂ = u·μ̃(h)^n·μf₃`
   (定数は `μ̃` で動かない——`hμF` と scalar tower)
2. `g_i^n = μf_i` を代入して `(g₁g₂/(μ̃(h)·g₃))^n = u`(定数)
3. 第 177 で `g₁g₂/(μ̃(h)·g₃) = ζ`(定数)
4. 第 181 の `aut_div_mul` に `z := μ̃(h)` を当てる——
   **`τ ∘ μ̃ = μ̃`** なので `z` は `τ` で不変

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `aut_div_mul_of_relation` | ★★★★★★★**双線型性の中核** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★★**双線型性の中核**——`f₁f₂ = u·h^n·f₃` から
`(τg₁/g₁)(τg₂/g₂) = τg₃/g₃`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 182(関係式)・第 181(`μ̃` と機構)・第 177(定数の `n` 乗根)の合流。 -/
theorem aut_div_mul_of_relation (h2 : IsUnit (2 : F)) {n : ℕ} (hn : 1 ≤ n)
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    (hcomp : (τ.toAlgHom.toRingHom).comp μ = μ)
    {f₁ f₂ f₃ : W.CoordinateRing} {hh : W.FunctionField} {u : F}
    (hhne : hh ≠ 0) (hu0 : u ≠ 0)
    (hrel : (algebraMap W.CoordinateRing W.FunctionField f₁)
        * (algebraMap W.CoordinateRing W.FunctionField f₂)
        = algebraMap F W.FunctionField u * hh ^ n
          * (algebraMap W.CoordinateRing W.FunctionField f₃))
    {g₁ g₂ g₃ : W.FunctionField} (hg₁ : g₁ ^ n = μ f₁) (hg₂ : g₂ ^ n = μ f₂)
    (hg₃ : g₃ ^ n = μ f₃) (hf₃ : μ f₃ ≠ 0) :
    (τ g₁ / g₁) * (τ g₂ / g₂) = τ g₃ / g₃ := by
  set z := muExt W hinj hh with hzdef
  have hzne : z ≠ 0 := by
    rw [hzdef]
    intro h0
    exact hhne ((muExt W hinj).injective (h0.trans (map_zero _).symm))
  have hg₃ne : g₃ ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hg₃
    exact hf₃ hg₃.symm
  have hconst : muExt W hinj (algebraMap F W.FunctionField u)
      = algebraMap F W.FunctionField u := by
    have h1 : (algebraMap F W.FunctionField u)
        = algebraMap W.CoordinateRing W.FunctionField (algebraMap F W.CoordinateRing u) :=
      IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField u
    rw [h1, muExt_algebraMap, hμF u]
    exact h1
  have hmap := congrArg (muExt W hinj) hrel
  rw [map_mul, map_mul, map_mul, map_pow, muExt_algebraMap, muExt_algebraMap,
    muExt_algebraMap, hconst, ← hzdef] at hmap
  have hzg : (g₁ * g₂ / (z * g₃)) ^ n = algebraMap F W.FunctionField u := by
    rw [div_pow, mul_pow, mul_pow, hg₁, hg₂, hg₃, hmap]
    field_simp
  obtain ⟨ζ, hζ0, hζ⟩ := const_of_pow_eq_const W h2 hn hu0 hzg
  have heq : g₁ * g₂ = algebraMap F W.FunctionField ζ * z * g₃ := by
    have hd := (div_eq_iff (mul_ne_zero hzne hg₃ne)).1 hζ
    rw [hd]; ring
  exact aut_div_mul W τ hζ0 hzne hg₃ne (aut_comp_muExt W τ hinj hcomp hh) heq

/-! ## ★出典の紐付け(`.src`) -/

def aut_div_mul_of_relation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性の中核)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
