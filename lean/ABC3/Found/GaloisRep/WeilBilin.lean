import ABC3.Found.GaloisRep.MulByNInj

/-!
# Galois (G5) 第 181 ブロック —— **★★★★★★双線型性の機構**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★双線型性は 1 つの等式に帰着する

    e_n(P₁+P₂, Q) = e_n(P₁,Q) · e_n(P₂,Q)

は、`g_i^n = μ f_{P_i}` として

    g₁ g₂ = ζ · z · g₃          (ζ は定数、z は τ_Q で不変)

の形の関係があれば出る(`aut_div_mul`):

    (τg₁/g₁)(τg₂/g₂) = (τg₁ τg₂)/(g₁g₂) = (ζ z τg₃)/(ζ z g₃) = τg₃/g₃

★`τ` は定数を固定するので `ζ` は消える。★★`z` が `τ` で不変なのが要点である。

### ★★★★★`z` は `μ̃(h)` の形になる

`XYIdeal(P₁)·XYIdeal(P₂)` と `XYIdeal(P₁+P₂)` は主イデアル `(h)` だけ違う
(mathlib の `mk_XYIdeal'_mul_mk_XYIdeal'`)。★`n` 乗して `μ` を当てると

    μf₁ · μf₂ = (定数) · μ̃(h)^n · μf₃

となり、`z := μ̃(h)` が取れる。★★そして **`τ ∘ μ̃ = μ̃`**(`aut_comp_muExt`)——
第 168 の `τ ∘ μ = μ` を関数体に延ばしたものである。

### ★★★★`μ` を関数体に延ばせるようになった

`μ̃ := IsFractionRing.lift hinj` は **`μ` の単射性**を要求する。
★第 180 でそれが定理になったので、延長が取れるようになった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `muExt` / `muExt_algebraMap` | ★★★★`μ` の関数体への延長 |
| `aut_comp_muExt` | ★★★★★★**`τ ∘ μ̃ = μ̃`** |
| `aut_div_mul` | ★★★★★★**双線型性の機構** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-! ## ★★★★`μ` の関数体への延長 -/

/-- ★★★★**`μ` を関数体に延ばす**(第 180 の単射性が要る)。 -/
noncomputable def muExt {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ) :
    W.FunctionField →+* W.FunctionField :=
  IsFractionRing.lift hinj

theorem muExt_algebraMap {μ : W.CoordinateRing →+* W.FunctionField}
    (hinj : Function.Injective μ) (r : W.CoordinateRing) :
    muExt W hinj (algebraMap W.CoordinateRing W.FunctionField r) = μ r :=
  IsFractionRing.lift_algebraMap hinj r

/-- ★★★★★★**`τ ∘ μ̃ = μ̃`**——第 168 を関数体に延ばしたもの。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem aut_comp_muExt (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hcomp : (τ.toAlgHom.toRingHom).comp μ = μ) (z : W.FunctionField) :
    τ (muExt W hinj z) = muExt W hinj z := by
  have hring : (τ.toAlgHom.toRingHom).comp (muExt W hinj) = muExt W hinj := by
    refine IsLocalization.ringHom_ext W.CoordinateRing⁰ ?_
    refine RingHom.ext (fun r => ?_)
    simp only [RingHom.comp_apply, muExt_algebraMap W hinj r]
    exact congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun)
      hcomp) r
  exact congrFun (congrArg (fun f => (f : W.FunctionField →+* W.FunctionField).toFun) hring) z

/-! ## ★★★★★★双線型性の機構 -/

/-- ★★★★★★**双線型性の機構**——`g₁ g₂ = ζ·z·g₃`(`ζ` は定数、`z` は `τ` で不変)なら
`(τg₁/g₁)(τg₂/g₂) = τg₃/g₃`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`τ` が定数を固定するので `ζ` は消える。★★`z` が `τ` で不変なのが要点。 -/
theorem aut_div_mul (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {g₁ g₂ g₃ z : W.FunctionField} {ζ : F} (hζ : ζ ≠ 0) (hz : z ≠ 0) (hg₃ : g₃ ≠ 0)
    (hτz : τ z = z)
    (heq : g₁ * g₂ = algebraMap F W.FunctionField ζ * z * g₃) :
    (τ g₁ / g₁) * (τ g₂ / g₂) = τ g₃ / g₃ := by
  have hζm : algebraMap F W.FunctionField ζ ≠ 0 := by simpa using hζ
  have hprod : g₁ * g₂ ≠ 0 := by
    rw [heq]
    exact mul_ne_zero (mul_ne_zero hζm hz) hg₃
  have h1 : g₁ ≠ 0 := fun h0 => hprod (by rw [h0, zero_mul])
  have h2 : g₂ ≠ 0 := fun h0 => hprod (by rw [h0, mul_zero])
  have hτ : τ g₁ * τ g₂ = algebraMap F W.FunctionField ζ * z * τ g₃ := by
    rw [← map_mul, heq, map_mul, map_mul, τ.commutes, hτz]
  rw [div_mul_div_comm, hτ, heq,
    mul_comm (algebraMap F W.FunctionField ζ * z) (τ g₃),
    mul_comm (algebraMap F W.FunctionField ζ * z) g₃,
    mul_div_mul_right _ _ (mul_ne_zero hζm hz)]

/-! ## ★出典の紐付け(`.src`) -/

def aut_comp_muExt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——μ の延長が平行移動で不変であること)",
    sectionId := "genell-thm-3-8" }

def aut_div_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性の機構)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
