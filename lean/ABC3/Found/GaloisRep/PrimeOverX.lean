import ABC3.Found.GaloisRep.UnramifiedStructure
import Mathlib.RingTheory.Ideal.GoingUp

/-!
# Galois (G5) 第 134 ブロック —— **★★★★★★極大イデアルの `x` 座標**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★零点定理は要らなかった

§9-442 では「極大イデアルが `XYIdeal` の形であること」に
**古典的零点定理**(mathlib に直接の形では無い)が要ると見積もっていた。
★**それは過大であった**——次の 3 つで足りる:

| 段 | 出所 |
|---|---|
| `P ≠ ⊥` ⟹ `P ∩ F[X] ≠ ⊥` | **mathlib**(`Ideal.IsIntegral.comap_ne_bot`) |
| `F[X]` の 0 でない素イデアルは `(X − c)` | `F[X]` が PID + `F` が代数閉 |
| `algebraMap(X − c) = x − c` | 第 116 ブロック |

★★整拡大(第 116)があれば、**零点定理を経由せずに** `x` 座標が取れる。

## ★★★★これが局所構造(第 132・133)の入口である

    P ≠ ⊥ 素   ⟹   ∃ c,  x − c ∈ P

★以後、`z² = Ψ₂Sq(x) ≡ Ψ₂Sq(c) (mod P)` から `z ≡ ±s` が出て、
第 132(`Ψ₂Sq(c) = 0`)と第 133(`Ψ₂Sq(c) ≠ 0`)の 2 つの場合に分かれる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_linear_mem_of_prime` | ★★★代数閉体上、`F[X]` の 0 でない素イデアルは `X − c` を含む |
| `algebraMap_X_sub_C` | ★`algebraMap (X − C c) = x − c` |
| `exists_genX_sub_mem` | ★★★★★★**`P ≠ ⊥` 素 ⟹ `∃ c, x − c ∈ P`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★**代数閉体上、`F[X]` の 0 でない素イデアルは `X − c` を含む**。

★`F[X]` は PID なので生成元 `q` が取れ、`q` は素元。
★★代数閉体上では既約多項式は 1 次(`IsAlgClosed.degree_eq_one_of_irreducible`)なので
`q = a(X − c)` であり、`a` は単元だから `X − c` も同じイデアルを生成する。 -/
theorem exists_linear_mem_of_prime (F : Type) [Field F] [IsAlgClosed F]
    (p : Ideal (Polynomial F)) (hp : p ≠ ⊥) [hpp : p.IsPrime] : ∃ c : F, X - C c ∈ p := by
  obtain ⟨q, hq⟩ := (IsPrincipalIdealRing.principal p).principal
  rw [show (Submodule.span (Polynomial F) {q}) = Ideal.span {q} from rfl] at hq
  have hq0 : q ≠ 0 := by
    rintro rfl
    exact hp (by rw [hq, Ideal.span_singleton_eq_bot.2 rfl])
  have hprime : Prime q := (Ideal.span_singleton_prime hq0).1 (hq ▸ hpp)
  have hdeg : q.degree = 1 := IsAlgClosed.degree_eq_one_of_irreducible F hprime.irreducible
  obtain ⟨a, b, ha, hq1⟩ : ∃ a b : F, a ≠ 0 ∧ q = C a * X + C b :=
    ⟨q.leadingCoeff, q.coeff 0, leadingCoeff_ne_zero.2 hq0, eq_X_add_C_of_degree_eq_one hdeg⟩
  refine ⟨-(b / a), ?_⟩
  rw [hq, Ideal.mem_span_singleton]
  refine ⟨C a⁻¹, ?_⟩
  rw [map_neg, sub_neg_eq_add, hq1, add_mul, mul_right_comm, ← map_mul,
    mul_inv_cancel₀ ha, map_one, one_mul, ← map_mul, div_eq_mul_inv]

/-- ★`algebraMap (X − C c) = x − c`——第 116 ブロックの `algebraMap` 記述から。 -/
theorem algebraMap_X_sub_C (W : WeierstrassCurve.Affine F) (c : F) :
    algebraMap (Polynomial F) W.CoordinateRing (X - C c)
      = genX W - algebraMap F W.CoordinateRing c := by
  rw [algebraMap_polynomial_coordinateRing, ← eval₂_genX, eval₂_sub, eval₂_X, eval₂_C]

/-- ★★★★★★**0 でない素イデアルは `x − c` を含む**(代数閉体上)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`F[W]` は `F[X]` 上整(第 116)なので `P ∩ F[X] ≠ ⊥`
(mathlib の `Ideal.IsIntegral.comap_ne_bot`)。
★★そこへ `F[X]` 側の補題を当てる。**古典的零点定理は要らない。** -/
theorem exists_genX_sub_mem [IsAlgClosed F] (W : WeierstrassCurve.Affine F)
    (P : Ideal W.CoordinateRing) (hP : P ≠ ⊥) [hPp : P.IsPrime] :
    ∃ c : F, genX W - algebraMap F W.CoordinateRing c ∈ P := by
  have hcomap : (P.comap (algebraMap (Polynomial F) W.CoordinateRing)) ≠ ⊥ :=
    Ideal.IsIntegral.comap_ne_bot (Polynomial F) hP
  haveI : (P.comap (algebraMap (Polynomial F) W.CoordinateRing)).IsPrime :=
    Ideal.comap_isPrime _ _
  obtain ⟨c, hc⟩ := exists_linear_mem_of_prime F _ hcomap
  exact ⟨c, by rw [← algebraMap_X_sub_C]; exact hc⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_genX_sub_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——0 でない素イデアルが x − c を含むこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
