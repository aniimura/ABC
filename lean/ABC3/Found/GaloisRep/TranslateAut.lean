import ABC3.Found.GaloisRep.Transcendence

/-!
# Galois (G5) 第 117 ブロック —— **★★★★★★平行移動は関数体に作用する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★葉 1 が閉じた(`Q` が 2 等分点でないとき)

第 116 ブロックで単射性が **`translateX` の超越性**に帰着した。★それを証明する。

## ★★★★★証明——`−Q` で評価する

`u := x − x₀`、`v := y − y₀` と置くと、加法公式から

    (x − x₀)² · X(τ_Q) = A,    A := v² + a₁·v·u − (a₂ + x + x₀)·u²

である(`uSq_mul_translateX`)。★ここで **`−Q` での評価** `ev : F[W] →+* F` を取る:

    ev(u) = x₀ − x₀ = 0
    ev(A) = ev(v)² = (negY(x₀,y₀) − y₀)²  ≠ 0     ← `Q` が 2 等分点でなければ

★★`p(X(τ_Q)) = 0` を仮定して分母を払うと `F[W]` の中で

    ∑_{i≤n} p_i · A^i · u^{2(n−i)} = 0

となる。★★★`ev` を当てると `u` の項がすべて消えて `p_n · ev(A)^n = 0` だけが残り、
`p_n ≠ 0`・`ev(A) ≠ 0` に矛盾する。

★★★★**Dedekind 環の理論も Zariski の補題も要らなかった**——1 点での評価だけで決まる。

## ★★★★足場(2026-08-20 実測)

| mathlib | 役割 |
|---|---|
| `AdjoinRoot.evalEval` | `F` 有理点での評価準同型 `F[W] →+* F` |
| `Affine.equation_neg` | `−Q` も曲線の点であること |
| `Polynomial.eval₂_eq_sum_range` | 評価を有限和に開く |
| `IsFractionRing.lift` | 単射な環準同型を分数体へ延ばす |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalPt` | ★★★`F` 有理点での評価準同型 |
| `uSq_mul_translateX` | ★★★★分母を払った形 |
| `translateX_transcendental` | ★★★★★★**`translateX` は超越的** |
| `translateHom_injective` | ★★★★★★**環準同型は単射** |
| `translateFieldHom` | ★★★★★★**関数体の自己準同型** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★`F` 有理点での評価 -/

/-- ★★★**`F` 有理点での評価準同型** `F[W] →+* F`。 -/
noncomputable def evalPt (W : WeierstrassCurve.Affine F) {a b : F} (h : W.Equation a b) :
    W.CoordinateRing →+* F :=
  AdjoinRoot.evalEval (x := a) (y := b) (p := W.polynomial) h

theorem evalPt_genX (W : WeierstrassCurve.Affine F) {a b : F} (h : W.Equation a b) :
    evalPt W h (genX W) = a := by
  rw [evalPt, genX, CoordinateRing.mk, AdjoinRoot.evalEval_mk]
  simp [Polynomial.evalEval]

theorem evalPt_genY (W : WeierstrassCurve.Affine F) {a b : F} (h : W.Equation a b) :
    evalPt W h (genY W) = b := by
  rw [evalPt, genY, CoordinateRing.mk, AdjoinRoot.evalEval_mk]
  simp [Polynomial.evalEval]

theorem evalPt_algebraMap (W : WeierstrassCurve.Affine F) {a b : F} (h : W.Equation a b) (c : F) :
    evalPt W h (algebraMap F W.CoordinateRing c) = c := by
  have hc : algebraMap F W.CoordinateRing c = CoordinateRing.mk W (Polynomial.C (Polynomial.C c)) :=
    rfl
  rw [hc, evalPt, CoordinateRing.mk, AdjoinRoot.evalEval_mk]
  simp [Polynomial.evalEval]

/-! ## ★★★★分母を払う -/

/-- ★`u = x − x₀`。 -/
noncomputable def uElt (W : WeierstrassCurve.Affine F) (x₀ : F) : W.CoordinateRing :=
  genX W - algebraMap F W.CoordinateRing x₀

/-- ★`v = y − y₀`。 -/
noncomputable def vElt (W : WeierstrassCurve.Affine F) (y₀ : F) : W.CoordinateRing :=
  genY W - algebraMap F W.CoordinateRing y₀

/-- ★★`A = v² + a₁vu − (a₂ + x + x₀)u²` —— `(x−x₀)²·X(τ_Q)` の分子。 -/
noncomputable def aElt (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) : W.CoordinateRing :=
  (vElt W y₀) ^ 2 + algebraMap F W.CoordinateRing W.a₁ * vElt W y₀ * uElt W x₀
    - (algebraMap F W.CoordinateRing W.a₂ + genX W + algebraMap F W.CoordinateRing x₀)
      * (uElt W x₀) ^ 2

theorem map_uElt (W : WeierstrassCurve.Affine F) (x₀ : F) :
    algebraMap W.CoordinateRing W.FunctionField (uElt W x₀)
      = coordX W - algebraMap F W.FunctionField x₀ := by
  rw [uElt, map_sub, ← IsScalarTower.algebraMap_apply]
  rfl

theorem map_vElt (W : WeierstrassCurve.Affine F) (y₀ : F) :
    algebraMap W.CoordinateRing W.FunctionField (vElt W y₀)
      = coordY W - algebraMap F W.FunctionField y₀ := by
  rw [vElt, map_sub, ← IsScalarTower.algebraMap_apply]
  rfl

theorem map_uElt_ne_zero (W : WeierstrassCurve.Affine F) (x₀ : F) :
    algebraMap W.CoordinateRing W.FunctionField (uElt W x₀) ≠ 0 := by
  rw [map_uElt, sub_ne_zero]
  exact coordX_ne_const W x₀

/-- ★★★★**`(x − x₀)²·X(τ_Q) = A`** —— 分母を払った形。 -/
theorem uSq_mul_translateX (W : WeierstrassCurve.Affine F) (x₀ y₀ : F) :
    (algebraMap W.CoordinateRing W.FunctionField (uElt W x₀)) ^ 2 * translateX W x₀ y₀
      = algebraMap W.CoordinateRing W.FunctionField (aElt W x₀ y₀) := by
  have hu := map_uElt_ne_zero W x₀
  rw [translateX, WeierstrassCurve.Affine.addX, slopeFF, aElt]
  simp only [map_sub, map_add, map_mul, map_pow, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂,
    map_uElt, map_vElt, ← IsScalarTower.algebraMap_apply]
  rw [map_uElt] at hu
  have hcX : algebraMap W.CoordinateRing W.FunctionField (genX W) = coordX W := rfl
  rw [hcX]
  field_simp
  ring

/-! ## ★★★★★`−Q` での評価 -/

theorem evalPt_uElt (W : WeierstrassCurve.Affine F) {x₀ y₀ : F}
    (h : W.Equation x₀ (W.negY x₀ y₀)) : evalPt W h (uElt W x₀) = 0 := by
  rw [uElt, map_sub, evalPt_genX, evalPt_algebraMap, sub_self]

theorem evalPt_vElt (W : WeierstrassCurve.Affine F) {x₀ y₀ : F}
    (h : W.Equation x₀ (W.negY x₀ y₀)) :
    evalPt W h (vElt W y₀) = W.negY x₀ y₀ - y₀ := by
  rw [vElt, map_sub, evalPt_genY, evalPt_algebraMap]

/-- ★★★★★**`ev_{−Q}(A) = (negY(x₀,y₀) − y₀)²`** —— `Q` が 2 等分点でなければ 0 でない。 -/
theorem evalPt_aElt (W : WeierstrassCurve.Affine F) {x₀ y₀ : F}
    (h : W.Equation x₀ (W.negY x₀ y₀)) :
    evalPt W h (aElt W x₀ y₀) = (W.negY x₀ y₀ - y₀) ^ 2 := by
  rw [aElt]
  simp only [map_sub, map_add, map_mul, map_pow, evalPt_uElt, evalPt_vElt, evalPt_algebraMap,
    evalPt_genX]
  ring

/-! ## ★★★★★★超越性 -/

/-- ★★★★★★**`translateX` は `F` 上超越的である**(`Q` が 2 等分点でないとき)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★分母を払って `−Q` で評価すると、`u` の項がすべて消えて先頭項だけが残る。 -/
theorem translateX_transcendental (W : WeierstrassCurve.Affine F)
    {x₀ y₀ : F} (hQ : W.Equation x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀)
    {p : Polynomial F} (hp : p ≠ 0) :
    Polynomial.eval₂ (algebraMap F W.FunctionField) (translateX W x₀ y₀) p ≠ 0 := by
  intro hzero
  set ι := algebraMap W.CoordinateRing W.FunctionField with hι
  set U := ι (uElt W x₀) with hU
  set A := ι (aElt W x₀ y₀) with hA
  set T := translateX W x₀ y₀ with hT
  set n := p.natDegree with hn
  have hUT : U ^ 2 * T = A := uSq_mul_translateX W x₀ y₀
  have hpow : ∀ i ≤ n, T ^ i * U ^ (2 * n) = A ^ i * U ^ (2 * (n - i)) := by
    intro i hi
    have h1 : U ^ (2 * i) * T ^ i = A ^ i := by
      rw [← hUT, mul_pow, ← pow_mul]
    have h2n : 2 * n = 2 * i + 2 * (n - i) := by omega
    rw [h2n, pow_add, ← mul_assoc, mul_comm (T ^ i), h1]
  have hsum := Polynomial.eval₂_eq_sum_range (algebraMap F W.FunctionField) (x := T) (p := p)
  rw [hsum] at hzero
  have hmul : (∑ i ∈ Finset.range (n + 1),
      algebraMap F W.FunctionField (p.coeff i) * T ^ i) * U ^ (2 * n) = 0 := by
    rw [hzero, zero_mul]
  rw [Finset.sum_mul] at hmul
  have hrw : ∀ i ∈ Finset.range (n + 1),
      algebraMap F W.FunctionField (p.coeff i) * T ^ i * U ^ (2 * n)
        = ι (algebraMap F W.CoordinateRing (p.coeff i) * aElt W x₀ y₀ ^ i
            * uElt W x₀ ^ (2 * (n - i))) := by
    intro i hi
    have hi' : i ≤ n := by
      have := Finset.mem_range.1 hi
      omega
    rw [mul_assoc, hpow i hi', hι]
    simp only [map_mul, map_pow, ← IsScalarTower.algebraMap_apply]
    ring
  rw [Finset.sum_congr rfl hrw, ← map_sum] at hmul
  have hinj : Function.Injective ι := IsFractionRing.injective _ _
  have hCR : (∑ i ∈ Finset.range (n + 1), algebraMap F W.CoordinateRing (p.coeff i)
      * aElt W x₀ y₀ ^ i * uElt W x₀ ^ (2 * (n - i))) = 0 := by
    apply hinj
    rw [hmul, map_zero]
  have hneg : W.Equation x₀ (W.negY x₀ y₀) := (WeierstrassCurve.Affine.equation_neg x₀ y₀).2 hQ
  have hev := congrArg (evalPt W hneg) hCR
  rw [map_zero, map_sum] at hev
  simp only [map_mul, map_pow, evalPt_uElt, evalPt_aElt, evalPt_algebraMap] at hev
  rw [Finset.sum_eq_single n] at hev
  · rw [Nat.sub_self] at hev
    simp only [Nat.mul_zero, pow_zero, mul_one] at hev
    have hlead : p.coeff n ≠ 0 := by
      rw [hn]
      exact Polynomial.leadingCoeff_ne_zero.2 hp
    have hA0 : (W.negY x₀ y₀ - y₀) ≠ 0 := sub_ne_zero.2 h2
    exact absurd hev (mul_ne_zero hlead (pow_ne_zero _ (pow_ne_zero _ hA0)))
  · intro i hi hne
    have hi' : i < n := by
      have := Finset.mem_range.1 hi
      omega
    have hz : 2 * (n - i) ≠ 0 := by omega
    rw [zero_pow hz, mul_zero]
  · intro h
    exact absurd (Finset.self_mem_range_succ n) h

/-! ## ★★★★★★単射性と関数体への延長 -/

/-- ★★★★★★**平行移動の環準同型は単射である**(`Q` が 2 等分点でないとき)。 -/
theorem translateHom_injective (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    Function.Injective (translateHom W hQ) :=
  translateHom_injective_of_transcendental W hQ
    (fun _ hp => translateX_transcendental W hQ.1 h2 hp)

/-- ★★★★★★**平行移動 `τ_Q` の関数体の自己準同型**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
noncomputable def translateFieldHom (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    W.FunctionField →+* W.FunctionField :=
  IsFractionRing.lift (translateHom_injective W hQ h2)

theorem translateFieldHom_coordX (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateFieldHom W hQ h2 (coordX W) = translateX W x₀ y₀ := by
  rw [translateFieldHom, coordX, IsFractionRing.lift_algebraMap, translateHom_genX]

theorem translateFieldHom_coordY (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (h2 : W.negY x₀ y₀ ≠ y₀) :
    translateFieldHom W hQ h2 (coordY W) = translateY W x₀ y₀ := by
  rw [translateFieldHom, coordY, IsFractionRing.lift_algebraMap, translateHom_genY]

/-! ## ★出典の紐付け(`.src`) -/

def translateFieldHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動が関数体に作用すること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
