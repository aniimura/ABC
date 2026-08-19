import ABC3.Found.GaloisRep.QDomain
import ABC3.Found.GaloisRep.TateXY
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point

/-!
# Galois (G6) 第 107 ブロック —— **★★★★★`X(u,q)`・`Y(u,q)` と反転則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★対 `(a, w)` で組む

第 106 ブロックの基本領域 `0 ≤ v(u) < v(q)` を取ると、
★`u` も `q/u` も付値が非負なので**どちらも `R` の元である**。
★★そこで `a ↦ u`、`w ↦ q/u`(`a·w = q`)という**対**を入力にする。
★★★すると両側和は次のように `R` の中だけで書ける:

    X = [f(a) + ∑_{m≥1} f(qᵐa)] + [f(w) + ∑_{m≥1} f(qᵐw)] − 2 s₁(q)
    Y = [g(a) + ∑_{m≥1} g(qᵐa)] − [f(w) + ∑_{m≥1} f(qᵐw)]
                                 − [g(w) + ∑_{m≥1} g(qᵐw)] + s₁(q)

★★★★`n ≤ −1` の項は `qⁿu = (q^{m−1}w)⁻¹` と書き直し、
第 105 ブロックの `f(1/t) = f(t)`・`g(t) + g(1/t) = −f(t)` を使って
**`w` の側の片側和に化けさせた**。

## ★★★★★★これで曲線の**反転則**が出る

`a` と `w` を入れ替えると `u ↦ q/u` であり、`Kˣ/q^ℤ` では `u ↦ u⁻¹` と同じ類である。
★上の式は入れ替えに対して

    X(q/u) = X(u),     Y(q/u) = −X(u) − Y(u)

を満たす。★★Tate 曲線は `a₁ = 1`, `a₃ = 0` なので mathlib の
`negY x y = −y − a₁x − a₃ = −y − x` と**一致する**:

    (X(q/u), Y(q/u)) = −(X(u), Y(u))     (`tateYpair_eq_negY`)

★★★★★★**級数と mathlib の曲線が初めて繋がった。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpair` / `tateYpair` | ★★★★★**`X(u,q)`・`Y(u,q)`** |
| `tateXpair_symm` | ★★★`X(q/u) = X(u)` |
| `tateYpair_swap` | ★★★★`Y(q/u) = −X(u) − Y(u)` |
| `tateYpair_eq_negY` | ★★★★★★**mathlib の `negY` と一致** |
| `exists_pair_of_normalized` | ★★★基本領域の代表元からは常に対が取れる |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★級数 -/

/-- ★★★★★**`X(u,q)`** —— `a` が `u`、`w` が `q/u` を表す(`a·w = q`)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`1 − a` が単元でないとき(`u ≡ 1 mod 𝔪`、すなわち原点の近く)は
`Ring.inverse` の既定値が入るので、この `R` 値の式は正しくない。
★★その場合は `K` の中で取り直すことになる——原点の近傍は形式群の側で扱う。 -/
noncomputable def tateXpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  (tateXterm a + tateXtail a q hq) + (tateXterm w + tateXtail w q hq)
    - 2 * evalAdic (sigmaSeries 1) q hq

/-- ★★★★★**`Y(u,q)`**。 -/
noncomputable def tateYpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  (tateYterm a + tateYtail a q hq) - (tateXterm w + tateXtail w q hq)
    - (tateYterm w + tateYtail w q hq) + evalAdic (sigmaSeries 1) q hq

/-! ## ★★★★反転則 -/

/-- ★★★**`X(q/u) = X(u)`**。 -/
theorem tateXpair_symm [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateXpair w a q hq = tateXpair a w q hq := by
  simp only [tateXpair]
  ring

/-- ★★★★**`Y(q/u) = −X(u) − Y(u)`**。 -/
theorem tateYpair_swap [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateYpair w a q hq = - tateXpair a w q hq - tateYpair a w q hq := by
  simp only [tateXpair, tateYpair]
  ring

/-- ★★★★★★**Tate 級数は曲線の反転則と整合する**
——`(X(q/u), Y(q/u)) = −(X(u), Y(u))`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★Tate 曲線は `a₁ = 1`, `a₃ = 0` なので `negY x y = −y − x` である。 -/
theorem tateYpair_eq_negY [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) :
    tateYpair w a q hq
      = (tateCurveAt q hq).toAffine.negY (tateXpair a w q hq) (tateYpair a w q hq) := by
  rw [WeierstrassCurve.Affine.negY]
  have h1 : (tateCurveAt q hq).a₁ = 1 := by
    simp [tateCurveAt, tateCurve, WeierstrassCurve.map]
  have h3 : (tateCurveAt q hq).a₃ = 0 := by
    simp [tateCurveAt, tateCurve, WeierstrassCurve.map]
  show tateYpair w a q hq
      = - tateYpair a w q hq - (tateCurveAt q hq).a₁ * tateXpair a w q hq
        - (tateCurveAt q hq).a₃
  rw [h1, h3, tateYpair_swap]
  ring

/-! ## ★★★基本領域との接続 -/

/-- ★★★**基本領域の代表元からは常に対 `(a, w)` が取れる**。

★`v(u) ≥ 0` なので `u` は `R` から来る。★★`v(q/u) = v(q) − v(u) > 0` なので `q/u` もそう。 -/
theorem exists_pair_of_normalized {K : Type} [Field K] [Algebra R K]
    (v : Kˣ →* Multiplicative ℤ) (Q u : Kˣ)
    (hu0 : 0 ≤ vAdd v u) (hu1 : vAdd v u < vAdd v Q)
    (hmem : ∀ x : Kˣ, 0 < vAdd v x → ∃ y ∈ I, algebraMap R K y = (x : K))
    (hmem0 : ∀ x : Kˣ, 0 ≤ vAdd v x → ∃ y : R, algebraMap R K y = (x : K)) :
    ∃ a w : R, algebraMap R K a = (u : K) ∧ w ∈ I ∧
      algebraMap R K w = ((Q * u⁻¹ : Kˣ) : K) ∧
      algebraMap R K (a * w) = (Q : K) := by
  obtain ⟨a, ha⟩ := hmem0 u hu0
  have hv : 0 < vAdd v (Q * u⁻¹) := by
    rw [vAdd_mul, vAdd_inv]
    omega
  obtain ⟨w, hwI, hw⟩ := hmem _ hv
  refine ⟨a, w, ha, hwI, hw, ?_⟩
  rw [map_mul, ha, hw]
  push_cast
  field_simp

/-! ## ★出典の紐付け(`.src`) -/

def tateXpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化の座標 X(u,q)・Y(u,q) と反転則)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
