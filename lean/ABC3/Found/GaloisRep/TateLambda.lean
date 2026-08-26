import ABC3.Found.GaloisRep.TateLinearize

/-!
# Galois (G6) 第 263 ブロック —— **★★★★★★★★母数座標と点の一意性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★母数座標 `Λ = Y/(X+Y)`

第 262 で `Y = a(X+Y)`(主要項)が出た。`X + Y` は単元なので

    Λ(a) := Y·(X+Y)⁻¹     ⟹     Λ(a) ≡ a  (mod I)

★**微分が 1 の座標**である。`X` や `Y` を単独で使うと微分が `f'(a) = (1+a)/(1−a)³` や
`g'(a) = a(2+a)/(1−a)⁴` になって、`a ≡ −1` や `a ≡ 0, −2` で単元でなくなるが、
`Λ` にはその問題がない。

## ★★★★★★★★点の一意性 —— 節点立方曲線の径数づけ

`I` を法とすると Tate 曲線は**節点立方曲線** `y² + xy = x³` になる。その滑らかな点は

    t ↦ (t/(1−t)², t²/(1−t)³)

で 1 対 1 に径数づけられる。★その `I` 進版が `curve_point_unique_of_rel`:

> `(x,y)`・`(x',y')` が曲線上、`x`・`x'` が単元、`y x' = y' x` なら `(x,y) = (x',y')`

### ★★★★★★証明は「`x` で消去」の 1 行

`y = μx`、`y' = μx'` と書けば(`μ := y/x`、`hrel` から `y' = μx'`)、
曲線の式を `x'²`・`x²` 倍して引くだけで

    (x − x')·(x²x'² − a₄xx' − a₆(x+x')) = 0

★★括弧は `x²x'²`(単元)+ `I` の元なので**単元**(`isUnit_add_mem`)。よって `x = x'`。
`μ²+μ` の消去も判別式も要らない——**`ring` と `linear_combination` 1 回**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateLambda` | ★★★★★★母数座標 `Y/(X+Y)` |
| `tateLambda_sub_mem` | ★★★★★★★`Λ(a) ≡ a` (mod `I`) |
| `eq_of_param_curve` | ★★★★★★★★径数づけられた点は `x` で決まる |
| `curve_point_unique_of_rel` | ★★★★★★★★**点の一意性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★母数座標 -/

/-- ★★★★★★**Tate 座標の「母数座標」** `Λ = Y/(X+Y)`。 -/
noncomputable def tateLambda [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) : R :=
  tateYpair a w q hq * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)

/-- ★★★★★★★**`Λ(a) ≡ a` (mod `I`)**——微分が 1 の座標。

★`X` や `Y` を単独で使うと微分が `a ≡ −1` や `a ≡ 0, −2` で単元でなくなるが、
`Λ` にはその問題がない。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateLambda_sub_mem [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (ha : IsUnit a) (hu : IsUnit (1 - a)) : tateLambda a w q hq - a ∈ I := by
  have hS : IsUnit (tateXpair a w q hq + tateYpair a w q hq) :=
    isUnit_tateXpair_add_tateYpair a w q hq hw ha hu
  have hkey : tateLambda a w q hq - a
      = (tateYpair a w q hq - a * (tateXpair a w q hq + tateYpair a w q hq))
        * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq) := by
    rw [tateLambda]
    have h2 : a = a * ((tateXpair a w q hq + tateYpair a w q hq)
        * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)) := by
      rw [Ring.mul_inverse_cancel _ hS, mul_one]
    calc tateYpair a w q hq * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq) - a
        = tateYpair a w q hq * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)
          - a * ((tateXpair a w q hq + tateYpair a w q hq)
            * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)) := by rw [← h2]
      _ = _ := by ring
  rw [hkey]
  exact Ideal.mul_mem_right _ _ (tateYpair_sub_mul_mem a w q hq hw hu)

/-! ## ★★★★★★★★点の一意性 -/

/-- ★★★★★★★★**径数づけられた曲線の点は `x` で決まる**——`x` が単元なら。

★曲線の式を `x'²`・`x²` 倍して引くと `(x − x')·(x²x'² − a₄xx' − a₆(x+x')) = 0`。 -/
theorem eq_of_param_curve [IsAdicComplete I R] {a₄ a₆ μ x x' : R}
    (h4 : a₄ ∈ I) (h6 : a₆ ∈ I)
    (he : (μ * x) ^ 2 + x * (μ * x) = x ^ 3 + a₄ * x + a₆)
    (he' : (μ * x') ^ 2 + x' * (μ * x') = x' ^ 3 + a₄ * x' + a₆)
    (hx : IsUnit x) (hx' : IsUnit x') : x = x' := by
  have hB : IsUnit (x ^ 2 * x' ^ 2 - a₄ * (x * x') - a₆ * (x + x')) := by
    have h1 : x ^ 2 * x' ^ 2 - a₄ * (x * x') - a₆ * (x + x')
        = x ^ 2 * x' ^ 2 + (-(a₄ * (x * x')) - a₆ * (x + x')) := by ring
    rw [h1]
    exact isUnit_add_mem ((hx.pow 2).mul (hx'.pow 2))
      (Ideal.sub_mem _ (neg_mem (Ideal.mul_mem_right _ _ h4)) (Ideal.mul_mem_right _ _ h6))
  have hz : (x - x') * (x ^ 2 * x' ^ 2 - a₄ * (x * x') - a₆ * (x + x')) = 0 := by
    linear_combination x ^ 2 * he' - x' ^ 2 * he
  exact sub_eq_zero.1 ((IsUnit.mul_left_eq_zero hB).1 hz)

/-- ★★★★★★★★**曲線の点は `y x' = y' x` と `x` の単元性で決まる**。

★`I` を法とすると Tate 曲線は節点立方曲線 `y² + xy = x³` になり、その滑らかな点は
`t ↦ (t/(1−t)², t²/(1−t)³)` で 1 対 1 に径数づけられる。その `I` 進版である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem curve_point_unique_of_rel [IsAdicComplete I R] {a₄ a₆ x y x' y' : R}
    (h4 : a₄ ∈ I) (h6 : a₆ ∈ I)
    (he : y ^ 2 + x * y = x ^ 3 + a₄ * x + a₆)
    (he' : y' ^ 2 + x' * y' = x' ^ 3 + a₄ * x' + a₆)
    (hx : IsUnit x) (hx' : IsUnit x') (hrel : y * x' = y' * x) :
    x = x' ∧ y = y' := by
  obtain ⟨ux, hux⟩ := hx
  have hxxi : x * ((ux⁻¹ : Rˣ) : R) = 1 := by
    rw [← hux, ← Units.val_mul, mul_inv_cancel, Units.val_one]
  set μ : R := y * ((ux⁻¹ : Rˣ) : R) with hμdef
  have hμx : μ * x = y := by
    rw [hμdef]
    calc y * ((ux⁻¹ : Rˣ) : R) * x = y * (x * ((ux⁻¹ : Rˣ) : R)) := by ring
      _ = y := by rw [hxxi, mul_one]
  have hμx' : μ * x' = y' := by
    rw [hμdef]
    calc y * ((ux⁻¹ : Rˣ) : R) * x' = (y * x') * ((ux⁻¹ : Rˣ) : R) := by ring
      _ = (y' * x) * ((ux⁻¹ : Rˣ) : R) := by rw [hrel]
      _ = y' * (x * ((ux⁻¹ : Rˣ) : R)) := by ring
      _ = y' := by rw [hxxi, mul_one]
  have hxx' := eq_of_param_curve (I := I) h4 h6 (by rw [hμx]; exact he)
    (by rw [hμx']; exact he') ⟨ux, hux⟩ hx'
  exact ⟨hxx', by rw [← hμx, ← hμx', hxx']⟩

/-! ## ★出典の紐付け(`.src`) -/

def tateLambda.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数座標 Y/(X+Y))",
    sectionId := "genell-def-3-3" }

def curve_point_unique_of_rel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——曲線の点の一意性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
