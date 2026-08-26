import ABC3.Found.GaloisRep.TateSurjInt

/-!
# Galois (G6) 第 287 ブロック —— **★★★★★★★★★★形式母数は `𝔪` の元(付値なし)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `1/x ∈ 𝔪`(すなわち `x` が `R` の外)で `(x,y)` が曲線上なら、
> **`z ∈ 𝔪` が在って `z·y = x`**(`tate_formal_param`)

★★★つまり形式群の座標 `x/y` が `𝔪` の元である。第 276 の `exists_tateZ_eq` に渡す
入力がこれで揃った。

## ★★★★★★★★★★付値が要らなかった

第 286 では「全射性は付値が要る——極の位数が要るから」と書いた。**これは外れである。**
`α := 1/x ∈ 𝔪` と置き、**`γ := y/x²` を見ると**

    γ² + α·γ = α·(1 + a₄α² + a₆α³) = α·u        (`u` は単元)

★★★★これは `γ` についての**モニックな 2 次式**である。したがって

| 段 | 使うもの | 結論 |
|---|---|---|
| 1 | `R` が整閉 | `γ ∈ R` |
| 2 | `γ² = αu − αγ ∈ 𝔪`、`𝔪` が素 | `γ ∈ 𝔪` |
| 3 | `α·u = γ(γ+α)` を `γ` で割る | `x/y = (γ+α)·u⁻¹ ∈ 𝔪` |

★★★★★**極の位数(`v(y) = (3/2)v(x)`)を数える代わりに、`y/x²` が整であることを使う。**
Newton 多角形の議論が、モニック多項式ひとつに化けた。
★★整閉性は「モニック 2 次式の根は `R` に入る」という形の仮定 `hquad` で受ける。

## ★★★なぜ `y/x²` なのか

`v(x) = −2n` なら `v(y) = −3n` なので `v(y/x²) = −3n + 4n = n > 0`。
★つまり `y/x²` は「ちょうど `𝔪` に落ちる」正規化である。`y/x` では `v = −n < 0`、
`y/x³` では `v = 3n` で情報が落ちる。★★**2 乗がただ一つの正しい重み**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_formal_param` | ★★★★★★★★★★**形式母数は `𝔪` の元** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**形式母数は `𝔪` の元**——付値を使わずに `x/y ∈ 𝔪` を出す。

★`γ := y/x²` がモニックな 2 次式 `γ² + αγ = αu` を満たすことがすべてである
(`α := 1/x`、`u := 1 + a₄α² + a₆α³` は単元)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_formal_param (hinj : Function.Injective (algebraMap R K)) (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (q : R) (hq : q ∈ I) (x y : K) (hx0 : x ≠ 0) (α₀ : R) (hα : α₀ ∈ I)
    (hαx : algebraMap R K α₀ = x⁻¹)
    (he : y ^ 2 + x * y = x ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * x
      + algebraMap R K ((tateCurveAt q hq).a₆)) :
    ∃ z₀ ∈ I, algebraMap R K z₀ * y = x := by
  have hA4 := tateCurveAt_a4_mem (I := I) q hq
  have hB6 := tateCurveAt_a6_mem (I := I) q hq
  set A := (tateCurveAt q hq).a₄ with hA
  set B := (tateCurveAt q hq).a₆ with hB
  have hax : algebraMap R K α₀ * x = 1 := by rw [hαx, inv_mul_cancel₀ hx0]
  have ha0 : algebraMap R K α₀ ≠ 0 := by
    intro h
    rw [h, zero_mul] at hax
    exact zero_ne_one hax
  have hquadγ : (y * algebraMap R K α₀ ^ 2) ^ 2
      + algebraMap R K α₀ * (y * algebraMap R K α₀ ^ 2)
      + algebraMap R K (-(α₀ + A * α₀ ^ 3 + B * α₀ ^ 4)) = 0 := by
    have h1 : algebraMap R K (-(α₀ + A * α₀ ^ 3 + B * α₀ ^ 4))
        = -(algebraMap R K α₀ + algebraMap R K A * algebraMap R K α₀ ^ 3
          + algebraMap R K B * algebraMap R K α₀ ^ 4) := by
      simp only [map_neg, map_add, map_mul, map_pow]
    rw [h1]
    linear_combination algebraMap R K α₀ ^ 4 * he
      + (-(algebraMap R K α₀ ^ 3 * y)
        + algebraMap R K α₀ * ((algebraMap R K α₀ * x) ^ 2 + algebraMap R K α₀ * x + 1)
        + algebraMap R K A * algebraMap R K α₀ ^ 3) * hax
  obtain ⟨γ₀, hγ₀⟩ := hquad _ α₀ (-(α₀ + A * α₀ ^ 3 + B * α₀ ^ 4)) hquadγ
  have hR : γ₀ ^ 2 + α₀ * γ₀ - (α₀ + A * α₀ ^ 3 + B * α₀ ^ 4) = 0 := by
    refine hinj ?_
    rw [map_zero, map_sub, map_add, map_pow, map_mul, hγ₀]
    have h2 : algebraMap R K (α₀ + A * α₀ ^ 3 + B * α₀ ^ 4)
        = -algebraMap R K (-(α₀ + A * α₀ ^ 3 + B * α₀ ^ 4)) := by rw [map_neg, neg_neg]
    rw [h2]
    linear_combination hquadγ
  set u₀ := 1 + A * α₀ ^ 2 + B * α₀ ^ 3 with hu
  have hu₀ : IsUnit u₀ := by
    rw [hu, show (1 : R) + A * α₀ ^ 2 + B * α₀ ^ 3 = 1 + (A * α₀ ^ 2 + B * α₀ ^ 3) by ring]
    refine isUnit_add_mem (I := I) isUnit_one (Ideal.add_mem _ ?_ ?_)
    · refine Ideal.mul_mem_left _ _ ?_
      rw [pow_succ]
      exact Ideal.mul_mem_left _ _ hα
    · refine Ideal.mul_mem_left _ _ ?_
      rw [pow_succ]
      exact Ideal.mul_mem_left _ _ hα
  have hR2 : γ₀ ^ 2 + α₀ * γ₀ = α₀ * u₀ := by rw [hu]; linear_combination hR
  have hγI : γ₀ ∈ I := by
    have hsq : γ₀ * γ₀ ∈ I := by
      have h : γ₀ * γ₀ = α₀ * u₀ - α₀ * γ₀ := by linear_combination hR2
      rw [h]
      exact Ideal.sub_mem _ (Ideal.mul_mem_right _ _ hα) (Ideal.mul_mem_right _ _ hα)
    rcases hI.mem_or_mem hsq with h | h <;> exact h
  set z₀ := (γ₀ + α₀) * Ring.inverse u₀ with hz
  have hz₀I : z₀ ∈ I := Ideal.mul_mem_right _ _ (Ideal.add_mem _ hγI hα)
  have hzγ₀ : z₀ * γ₀ = α₀ := by
    have hiu : u₀ * Ring.inverse u₀ = 1 := Ring.mul_inverse_cancel _ hu₀
    rw [hz]
    linear_combination Ring.inverse u₀ * hR2 + α₀ * hiu
  refine ⟨z₀, hz₀I, ?_⟩
  have hzγ : algebraMap R K z₀ * (y * algebraMap R K α₀ ^ 2) = algebraMap R K α₀ := by
    rw [← hγ₀, ← map_mul, hzγ₀]
  have h4 : algebraMap R K z₀ * y * (algebraMap R K α₀ * x) ^ 2
      = x * (algebraMap R K α₀ * x) := by
    linear_combination x ^ 2 * hzγ
  rw [hax] at h4
  simpa using h4

/-! ## ★出典の紐付け(`.src`) -/

def tate_formal_param.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——形式母数は m の元)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
