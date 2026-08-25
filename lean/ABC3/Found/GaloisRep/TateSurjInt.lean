import ABC3.Found.GaloisRep.TatePhiInj

/-!
# Galois (G6) 第 286 ブロック —— **★★★★★★★★整な点についての全射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `x, y ∈ R` が曲線の式を満たすなら、対 `(a,w)` が在って `X(a,w) = x`、`Y(a,w) = y`
> (`tate_surjective_integral`)

★★★局所環の二分法で `x` が単元か `I` の元かに分け、第 266(単元の領域)と
第 271(環帯)を貼り合わせるだけである。

## ★★★★★`x ∈ I` なら `y ∈ I` は式から出る

    y² = (x³ + a₄x + a₆) − xy

★`a₄, a₆ ∈ I`(第 212)、`x ∈ I` なので右辺は `I` の元。したがって `y² ∈ I`、
**`I` が素なら `y ∈ I`**(`tate_y_mem_of_x_mem`)。
★★環帯の全射性(第 271)は `x, y` の両方が `I` に入ることを要求するので、この一行が要る。

## ★★残っているもの——原点近傍の全射性

`x` が `R` の像に**入らない**点(原点近傍)については、まだ全射性が無い。要るのは
**付値の評価**である:

    v(x) < 0  ⟹  v(y) = (3/2)·v(x)  ⟹  v(−x/y) = −v(x)/2 > 0

★これは `y² + xy = x³ + a₄x + a₆` の Newton 多角形の議論で、
`v(x³) = 3v(x)` が右辺で厳密に最小になることから出る。
★★★`z = −x/y ∈ 𝔪` が言えれば第 276 の `exists_tateZ_eq` で母数が取れ、
`s = −1/y` の一意性(縮小)で点が一致する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_y_mem_of_x_mem` | ★★★★★`x ∈ I` なら `y ∈ I` |
| `tate_surjective_integral` | ★★★★★★★★**整な点についての全射性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★**曲線上の `R` 点では `x ∈ I` なら `y ∈ I`**——`a₄, a₆ ∈ I` だから。 -/
theorem tate_y_mem_of_x_mem (hI : I.IsPrime) (x y q : R) (hq : q ∈ I) (hxm : x ∈ I)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    y ∈ I := by
  have h4 := tateCurveAt_a4_mem (I := I) q hq
  have h6 := tateCurveAt_a6_mem (I := I) q hq
  have hsq : y * y ∈ I := by
    have hkey : y * y
        = (x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) - x * y := by
      linear_combination he
    rw [hkey]
    refine Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _ ?_ ?_) h6) ?_
    · rw [pow_succ]
      exact Ideal.mul_mem_left _ _ hxm
    · exact Ideal.mul_mem_right _ _ h4
    · exact Ideal.mul_mem_right _ _ hxm
  rcases hI.mem_or_mem hsq with h | h <;> exact h

/-- ★★★★★★★★**整な点についての全射性**——単元の領域と環帯を貼り合わせた。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_surjective_integral (hloc : ∀ x : R, IsUnit x ∨ x ∈ I) (hI : I.IsPrime)
    (x y q : R) (hq : q ∈ I)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    ∃ a w : R, a * w = q ∧ w ∈ I ∧
      tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  rcases hloc x with hxu | hxm
  · obtain ⟨a, w, _, hwI, haw, hX, hY⟩ := tate_surjective x y q hq hxu he
    exact ⟨a, w, haw, hwI, hX, hY⟩
  · obtain ⟨a, _, w, hwI, haw, hX, hY⟩ :=
      tate_surjective_annulus x y q hq hxm (tate_y_mem_of_x_mem hI x y q hq hxm he) he
    exact ⟨a, w, haw, hwI, hX, hY⟩

/-! ## ★出典の紐付け(`.src`) -/

def tate_surjective_integral.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——整な点についての全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
