import ABC3.Found.GenEll.ArchCompat
import ABC3.Found.GenEll.ArchADivBase
import ABC3.Found.GenEll.HeightBaseChange

/-!
# [GenEll] Definition 1.2, (i) —— **高さは定義体に依らない(仮定を discharge した形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★本日の到達点をひとつの定理にする

`HeightBaseChange.lean` の `htArith_baseChange` は 3 つの仮定を受けていた。
★本日そのすべてを定理にした:

| 仮定 | 定理にした場所 |
|---|---|
| `hfin`(有限素点側) | `PullbackNatural.lean` |
| `harch`(アルキメデス側) | `ArchADivBase.lean` + `ArchCompat.lean` |
| `hlies` | `FinitePlaceRel.lean`(instance) |

★★★**残る仮定は 2 つだけで、どちらも原文自身のものである:**

- `IsConjInvariant D.green` —— 原文の「計量は `ι_X` と両立する」(`Definition 1.1`)
- `x_F^* D ≠ 0` —— 「`x` が `D` を通らない」

## ★★2 つ目の仮定は因子表示の境界である

原文の `ht_L̄` は**可逆層**の引き戻しなので条件が要らない。
★因子で表すと「通らない」が要るのは**表示の側の制約**である。

★★★**原文が `Proposition 1.6` を `U(ℚ̄) = X(ℚ̄) \ D` の上で述べているのは、
まさにこの境界に沿っている**——本定理はそこでそのまま使える。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
  [Algebra F K]

/-! ## ★★★仮定を discharge した底変換不変性 -/

open scoped Classical in
/-- ★★★**[GenEll] Definition 1.2, (i)** —— 高さは定義体の取り方に依らない
(**仮定を discharge した形**)。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**残る仮定は 2 つだけで、どちらも原文自身のものである:**

- `hg : IsConjInvariant D.green` —— 原文の「計量は `ι_X` と両立する」
- `hJ : x_F^* D ≠ 0` —— 「`x` が `D` を通らない」

★★これで `ht` が `X(ℚ̄) \ D` の上で well-defined になる。
★原文が `Proposition 1.6` を `U(ℚ̄)` の上で述べているのと同じ範囲である。 -/
theorem htArith_baseChange_natural {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (hg : IsConjInvariant D.green)
    (hJ : pullbackIdeal F D.divisor xF ≠ 0) :
    htArith K D
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = htArith F D xF := by
  refine htArith_baseChange F K D xF
    (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) hJ ?_ ?_ ?_
  · -- ★有限素点側 —— `PullbackNatural.lean`
    exact pullbackIdeal_specMap F K D.divisor xF _
  · -- ★アルキメデス側 —— `ArchADivBase.lean` + `ArchCompat.lean`
    exact archADiv_baseChange F K D.green hg xF _
      (fun w => archRingHom_compat F K w) (pullbackADiv F D xF) rfl
  · -- ★`LiesOver` は instance
    exact fun W => inferInstance

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2` 全体には
`X(ℚ̄)` の型そのものの構成(数体についての colimit)と、
`D` を通らない代表を選ぶ移動補題が残っている。 -/

def htArith_baseChange_natural.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(高さの底変換不変性——U(ℚ̄) の上で)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
