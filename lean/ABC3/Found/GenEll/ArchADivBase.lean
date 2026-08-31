import ABC3.Found.GenEll.ArchBaseChange
import ABC3.Found.GenEll.IdealADivBase

/-!
# [GenEll] Definition 1.2, (i) —— **アルキメデス側を `ADiv` に繋ぐ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★段 2 —— 直前の訂正がこれを可能にした

`archADiv` の係数に**局所次数 `mult` の重みが抜けていた**ことを直前に訂正した。
★★★**その訂正がそのまま段 2 を通す。**

    `baseChangeArc F K a w = (mult w / mult v) · a.arc v`   (`v = w|_F`)

- 訂正**前**: `archADiv K g x_K w = g(x_w) = g(x_v)` で、
  `(mult w / mult v)` の分だけ**合わなかった**
- 訂正**後**: `archADiv K g x_K w = mult w · g(x_w) = mult w · g(x_v)`
  であり、`baseChangeArc` の右辺も
  `(mult w / mult v) · (mult v · g(x_v)) = mult w · g(x_v)` ★★**一致する**

★★★**定義の誤りを直したら、次の段がそのまま通った。**

## ★仮定

- `IsConjInvariant g` —— 原文の「計量は `ι_X` と両立する」
- 埋め込みの両立(共役を除く)—— `ArchBaseChange.lean` の場合分けを通す
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
  [Algebra F K]

/-! ## ★局所次数は正 -/

/-- ★**局所次数は 0 でない**(`mult` は 1 か 2)。 -/
theorem mult_ne_zero (v : InfinitePlace F) : (InfinitePlace.mult v : ℝ) ≠ 0 := by
  have : 0 < InfinitePlace.mult v := InfinitePlace.mult_pos
  positivity

/-! ## ★★★アルキメデス側の底変換 -/

open scoped Classical in
/-- ★★★**Green 関数の引き戻しは底変換と両立する**(アルキメデス側)。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `archADiv K g x_K = baseChangeArc F K (archADiv F g x_F)`

★★★**これが `Definition 1.2, (i)` の段 2 である。**

★機構は `green_archPoint_baseChange`(共役の場合分けを `IsConjInvariant` が吸収)と
`mult v ≠ 0` による約分だけ。 -/
theorem archADiv_baseChange {X : Scheme.{0}} (g : GreenFn X)
    (hg : IsConjInvariant g) (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K))
    (hcompat : ∀ w : InfinitePlace K,
      (archRingHom K w).comp φ.hom
          = archRingHom F (w.comap (algebraMap F K))
        ∨ (archRingHom K w).comp φ.hom
          = (starRingEnd ℂ).comp (archRingHom F (w.comap (algebraMap F K))))
    (a : ADiv F) (ha : a.arc = archADiv F g xF) :
    archADiv K g (Spec.map φ ≫ xF) = baseChangeArc F K a := by
  ext w
  set v := w.comap (algebraMap F K) with hv
  have hpt : g (archPoint (Spec.map φ ≫ xF) w) = g (archPoint xF v) :=
    green_archPoint_baseChange F K v w xF φ g hg (hcompat w)
  rw [archADiv_apply, hpt]
  show _ = ((InfinitePlace.mult w : ℝ) / (InfinitePlace.mult v : ℝ)) * a.arc v
  rw [ha, archADiv_apply]
  field_simp

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
`degNormalized_baseChange` との接続と `X(ℚ̄)` の型は `AlgPointClass.lean`(§9-744)で入った。 -/

def archADiv_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(アルキメデス側を ADiv に繋ぐ段のみ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
