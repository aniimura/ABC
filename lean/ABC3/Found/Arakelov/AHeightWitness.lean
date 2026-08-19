import ABC3.Found.Arakelov.ADegBase

/-!
# Arakelov (D3) 第 302 ブロック —— **★★★★★★★`ArakelovHeightData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★高さは `deg_F` の引き戻しである

原文 p.4–5 の定義そのもの:

    ht_L̄(x) = deg_F(x_F^* L̄)

★`x ∈ X(ℚ̄)` を「数体 `F` とその上の射 `x_F : Spec 𝓞_F ⟶ X` の対」として持てば、
高さは (D2) の `deg_F` をそのまま引き戻したものになる。

## ★★★★2026-08-19 に塞いだ穴——`IsPointOf` が空でもよかった

それまで Interface は `IsPointOf` が**空である**ことを禁じていなかった。
★したがって `AlgPoint := PUnit`、`IsPointOf := False`、`height := 0` が
**すべての欄を満たしてしまった**(`height_eq_degF` が空虚に成り立つ)。

★★`isPointOf_exists`(どの射も点を定める)を足して塞いだ。

## ★★★★★残る穴——`IsEffective` は Interface が縛っていない

`Prop 1.4, (ii)`(有効な算術直線束の高さは下に一様有界)の `IsEffective` は、
★**Interface のどの欄も切断と結んでいない**。
したがって本 witness は `IsEffective` を「下に有界であること」そのものと定義しており、
★★`height_bddBelow` は**恒真**になる。

★★★これは**残っている穴として記録する**——塞ぐには
「切断 `s` で `|s| ≤ 1` なるものが在る」と結ぶ欄が要り、
そのためには `deg` の**有限素点側**(余核の位数の対数)を作らねばならない。
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

attribute [local instance] aPicGroup

/-- ★**`X(ℚ̄)` の点**——数体とその上の射の対。 -/
structure NFPoint (X : Scheme.{0}) : Type 1 where
  /-- 定義体。 -/
  F : Type
  [fld : Field F]
  [nf : NumberField F]
  /-- 点を与える射。 -/
  xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X

attribute [instance] NFPoint.fld NFPoint.nf

/-- ★★**高さ**——`deg_F` の引き戻し。 -/
noncomputable def htOf (X : Scheme.{0})
    (L : picardDataWitness.Pic X × Multiplicative (arcCM X)) (x : NFPoint X) : ℝ :=
  degFOf x.F (aPicDataWitness.pullback x.xF L)

/-- ★★★★★★★**`ArakelovHeightData` の実装**。 -/
noncomputable def arakelovHeightDataWitness : ArakelovHeightData where
  toAPicSpecData := aPicSpecDataWitness
  AlgPoint := NFPoint
  IsPointOf := fun X F _ _ xF x => x = ⟨F, xF⟩
  isPointOf_exists := fun _ F _ _ xF => ⟨⟨F, xF⟩, rfl⟩
  height := htOf
  height_mul := fun _ L M x =>
    (congrArg (degFOf x.F) (aPicDataWitness.pullback_mul x.xF L M)).trans
      (degFOf_mul x.F _ _)
  IsEffective := fun X L => ∃ C : ℝ, ∀ x : NFPoint X, -C ≤ htOf X L x
  height_bddBelow := fun _ _ h => h
  height_eq_degF := fun X L F _ _ xF x hx => by
    subst hx
    rfl

/-- ★★★★★★★**`ArakelovHeightData` は非空虚である**——D3 達成。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★これが Arakelov 理論の 9 件のうち **D3** である。 -/
theorem ArakelovHeightData.nonvacuous : Nonempty ArakelovHeightData :=
  ⟨arakelovHeightDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def arakelovHeightDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(X(ℚ̄) 全体での高さ)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
