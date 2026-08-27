/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightClass
import ABC3.Found.GenEll.BaseChange

/-!
# [GenEll] 垂直なひねりは高さを定数だけ動かす（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

## ★★★★★★★★原文が 1 行で畳んでいる段

`Proposition 1.4, (ii)` の証明の途中に、上の 1 文がある。
続く式は `ht_{L⊗M} ≈_{X(ℚ̄)} ht_L` である。

★★これは **`Remark 1.4.1` を「2 つの `ℤ`-モデルにまたがる形」で言うときの鍵**でもある:
2 つのモデルの算術直線束は `ℚ` 上一致していても、
`Σ` の上（＝**垂直**な方向）でずれ得る。そのずれが `Spec ℤ` から来るなら、
高さは**定数だけ**動く。

## ★★★★★★実際には `≈` ではなく**等式**である

`M` が `Spec ℤ` から来るなら、点 `x : Spec 𝓞_F → X` に沿った引き戻しは

    x^*(f^*M) = (f ∘ x)^* M = (Spec 𝓞_F → Spec ℤ)^* M

であり、これは `M` の `ADiv ℚ` としての**底変換**そのものである。
★★`degNormalized` は底変換で不変（`BaseChange.lean` の `degNormalized_baseChange`）だから、

    ht_{D}(x) − ht_{E}(x) = degNormalized(N)   （**点にも定義体にも依らない**）

★★★原文の `≈` は、この**等式**の帰結である。定数は `|deg_ℚ(N)|`。

## ★逸脱（明示）

★`f^*` を `ArithCartier` の水準で構成する代わりに、
**「引き戻した算術因子の差が固定の `ADiv ℚ` の底変換である」**を仮説に置いた。
★★`f^*M` を作れば仮説は自動で満たされるので、制限にはならない
——構成の側は `Definition 1.1` の担当である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

variable (F : Type) [Field F] [NumberField F]

/-! ## ★正規化次数は差を保つ -/

/-- ★**正規化次数は差を保つ**（`degNormalized_add` から）。 -/
theorem degNormalized_sub (a b : ADiv F) :
    degNormalized (a - b) = degNormalized a - degNormalized b := by
  have h := degNormalized_add (a - b) b
  rw [sub_add_cancel] at h
  linarith

/-! ## ★★★★★★★★垂直なひねりは高さを定数だけ動かす -/

/-- ★★★★★★★★**`Spec ℤ` から来る分だけ違うなら、高さの差は定数**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★原文は続けて `ht_{L⊗M} ≈ ht_L` と書く。
★★実際には**等式**で、差は `deg_ℚ(N)`——**点にも定義体にも依らない**。

★★★機構は `degNormalized_baseChange`（`BaseChange.lean`）1 本である。 -/
theorem htArith_sub_eq_of_baseChange {X : Scheme.{0}} (D E : ArithCartier X)
    (N : ADiv ℚ) (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackADiv F D xF - pullbackADiv F E xF = baseChange ℚ F N) :
    htArith F D xF - htArith F E xF = degNormalized N := by
  rw [htArith, htArith, ← degNormalized_sub, h, degNormalized_baseChange]

/-- ★★★★★★★★**原文の `≈` そのもの** —— 定数は `|deg_ℚ(N)|`。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★★定数が**点にも定義体にも依らない**ことが本質である
——`N` は `X` の上ではなく `Spec ℤ` の上にあるからである。 -/
theorem htArith_bdeq_of_baseChange {X : Scheme.{0}} (D E : ArithCartier X)
    (N : ADiv ℚ)
    (h : ∀ xF : specRingOfIntegers F ⟶ X,
      pullbackADiv F D xF - pullbackADiv F E xF = baseChange ℚ F N) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E xF) :=
  ⟨|degNormalized N|, fun xF => by
    show |htArith F D xF - htArith F E xF| ≤ |degNormalized N|
    rw [htArith_sub_eq_of_baseChange F D E N xF (h xF)]⟩

/-- ★★★★★★★★★**2 つの `ℤ`-モデルにまたがる形**（点の対応 `ePt` を挟む）。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★`Remark 1.4.1`（`Remark141.lean`）が消費する形。
★★`D` と `E` は**別のスキームの上**にあってよい——
`ℚ` 上対応していれば差は垂直で、`Spec ℤ` から来る分だけである。 -/
theorem htArith_bdeq_of_baseChange' {X X' : Scheme.{0}}
    (D : ArithCartier X) (E : ArithCartier X')
    (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'))
    (N : ADiv ℚ)
    (h : ∀ xF, pullbackADiv F D xF - pullbackADiv F E (ePt xF) = baseChange ℚ F N) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF)) :=
  ⟨|degNormalized N|, fun xF => by
    show |htArith F D xF - htArith F E (ePt xF)| ≤ |degNormalized N|
    rw [htArith, htArith, ← degNormalized_sub, h xF, degNormalized_baseChange]⟩

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (ii)` の証明の中の 1 段である。命題全体ではない。 -/

def htArith_sub_eq_of_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——Spec ℤ から来る分は高さを定数だけ動かす)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——ht_{L⊗M} ≈ ht_L)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_baseChange'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——2 つの ℤ-モデルにまたがる形)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degNormalized_baseChange(正規化次数の底変換不変性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.degNormalized_baseChange") 4,
    .implicitStep
      ("★★★原文は ≈ と書くが、実際には**等式**である" ++
       "——差はちょうど deg_ℚ(N) で、点にも定義体にも依らない") 6,
    .implicitStep
      ("★逸脱: f^* を ArithCartier の水準で構成する代わりに、" ++
       "「引き戻した差が固定の ADiv ℚ の底変換である」を仮説に置いた。" ++
       "★★f^*M を作れば仮説は自動で満たされるので制限にはならない") 6 ]

end ABC3.Found.GenEll
