import ABC3.Found.FrdI.Prop44Gp
import ABC3.Found.FrdI.Prop33Coa
import ABC3.Found.FrdI.Def24
import ABC3.Found.FrdI.Def27
import ABC3.Found.FrdI.Prop113

/-!
# [FrdI] Definition 4.5 —— 実装できる 2 条(と、待っている 2 条)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.86。

原文 (FrdI p.86):
> Definition 4.5.

## ★★4 条の依存(測定、2026-08-17)

| 条 | 内容 | 依存 | 状況 |
|---|---|---|---|
| (i) | `birationally Frobenius-normalized` / 同型 / `model type` | ★**`𝒞^birat` の pre-Frobenioid 構造**(取得済み) | ★**本ファイルで実装** |
| (ii) | `strictly rational` / `rational` / 同型 | ★★**`Φ^birat`**(`Proposition 4.4, (iii)`) | 待ち |
| (iii) | `rationally standard type` | (ii) ＋ `(𝒞^un-tr)^birat` ＋ `Frobenius-compact` | 待ち |
| (iv) | `Div-slim` | ★**`Φ` の関手性だけ** | ★**本ファイルで実装** |

★★**(iv) はこの命題群の中で唯一「`𝒞` を見ない」条件**である ——
`𝒟` と `Φ` だけで決まる。★`Corollary 4.11` / `Corollary 5.4` / `Corollary 5.7` が
いずれも仮定として引くので、**先に置いておく価値がある**。

## ★(iv) の Lean での書き方

原文 (FrdI p.86):
> Aut(DA →D) →Aut(DA →Mon)

★我々の `MonoidOn` は `Φ.functor : 𝒟ᵒᵖ ⥤ AddCommMon` という**反変**の形なので、
`𝒟_A → Mon` は **`(𝒟_A)ᵒᵖ ⥤ AddCommMon`**、すなわち
`(Over.forget A).op ⋙ Φ.functor` である。
★`Aut` の間の写像は `NatIso.op` と `isoWhiskerRight` の合成で書ける。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★(iv) —— `Div-slim` -/

variable (Φ) in
/-- ★**`Φ` をスライス `𝒟_A` の上へ引いた関手** —— 原文の `𝒟_A → Mon`。

★`Φ` は反変なので、行き先は `(𝒟_A)ᵒᵖ ⥤ AddCommMon` である。 -/
def overPhi (A : D) : (Over A)ᵒᵖ ⥤ AddCommMonCat.{w} :=
  (Over.forget A).op ⋙ Φ.functor

variable (Φ) in
/-- ★★**`Aut(𝒟_A → 𝒟) → Aut(𝒟_A → Mon)`** —— 原文の「induced by composition with
the functor `Φ`」。 -/
def overPhiAut (A : D) (η : Aut (Over.forget A)) : Aut (overPhi Φ A) :=
  Functor.isoWhiskerRight (NatIso.op η) Φ.functor

variable (Φ) in
/-- ★★★**[FrdI] Definition 4.5, (iv)** —— `𝒟` が `Φ` に関して **Div-slim**。

原文 (FrdI p.86):
> Aut(DA →D) →Aut(DA →Mon)

★★**`𝒞` を一切見ない条件**である —— `𝒟` と `Φ` だけで決まる。 -/
def IsDivSlim : Prop := ∀ A : D, Function.Injective (overPhiAut Φ A)

variable (Φ) in
/-- ★★**slim なら Div-slim**(原文の「Thus, if `D` is slim, then it is Div-slim.」)。

原文 (FrdI p.86):
> Div-slim [relative to Φ] if, for every A ∈Ob(D),

★`Aut(𝒟_A → 𝒟)` が自明なら、そこからの写像は何であれ単射である。 -/
theorem isDivSlim_of_isSlim (hslim : IsSlimCat D) : IsDivSlim Φ := by
  intro A η₁ η₂ _
  rw [hslim A η₁, hslim A η₂]

/-! ## ★(i) —— `birationally Frobenius-normalized` -/

variable (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-- ★★★**[FrdI] Definition 4.5, (i)** —— `𝒞` の対象が
**birationally Frobenius-normalized** であること。

原文 (FrdI p.86):
> (i) We shall say that an object of C is birationally Frobenius-normalized if

★★`𝒞^birat` の pre-Frobenioid 構造(`biratPre`)の上で
`Frobenius-normalized` であることをそのまま言う。 -/
def IsBirationallyFrobeniusNormalized (A : C) : Prop :=
  IsFrobeniusNormalized (biratPre P G) ((toBiratCat P G).obj A)

variable (C) in
/-- ★★**[FrdI] Definition 4.5, (i)** —— `𝒞` が
**birationally Frobenius-normalized 型**であること。 -/
def IsOfBirationallyFrobeniusNormalizedType : Prop :=
  ∀ A : C, IsBirationallyFrobeniusNormalized P G A

variable (C) in
/-- ★★★**[FrdI] Definition 4.5, (i)** —— `𝒞` が **model 型**であること。

原文 (FrdI p.86):
> normalized type. If C is of pre-model and birationally Frobenius-

★`pre-model 型`(`Definition 2.7, (iii)`)かつ
`birationally Frobenius-normalized 型`。 -/
def IsOfModelType [IsConnected D] : Prop :=
  IsPreModelType P ∧ IsOfBirationallyFrobeniusNormalizedType C P G

/-! ## ★(ii) —— `strictly rational` / `rational`

原文 (FrdI p.86):
> (ii) Suppose that Φ is perf-factorial; A ∈Ob(C). Then we shall say that A

★★`Φ^birat(A)` は `Prop44Gp.lean` の `phiBiratAt`(`𝒪^×(A^birat)` の `Div^gp` 像)。
★`Supp` は `Definition 2.4, (i), (d)`(`Def24.lean`)のもので、
**`M^pf` の元に対して**定義されているので、`M` の元は `Pf.mk a 1` で持ち上げる。 -/

open scoped NNReal

/-- ★**単系の元の台** —— `Definition 2.4, (i), (d)` の `Supp` を `M` の元に当てたもの。 -/
noncomputable def SuppElt {M : Type w} [AddCommMonoid M]
    (ι : Prime M → Pf M → ℝ≥0) (a : M) : Set (Prime M) :=
  Supp (factorMap ι (Pf.mk a 1))

/-- ★★★**[FrdI] Definition 4.5, (ii)** —— `A` が **strictly rational**。

原文 (FrdI p.86):
> is strictly rational if, for every prime p ∈Prime(Φ(A)), there exists an element

★各素点 `p` について、`Φ^birat(A)` の元 `a − b`(`a, b ∈ Φ(A)`)で
`p ∈ Supp(a)`、`p ∉ Supp(b)` となるものが取れること。 -/
def IsStrictlyRational (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0) (A : C) : Prop :=
  ∀ p : Prime (Φ.val (P.toElem.obj A).base),
    ∃ a b : Φ.val (P.toElem.obj A).base,
      toGp _ a - toGp _ b ∈ phiBiratAt P G (show BiratCat P G from A) ∧
      p ∈ SuppElt (ι _) a ∧ p ∉ SuppElt (ι _) b

/-- ★★★**[FrdI] Definition 4.5, (ii)** —— `A` が **rational**。

原文 (FrdI p.86):
> Definition 2.4, (i), (d)]. We shall say that A is rational if there exists a pull-back

★strictly rational な対象からの **pull-back 射**が入ること。 -/
def IsRational (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0) (A : C) : Prop :=
  ∃ (B : C) (φ : B ⟶ A), IsPullBack P φ ∧ IsStrictlyRational P G ι B

variable (C) in
/-- ★★**[FrdI] Definition 4.5, (ii)** —— `𝒞` が **rational 型**。 -/
def IsOfRationalType (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0) : Prop :=
  ∀ A : C, IsRational P G ι A

variable (C) in
/-- ★★**[FrdI] Definition 4.5, (ii)** —— `𝒞` が **strictly rational 型**。 -/
def IsOfStrictlyRationalType (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0) : Prop :=
  ∀ A : C, IsStrictlyRational P G ι A

/-! ## ★(iii) —— `rationally standard type`

原文 (FrdI p.86):
> (iii) We shall say that C is of rationally standard type if the following conditions

★★**`(𝒞^un-tr)^birat` が書けるようになった**のが要点である ——
`𝒞^un-tr` が Frobenioid であること(`unTr_frobenioid`、`Prop33Coa.lean`)と、
Frobenioid の birationalization(`biratPre`)が揃ったので、
**そのまま合成できる**。 -/

/-- ★`(𝒞^un-tr)^birat` の pre-Frobenioid 構造。 -/
noncomputable def unTrBiratPre (Fc : FrobenioidCore P) (G' : Frobenioid P) :=
  biratPre (unTrPre P Fc) (unTr_frobenioid P Fc G')

/-- ★★★**[FrdI] Definition 4.5, (iii)** —— `𝒞` が **rationally standard 型**。

原文 (FrdI p.86):
> (iii) We shall say that C is of rationally standard type if the following conditions

| 原文 | フィールド |
|---|---|
| (a) birationally Frobenius-normalized 型 | `biratFrobNormalized` |
| (a) rational 型 | `rational` |
| (a) standard 型 | `standard` |
| (b) `(𝒞^un-tr)^birat` が Frobenius-compact 対象を持つ | `unTrBiratCompact` | -/
structure IsOfRationallyStandardType
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0) : Prop where
  /-- **(a)** birationally Frobenius-normalized 型。 -/
  biratFrobNormalized : IsOfBirationallyFrobeniusNormalizedType C P G
  /-- **(a)** rational 型。 -/
  rational : IsOfRationalType C P G ι
  /-- **(a)** standard 型。 -/
  standard : IsOfStandardType D C P G.core
  /-- **(b)** `(𝒞^un-tr)^birat` が Frobenius-compact 対象を持つ。 -/
  unTrBiratCompact : ∃ X : BiratCat (unTrPre P G.core) (unTr_frobenioid P G.core G),
    IsFrobeniusCompact (unTrBiratPre P G.core G) X

/-! ## ★`Definition 4.5` の 4 条がすべて実装された

| 条 | 実装 |
|---|---|
| (i) | `IsBirationallyFrobeniusNormalized` / `IsOfBirationallyFrobeniusNormalizedType` / `IsOfModelType` |
| (ii) | `IsStrictlyRational` / `IsRational` / `IsOfRationalType` / `IsOfStrictlyRationalType` |
| (iii) | `IsOfRationallyStandardType` |
| (iv) | `IsDivSlim`(＋ `isDivSlim_of_isSlim`) |
-/

/-- ★★★★**[FrdI] Definition 4.5** —— 4 条がすべて実装された。 -/
def IsDivSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86, item := "Definition 4.5",
    sectionId := "frdi-def-4-5" }

end ABC3.Found.FrdI
