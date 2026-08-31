/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AlgPoint
import ABC3.Meta.Claim

/-!
# **`X(ℚ̄)` の構成と、その上の高さ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★原文が言っていること

`Definition 1.1 (ii)` の末尾で原文は

    `x ∈ X(ℚ̄) = ⋃_{ℚ̄ ⊇ F, [F:ℚ]<∞} X(F)`,   `ht_M̄(x) ≝ deg_F(x_F^* M̄)`

と書き、続けて「`x_F` は **`x` を与える任意の**射」と言う。
★したがって要るのは

| 段 | 内容 | 場所 |
|---|---|---|
| 1 | 値が定義体の取り方に依らないこと | `HeightInvariant.lean`（`htArith_baseChange_natural`） |
| 2 | 型の言葉での同上 | `AlgPoint.lean`（`htOff_baseChange`） |
| 3 | **`⋃_F X(F)` そのものの構成** | ★★**本ファイル** |

★★段 1・2 は済んでいた。`AlgPoint.lean` は自身の docstring で
「残るのは商そのものの構成(設計)だけである」と記録していた——**それが本ファイルである**。

## ★★★★★★★★★★合同関係を「生成関係」で済ませる

素朴には「`(F, x) ≈ (K, y)` ⟺ 共通の拡大 `L` の上で一致」と定義したくなるが、
★**推移律に合成体（compositum）が要る**——`L₁ ⊇ F, K` と `L₂ ⊇ K, M` から
`L₃ ⊇ L₁, L₂` を作る段である。`IntermediateField` に埋め込み直せば
`IntermediateField.finiteDimensional_sup` で作れるが、
本プロジェクトの `AlgPointOff` は定義体を**抽象な `Type`** として持つので、
埋め込み直しのほうが高くつく。

★★★そこで **`Quot`（`Quotient` ではない）** を使う:

    `Quot r` は **`r` が生成する同値関係**による商である

集合の圏の colimit `colim_F X(F)` は、直和を「底変換で移り合う」関係が
**生成する**同値関係で割ったものそのものなので、これが求める型である。
★★★★合成体は**一度も要らない**。

## ★★★★★★★★降ろす条件

`Quot.lift` が要求するのは「`r` で結ばれる 2 点で値が等しいこと」だけであり、
それはちょうど `htOff_baseChange`（段 2）である。

## ★2 つの型

| 型 | 意味 |
|---|---|
| `AlgPointAnyClass X` | 原文の **`X(ℚ̄)`** そのもの |
| `AlgPointClass X D` | **`U_X(ℚ̄) = X(ℚ̄) \ D`**——因子表示の高さが定義できる範囲 |

★★因子表示では「`x` が `D` を通らない」が要る（`HeightInvariant.lean` の測定）ので、
高さが降りるのは後者の上である。★★★原文が `Proposition 1.6` を
`U_X(ℚ̄)` の上で述べているのと同じ範囲である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★★★★`X(ℚ̄)` —— 因子の条件なし -/

/-- ★★**定義体つきの代数的点**（因子の条件なし）。 -/
structure AlgPointAny (X : Scheme.{0}) where
  /-- 定義体 -/
  fld : Type
  instField : Field fld
  instNF : @NumberField fld instField
  /-- `x_F : Spec 𝓞_F ⟶ X` -/
  map : @specRingOfIntegers fld instField ⟶ X

/-- ★定義体と射から点を作る（依存インスタンスを避けた形）。 -/
def algPointAny (F : Type) [instF : Field F] [instN : NumberField F] {X : Scheme.{0}}
    (xF : specRingOfIntegers F ⟶ X) : AlgPointAny X :=
  ⟨F, instF, instN, xF⟩

/-- ★★★**「同じ点」の生成関係**——定義体を上げると移り合うこと。

★これが**生成する**同値関係で割るので、推移律（＝合成体）は要らない。 -/
inductive AlgPointAnyRel {X : Scheme.{0}} : AlgPointAny X → AlgPointAny X → Prop
  | base (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
      (xF : specRingOfIntegers F ⟶ X) :
      AlgPointAnyRel (algPointAny F xF)
        (algPointAny K (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF))

/-- ★★★★★★★★★★**`X(ℚ̄) = ⋃_{F} X(F)`**——数体についての colimit。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★集合の圏の colimit は「直和を、移り合う関係が**生成する**同値関係で割ったもの」なので、
`Quot` がそのまま求める型を与える。 -/
def AlgPointAnyClass (X : Scheme.{0}) : Type 1 := Quot (@AlgPointAnyRel X)

/-- ★定義体つきの点から、その定める `X(ℚ̄)` の点。 -/
def AlgPointAnyClass.mk {X : Scheme.{0}} (p : AlgPointAny X) : AlgPointAnyClass X := Quot.mk _ p

/-- ★★★**定義体を上げても同じ点である**。 -/
theorem AlgPointAnyClass.sound_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] {X : Scheme.{0}}
    (xF : specRingOfIntegers F ⟶ X) :
    AlgPointAnyClass.mk
        (algPointAny K (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF))
      = AlgPointAnyClass.mk (algPointAny F xF) :=
  (Quot.sound (AlgPointAnyRel.base F K xF)).symm

theorem AlgPointAnyClass.exists_rep {X : Scheme.{0}} (x : AlgPointAnyClass X) :
    ∃ p : AlgPointAny X, AlgPointAnyClass.mk p = x := Quot.exists_rep x

@[elab_as_elim] theorem AlgPointAnyClass.ind {X : Scheme.{0}}
    {motive : AlgPointAnyClass X → Prop} (h : ∀ p : AlgPointAny X, motive (mk p))
    (x : AlgPointAnyClass X) : motive x := Quot.ind h x

/-! ## ★★★★★★`U_X(ℚ̄)` —— 因子を通らない範囲 -/

/-- ★★★**「同じ点」の生成関係**（`D` を通らない点の側）。 -/
inductive AlgPointRel {X : Scheme.{0}} {D : X.IdealSheafData} :
    AlgPointOff X D → AlgPointOff X D → Prop
  | base (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
      (xF : specRingOfIntegers F ⟶ X)
      (hF : pullbackIdeal F D xF ≠ 0)
      (hK : pullbackIdeal K D
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
      AlgPointRel (algPointOff F xF hF)
        (algPointOff K (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hK)

/-- ★★★★★★★★**`U_X(ℚ̄) = X(ℚ̄) \ D`**——数体についての colimit。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x. -/
def AlgPointClass (X : Scheme.{0}) (D : X.IdealSheafData) : Type 1 :=
  Quot (@AlgPointRel X D)

/-- ★定義体つきの点から、その定める `U_X(ℚ̄)` の点。 -/
def AlgPointClass.mk {X : Scheme.{0}} {D : X.IdealSheafData} (p : AlgPointOff X D) :
    AlgPointClass X D := Quot.mk _ p

/-- ★★★**定義体を上げても同じ点である**。 -/
theorem AlgPointClass.sound_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] {X : Scheme.{0}} {D : X.IdealSheafData}
    (xF : specRingOfIntegers F ⟶ X) (hF : pullbackIdeal F D xF ≠ 0)
    (hK : pullbackIdeal K D
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
    AlgPointClass.mk
        (algPointOff K (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hK)
      = AlgPointClass.mk (algPointOff F xF hF) :=
  (Quot.sound (AlgPointRel.base F K xF hF hK)).symm

theorem AlgPointClass.exists_rep {X : Scheme.{0}} {D : X.IdealSheafData}
    (x : AlgPointClass X D) : ∃ p : AlgPointOff X D, AlgPointClass.mk p = x :=
  Quot.exists_rep x

@[elab_as_elim] theorem AlgPointClass.ind {X : Scheme.{0}} {D : X.IdealSheafData}
    {motive : AlgPointClass X D → Prop} (h : ∀ p : AlgPointOff X D, motive (mk p))
    (x : AlgPointClass X D) : motive x := Quot.ind h x

/-! ## ★★`U_X(ℚ̄) ⊆ X(ℚ̄)` -/

/-- ★因子の条件を忘れる。 -/
def AlgPointOff.forget {X : Scheme.{0}} {D : X.IdealSheafData} (p : AlgPointOff X D) :
    AlgPointAny X := ⟨p.fld, p.instField, p.instNF, p.map⟩

/-- ★★**`U_X(ℚ̄) → X(ℚ̄)`**。 -/
def forgetOff {X : Scheme.{0}} {D : X.IdealSheafData} :
    AlgPointClass X D → AlgPointAnyClass X :=
  Quot.lift (fun p => AlgPointAnyClass.mk (AlgPointOff.forget p)) (by
    rintro _ _ ⟨F, K, xF, hF, hK⟩
    exact (AlgPointAnyClass.sound_baseChange F K xF).symm)

@[simp] theorem forgetOff_mk {X : Scheme.{0}} {D : X.IdealSheafData} (p : AlgPointOff X D) :
    forgetOff (AlgPointClass.mk p) = AlgPointAnyClass.mk (AlgPointOff.forget p) := rfl

/-! ## ★★★★★★★★★★高さが `U_X(ℚ̄)` の上に降りること -/

/-- ★★★★★★★★★★**高さは `U_X(ℚ̄)` の上の関数である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★これが原文の「`x_F` は `x` を与える**任意の**射」に当たる段である
——値が代表元に依らないことを `Quot.lift` が型で保証する。

★★降ろす条件は `htOff_baseChange`（`AlgPoint.lean`）そのものであり、
その仮定は原文自身のもの 2 つ（`ι_X` 両立と「`D` を通らない」）だけである。 -/
noncomputable def htClass {X : Scheme.{0}} (D : ArithCartier X)
    (hg : IsConjInvariant D.green) : AlgPointClass X D.divisor → ℝ :=
  Quot.lift (htOff D) (by
    rintro _ _ ⟨F, K, xF, hF, hK⟩
    exact (htOff_baseChange F K D xF hg hF hK).symm)

@[simp] theorem htClass_mk {X : Scheme.{0}} (D : ArithCartier X)
    (hg : IsConjInvariant D.green) (p : AlgPointOff X D.divisor) :
    htClass D hg (AlgPointClass.mk p) = htOff D p := rfl

/-- ★★**定義体と射を与えたときの値**——原文の `ht_M̄(x) = deg_F(x_F^* M̄)` の形。 -/
@[simp] theorem htClass_algPointOff (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (D : ArithCartier X) (hg : IsConjInvariant D.green)
    (xF : specRingOfIntegers F ⟶ X) (h : pullbackIdeal F D.divisor xF ≠ 0) :
    htClass D hg (AlgPointClass.mk (algPointOff F xF h)) = htArith F D xF := rfl

/-! ### ★出典の紐付け(`.src`) -/

def AlgPointAnyClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(X(ℚ̄) = ⋃_F X(F) —— 数体についての colimit)",
    sectionId := "genell-def-1-1-ii" }

def AlgPointClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(U_X(ℚ̄) = X(ℚ̄) \\ D —— 因子表示で高さが定義できる範囲)",
    sectionId := "genell-def-1-1-ii" }

def htClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(高さが U_X(ℚ̄) の上の関数に降りること)",
    sectionId := "genell-def-1-1-ii" }

def htClass.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htOff_baseChange(高さは定義体の取り方に依らない、型の言葉で)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htOff_baseChange") 4,
    .citation "[ABC3]" "htArith_baseChange_natural(技術的仮定なしの底変換不変性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_baseChange_natural") 4,
    .implicitStep
      ("★★★合同関係に**合成体は要らない**——`Quot`(`Quotient` ではない)は " ++
       "与えた関係が**生成する**同値関係による商であり、" ++
       "集合の圏の colimit はまさにそれだからである。" ++
       "素朴に「共通の拡大で一致」と定義すると推移律に合成体が要る") 4,
    .implicitStep
      ("★因子表示では「x が D を通らない」が要るので、高さが降りるのは " ++
       "U_X(ℚ̄) の上である。原文が Proposition 1.6 を U_X(ℚ̄) の上で述べているのと同じ範囲") 4,
    .implicitStep
      ("★★残っている段: 計量表示の側の ht_M̄(x) = deg_F(x_F^* M̄)。" ++
       "引き戻し x_F^* M̄ は §9-742、その乗法性は §9-743 で入ったが、" ++
       "deg_F を APicM (Spec 𝓞_F) の上で読むには " ++
       "ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)([Szp] Prop 1.1、原文自身の引用)が要る") 4 ]

end ABC3.Found.GenEll
