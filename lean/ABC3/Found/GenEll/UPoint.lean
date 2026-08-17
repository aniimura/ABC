import ABC3.Found.GenEll.AlgPoint

/-!
# [GenEll] Definition 1.2, (i) —— **`U_X(ℚ̄)` を商として作り、高さを降ろす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★最終段 —— 商を作る

`AlgPoint.lean` は「定義体つきの点」の型と、その上の高さを作った。
★同じ点が**複数の定義体で表される**ので、それを同一視する必要がある。

★★★**本ファイルで商を作り、高さを降ろす。**

    `ht_M̄ : U_X(ℚ̄) → ℝ`

## ★★降りることは証明済み

`AlgPoint.lean` の `htOff_baseChange` が
「定義体を上げても値は変わらない」を与えている。
★★`Quot.lift` はそれをそのまま受ける。

## ★仮定は原文自身のもの 2 つだけ

- `IsConjInvariant D.green` —— 原文の「計量は `ι_X` と両立する」(`Definition 1.1`)
- 各点の `x_F^* D ≠ 0` —— 型 `AlgPointOff` に組み込まれている

★★★**これで `Definition 1.2, (i)` が因子表示の範囲で完成する。**
(ii)(BD-class)は `BDClass.lean` で取得済である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★同一視の関係 -/

/-- ★★**底変換による同一視の生成関係**。

★`x_F` と、それを `K ⊇ F` へ底変換した `x_K` を同一視する。 -/
inductive BaseChangeRel {X : Scheme.{0}} {D : X.IdealSheafData} :
    AlgPointOff X D → AlgPointOff X D → Prop
  | up (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
      (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D xF ≠ 0)
      (hJ' : pullbackIdeal K D
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
      BaseChangeRel (algPointOff F xF hJ)
        (algPointOff K
          (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hJ')

/-! ## ★★★`U_X(ℚ̄)` -/

/-- ★★★**`U_X(ℚ̄) = X(ℚ̄) \ D`** —— 定義体の取り方で同一視した点の型。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★原文の `U_X(ℚ̄) ⊆ X(ℚ̄)` に当たる。
★★因子表示なので「`D` を通らない」が型に組み込まれている
(`Proposition 1.6` が述べられる範囲と同じ)。 -/
def UPoint (X : Scheme.{0}) (D : X.IdealSheafData) : Type 1 :=
  Quot (@BaseChangeRel X D)

/-- ★点から商の元を作る。 -/
def UPoint.mk {X : Scheme.{0}} {D : X.IdealSheafData} (p : AlgPointOff X D) :
    UPoint X D :=
  Quot.mk _ p

/-! ## ★★★高さを商の上に降ろす -/

/-- ★★★**[GenEll] Definition 1.2, (i)** —— `ht_M̄ : U_X(ℚ̄) → ℝ`。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これが原文の高さ関数そのものである**(因子表示の範囲で)。

★降りることは `htOff_baseChange`(`AlgPoint.lean`)が保証している。
★★仮定は原文自身のもの——`IsConjInvariant D.green`
(原文の「計量は `ι_X` と両立する」)だけである。 -/
noncomputable def htU {X : Scheme.{0}} (D : ArithCartier X)
    (hg : IsConjInvariant D.green) (p : UPoint X D.divisor) : ℝ :=
  Quot.lift (htOff D) (by
    intro a b h
    cases h
    exact (htOff_baseChange _ _ D _ hg _ _).symm) p

@[simp] theorem htU_mk {X : Scheme.{0}} (D : ArithCartier X)
    (hg : IsConjInvariant D.green) (p : AlgPointOff X D.divisor) :
    htU D hg (UPoint.mk p) = htOff D p := rfl

/-- ★★**同じ点の別の定義体による表示は、商では同じ元である**。 -/
theorem UPoint.mk_baseChange {X : Scheme.{0}} {D : X.IdealSheafData}
    (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
    (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D xF ≠ 0)
    (hJ' : pullbackIdeal K D
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
    UPoint.mk (algPointOff F xF hJ)
      = UPoint.mk (algPointOff K
          (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hJ') :=
  Quot.sound (BaseChangeRel.up F K xF hJ hJ')

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `Definition 1.2, (i)` は `X(ℚ̄)` **全体**で定めており、
因子表示で全域化するには移動補題(`D` を動かして `x` を避ける)が要る。
★★**原文が実際に使う `U_X(ℚ̄)` の範囲では完成している。** -/

def htU.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(U_X(ℚ̄) 上の高さ関数——全域化の移動補題は含まない)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
