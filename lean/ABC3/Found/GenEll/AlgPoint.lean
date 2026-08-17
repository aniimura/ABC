import ABC3.Found.GenEll.HeightInvariant

/-!
# [GenEll] Definition 1.2, (i) —— **`U_X(ℚ̄)` の型と、その上の高さ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★型を作る

原文の高さは **`ht_M̄ : X(ℚ̄) → ℝ`** である。
★`X(ℚ̄)` の点は「定義体 `F` と `x_F : Spec 𝓞_F → X` の組」で表す。

★★因子表示では「`x` が `D` を通らない」が要る(`HeightInvariant.lean` の測定)ので、
本ファイルが作るのは **`U_X(ℚ̄) = X(ℚ̄) \ D`** に当たる型である。

★★★**原文が `Proposition 1.6` を `U_X(ℚ̄)` の上で述べているのと同じ範囲**である。

## ★★底変換で値が変わらないことは証明済み

`HeightInvariant.lean` の `htArith_baseChange_natural` が
**技術的仮定なしで**それを与える。
★したがって本ファイルの `htOff` は、**同じ点の別の定義体による表示**でも
同じ値を返す(`htOff_baseChange`)。

## ★残る設計

「同じ点」を商として実装する(数体についての colimit)。
★★**降りることは本ファイルの `htOff_baseChange` が保証している**——
残るのは商の作り方だけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★型 -/

/-- ★★**`D` を通らない代数的点**(定義体つき)。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★原文の `x ∈ U_X(F) ⊆ U_X(ℚ̄)` に当たる。
★★`off` が「`x` が `D` を通らない」である。 -/
structure AlgPointOff (X : Scheme.{0}) (D : X.IdealSheafData) where
  /-- 定義体 -/
  fld : Type
  instField : Field fld
  instNF : @NumberField fld instField
  /-- `x_F : Spec 𝓞_F ⟶ X` -/
  map : @specRingOfIntegers fld instField ⟶ X
  /-- `x` は `D` を通らない -/
  off : @pullbackIdeal fld instField X D map ≠ 0

/-! ## ★★★その上の高さ -/

/-- ★★★**`U_X(ℚ̄)` 上の高さ**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★これが `Definition 1.2, (i)` の関数そのものである(因子表示の範囲で)。 -/
noncomputable def htOff {X : Scheme.{0}} (D : ArithCartier X)
    (p : AlgPointOff X D.divisor) : ℝ :=
  letI := p.instField
  letI := p.instNF
  htArith p.fld D p.map

/-- ★**定義体と射から点を作る**(依存インスタンスを避けた形)。 -/
def algPointOff (F : Type) [instF : Field F] [instN : NumberField F] {X : Scheme.{0}}
    {D : X.IdealSheafData} (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D xF ≠ 0) : AlgPointOff X D where
  fld := F
  instField := instF
  instNF := instN
  map := xF
  off := h

@[simp] theorem htOff_algPointOff (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} (D : ArithCartier X) (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D.divisor xF ≠ 0) :
    htOff D (algPointOff F xF h) = htArith F D xF := rfl

/-! ## ★★★定義体を上げても値は変わらない -/

/-- ★★★**高さは定義体の取り方に依らない**(型の言葉で)。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これで `ht` を商の上に降ろす準備が整った**——
残るのは商そのものの構成(設計)だけである。

★仮定は原文自身のもの 2 つ(`ι_X` 両立と「`D` を通らない」)だけである。 -/
theorem htOff_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] {X : Scheme.{0}}
    (D : ArithCartier X) (xF : specRingOfIntegers F ⟶ X)
    (hg : IsConjInvariant D.green)
    (hJ : pullbackIdeal F D.divisor xF ≠ 0)
    (hJ' : pullbackIdeal K D.divisor
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) ≠ 0) :
    htOff D (algPointOff K
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) hJ')
      = htOff D (algPointOff F xF hJ) := by
  rw [htOff_algPointOff, htOff_algPointOff]
  exact htArith_baseChange_natural F K D xF hg hJ

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2` 全体には
「同じ点」を商として実装する段(数体についての colimit)が残っている。 -/

def htOff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(U_X(ℚ̄) 上の高さ——商の構成は含まない)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
