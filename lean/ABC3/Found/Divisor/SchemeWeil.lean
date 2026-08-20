/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.Divisor.HeightOneDVR
import Mathlib.AlgebraicGeometry.FunctionField
import Mathlib.AlgebraicGeometry.Noetherian
import Mathlib.RingTheory.KrullDimension.LocalRing
import Mathlib.RingTheory.DiscreteValuationRing.TFAE

/-!
# スキーム上の Weil 因子 —— 余次元 1 の点と茎の DVR 性(鎖 `weil`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> normal [geometrically integral] variety over a field k; K the function field of V ;

## ★★★見立ての訂正(2026-08-20)

`Skeleton/Divisor/SchemeWeil.lean` は「茎が DVR であることには**アフィン開との
同一視(`affine-compat`)が要る**」と書いていた。★**それは誤りである**。

`IsDiscreteValuationRing.TFAE` の第 4 項
「整閉 ∧ 非零素イデアルがちょうど 1 つ」は**局所 Noether 整域**についての条件で、
茎はそのすべてを満たす:

| 条件 | 根拠 |
|---|---|
| 局所 | `X.presheaf.stalk x` は局所環(mathlib の instance) |
| Noether | `X` が Noether スキームなら茎も Noether(mathlib の instance) |
| 整域 | `X` が整スキームなら茎も整域(mathlib の instance) |
| 整閉 | `IsNormalScheme X` の定義そのもの |
| 非零素が 1 つ | `ringKrullDim = 1` ＋ 局所 —— 非零素はすべて極大 |

★★したがって**アフィン開への降下は要らない**。`affine-compat` が要るのは
「余次元 1 の点 ↔ 高さ 1 の素イデアル」の**対応を使う**場面だけである。

## ★本ファイルで閉じること

| 節点 | 定理 |
|---|---|
| `weil:normal-scheme` | `IsNormalScheme` |
| `weil:codim1-pt` | `IsCodimOnePt` / `PrimeDivisorPt` / `WeilDiv` |
| `weil:stalk-dvr` | `isDiscreteValuationRing_stalk_of_codimOne` |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory

universe u

/-! ## ★1. 環の層 —— 局所 Noether 整域で次元 1 かつ整閉なら DVR -/

/-- ★★★★**局所 Noether 整域が Krull 次元 1 かつ整閉なら DVR**。

★`IsDiscreteValuationRing.TFAE` の第 4 項「整閉 ∧ 非零素がちょうど 1 つ」に落とす。
★非零素が 1 つであることは **`ringKrullDim = 1` ＋ 局所**から出る ——
次元 `≤ 1` なら非零素はすべて極大(`Ring.krullDimLE_one_iff_of_noZeroDivisors`)、
局所なら極大イデアルは 1 つだけ。 -/
theorem isDiscreteValuationRing_of_krullDim_one {R : Type*} [CommRing R] [IsDomain R]
    [IsNoetherianRing R] [IsLocalRing R] [IsIntegrallyClosed R]
    (h : ringKrullDim R = 1) : IsDiscreteValuationRing R := by
  have hnf : ¬ IsField R := fun h' => zero_ne_one ((ringKrullDim_eq_zero_of_isField h') ▸ h)
  have hunique : ∃! P : Ideal R, P ≠ ⊥ ∧ P.IsPrime := by
    refine ⟨IsLocalRing.maximalIdeal R,
      ⟨IsLocalRing.isField_iff_maximalIdeal_eq.not.mp hnf, inferInstance⟩, ?_⟩
    rintro P ⟨hP0, hPp⟩
    have hdimle : Ring.KrullDimLE 1 R := Ring.krullDimLE_iff.mpr (by simp [h])
    have hmax : P.IsMaximal :=
      Ring.krullDimLE_one_iff_of_noZeroDivisors.mp hdimle P hP0 hPp
    exact IsLocalRing.eq_maximalIdeal hmax
  have key : IsIntegrallyClosed R ∧ ∃! P : Ideal R, P ≠ ⊥ ∧ P.IsPrime :=
    ⟨inferInstance, hunique⟩
  exact ((IsDiscreteValuationRing.TFAE R hnf).out 3 0).mp key

/-! ## ★2. 余次元 1 の点(鎖 `weil` の `codim1-pt`) -/

/-- ★**余次元 1 の点** —— 局所環の Krull 次元が 1。

★整スキームでは「余次元 1 の既約閉部分多様体」と生成点で 1 対 1 に対応する。
原文の `prime divisor` はこの意味である。 -/
def IsCodimOnePt (X : Scheme.{u}) (x : X) : Prop :=
  ringKrullDim (X.presheaf.stalk x) = 1

/-- ★**素因子の型** —— 余次元 1 の点。 -/
def PrimeDivisorPt (X : Scheme.{u}) : Type u := {x : X // IsCodimOnePt X x}

instance (X : Scheme.{u}) : CoeOut (PrimeDivisorPt X) X := ⟨Subtype.val⟩

/-- ★**Weil 因子の群** —— 素因子で生成される自由アーベル群。 -/
abbrev WeilDiv (X : Scheme.{u}) : Type u := PrimeDivisorPt X →₀ ℤ

/-! ## ★3. 正規スキーム(鎖 `weil` の `normal-scheme`)

★mathlib に `Scheme.IsNormal` に相当する述語は**無い**(grep 0 件、2026-08-20)。
原文の「normal variety」を写すために自分で置く。 -/

/-- ★**正規スキーム** —— 各茎が整閉。

★原文の `normal variety` はこれに `IsIntegral`(整)と proper を合わせたもの。 -/
def IsNormalScheme (X : Scheme.{u}) : Prop :=
  ∀ x : X, IsIntegrallyClosed (X.presheaf.stalk x)

/-! ## ★4. 余次元 1 の茎は DVR(鎖 `weil` の `stalk-dvr`) -/

/-- ★★★★★**正規 Noether 整スキームの余次元 1 の点の茎は DVR**。

★★**アフィン開への降下は要らない** —— 茎はそれ自身が局所 Noether 整域だから、
`isDiscreteValuationRing_of_krullDim_one` をそのまま当てられる。 -/
theorem isDiscreteValuationRing_stalk_of_codimOne (X : Scheme.{u}) [IsIntegral X]
    [IsLocallyNoetherian X] (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    IsDiscreteValuationRing (X.presheaf.stalk x.1) := by
  haveI := hnorm x.1
  exact isDiscreteValuationRing_of_krullDim_one x.2

/-- ★★茎は関数体の中で分数環になる(`X` が整かつ局所 Noether のとき)。

★これで `ord_x` を関数体の上で定義する道具が揃う。 -/
theorem isFractionRing_stalk (X : Scheme.{u}) [IsIntegral X] [IsLocallyNoetherian X] (x : X) :
    IsFractionRing (X.presheaf.stalk x) X.functionField := inferInstance

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.1` の「正規多様体」の述語。 -/
def IsNormalScheme.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 正規スキームの述語",
    sectionId := "frdi-example-6-1" }

/-- ★★locator —— `Example 6.1` の素因子(余次元 1 の点)と Weil 因子の群。 -/
def PrimeDivisorPt.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 余次元 1 の点と Weil 因子の群",
    sectionId := "frdi-example-6-1" }

/-- ★★★locator —— `Example 6.1` の「余次元 1 の茎は DVR」。 -/
def isDiscreteValuationRing_stalk_of_codimOne.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 正規 Noether 整スキームの余次元 1 の茎は DVR",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
