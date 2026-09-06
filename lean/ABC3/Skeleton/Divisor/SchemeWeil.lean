/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.Divisor.SchemeWeil
import ABC3.Found.Divisor.SchemeDivFinite
import Mathlib.AlgebraicGeometry.FunctionField
import Mathlib.RingTheory.KrullDimension.Basic
import Mathlib.RingTheory.IntegralClosure.IntegrallyClosed
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import Mathlib.AlgebraicGeometry.Noetherian

/-!
# 因子論の第 1 層 —— スキーム上の Weil 因子(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109(`Example 6.1`)。

原文 (FrdI p.109):
> normal [geometrically integral] variety over a field k; K the function field of V ;

原文 (FrdI p.109):
> DK a set of Q-Cartier prime divisors on V . The connected objects of the Galois category

## ★★このファイルの位置づけ

★2026-08-20 まで、`Example 6.1` は依存グラフで**1 節点に畳んだまま**だった。
CLAUDE.md の「進め方」——スケルトンで依存グラフを作ってから葉を形式化する——に従い、
**因子論を 4 層に割った**(`ResearchPaper/frdi-decomposition.json` の鎖
`weil` / `cartier` / `normalize` / `arith`)。本ファイルは第 1 層である。

★**`sorry` は「正しい状態」である** —— `Skeleton/` は statement 専用トラックだからである。

## ★★在庫測定(2026-08-20、`lean/.lake/packages/mathlib`)

| 要るもの | mathlib | 判定 |
|---|---|---|
| `ringKrullDim` | `RingTheory/KrullDimension/Basic.lean` | ★**ある** |
| `IsIntegrallyClosed` | `RingTheory/IntegralClosure/IntegrallyClosed.lean` | ★**ある** |
| `AlgebraicGeometry.IsNoetherian` | `AlgebraicGeometry/Noetherian.lean` | ★**ある**(locally Noetherian ＋ quasi-compact) |
| `Scheme.functionField` | `AlgebraicGeometry/FunctionField.lean` | ★**ある**(整スキームで体) |
| `IsDiscreteValuationRing.TFAE` | `RingTheory/DiscreteValuationRing/TFAE.lean` | ★**ある** |
| `WeilDivisor` / `PrimeDivisor` / `CartierDivisor` | 全体を grep して **0 件** | ★**無い** |

★★**したがって「因子論は一分野ぶん」という見立ては正しくない。**
環の層(`Found/Divisor/HeightOneDVR.lean`)は**すでに閉じている** ——
`isDiscreteValuationRing_atPrime_of_minimal`(正規 Noether 整域の高さ 1 の素での
局所化は DVR)から `ordZ`・台の有限性・`divOfFrac` まで。
★本ファイルが受け持つのは**スキーム層への橋**だけである。
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta

universe u

/-! ## ★1. 余次元 1 の点(鎖 `weil` の `codim1-pt`) -/

/-- ★**余次元 1 の点** —— 局所環の Krull 次元が 1。

★整スキームでは「余次元 1 の既約閉部分多様体」と生成点で 1 対 1 に対応する。
原文の `prime divisor` はこの意味である。 -/
def IsCodimOnePt (X : Scheme.{u}) (x : X) : Prop :=
  ringKrullDim (X.presheaf.stalk x) = 1

/-- ★**素因子の型** —— 余次元 1 の点。 -/
def PrimeDivisorPt (X : Scheme.{u}) : Type u := {x : X // IsCodimOnePt X x}

/-- ★**Weil 因子の群** —— 素因子で生成される自由アーベル群。 -/
abbrev WeilDiv (X : Scheme.{u}) : Type u := PrimeDivisorPt X →₀ ℤ

/-! ## ★2. 正規スキーム(鎖 `weil` の `normal-scheme`)

★mathlib に `Scheme.IsNormal` に相当する述語は**無い**(grep 0 件、2026-08-20)。
原文の「normal variety」を写すために自分で置く。 -/

/-- ★**正規スキーム** —— 各茎が整閉。

★原文の `normal variety` はこれに `IsIntegral`(整)と proper を合わせたもの。 -/
def IsNormalScheme (X : Scheme.{u}) : Prop :=
  ∀ x : X, IsIntegrallyClosed (X.presheaf.stalk x)

/-! ## ★3. 余次元 1 の茎は DVR(鎖 `weil` の `stalk-dvr`)

★★**環の層は閉じている** —— `Found/Divisor/HeightOneDVR.lean` の
`isDiscreteValuationRing_atPrime_of_minimal`。
`IsDiscreteValuationRing.TFAE` の第 4 項「整閉 ∧ 非零素イデアルがちょうど 1 つ」に
`isIntegrallyClosed_of_isLocalization` と `IsLocalization.AtPrime.orderIsoOfPrime` を
合わせると**任意次元**で出る(mathlib は Dedekind の場合しか持っていない)。

## ★★記録の訂正(2026-09-06)

★2026-08-20 の本ファイルは「ここに残るのは『茎 = アフィン開の局所化』という**橋**だけ」と
書いていた。★**それは過大な見積もりだった** —— `Found/Divisor/SchemeWeil.lean` 冒頭の
「見立ての訂正」のとおり、**アフィン開への降下は要らない**。茎はそれ自身が局所 Noether 整域
なので `isDiscreteValuationRing_of_krullDim_one` をそのまま当てられる。
★アフィン開との両立(鎖 `weil` の `affine-compat`)が要るのは、
**素因子を環の言葉で数える**とき(`div(f)` の台の有限性)だけである。 -/

/-- ★★**正規 Noether 整スキームの余次元 1 の点の茎は DVR**。

★★★2026-09-06 に `Found` へ配線した(同名定理 `ABC3.Found.Divisor.*`)。
★**逸脱は無い** —— `Found` 側は `[IsLocallyNoetherian X]` しか要求しないのに対し、
本 statement は `[AlgebraicGeometry.IsNoetherian X]`(= 局所 Noether ＋ 準コンパクト)と
`[IsDomain (X.presheaf.stalk x.1)]` を**余分に**持っている。
余分な仮定は使わないまま、そのまま適用できる。 -/
theorem isDiscreteValuationRing_stalk_of_codimOne
    (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (_hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    [IsDomain (X.presheaf.stalk x.1)] :
    IsDiscreteValuationRing (X.presheaf.stalk x.1) :=
  ABC3.Found.Divisor.isDiscreteValuationRing_stalk_of_codimOne X _hnorm x

/-! ## ★4. `ord` と `div(f)`(鎖 `weil` の `ord-pt` / `div-finite` / `weil-group`)

## ★★★2026-09-06(D8 採用): `Found` の実物へ配線した

★**スキーム層の実物は `Found/` に sorry ゼロで在る**:
`Found/Divisor/SchemeWeilOrd.lean` の `ordPt` / `ordPt_mul`、
`Found/Divisor/SchemeDivFinite.lean` の `finite_ordPt_ne_zero` / `weilDivOfFn` /
`weilDivOfFn_mul`。★環の層(`ordAtPrime` ほか)を経由する必要はもう無い。
本節の 5 宣言はこの 5 本にそのまま対応する。

## ★★逸脱の記録(2026-09-06、D8 採用)

**2026-09-06(D8 採用): 原典の proper normal variety の正規性を復元した。
Skeleton が落としていた仮定を戻すもので、原典からの逸脱ではない。
`[IsLocallyNoetherian X]` は `ordPt` の要求。`[CompactSpace X]` は `IsNoetherian` から出るので
足していない。**

## ★★零写像で埋める道を採らなかった理由

★正規性が無ければ余次元 1 の茎は DVR ではなく `ord` は定義できない。
`ordAtDiv` を零写像で埋めると本節が丸ごと自明になることは
`Check/FrdI/Ex61OrdDegenerate.lean` が固定している。
★`ordAtDiv` を `if IsLocallyNoetherian X ∧ IsNormalScheme X then ordPt … else 0` の
`dite` で書けば statement を変えずに閉じるが、それは
**非正規スキームの枝を零写像で埋める**ことに他ならず、8 例目の退化そのものである。**採らなかった。**
★退化を排除する錨は `Found/Divisor/SchemeWeilOrd.lean` の
`exists_ordPt_eq_one` / `not_forall_ordPt_eq_zero` が供給する。 -/

/-- ★**`ord_x : K(X) → ℤ`** —— 余次元 1 の点での位数。

★★★2026-09-06 に `Found` へ配線した(`ABC3.Found.Divisor.ordPt`)。

★**逸脱の記録**: 2026-09-06(D8 採用): 原典の proper normal variety の正規性を復元した。
Skeleton が落としていた仮定を戻すもので、原典からの逸脱ではない。
`[IsLocallyNoetherian X]` は `ordPt` の要求。`[CompactSpace X]` は `IsNoetherian` から出るので
足していない。 -/
noncomputable def ordAtDiv (X : Scheme.{u}) [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (f : X.functionField) : ℤ :=
  ABC3.Found.Divisor.ordPt X hnorm x f

/-- ★**`ord` は乗法的**。

★★★2026-09-06 に `Found` へ配線した(`ABC3.Found.Divisor.ordPt_mul`)。
仮定の追加は `ordAtDiv` と同じ(正規性の復元 ＋ `ordPt` が要求する局所 Noether 性)。 -/
theorem ordAtDiv_mul (X : Scheme.{u}) [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (f g : X.functionField) (hf : f ≠ 0) (hg : g ≠ 0) :
    ordAtDiv X hnorm x (f * g) = ordAtDiv X hnorm x f + ordAtDiv X hnorm x g :=
  ABC3.Found.Divisor.ordPt_mul hnorm x hf hg

/-- ★★**`div(f)` の台は有限** —— これが `Finsupp` を作る本体。

★★★2026-09-06 に `Found` へ配線した(`ABC3.Found.Divisor.finite_ordPt_ne_zero`)。
足したのは `hnorm` **だけ**である ——
`[CompactSpace X]`(準コンパクト性でアフィン被覆を有限に取る段)は
`[AlgebraicGeometry.IsNoetherian X]` から instance で出るので足していない。 -/
theorem finite_support_ordAtDiv (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (hnorm : IsNormalScheme X) (f : X.functionField) (hf : f ≠ 0) :
    {x : PrimeDivisorPt X | ordAtDiv X hnorm x f ≠ 0}.Finite :=
  ABC3.Found.Divisor.finite_ordPt_ne_zero hnorm hf

/-- ★**有理函数の因子** `div(f) = (f)_0 − (f)_∞`。

★★★2026-09-06 に `Found` へ配線した(`ABC3.Found.Divisor.weilDivOfFn`)。
足したのは `hnorm` だけである。 -/
noncomputable def divOfFn (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (hnorm : IsNormalScheme X) (f : (X.functionField)ˣ) : WeilDiv X :=
  ABC3.Found.Divisor.weilDivOfFn hnorm (Units.ne_zero f)

/-- ★★**`div : K(X)^× → WeilDiv X` は群準同型**。

★★★2026-09-06 に `Found` へ配線した(`ABC3.Found.Divisor.weilDivOfFn_mul`)。
足したのは `hnorm` だけである。 -/
theorem divOfFn_mul (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (hnorm : IsNormalScheme X) (f g : (X.functionField)ˣ) :
    divOfFn X hnorm (f * g) = divOfFn X hnorm f + divOfFn X hnorm g :=
  ABC3.Found.Divisor.weilDivOfFn_mul hnorm (Units.ne_zero f) (Units.ne_zero g)

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def IsCodimOnePt.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — prime divisor(余次元 1 の点)",
    sectionId := "frdi-example-6-1" }

def PrimeDivisorPt.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — prime divisor(素因子の型)",
    sectionId := "frdi-example-6-1" }

def WeilDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — ℤ[D_L](Weil 因子の群)",
    sectionId := "frdi-example-6-1" }

def IsNormalScheme.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — proper normal variety(正規性)",
    sectionId := "frdi-example-6-1" }

def isDiscreteValuationRing_stalk_of_codimOne.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — 余次元 1 の茎は DVR",
    sectionId := "frdi-example-6-1" }

/-- ★★★2026-09-06 に埋まった。**この `.needs` は過大だった** ——
「茎とアフィン開の局所化の同一視」は**要らなかった**ので落とした
(根拠は `Found/Divisor/SchemeWeil.lean` 冒頭の「見立ての訂正」)。
実際に使ったのは `Found` の同名定理 1 本だけである。 -/
def isDiscreteValuationRing_stalk_of_codimOne.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isDiscreteValuationRing_stalk_of_codimOne(茎は局所 Noether 整域で整閉、次元 1)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isDiscreteValuationRing_stalk_of_codimOne") 109,
    .citation "[mathlib]" "IsDiscreteValuationRing.TFAE(整閉 ∧ 非零素がちょうど 1 つ ⟹ DVR)"
      (.inMathlib "IsDiscreteValuationRing.TFAE") 109,
    .implicitStep
      "★原文は「normal variety」と書くだけで、余次元 1 の局所環が DVR であることは前提として使う" 109 ]

def ordAtDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — ord(有理函数の位数)",
    sectionId := "frdi-example-6-1" }

def ordAtDiv_mul.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — ord の乗法性",
    sectionId := "frdi-example-6-1" }

def ordAtDiv_mul.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ordAtPrime_mul(環の層での乗法性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ordAtPrime_mul") 109,
    .derivation "DVR の加法付値が乗法を和に送ること" 109 ]

def finite_support_ordAtDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — div(f) の台の有限性",
    sectionId := "frdi-example-6-1" }

def finite_support_ordAtDiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "finite_support_frac(環の層での台の有限性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.finite_support_frac") 109,
    .derivation "準コンパクト性で有限アフィン被覆を取り、各アフィンで環の層の有限性を当てる" 109,
    .implicitStep
      "★原文は「the divisor of poles」「the divisor of zeroes」と書くだけで、有限性は前提として使う" 109 ]

def divOfFn.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — div(f) = (f)_0 − (f)_∞",
    sectionId := "frdi-example-6-1" }

def divOfFn_mul.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — B(L) → Φ(L)^gp が準同型",
    sectionId := "frdi-example-6-1" }

def divOfFn_mul.needs : List ProofObligation :=
  [ .citation "[ABC3]" "divOfFrac_mul(環の層での準同型性)"
      (.inProject "ABC3" "ABC3.Found.Divisor.divOfFrac_mul") 109,
    .derivation "各素因子成分で ordAtDiv_mul を当てる" 109 ]

end ABC3.Skeleton.Divisor
