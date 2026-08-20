/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem

/-!
# 因子論の第 4 層 —— 算術因子(`Example 6.3`)(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.112–114。

原文 (FrdI p.112):
> Example 6.3. A Frobenioid of Arithmetic Origin. Let F be a number

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group

## ★★幾何側とは**別の**因子論である

原文は `Example 6.1`(幾何)と `Example 6.3`(算術)を並べる。
底の圏はどちらも `B(G)⁰` だが、因子の単系が違う:

| | 素点の集合 | `ord(O_v^▷)` | `ord(F_v)` |
|---|---|---|---|
| 非アルキメデス | `HeightOneSpectrum (𝓞 L)` | `ℕ` | `ℤ` |
| アルキメデス | `InfinitePlace L` | `ℝ≥0` | `ℝ` |

★★**アルキメデス成分が `ℝ≥0` である**ことが幾何側との唯一の構造的な違いで、
`Found/Divisor/CartierPerfFactorial.lean` の `not_isLambdaMonoprime_real`
(「`M_p ≃+ ℕ` なので `ℝ≥0` にはならない」)が**そのままでは使えない**所である。

## ★★在庫測定(2026-08-20)—— 算術側は**良い**

| 要るもの | mathlib | 判定 |
|---|---|---|
| 有限素点 `NumberField.FinitePlace` | `NumberField/Completion/FinitePlace.lean` | ★**ある** |
| 無限素点 `NumberField.InfinitePlace` | `NumberField/InfinitePlace/Basic.lean` | ★**ある** |
| **積公式** `NumberField.prod_abs_eq_one` | `NumberField/ProductFormula.lean` | ★★**ある** |
| **Dirichlet 単数定理** | `NumberField/Units/DirichletTheorem.lean` | ★★**ある** |
| Tchebotarev 密度定理 | 全体を `Chebotarev\|Tchebotarev` で grep して **0 件** | ★**無い** |

★★**したがって `Example 6.3` を止めている在庫不足は無い。**
`Theorem 6.4, (iv)` だけが Tchebotarev を要求する。

## ★★逸脱(記録)

★原文は `ord(O_v^▷) := O_v^▷/O_v^×` と**商として**定義し、そのうえで
「`≅ ℕ`(非アルキメデス)/ `≅ ℝ≥0`(アルキメデス)」と述べる。
★**本ファイルは同一視した先(`ℕ` / `ℝ≥0`)を定義に採る**。
同一視そのものは鎖 `arith` の節点 `ord-mon`(下の `ordMon_nonarch_equiv` /
`ordMon_arch_equiv`)であり、**まだ証明していない**。
-/

namespace ABC3.Skeleton.Divisor

open ABC3.Meta NumberField

universe u

variable (L : Type u) [Field L] [NumberField L]

/-! ## ★1. 素点と `ord`(鎖 `arith` の `places` / `ord-mon`) -/

/-- ★**`V(L)`** —— 有限素点と無限素点の合併。 -/
abbrev Places : Type u := FinitePlace L ⊕ InfinitePlace L

/-- ★★**`ord(O_v^▷) ≅ ℕ`(非アルキメデス)**。

★原文の定義は `O_v^▷/O_v^×`。同一視は付値が全射 `ℤ` であることから。 -/
theorem ordMon_nonarch_equiv (v : FinitePlace L) :
    Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → ℕ) := by
  sorry

/-- ★★**`ord(O_v^▷) ≅ ℝ≥0`(アルキメデス)**。

★アルキメデス素点では `O_v^▷ = {|x|_v ≤ 1}` で、`ord` は `−log|x|_v`。 -/
theorem ordMon_arch_equiv (v : InfinitePlace L) :
    Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → NNReal) := by
  sorry

/-! ## ★2. 有効算術因子とその群化(鎖 `arith` の `arith-phi`)

★同一視した先を定義に採る(上の「逸脱」を見よ)。 -/

/-- ★★**有効算術因子** `Φ(L) = ⊕_{v ∈ V(L)} ord(O_v^▷)`。

★有限素点成分は `ℕ`、無限素点成分は `ℝ≥0`。直和なので**台は有限**。 -/
abbrev ArithPhi : Type u :=
  (FinitePlace L →₀ ℕ) × (InfinitePlace L →₀ NNReal)

/-- ★**算術因子** `Φ(L)^gp = ⊕_{v ∈ V(L)} ord(F_v)`。 -/
abbrev ArithPhiGp : Type u :=
  (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)

/-! ## ★3. `B(L) = L^× → Φ(L)^gp`(鎖 `arith` の `arith-B`) -/

/-- ★★**`B(L) = L^× → Φ(L)^gp`** —— 各素点での位数を並べる。

原文 (FrdI p.113):
> [given by mapping an element f ∈F × to the images of f in the various factors

★「all but a finite number of which are zero」が `Finsupp` の中身である。 -/
noncomputable def arithDivOfElt (_f : Lˣ) : ArithPhiGp L :=
  sorry

/-- ★**`B(L) → Φ(L)^gp` は群準同型**。 -/
theorem arithDivOfElt_mul (f g : Lˣ) :
    arithDivOfElt L (f * g) = arithDivOfElt L f + arithDivOfElt L g := by
  sorry

/-! ## ★4. 単系としての性質(鎖 `arith` の `arith-perffactorial`)

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

★原文は「one verifies immediately」で畳む。★幾何側(`Example 6.1`)と**同じ形**だが、
アルキメデス成分が `ℝ≥0` なので `M_p ≃+ ℕ` の段は別に要る。 -/

/-- ★★**`Prime(Φ(L)) ≃ V(L)`**。 -/
theorem prime_arithPhi_equiv_places :
    Nonempty (Places L → Places L) := by
  sorry

/-- ★★**`Φ(L)` の元の台はちょうど `V(L)` の有限部分集合**。 -/
theorem support_arithPhi_eq_finite (T : Finset (FinitePlace L)) :
    ∃ a : ArithPhi L, a.1.support = T := by
  sorry

/-! ## ★5. 算術次数(鎖 `arith` の `deg-arith` / `deg-vanish`)

原文 (FrdI p.113):
> degarith

★アルキメデス成分は `−[F_v:ℝ]·log λ`、非アルキメデス成分は `log #(O_v/(λ))`。 -/

/-- ★**算術次数** `deg^arith_L : Φ(L)^gp → ℝ`。 -/
noncomputable def degArith (_x : ArithPhiGp L) : ℝ :=
  sorry

/-- ★`deg^arith` は加法的。 -/
theorem degArith_add (x y : ArithPhiGp L) :
    degArith L (x + y) = degArith L x + degArith L y := by
  sorry

/-- ★★★**`deg^arith` は `B(L)` の像で消える** —— これが**積公式**である。

原文 (FrdI p.114):
> vanishes on the image of B(L).

★★mathlib の `NumberField.prod_abs_eq_one` の対数を取ったものに他ならない。 -/
theorem degArith_arithDivOfElt (f : Lˣ) : degArith L (arithDivOfElt L f) = 0 := by
  sorry

/-! ## ★6. 単元と Dirichlet(鎖 `arith` の `mu-units` / `dirichlet`)

原文 (FrdI p.113):
> to Spec(L) ∈Ob(B(G)0), we have -/

/-- ★★**`O^×(A) = O^▷(A) = μ(L)`** —— 全素点で位数 0 の元は 1 の冪根。 -/
theorem units_eq_roots_of_unity (f : Lˣ) (_h : arithDivOfElt L f = 0) :
    ∃ n : ℕ, 0 < n ∧ (f : L) ^ n = 1 := by
  sorry

/-- ★★★**`δ_A : Pic_Φ(A) ≅ ℝ`** —— `Theorem 6.4, (i)` の最後。

原文 (FrdI p.115):
> support whose image under degarith

★★原文は「an immediate consequence of the well-known **Dirichlet unit theorem**」と書く。
mathlib に在る(`NumberField.Units.dirichletUnitTheorem`)。 -/
theorem degArith_surjective_and_kernel_eq_image :
    Function.Surjective (degArith L) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def Places.src : Source :=
  { paper := "FrdI", pdfPage := 112, item := "Example 6.3 — V(F)(素点の集合)",
    sectionId := "frdi-example-6-3" }

def ordMon_nonarch_equiv.src : Source :=
  { paper := "FrdI", pdfPage := 112, item := "Example 6.3 — ord(O_v^▷) ≅ ℤ≥0(非アルキメデス)",
    sectionId := "frdi-example-6-3" }

def ordMon_nonarch_equiv.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsDedekindDomain.HeightOneSpectrum.valuation(付値)"
      (.inMathlib "IsDedekindDomain.HeightOneSpectrum.valuation") 112,
    .derivation "付値が ℤ へ全射であることから O_v^▷/O_v^× ≅ ℕ" 112 ]

def ordMon_arch_equiv.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — ord(O_v^▷) ≅ ℝ≥0(アルキメデス)",
    sectionId := "frdi-example-6-3" }

def ordMon_arch_equiv.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.InfinitePlace(無限素点)"
      (.inMathlib "NumberField.InfinitePlace") 113,
    .derivation "`−log|x|_v` が `O_v^▷/O_v^× ≅ ℝ≥0` を与える" 113 ]

def ArithPhi.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — Φ(F)(有効算術因子)",
    sectionId := "frdi-example-6-3" }

def ArithPhiGp.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — Φ(F)^gp(算術因子)",
    sectionId := "frdi-example-6-3" }

def arithDivOfElt.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — B(F) = F^× → Φ(F)^gp",
    sectionId := "frdi-example-6-3" }

def arithDivOfElt_mul.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — B(F) → Φ(F)^gp は準同型",
    sectionId := "frdi-example-6-3" }

def arithDivOfElt_mul.needs : List ProofObligation :=
  [ .derivation "各素点成分で付値の乗法性を使う" 113,
    .implicitStep "★台の有限性(「all but a finite number of which are zero」)は原文が畳んでいる" 113 ]

def prime_arithPhi_equiv_places.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — Prime(Φ(L)) ≃ V(L)",
    sectionId := "frdi-example-6-3" }

def prime_arithPhi_equiv_places.needs : List ProofObligation :=
  [ .citation "[ABC3]" "effSubPrimeEquiv(幾何側の同じ形)"
      (.inProject "ABC3" "ABC3.Found.FrdI.effSubPrimeEquiv") 113,
    .derivation "primary 元の台が 1 点であること(幾何側の `mprec_effSub_iff` の算術版)" 113,
    .implicitStep "★原文は「one verifies immediately」で畳む" 113 ]

def support_arithPhi_eq_finite.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — 台はちょうど V(L) の有限部分集合",
    sectionId := "frdi-example-6-3" }

def support_arithPhi_eq_finite.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_effSub_support_eq(幾何側の同じ形)"
      (.inProject "ABC3" "ABC3.Found.FrdI.exists_effSub_support_eq") 113,
    .implicitStep "★原文は「one verifies immediately」で畳む" 113 ]

def degArith.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — deg^arith_L : Φ(L)^gp → ℝ",
    sectionId := "frdi-example-6-3" }

def degArith_add.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — deg^arith は加法的",
    sectionId := "frdi-example-6-3" }

def degArith_add.needs : List ProofObligation :=
  [ .derivation "成分ごとに定義するので加法性は自明" 113 ]

def degArith_arithDivOfElt.src : Source :=
  { paper := "FrdI", pdfPage := 114, item := "Example 6.3 — deg^arith は B(L) の像で消える",
    sectionId := "frdi-example-6-3" }

/-- ★★これは**積公式**そのもので、mathlib に在る。 -/
def degArith_arithDivOfElt.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.prod_abs_eq_one(積公式: ∏_v |x|_v = 1)"
      (.inMathlib "NumberField.prod_abs_eq_one") 114,
    .derivation "積公式の両辺の対数を取ると deg^arith(div(f)) = 0 になる" 114,
    .implicitStep "★原文は「one verifies immediately that deg^arith vanishes on the image of B(L)」で畳む" 114 ]

def units_eq_roots_of_unity.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — O^×(A) = O^▷(A) = μ(L)",
    sectionId := "frdi-example-6-3" }

def units_eq_roots_of_unity.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.Units.torsion(単数群の捻れ部分 = 1 の冪根)"
      (.inMathlib "NumberField.Units.torsion") 113,
    .derivation "全素点で絶対値 1 の代数的数は 1 の冪根(Kronecker)" 113 ]

def degArith_surjective_and_kernel_eq_image.src : Source :=
  { paper := "FrdI", pdfPage := 114, item := "Theorem 6.4, (i) — δ_A : Pic_Φ(A) ≅ ℝ",
    sectionId := "frdi-thm-6-4" }

/-- ★★原文が「well-known」で畳んだ所。mathlib に在る。 -/
def degArith_surjective_and_kernel_eq_image.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.Units.dirichletUnitTheorem(Dirichlet 単数定理)"
      (.inMathlib "NumberField.Units.dirichletUnitTheorem") 114,
    .derivation "Φ^birat(L) ⊗ ℝ の像が deg^arith の核に一致する(単数定理の階数の主張)" 114,
    .implicitStep "★原文は「an immediate consequence of the well-known Dirichlet unit theorem」で畳む" 114 ]

end ABC3.Skeleton.Divisor
