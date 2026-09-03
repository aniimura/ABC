/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem

/-!
# ArithDivisor —— `[FrdI] Example 6.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
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

end ABC3.Skeleton.Divisor
