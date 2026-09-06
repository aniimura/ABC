/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.Divisor.ArithOrd
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem

/-!
# ArithDivisor —— `[FrdI] Example 6.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。

## ★★★2026-09-06(判断 D10、第 1 波): `Found` の実物へ配線した

★**statement は 1 字も変えていない。** 変えたのは `sorry` だった本体を入れたことだけである。

`arithDivOfElt` と `degArith` は **`sorry` を本体に持つ `def`** だった。
`sorry` 本体の `def` は**定数関数に展開される**ので、下流(`Theorem64.lean`)の主張は
証明も反証もできる状態にあった(判断 D21 の一般則)。
★実際 `Function.Surjective (degArith L)` からは `(0 : ℝ) = 1` が出る、と測定されていた。
**本体を入れることでその危険は消える。**

★★**型が既に一致していた** —— `Found/Divisor/ArithOrd.lean` の
`arithDivGroupEquiv : ((FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)) ≃+ arithDivGroup L`
の定義域は本ファイルの `ArithPhiGp L` と**リテラルに同一**(宇宙も一致)なので、
型の詰め替えは要らなかった。`ArithOrd.lean` は `Found` と mathlib しか import しないので
循環もしない。

★★**本体を入れる方向に偽は無い**(5 件確認済み) ——
`arithDivOfElt_mul` / `degArith_add` / `degArith_arithDivOfElt` / `units_eq_roots_of_unity` /
`Theorem64.lean::degArith_surjective_and_kernel_eq_image` は**全部真になる**。
★これは判断 D16 の「核の条を足すと偽になる」とは**別物**である
(あちらは statement を強める向き、こちらは本体を入れる向き)。

## ★★残っている 4 つの `sorry` は「零写像で埋めない」ために残している

`ordMon_nonarch_equiv` / `ordMon_arch_equiv` / `prime_arithPhi_equiv_places` /
`support_arithPhi_eq_finite` は**第 2 波(statement 変更)待ち**である。
★前 3 つは現在の型のままなら `fun _ => 0` や `id` で閉じるが、それは
`Check/FrdI/Ex63DegDegenerate.lean` が固定した 9 例目の退化そのものなので**採らない**。
★★とくに `ordMon_arch_equiv` は**素朴に強めると偽になる** ——
`L` は可算なので `{|x|_v : x ∈ L^×}` は `ℝ_{>0}` の可算部分群にすぎず、
`≅ ℝ≥0` は成り立たない(`Check/FrdI/Ex63DegDegenerate.lean` の
`ordMon_arch_equiv_strengthened_false` が証明済み)。原文の `O_v^▷` は
**完備化 `F_v` の中**の対象であり、修復は「定義域を `v.Completion` へ移す」ことと
**一体でしか動かせない**。
-/

namespace ABC3.Skeleton.Divisor

open ABC3.Meta NumberField
universe u
variable (L : Type u) [Field L] [NumberField L]

/-! ## ★1. 素点と `ord`(鎖 `arith` の `places` / `ord-mon`) -/

/-- ★**`V(L)`** —— 有限素点と無限素点の合併。 -/
abbrev Places : Type u := FinitePlace L ⊕ InfinitePlace L

/-- ★★**`ord(O_v^▷) ≅ ℕ`(非アルキメデス)**。

★原文の定義は `O_v^▷/O_v^×`。同一視は付値が全射 `ℤ` であることから。

## ★★第 2 波(statement 変更)待ち —— 2026-09-06、判断 D10

★現在の型 `Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → ℕ)` は `fun _ => 0` で閉じる。
**それは 9 例目の退化そのもの**(`Check/FrdI/Ex63DegDegenerate.lean`)なので採らない。
★実物は `Found/Divisor/ArithDivisor.lean` の `ordFin` と
`Found/Divisor/ArithOrd.lean` の `arithOfParts_mem_arithEff_iff` に sorry 0 で在る。
配線には statement を「単系 `O_v^▷/O_v^×` の `ℕ` との同型」へ言い直すことが要る。 -/
theorem ordMon_nonarch_equiv (v : FinitePlace L) :
    Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → ℕ) := by
  sorry

/-- ★★**`ord(O_v^▷) ≅ ℝ≥0`(アルキメデス)**。

★アルキメデス素点では `O_v^▷ = {|x|_v ≤ 1}` で、`ord` は `−log|x|_v`。

## ★★★第 2 波(statement 変更)待ち —— 2026-09-06、判断 D10

★★**素朴な修復は偽を作る。** `L` の上のまま `≅ ℝ≥0` に強めると
`L` は可算なので `{|x|_v : x ∈ L^×}` は `ℝ_{>0}` の可算部分群にすぎず**偽**になる
(`Check/FrdI/Ex63DegDegenerate.lean::ordMon_arch_equiv_strengthened_false`)。
原文の `O_v^▷` は**完備化 `F_v` の中**の対象である。
★`Found/Divisor/ArithOrd.lean` はここを正しく扱っている ——
`ordArch (v : InfinitePlace L) (x : v.Completion) : ℝ` と `ordArch_surjective`。
**修復は「定義域を `v.Completion` へ移す」ことと一体でしか動かせない。** -/
theorem ordMon_arch_equiv (v : InfinitePlace L) :
    Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → NNReal) := by
  sorry

/-! ## ★2. 有効算術因子とその群化(鎖 `arith` の `arith-phi`)

★同一視した先を定義に採る(上の「逸脱」を見よ)。 -/

/-- ★★**有効算術因子** `Φ(L) = ⊕_{v ∈ V(L)} ord(O_v^▷)`。

★有限素点成分は `ℕ`、無限素点成分は `ℝ≥0`。直和なので**台は有限**。 -/
abbrev ArithPhi : Type u :=
  (FinitePlace L →₀ ℕ) × (InfinitePlace L →₀ NNReal)

/-- ★**算術因子** `Φ(L)^gp = ⊕_{v ∈ V(L)} ord(F_v)`。

★★`Found/Divisor/ArithOrd.lean` の `arithDivGroupEquiv` の定義域と**リテラルに同一**である
(2026-09-06、D10 の測定)。 -/
abbrev ArithPhiGp : Type u :=
  (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)

/-! ## ★3. `B(L) = L^× → Φ(L)^gp`(鎖 `arith` の `arith-B`) -/

/-- ★★**`B(L) = L^× → Φ(L)^gp`** —— 各素点での位数を並べる。

原文 (FrdI p.113):
> [given by mapping an element f ∈F × to the images of f in the various factors

★「all but a finite number of which are zero」が `Finsupp` の中身である。

★★★2026-09-06(D10、第 1 波)に `Found` へ配線した。**statement は変えていない。**
本体は `arithDiv f`(各素点成分が `−log|f|_v`)を `Φ(L)^gp = ⊕_v ord(F_v)` の
成分表示へ戻したものである(`arithDivGroupEquiv` の逆)。 -/
noncomputable def arithDivOfElt (f : Lˣ) : ArithPhiGp L :=
  (ABC3.Found.Divisor.arithDivGroupEquiv (L := L)).symm
    ⟨ABC3.Found.Divisor.arithDiv f, ABC3.Found.Divisor.arithDiv_mem_arithDivGroup f⟩

/-- ★**`B(L) → Φ(L)^gp` は群準同型**。

★★★2026-09-06(D10)に埋まった。`Found` の `arithDiv_mul` に
`arithDivGroupEquiv.symm` の加法性を合わせるだけである。 -/
theorem arithDivOfElt_mul (f g : Lˣ) :
    arithDivOfElt L (f * g) = arithDivOfElt L f + arithDivOfElt L g := by
  unfold arithDivOfElt
  rw [← map_add]
  congr 1
  exact Subtype.ext (ABC3.Found.Divisor.arithDiv_mul f g)

/-- ★成分表示へ戻す前の姿(`Found` の `arithDiv`)へ戻る。 -/
theorem arithOfParts_arithDivOfElt (f : Lˣ) :
    ABC3.Found.Divisor.arithOfParts (L := L) (arithDivOfElt L f)
      = ABC3.Found.Divisor.arithDiv f :=
  congrArg Subtype.val ((ABC3.Found.Divisor.arithDivGroupEquiv (L := L)).apply_symm_apply
    ⟨ABC3.Found.Divisor.arithDiv f, ABC3.Found.Divisor.arithDiv_mem_arithDivGroup f⟩)

/-! ## ★4. 単系としての性質(鎖 `arith` の `arith-perffactorial`)

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

★原文は「one verifies immediately」で畳む。★幾何側(`Example 6.1`)と**同じ形**だが、
アルキメデス成分が `ℝ≥0` なので `M_p ≃+ ℕ` の段は別に要る。 -/

/-- ★★**`Prime(Φ(L)) ≃ V(L)`**。

## ★★第 2 波(statement 変更)待ち —— 2026-09-06、判断 D10

★現在の型 `Nonempty (Places L → Places L)` は `id` で閉じるが、
**「Prime」が型のどこにも現れない**ので中身がない(9 例目の退化)。
★実物は `Found/Divisor/ArithDivisor*.lean` の `arithPrimeEquiv` に sorry 0 で在る。 -/
theorem prime_arithPhi_equiv_places :
    Nonempty (Places L → Places L) := by
  sorry

/-- ★★**`Φ(L)` の元の台はちょうど `V(L)` の有限部分集合**。

## ★第 2 波(statement 変更)待ち —— 2026-09-06、判断 D10

★現在の型は有限素点側にしか量化していない。実物
(`Found` の `exists_arithEff_support_eq`)は `ArithPlace L` 全体(無限素点込み)で
量化しているので、statement の言い直しと一体で配線する。 -/
theorem support_arithPhi_eq_finite (T : Finset (FinitePlace L)) :
    ∃ a : ArithPhi L, a.1.support = T := by
  sorry

/-! ## ★5. 算術次数(鎖 `arith` の `deg-arith` / `deg-vanish`)

原文 (FrdI p.113):
> degarith

★アルキメデス成分は `−[F_v:ℝ]·log λ`、非アルキメデス成分は `log #(O_v/(λ))`。 -/

/-- ★**算術次数** `deg^arith_L : Φ(L)^gp → ℝ`。

★★★2026-09-06(D10、第 1 波)に `Found` へ配線した。**statement は変えていない。**
本体は `arithOfParts`(成分表示を `⊕_v ord(F_v)` へ戻す)に `Found` の
`arithDegree`(成分の総和)を合成したものである。
★非アルキメデス成分に `log(N v)` の重みが入るのは `arithOfParts` の中(`scaleFin`)であり、
原文の「非アルキメデス成分は `log #(O_v/(λ))`」がそこに対応する。 -/
noncomputable def degArith (x : ArithPhiGp L) : ℝ :=
  ABC3.Found.Divisor.arithDegree (ABC3.Found.Divisor.arithOfParts (L := L) x)

/-- ★`deg^arith` は加法的。

★★★2026-09-06(D10)に埋まった。`arithOfParts` と `arithDegree` が
どちらも `AddMonoidHom` なので `map_add` を 2 回当てるだけである。 -/
theorem degArith_add (x y : ArithPhiGp L) :
    degArith L (x + y) = degArith L x + degArith L y := by
  simp only [degArith, map_add]

/-- ★★★**`deg^arith` は `B(L)` の像で消える** —— これが**積公式**である。

原文 (FrdI p.114):
> vanishes on the image of B(L).

★★mathlib の `NumberField.prod_abs_eq_one` の対数を取ったものに他ならない。

★★★2026-09-06(D10)に埋まった。`Found` の `arithDegree_arithDiv`(積公式、sorry 0)へ
`arithOfParts_arithDivOfElt` で乗り換えるだけである。 -/
theorem degArith_arithDivOfElt (f : Lˣ) : degArith L (arithDivOfElt L f) = 0 := by
  rw [degArith, arithOfParts_arithDivOfElt, ABC3.Found.Divisor.arithDegree_arithDiv]

/-- ★無限素点 1 つに値 `r` を置いた算術因子は `⊕_v ord(F_v)` の中で `single` になる。

★`Theorem64.lean` の全射性の錨。 -/
theorem arithOfParts_single_inr (v : InfinitePlace L) (r : ℝ) :
    ABC3.Found.Divisor.arithOfParts (L := L) (0, Finsupp.single v r)
      = Finsupp.single (Sum.inr v) r := by
  classical
  ext x
  cases x with
  | inl w =>
      rw [ABC3.Found.Divisor.arithOfParts_inl]
      simp
  | inr w =>
      rw [ABC3.Found.Divisor.arithOfParts_inr]
      show (Finsupp.single v r) w = (Finsupp.single (Sum.inr v) r) (Sum.inr w)
      simp [Finsupp.single_apply]

/-- ★★**錨** —— `deg^arith` は零写像ではない。無限素点 `v` に `r` を置けば値は `r`。

★`degArith` が `sorry` 本体の `def` だったころは
`degArith L x = degArith L y` が `rfl` で通っていた(判断 D21)。
**この 1 本がその退化を排除する。** -/
theorem degArith_single_infinite (v : InfinitePlace L) (r : ℝ) :
    degArith L (0, Finsupp.single v r) = r := by
  rw [degArith, arithOfParts_single_inr]
  simp [ABC3.Found.Divisor.arithDegree]

/-! ## ★6. 単元と Dirichlet(鎖 `arith` の `mu-units` / `dirichlet`)

原文 (FrdI p.113):
> to Spec(L) ∈Ob(B(G)0), we have -/

/-- ★★**`O^×(A) = O^▷(A) = μ(L)`** —— 全素点で位数 0 の元は 1 の冪根。

★★★2026-09-06(D10、第 2 波(段階 2))に `Found` へ配線した。
**Kronecker の定理は新規に書く必要が無かった** —— `Found/Divisor/ArithDivisor.lean` の
`exists_pow_eq_one_of_arithDiv_eq_zero`(sorry 0)が実物である。
中身は「非アルキメデス側で `x ∈ (𝓞 L)ˣ` を出し
(`IsDedekindDomain.HeightOneSpectrum.mem_integers_of_valuation_le_one`)、
アルキメデス側と `NumberField.Units.mem_torsion` で有限位数にする」。 -/
theorem units_eq_roots_of_unity (f : Lˣ) (h : arithDivOfElt L f = 0) :
    ∃ n : ℕ, 0 < n ∧ (f : L) ^ n = 1 := by
  refine ABC3.Found.Divisor.exists_pow_eq_one_of_arithDiv_eq_zero f ?_
  have h1 := congrArg (ABC3.Found.Divisor.arithOfParts (L := L)) h
  rw [arithOfParts_arithDivOfElt, map_zero] at h1
  exact h1

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

/-- ★★2026-09-06(D10)に新設。`arithDivOfElt` には `.needs` が無かった
(★**退化のあった 2 つ**にだけ `.needs` が無いという偏りだった)。 -/
def arithDivOfElt.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithDiv(Found 側の本体。各素点成分が −log|f|_v、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDiv") 113,
    .citation "[ABC3]" "arithDiv_mem_arithDivGroup(非アルキメデス成分が log(N v) の整数倍であること)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDiv_mem_arithDivGroup") 113,
    .citation "[ABC3]" "arithDivGroupEquiv(Φ(F)^gp = ⊕_v ord(F_v) の成分表示、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDivGroupEquiv") 113,
    .derivation "各素点で f の位数を取り、非アルキメデス成分は log(N v) で割って整数に戻す" 113 ]

def arithOfParts_arithDivOfElt.src : Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 成分表示から ⊕_v ord(F_v) へ戻す(Skeleton ↔ Found の橋)",
    sectionId := "frdi-example-6-3" }

def arithOfParts_arithDivOfElt.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithDivGroupEquiv(Found 側。成分表示との同型、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDivGroupEquiv") 113,
    .derivation "`AddEquiv.apply_symm_apply` に `Subtype.val` を当てるだけ" 113 ]

def arithDivOfElt_mul.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — B(F) → Φ(F)^gp は準同型",
    sectionId := "frdi-example-6-3" }

def arithDivOfElt_mul.needs : List ProofObligation :=
  [ .derivation "各素点成分で付値の乗法性を使う" 113,
    .citation "[ABC3]" "arithDiv_mul(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDiv_mul") 113,
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

/-- ★★2026-09-06(D10)に新設。`degArith` にも `.needs` が無かった。
★原文は 2 場合の式(アルキメデスは `−[F_v:ℝ]·log λ`、非アルキメデスは `log #(O_v/(λ))`)を
明示しているので、下界に写す。 -/
def degArith.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithDegree(Found 側の本体。成分の総和、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDegree") 113,
    .citation "[ABC3]" "arithOfParts(成分表示を ⊕_v ord(F_v) へ戻す。非アルキメデス成分に log(N v) の重みが入る)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithOfParts") 113,
    .derivation "アルキメデス成分は −[F_v:ℝ]·log λ、非アルキメデス成分は log #(O_v/(λ)) を足し合わせる" 113 ]

def arithOfParts_single_inr.src : Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 無限素点 1 つに値を置いた算術因子(⊕_v ord(F_v) の中では single)",
    sectionId := "frdi-example-6-3" }

def arithOfParts_single_inr.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithOfParts_inr(Found 側。無限素点成分はそのまま入る)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithOfParts_inr") 113,
    .citation "[ABC3]" "arithOfParts_inl(Found 側。有限素点成分は log(N v) の重みつき)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithOfParts_inl") 113,
    .derivation "有限素点成分は 0、無限素点成分は `Finsupp.single` の成分ごとの比較" 113 ]

def degArith_single_infinite.src : Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — deg^arith の錨(無限素点 1 つに r を置けば deg = r)",
    sectionId := "frdi-example-6-3" }

/-- ★★**退化を排除する錨**である。`degArith` が `sorry` 本体の `def` だったころは
`degArith L x = degArith L y` が `rfl` で通っていた(判断 D21)。 -/
def degArith_single_infinite.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithDegree(Found 側。成分の総和、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDegree") 113,
    .derivation "`arithDegree (Finsupp.single a r) = r`(`Finsupp.liftAddHom` の single での値)" 113 ]

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
    .citation "[ABC3]" "arithDegree_arithDiv(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDegree_arithDiv") 114,
    .derivation "積公式の両辺の対数を取ると deg^arith(div(f)) = 0 になる" 114,
    .implicitStep "★原文は「one verifies immediately that deg^arith vanishes on the image of B(L)」で畳む" 114 ]

def units_eq_roots_of_unity.src : Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — O^×(A) = O^▷(A) = μ(L)",
    sectionId := "frdi-example-6-3" }

def units_eq_roots_of_unity.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.Units.torsion(単数群の捻れ部分 = 1 の冪根)"
      (.inMathlib "NumberField.Units.torsion") 113,
    .citation "[ABC3]" "exists_pow_eq_one_of_arithDiv_eq_zero(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_pow_eq_one_of_arithDiv_eq_zero") 113,
    .derivation "全素点で絶対値 1 の代数的数は 1 の冪根(Kronecker)" 113 ]

end ABC3.Skeleton.Divisor
