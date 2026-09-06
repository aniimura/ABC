/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem
import Mathlib.NumberTheory.Real.Irrational
import Mathlib.FieldTheory.KummerPolynomial
import Mathlib.Analysis.Real.Cardinality

/-!
# [FrdI] Theorem 6.4, (i) 末尾——名前が約束した半分が型に無く、素直に足すと偽になる

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114–115
(`Theorem 6.4, (i)` の最後の一文と、その証明)。

原文 (FrdI p.114):
> determines an isomorphism of groups

原文 (FrdI p.115):
> to show that the surjection

原文 (FrdI p.115):
> it suffices to verify that the image of

原文 (FrdI p.115):
> with finite support whose image under

原文 (FrdI p.115):
> But this is an immediate consequence of the well-known Dirichlet unit theorem of classical number theory

## ★★検査の対象

`Skeleton/Divisor/ArithDivisor/Theorem64.lean` の宣言(2026-09-06 時点):

```
theorem degArith_surjective_and_kernel_eq_image :
    Function.Surjective (degArith L)
```

★名前は「全射 **かつ** 核 = 像」を約束しているのに、型は全射性しか言っていない。
核の等式は needs 欄の derivation に日本語で書かれているだけで、型のどこにも無い。

## ★★★退化 その 1——落ちたのは、原文が証明していた側である

原文 p.115 の証明は「δ_A は全射である」を**既知として**始まり、
「同型であることを示すには、Φ^birat(L) ⊗_Z R = (L^×) ⊗_Z R の像が
(Φ^rlf-factor)^gp(L) の元のうち有限台をもち deg^arith が 0 のもの全体に一致することを
確かめればよい」と述べ、そこで **Dirichlet 単数定理**を引く。
★つまり原文が畳んだ「well-known」は**核の等式のほう**に掛かっており、
全射性は原文にとって出発点である。Skeleton は易しい半分だけを型に残した。

`surjectivity_alone_is_the_easy_half` は、`div ≡ 0`・`deg = (無限素点 v₀ での係数)` という組が

* `Example 6.3` の課す条件(準同型性・加法性・積公式 `deg (div f) = 0`)
* `Theorem 6.4` の全射性

を**すべて**満たすことを示す。素点も対数も積公式の中身も 1 ミリも使っていない。
そして同じ組は「核 = 像」を**全く**満たさない。

## ★★★★退化 その 2——名前どおりに核の条を足すと**偽**になる

Skeleton の `ArithPhiGp L = (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)` は
**実現化していない** Φ^gp である。この型のまま「`deg` の核 = `L^×` の像」を足すと
偽になる(`arithPic_ker_not_principal_subgroup`)。理屈は 2 行である:

* `deg` がアルキメデス成分で ℝ-斉次なら(`Found` の `arithDegreeLin` はそうである)、
  無限素点が 2 つあるとき、核は**アルキメデス方向の直線を丸ごと含む**——連続体濃度。
* `L^×` は可算(数体は ℚ 上有限次元)なので、像は可算。

したがって「核 ⊆ 像」はありえない。★逆の包含「像 ⊆ 核」は**積公式そのもの**で、
`range_div_subset_ker` として無条件に真である。偽になるのは、
「核 = 像」が付け加える側**だけ**である。

★具体例は `L = ℚ(√2)`(`arithPic_ker_not_principal_subgroup_qsqrt2`)。
無限素点が 2 つ以上あること(実素点が 2 個)を `one_lt_card_infinitePlace_Qsqrt2` で示した。

## ★★★★★なぜ実現化すると直るのか

原文は **2 か所**で実現化している:

1. 台となる群が (Φ^rlf-factor)^gp(L)——「有限台の**実**係数関数」である。
2. 割る側が Φ^birat(L) ⊗_Z R = (L^×) ⊗_Z R——主因子の **ℝ-span** である。

★Skeleton の型では 2 のほうを**書くことすらできない**。`Submodule.span ℝ` は
台となる群が ℝ-加群であることを要求するが、`FinitePlace L →₀ ℤ` は ℝ-加群ではないからである。
`Found/Divisor/ArithPicR.lean` が `ArithPlace L →₀ ℝ`(全成分が実係数)に移り、
`principalSpan L = Submodule.span ℝ (Set.range arithDiv)` を取っているのは、
まさにこの 2 点を同時に直すためである。そこでは `principalSpan_eq_ker`
(Dirichlet 単数定理 + 類数有限、sorry 0)が成り立ち、
`Found/Divisor/Ex63RlfPic.lean` の `rlfDeltaA : (Φ^rlf(A)^gp ⧸ Φ^birat) ≃+ ℝ` が
原文の結論そのものとして出ている。

★★**対比が要点である**——同じ「核 = 像」が、実現化すれば真、非実現化なら偽になる。
退化の修復は「名前どおりに条を足す」ではなく「型を実現化する」でなければならない。

## ★★★★★★何が残るのか(円)

`L = ℚ(√2)` では類数 1・単数の階数 1 で、次数 0 のアルキメデス因子の空間は 1 次元、
単数の格子はその中の階数 1 の格子だから、商はちょうど円 ℝ/ℤ である。
★本ファイルが形式化したのはその**濃度版**である(`uncountably_many_classes_remain`:
直線上の点のうち像に入るものは可算、入らないものは非可算)。
ℝ/ℤ との同型そのものは Dirichlet と類数有限を要するのでここでは主張していない
——それは `Found` 側(`ArithPicR.lean`)の仕事である。

## ★★反例は無限素点が 2 つ要る(ℚ では起きない)

`arch_kernel_trivial_of_subsingleton`——無限素点が 1 つしかなく `deg` がそこで消えないなら、
アルキメデス方向の核は 0 だけである。★実際 `L = ℚ` では類数 1・単数の階数 0 で、
「核 = 像」は非実現化のままでも真になる。反例が生きるのは無限素点が 2 つ以上のとき、
すなわち**単数の階数が 1 以上**のとき(`two_infinite_places_iff_units_rank_pos`)——
つまり原文が Dirichlet を引く必要が生じるのと**同じ条件**である。

## ★逸脱の記録

★`arithPic_ker_not_principal_subgroup` は Skeleton に無い仮定を 1 つ足している:
`deg` がアルキメデス成分で ℝ-斉次であること(`hsmul`)。
Skeleton の `degArith` / `arithDivOfElt` は `sorry` 本体の `def`(不透明定数)なので
それ自身については何も証明できない(9 例目と同じ事情)。そこで課されている条件だけを
`Thm64Spec` に括り出し、「その条件を満たす組」を走らせている。
★`hsmul` は `Found` の `arithDegreeLin : (ArithPlace L →₀ ℝ) →ₗ[ℝ] ℝ` が満たす性質なので、
原典に忠実な側の仮定である。
★これを外して加法性だけにすると反証はできない(ℝ を ℚ 上のベクトル空間と見た
Hamel 基底を使えば単射な加法写像が作れてしまう)。この点は形式化していない見立てである。

## ★位置づけ

「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例の **11 個目**
(`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・`Cor33Degenerate`・
`Prop22Degenerate`・`Cor13Degenerate`・`Prop12Degenerate`・`Ex61OrdDegenerate`・
`Ex63DegDegenerate`・`Prop12ForallRD`・本件)。
9 例目(`Ex63DegDegenerate`)で出た新種「素朴な修復が偽の主張を作る」の 2 例目であり、
10 例目(`Prop12ForallRD`)の「修復が強すぎる」とは別の形(こちらは端的に偽)である。

★本ファイルは `Skeleton` を import しない。検査時点(2026-09-06)の statement を
`ABC3/Skeleton/Divisor/ArithDivisor/Theorem64.lean` と
`ABC3/Skeleton/Divisor/ArithDivisor/Example63.lean` から写し取って固定している。
-/

namespace ABC3.Check.FrdI

open NumberField

universe u

/-! ## ★0. `Skeleton` から写し取った型(2026-09-06 時点)

★`degArith` / `arithDivOfElt` は `sorry` 本体の `def` なので、それ自身については
何も証明できない。そこで**課されている条件だけ**を構造体に括り出し、
「その条件を満たす組」を走らせる形に書き直す(`Ex63DegDegenerate` の
`ArithDivSpecOld8` と同じ手口)。 -/

/-- `Skeleton.Divisor.ArithPhiGp` の写し —— `Φ(L)^gp`。★**実現化していない**:
有限素点成分は `ℤ` 係数のままである。 -/
abbrev Thm64Gp (L : Type u) [Field L] [NumberField L] : Type u :=
  (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)

/-- **`Example 6.3` が `B(L) → Φ(L)^gp` と `deg^arith` に課している条件**を括り出したもの。

`div` / `div_mul` は `arithDivOfElt` / `arithDivOfElt_mul`、
`deg` / `deg_add` / `deg_div` は `degArith` / `degArith_add` / `degArith_arithDivOfElt`
(= 積公式)に対応する。★`Theorem 6.4` が足すのは `Function.Surjective deg` だけである。 -/
structure Thm64Spec (L : Type u) [Field L] [NumberField L] where
  /-- `Skeleton.Divisor.arithDivOfElt` —— `B(L) = L^× → Φ(L)^gp`。 -/
  div : Lˣ → Thm64Gp L
  /-- `Skeleton.Divisor.arithDivOfElt_mul`。 -/
  div_mul : ∀ f g : Lˣ, div (f * g) = div f + div g
  /-- `Skeleton.Divisor.degArith` —— `deg^arith_L : Φ(L)^gp → ℝ`。 -/
  deg : Thm64Gp L → ℝ
  /-- `Skeleton.Divisor.degArith_add`。 -/
  deg_add : ∀ x y : Thm64Gp L, deg (x + y) = deg x + deg y
  /-- `Skeleton.Divisor.degArith_arithDivOfElt` —— ★これが**積公式**である。 -/
  deg_div : ∀ f : Lˣ, deg (div f) = 0

/-- ★★**「像 ⊆ 核」は無条件に真**——これが積公式(`Example 6.3` の #9)である。

★したがって「核 = 像」が付け加えているのは**逆向きの包含だけ**であり、
偽になるのもそちらだけである。 -/
theorem range_div_subset_ker {L : Type u} [Field L] [NumberField L] (S : Thm64Spec L) :
    Set.range S.div ⊆ {x : Thm64Gp L | S.deg x = 0} := by
  rintro x ⟨f, rfl⟩
  exact S.deg_div f

/-! ## ★1. 退化 その 1——全射性は易しい半分である -/

/-- **無限素点 `v₀` の係数を読むだけの組**。`div ≡ 0`、`deg (a, r) = r v₀`。

★素点の対数も、`ord` も、積公式の中身も使っていない。 -/
def evalSpec (L : Type u) [Field L] [NumberField L] (v₀ : InfinitePlace L) : Thm64Spec L where
  div _ := 0
  div_mul _ _ := by simp
  deg x := x.2 v₀
  deg_add _ _ := rfl
  deg_div _ := by simp

/-- ★**`Theorem 6.4` の statement(全射性)は、係数を 1 つ読むだけで満たされる**。

`Finsupp.single v₀ r` を取ればよい。★アルキメデス成分が実係数であることだけが効いている。 -/
theorem evalSpec_surjective (L : Type u) [Field L] [NumberField L] (v₀ : InfinitePlace L) :
    Function.Surjective (evalSpec L v₀).deg := by
  intro r
  exact ⟨(0, Finsupp.single v₀ r), by simp [evalSpec]⟩

/-- ★★★★★**全射性だけでは、宣言名が約束した内容にならない**(退化 その 1)。

`Example 6.3` の条件をすべて満たし、`Theorem 6.4` の全射性も満たすのに、
主因子の像が `{0}` しかない組が存在する。★そのような組では
「核 ⊆ 像」は成り立たない(`Finsupp.single w 1` が核にあって像に無い)。

★原文 p.115 の証明は δ_A の全射性を既知として始め、Dirichlet 単数定理を
**核の等式**のために引く。Skeleton は易しい半分だけを型に残した。 -/
theorem surjectivity_alone_is_the_easy_half (L : Type u) [Field L] [NumberField L]
    {v w : InfinitePlace L} (hvw : v ≠ w) :
    ∃ S : Thm64Spec L,
      Function.Surjective S.deg ∧ (∀ f : Lˣ, S.div f = 0) ∧
        ¬ ({x : Thm64Gp L | S.deg x = 0} ⊆ Set.range S.div) := by
  refine ⟨evalSpec L v, evalSpec_surjective L v, fun _ => rfl, ?_⟩
  intro hsub
  have hx : ((0 : FinitePlace L →₀ ℤ), Finsupp.single w (1 : ℝ))
      ∈ {x : Thm64Gp L | (evalSpec L v).deg x = 0} := by
    show Finsupp.single w (1 : ℝ) v = 0
    simp [hvw]
  obtain ⟨f, hf⟩ := hsub hx
  have hf0 : ((0 : FinitePlace L →₀ ℤ), (0 : InfinitePlace L →₀ ℝ))
      = ((0 : FinitePlace L →₀ ℤ), Finsupp.single w (1 : ℝ)) := hf
  have h1 : Finsupp.single w (1 : ℝ) = 0 := (congrArg Prod.snd hf0).symm
  exact one_ne_zero (Finsupp.single_eq_zero.mp h1)

/-! ## ★2. ★★本命——非実現化の型で「核 = 主因子の像」を足すと偽になる -/

/-- ★数体の乗法群は可算である(`L` が `ℚ` 上有限次元だから)。 -/
theorem countable_units_of_numberField (L : Type u) [Field L] [NumberField L] :
    Countable Lˣ := by
  have b := Module.Free.chooseBasis ℚ L
  haveI : Countable L := (b.equivFun.toEquiv).countable_iff.mpr inferInstance
  exact Function.Injective.countable (f := (Units.val : Lˣ → L)) Units.val_injective

/-- ★★**アルキメデス方向の直線のうち、主因子の像に入る点は可算個しかない**。

`div f = (0, t • e)` なら `t` は `(div f) u / e u`(`e u ≠ 0` なる座標 `u`)で決まるので、
そのような `t` の全体は可算集合 `Lˣ` の像に含まれる。 -/
theorem countable_line_in_range (L : Type u) [Field L] [NumberField L]
    (div : Lˣ → Thm64Gp L) {e : InfinitePlace L →₀ ℝ} (he : e ≠ 0) :
    {t : ℝ | ∃ f : Lˣ, div f = (0, t • e)}.Countable := by
  classical
  obtain ⟨u, hu⟩ := Finsupp.ne_iff.mp he
  simp only [Finsupp.coe_zero, Pi.zero_apply] at hu
  haveI := countable_units_of_numberField L
  refine Set.Countable.mono ?_ (Set.countable_range (fun f : Lˣ => (div f).2 u / e u))
  rintro t ⟨f, hf⟩
  refine ⟨f, ?_⟩
  show (div f).2 u / e u = t
  rw [hf]
  have h2 : ((0 : FinitePlace L →₀ ℤ), t • e).2 u = t * e u := by simp
  rw [h2, mul_div_assoc, div_self hu, mul_one]

/-- ★★★★**核にアルキメデス方向の非零ベクトルが 1 本でもあれば、「核 ⊆ 像」は偽**。

`deg (0, t • e) = t * deg (0, e) = 0` なので直線 `{(0, t • e) | t ∈ ℝ}` は丸ごと核に入る。
これが像に含まれるなら ℝ が可算になってしまう。 -/
theorem no_nonzero_arch_kernel (L : Type u) [Field L] [NumberField L]
    (div : Lˣ → Thm64Gp L) (deg : Thm64Gp L → ℝ)
    (hsmul : ∀ (c : ℝ) (r : InfinitePlace L →₀ ℝ), deg (0, c • r) = c * deg (0, r))
    (hker : ∀ x : Thm64Gp L, deg x = 0 → ∃ f : Lˣ, div f = x)
    {e : InfinitePlace L →₀ ℝ} (he : e ≠ 0) (h0 : deg (0, e) = 0) : False := by
  refine Cardinal.not_countable_real (Set.Countable.mono ?_ (countable_line_in_range L div he))
  intro t _
  exact hker _ (by rw [hsmul, h0, mul_zero])

/-- ★★★**商には非可算個の類が残る**——`L = ℚ(√2)` で残る「円」の濃度版。

直線 `{(0, t • e)}` の点のうち主因子の像に入るものは可算、入らないものは**非可算**。
★`ℝ/ℤ`(円)との同型そのものは Dirichlet 単数定理と類数有限を要するので、
ここでは主張していない(`Found/Divisor/ArithPicR.lean` の仕事)。 -/
theorem uncountably_many_classes_remain (L : Type u) [Field L] [NumberField L]
    (div : Lˣ → Thm64Gp L) {e : InfinitePlace L →₀ ℝ} (he : e ≠ 0) :
    ¬ ({t : ℝ | ∃ f : Lˣ, div f = (0, t • e)}ᶜ).Countable := by
  intro h
  have hu := (countable_line_in_range L div he).union h
  rw [Set.union_compl_self] at hu
  exact Cardinal.not_countable_real hu

/-- ★★★★★★★**非実現化の `Φ^gp` で「核 = 主因子の像」を足すと偽になる**(退化 その 2)。

無限素点が 2 つあれば、`deg` をアルキメデス成分で ℝ-斉次とするどんな組についても
「核 ⊆ 主因子の像」は成り立たない。

証明: `a = deg (0, δ_v)`、`b = deg (0, δ_w)` として

* `a = 0` なら `e = δ_v` が核の非零元、
* `a ≠ 0` なら `e = b·δ_v − a·δ_w` が核の非零元(`e w = −a ≠ 0`)

を取り、`no_nonzero_arch_kernel` に渡す。★核は連続体濃度、像は可算。

★★「名前どおりに核の条を足す」修復は `False` を作り込む。正しい直し方は
台となる群を実現化することで、それが `Found/Divisor/ArithPicR.lean` の
`principalSpan_eq_ker` と `Found/Divisor/Ex63RlfPic.lean` の `rlfDeltaA` である。 -/
theorem arithPic_ker_not_principal_subgroup (L : Type u) [Field L] [NumberField L]
    {v w : InfinitePlace L} (hvw : v ≠ w) (S : Thm64Spec L)
    (hsmul : ∀ (c : ℝ) (r : InfinitePlace L →₀ ℝ), S.deg (0, c • r) = c * S.deg (0, r)) :
    ¬ ({x : Thm64Gp L | S.deg x = 0} ⊆ Set.range S.div) := by
  intro hsub
  have hker : ∀ x : Thm64Gp L, S.deg x = 0 → ∃ f : Lˣ, S.div f = x := fun x hx => hsub hx
  by_cases hv : S.deg (0, Finsupp.single v (1 : ℝ)) = 0
  · exact no_nonzero_arch_kernel L S.div S.deg hsmul hker
      (Finsupp.single_ne_zero.mpr one_ne_zero) hv
  · set a := S.deg ((0 : FinitePlace L →₀ ℤ), Finsupp.single v (1 : ℝ)) with ha
    set b := S.deg ((0 : FinitePlace L →₀ ℤ), Finsupp.single w (1 : ℝ)) with hb
    refine no_nonzero_arch_kernel L S.div S.deg hsmul hker
      (e := b • Finsupp.single v (1 : ℝ) - a • Finsupp.single w (1 : ℝ)) ?_ ?_
    · intro hzero
      have hw : (b • Finsupp.single v (1 : ℝ) - a • Finsupp.single w (1 : ℝ)) w = -a := by
        simp [Ne.symm hvw]
      rw [hzero] at hw
      simp at hw
      exact hv (by linarith)
    · have hsplit : ((0 : FinitePlace L →₀ ℤ),
          b • Finsupp.single v (1 : ℝ) - a • Finsupp.single w (1 : ℝ))
          = ((0 : FinitePlace L →₀ ℤ), b • Finsupp.single v (1 : ℝ))
            + ((0 : FinitePlace L →₀ ℤ), (-a) • Finsupp.single w (1 : ℝ)) := by
        simp [sub_eq_add_neg, neg_smul]
      rw [hsplit, S.deg_add, hsmul, hsmul, ← ha, ← hb]
      ring

/-- ★**集合の等式としての形** —— 核と主因子の像は一致しない。 -/
theorem arithPic_ker_ne_principal_subgroup (L : Type u) [Field L] [NumberField L]
    {v w : InfinitePlace L} (hvw : v ≠ w) (S : Thm64Spec L)
    (hsmul : ∀ (c : ℝ) (r : InfinitePlace L →₀ ℝ), S.deg (0, c • r) = c * S.deg (0, r)) :
    {x : Thm64Gp L | S.deg x = 0} ≠ Set.range S.div :=
  fun h => arithPic_ker_not_principal_subgroup L hvw S hsmul h.subset

/-! ## ★3. 反例には無限素点が 2 つ要る(= 単数の階数が 1 以上) -/

/-- ★★**無限素点が 1 つしかなければ、アルキメデス方向の核は 0 だけ**。

★したがって上の反例は `L = ℚ` では作れない。実際 `ℚ` は類数 1・単数の階数 0 で、
「核 = 像」は非実現化のままでも真になる(その事実自体はここでは形式化していない)。

★★ただし「無限素点が 1 つなら安全」ではない。虚二次体(例: `ℚ(√−5)`)は無限素点が
1 つだが類数が 1 でないので、有限素点側で「核 = 像」は破れる。★本ファイルが押さえたのは
**アルキメデス方向の破れ**だけで、類数による破れは別の水源である(ここでは形式化していない)。 -/
theorem arch_kernel_trivial_of_subsingleton (L : Type u) [Field L] [NumberField L]
    [Subsingleton (InfinitePlace L)] (v : InfinitePlace L) (deg : Thm64Gp L → ℝ)
    (hsmul : ∀ (c : ℝ) (r : InfinitePlace L →₀ ℝ), deg (0, c • r) = c * deg (0, r))
    (hv : deg (0, Finsupp.single v (1 : ℝ)) ≠ 0)
    {r : InfinitePlace L →₀ ℝ} (hr : deg (0, r) = 0) : r = 0 := by
  have hrv : r = (r v) • Finsupp.single v (1 : ℝ) := by
    refine Finsupp.ext fun u => ?_
    rw [Subsingleton.elim u v]
    simp
  rw [hrv, hsmul] at hr
  have hz : r v = 0 := by
    rcases mul_eq_zero.mp hr with h | h
    · exact h
    · exact absurd h hv
  refine Finsupp.ext fun u => ?_
  rw [Subsingleton.elim u v]
  simpa using hz

/-- ★★★**反例が生きる条件は「単数の階数が 1 以上」と同値**である。

`NumberField.Units.rank L = #V_∞(L) − 1` だから。★原文が Dirichlet 単数定理を
引く必要が生じるのと**同じ条件**で、非実現化の型は壊れる。 -/
theorem two_infinite_places_iff_units_rank_pos (L : Type u) [Field L] [NumberField L] :
    0 < NumberField.Units.rank L ↔ 1 < Fintype.card (InfinitePlace L) := by
  have hr : NumberField.Units.rank L = Fintype.card (InfinitePlace L) - 1 := rfl
  omega

/-! ## ★4. `L = ℚ(√2)` —— 反例が実在することの確認

★`math-planner` の見立て(`ℚ(√2)` で商に円が残る)を、無限素点が 2 つあることの
確認まで実行する。円そのものではなく**非可算性**が形式化された内容である。 -/

open Polynomial IntermediateField in
/-- `√2` は `ℚ` 上整である(`X^2 − 2` の根)。 -/
theorem sqrt2_isIntegral : IsIntegral ℚ (Real.sqrt 2) := by
  refine ⟨X ^ 2 - C 2, monic_X_pow_sub_C _ two_ne_zero, ?_⟩
  have h : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp [h]

open IntermediateField in
/-- **実二次体 `ℚ(√2)`** —— `ℝ` の中の `ℚ` 上の中間体として取る。 -/
noncomputable abbrev Qsqrt2 : IntermediateField ℚ ℝ := ℚ⟮Real.sqrt 2⟯

instance : FiniteDimensional ℚ Qsqrt2 :=
  IntermediateField.adjoin.finiteDimensional sqrt2_isIntegral

instance : NumberField Qsqrt2 := ⟨⟩

/-- `2` は有理数の平方ではない(`√2` の無理性)。 -/
theorem rat_sq_ne_two (b : ℚ) : b ^ 2 ≠ 2 := by
  intro hb
  have hb' : ((b : ℝ)) ^ 2 = 2 := by exact_mod_cast congrArg (fun q : ℚ => (q : ℝ)) hb
  refine irrational_sqrt_two ⟨|b|, ?_⟩
  rw [← hb', Real.sqrt_sq_eq_abs]
  push_cast
  ring

open Polynomial in
/-- `√2` の最小多項式は `X^2 − 2`(Eisenstein/Kummer 型の既約性)。 -/
theorem minpoly_sqrt2 : minpoly ℚ (Real.sqrt 2) = X ^ 2 - C 2 := by
  refine (minpoly.eq_of_irreducible_of_monic
    (X_pow_sub_C_irreducible_of_prime Nat.prime_two rat_sq_ne_two) ?_
    (monic_X_pow_sub_C _ two_ne_zero)).symm
  have h : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp [h]

open Polynomial in
/-- `[ℚ(√2) : ℚ] = 2`。 -/
theorem finrank_Qsqrt2 : Module.finrank ℚ Qsqrt2 = 2 := by
  rw [IntermediateField.adjoin.finrank sqrt2_isIntegral, minpoly_sqrt2, natDegree_X_pow_sub_C]

/-- `ℚ(√2) ⊆ ℝ ⊆ ℂ` という埋め込み。 -/
noncomputable def qsqrt2Embed : Qsqrt2 →+* ℂ :=
  Complex.ofRealHom.comp (IntermediateField.val Qsqrt2).toRingHom

/-- その埋め込みは実である(複素共役で不変)。 -/
theorem qsqrt2Embed_isReal : NumberField.ComplexEmbedding.IsReal qsqrt2Embed := by
  refine ComplexEmbedding.isReal_iff.mpr ?_
  ext x
  simp [qsqrt2Embed, ComplexEmbedding.conjugate]

/-- ★★**`ℚ(√2)` の無限素点は 2 つ以上ある**。

実埋め込みが 1 つあるので実素点の個数は 1 以上、`#実 + 2·#複素 = [K:ℚ] = 2` から
実素点はちょうど 2 個・複素素点は 0 個になる。 -/
theorem one_lt_card_infinitePlace_Qsqrt2 : 1 < Fintype.card (InfinitePlace Qsqrt2) := by
  classical
  have h1 : 0 < InfinitePlace.nrRealPlaces Qsqrt2 := by
    rw [← InfinitePlace.card_real_embeddings]
    exact Fintype.card_pos_iff.mpr ⟨⟨qsqrt2Embed, qsqrt2Embed_isReal⟩⟩
  have h2 := InfinitePlace.card_add_two_mul_card_eq_rank Qsqrt2
  rw [finrank_Qsqrt2] at h2
  have h3 := InfinitePlace.card_eq_nrRealPlaces_add_nrComplexPlaces Qsqrt2
  omega

/-- ★★★★★★★★**`L = ℚ(√2)` での反例**——`math-planner` の見立てのとおりである。

非実現化の `Φ^gp` の上では、`Example 6.3` の条件を満たし `deg` がアルキメデス成分で
ℝ-斉次などんな組についても、「核 ⊆ 主因子の像」は成り立たない。 -/
theorem arithPic_ker_not_principal_subgroup_qsqrt2 (S : Thm64Spec Qsqrt2)
    (hsmul : ∀ (c : ℝ) (r : InfinitePlace Qsqrt2 →₀ ℝ), S.deg (0, c • r) = c * S.deg (0, r)) :
    ¬ ({x : Thm64Gp Qsqrt2 | S.deg x = 0} ⊆ Set.range S.div) := by
  obtain ⟨v, w, hvw⟩ := Fintype.exists_pair_of_one_lt_card one_lt_card_infinitePlace_Qsqrt2
  exact arithPic_ker_not_principal_subgroup Qsqrt2 hvw S hsmul

#print axioms surjectivity_alone_is_the_easy_half
#print axioms arithPic_ker_not_principal_subgroup
#print axioms arithPic_ker_not_principal_subgroup_qsqrt2

end ABC3.Check.FrdI
