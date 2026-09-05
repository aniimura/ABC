/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem
import Mathlib.Analysis.Real.Cardinality

/-!
# [FrdI] Example 6.3 の算術因子論——現行 statement は 10 件中 8 件が零写像・恒等写像で潰れる

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.112–114(`Example 6.3`)。

原文 (FrdI p.112):
> for the multiplicative monoid of elements of valuation ≤1; μ(F) ⊆O×

原文 (FrdI p.113):
> [given by mapping an element f ∈F × to the images of f in the various factors

原文 (FrdI p.114):
> vanishes on the image of B(L).

## ★★なぜこの検査が要るか

`Skeleton/Divisor/ArithDivisor/Example63.lean` には `sorry` が **10 個**ある
(`def` 2 個 = `arithDivOfElt` / `degArith`、`theorem` 8 個)。
★このうち **8 個は、零写像と恒等写像だけで「正直に」埋まってしまう**
——つまり `ord`・`deg` について何も言っていない。

| # | 宣言 | 何で閉じるか |
|---|---|---|
| 1 | `ordMon_nonarch_equiv` | `⟨fun _ => 0⟩` |
| 2 | `ordMon_arch_equiv` | `⟨fun _ => 0⟩` |
| 3 | `arithDivOfElt` | `(0, 0)` |
| 4 | `arithDivOfElt_mul` | `0 = 0 + 0` |
| 5 | `prime_arithPhi_equiv_places` | `⟨id⟩` |
| 7 | `degArith` | `0` |
| 8 | `degArith_add` | `0 = 0 + 0` |
| 9 | `degArith_arithDivOfElt`(積公式) | `0 = 0` |

残る 2 個は:

| # | 宣言 | 状況 |
|---|---|---|
| 6 | `support_arithPhi_eq_finite` | 零写像では閉じない。ただし `Finsupp` の一般論(指示関数)だけで閉じ、`ord` の内容は要らない(`exists_arithPhi_support_eq_contentless`) |
| 10 | `units_eq_roots_of_unity` | ★零写像を採ると **偽になる**。唯一の錨(`zeroArithDiv_not_units_eq_roots_of_unity`) |

★★したがって全体としては「正直に零で埋める」ことはできない。しかしそれは
**#10 という 1 本の錨のおかげ**であって、`ord`・`deg` の内容が statement に
書かれているからではない。#10 を落とせば `Example 6.3` の形式化は
**まるごと空になる**。

## ★★★#1・#2・#5 は永久に無内容である

この 3 つは `L` や `v` を何に取っても真である:

* `ordMon_nonarch_equiv : Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → ℕ)`
  ——`Nonempty (A → ℕ)` はどんな `A` でも真(`fun _ => 0`)
* `ordMon_arch_equiv : Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → NNReal)`
  ——同上
* `prime_arithPhi_equiv_places : Nonempty (Places L → Places L)`
  ——`Nonempty (X → X)` はどんな `X` でも真(`id`)

★とくに問題なのは語彙である。`ordMon_*` の statement に `ord` の語は
一度も現れない(`≅` を単なる関数 `→` に落としているので、`ord` を
特徴づける乗法性・全射性が消えている)。`prime_arithPhi_equiv_places` に至っては
「`Prime`」が型のどこにも無い——`Places L → Places L` は `Prime(Φ(L))` の
話をしていない。名前だけが原文を指している。

これらは `Skeleton` の他の `sorry` と違い、本体を書けば埋まる類ではない。
`sorry` を消しても意味が入らないので、statement の側を書き直すしかない。
`Found/Divisor/ArithOrd.lean` の `ordArch` / `ordFin` が実物である。

## ★★★★もう 1 つの発見——「素朴に退化を直すと偽になる」(新種)

上の #2 を「素直に強める」と、`Nonempty (… → NNReal)` は
`{x : L // x ≠ 0 ∧ v x ≤ 1}` を単元で割ったものが `ℝ≥0` と同型、になる。
★この強めた statement は偽である(`ordMon_arch_equiv_strengthened_false`)。

理屈は可算性である。`L` は数体なので `ℚ` 上有限次元、したがって可算。
その部分集合の商もまた可算だが、`ℝ≥0` は非可算。よってどんな商を取っても
同型にはならない。

★原文の `O_v^▷` は完備化 `F_v` の中の対象であって `L` の部分集合ではない
(原文 p.112 の「`ord(F_v) ≅ ℝ`(`v` がアルキメデス)」はそのために出る)。
`Skeleton` が `{x : L // …}` と書いたのが読み違いで、`Found/Divisor/ArithOrd.lean` は
そこを正しく `ordArch (v : InfinitePlace L) (x : v.Completion) : ℝ` と扱っている。

★★これは既知 8 例とは別種の危険である。既知の 8 例は
「落とした条件は、主張を偽にするか自明にするかのどちらかになる」だったが、
本件はその先——自明な形を「素直に」強めると、今度は偽が作り込まれる。
弱すぎる statement を強めるとき、強める方向を間違えると `False` を証明できる
`Skeleton` ができあがる。★退化の修復は「強めればよい」ではない。

## ★この検査の位置づけ

「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例の **9 つ目**
(`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・`Cor33Degenerate`・
`Prop22Degenerate`・`Cor13Degenerate`・`Prop12Degenerate`・`Ex61OrdDegenerate`・本件。
`VLocFalse` は枠そのものの反証なのでこの通番には入れていない)。

★本ファイルは `Skeleton/` を import しない。`Skeleton` 側は本体が配線して
statement が変わる予定なので、検査時点(2026-09-06)の statement を写し取って
固定している。写し元は `ABC3/Skeleton/Divisor/ArithDivisor/Example63.lean`。
-/

namespace ABC3.Check.FrdI

open NumberField

universe u

/-! ## ★0. `Skeleton/Divisor/ArithDivisor/Example63.lean` から写し取った型(2026-09-06 時点) -/

/-- `Skeleton.Divisor.Places` の写し —— `V(L)`。 -/
abbrev PlacesOld (L : Type u) [Field L] [NumberField L] : Type u :=
  FinitePlace L ⊕ InfinitePlace L

/-- `Skeleton.Divisor.ArithPhi` の写し —— `Φ(L)`。 -/
abbrev ArithPhiOld (L : Type u) [Field L] [NumberField L] : Type u :=
  (FinitePlace L →₀ ℕ) × (InfinitePlace L →₀ NNReal)

/-- `Skeleton.Divisor.ArithPhiGp` の写し —— `Φ(L)^gp`。 -/
abbrev ArithPhiGpOld (L : Type u) [Field L] [NumberField L] : Type u :=
  (FinitePlace L →₀ ℤ) × (InfinitePlace L →₀ ℝ)

/-! ## ★1. 旧形——`Example63.lean` が課している条件のうち、退化する 8 件

★`arithDivOfElt` / `degArith` は `sorry` 本体の `def` なので、それ自身については
何も証明できない。そこで課されている条件だけを構造体に括り出し、
「その条件を満たす組」を走らせる形に書き直す(`Ex61OrdDegenerate` の
`WeilOrdSpecOld` と同じ手口)。 -/

/-- **旧形(退化する 8 件)** —— `Example63.lean` の #1〜#5・#7〜#9 を逐語で写したもの。
★#6(`support_arithPhi_eq_finite`)と #10(`units_eq_roots_of_unity`)は
入れていない。理由は下の §3・§4 を見よ。 -/
structure ArithDivSpecOld8 (L : Type u) [Field L] [NumberField L] where
  /-- #1 `Skeleton.Divisor.ordMon_nonarch_equiv`。 -/
  ordMonNonarch : ∀ v : FinitePlace L, Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → ℕ)
  /-- #2 `Skeleton.Divisor.ordMon_arch_equiv`。 -/
  ordMonArch : ∀ v : InfinitePlace L, Nonempty ({x : L // x ≠ 0 ∧ v x ≤ 1} → NNReal)
  /-- #3 `Skeleton.Divisor.arithDivOfElt` —— `B(L) = L^× → Φ(L)^gp`。 -/
  div : Lˣ → ArithPhiGpOld L
  /-- #4 `Skeleton.Divisor.arithDivOfElt_mul`。 -/
  div_mul : ∀ f g : Lˣ, div (f * g) = div f + div g
  /-- #5 `Skeleton.Divisor.prime_arithPhi_equiv_places`。 -/
  primeEquiv : Nonempty (PlacesOld L → PlacesOld L)
  /-- #7 `Skeleton.Divisor.degArith` —— `deg^arith_L : Φ(L)^gp → ℝ`。 -/
  deg : ArithPhiGpOld L → ℝ
  /-- #8 `Skeleton.Divisor.degArith_add`。 -/
  deg_add : ∀ x y : ArithPhiGpOld L, deg (x + y) = deg x + deg y
  /-- #9 `Skeleton.Divisor.degArith_arithDivOfElt` —— ★これが積公式である。 -/
  deg_div : ∀ f : Lˣ, deg (div f) = 0

/-- ★★**零写像・恒等写像は旧形 8 件をすべて充足する**。

`⟨fun _ => 0⟩`・`(0, 0)`・`⟨id⟩`・`0` を並べるだけである。★とくに #9 は
`Example 6.3` の中で唯一の深い主張(積公式 `∏_v |x|_v = 1`)なのに、
`deg ≡ 0` のもとでは `0 = 0` に潰れる。 -/
def zeroArithDiv (L : Type u) [Field L] [NumberField L] : ArithDivSpecOld8 L where
  ordMonNonarch _ := ⟨fun _ => 0⟩
  ordMonArch _ := ⟨fun _ => 0⟩
  div _ := 0
  div_mul _ _ := by simp
  primeEquiv := ⟨id⟩
  deg _ := 0
  deg_add _ _ := by simp
  deg_div _ := rfl

@[simp] theorem zeroArithDiv_div (L : Type u) [Field L] [NumberField L] (f : Lˣ) :
    (zeroArithDiv L).div f = 0 := rfl

@[simp] theorem zeroArithDiv_deg (L : Type u) [Field L] [NumberField L] (x : ArithPhiGpOld L) :
    (zeroArithDiv L).deg x = 0 := rfl

/-- ★★★★★★**[FrdI] Example 6.3 の現行 statement は 10 件中 8 件が零写像で充足される**。

`Example63.lean` の `sorry` のうち #1〜#5・#7〜#9 は、数学的内容ゼロで
いっぺんに埋まる。 -/
theorem example_6_3_eight_of_ten_satisfied_by_zero (L : Type u) [Field L] [NumberField L] :
    ∃ S : ArithDivSpecOld8 L,
      (∀ f : Lˣ, S.div f = 0) ∧ (∀ x : ArithPhiGpOld L, S.deg x = 0) :=
  ⟨zeroArithDiv L, fun _ => rfl, fun _ => rfl⟩

/-! ## ★2. #1・#2・#5 は「永久に」無内容である

★これらは `L`・`v` に何を入れても真になる。つまり `Example 6.3` の本文を
どう読み替えても、この 3 つの statement からは `ord` も `Prime` も出てこない。 -/

/-- ★**#1 の statement は任意の型で真** —— `ord` の語はどこにも効いていない。 -/
theorem ordMon_nonarch_statement_vacuous (A : Type*) : Nonempty (A → ℕ) := ⟨fun _ => 0⟩

/-- ★**#2 の statement は任意の型で真** —— 同上。 -/
theorem ordMon_arch_statement_vacuous (A : Type*) : Nonempty (A → NNReal) := ⟨fun _ => 0⟩

/-- ★★**#5 の statement は任意の型で真** —— `Nonempty (X → X)` は `id` で埋まる。
★`prime_arithPhi_equiv_places` の型には「`Prime`」が一度も現れない。 -/
theorem prime_arithPhi_statement_vacuous (X : Type*) : Nonempty (X → X) := ⟨id⟩

/-! ## ★3. #6 は零写像では閉じないが、`ord` の内容も要らない

`support_arithPhi_eq_finite : ∀ T, ∃ a : ArithPhi L, a.1.support = T` は
`a = 0` では閉じない(`T` が空でなければ台が合わない)。しかし
`T` の指示関数を取れば閉じる。★`Finsupp` の一般論だけで、`ord` も素点も
使わない。 -/

/-- ★**#6 は `Finsupp` の一般論だけで閉じる** —— `T` の指示関数を取ればよい。
`ord` の内容は 1 ミリも要らない。 -/
theorem exists_arithPhi_support_eq_contentless (L : Type u) [Field L] [NumberField L]
    (T : Finset (FinitePlace L)) : ∃ a : ArithPhiOld L, a.1.support = T := by
  classical
  refine ⟨(⟨T, fun t => if t ∈ T then 1 else 0, fun t => ?_⟩, 0), rfl⟩
  by_cases h : t ∈ T <;> simp [h]

/-! ## ★4. #10 が唯一の錨——零写像を採ると偽になる

`units_eq_roots_of_unity : arithDivOfElt L f = 0 → ∃ n > 0, f^n = 1` は、
`arithDivOfElt ≡ 0` のもとでは「`L^×` の元がすべて 1 の冪根」を意味する。
`2 ∈ L^×` は 1 の冪根ではないので、これは偽である。 -/

/-- ★`2` は 1 の冪根ではない(標数 0)。 -/
theorem two_pow_ne_one (L : Type u) [Field L] [NumberField L] {n : ℕ} (hn : 0 < n) :
    (2 : L) ^ n ≠ 1 := by
  intro h
  have hcast : ((2 ^ n : ℕ) : L) = ((1 : ℕ) : L) := by push_cast; simpa using h
  have h2 : (2 : ℕ) ^ n = 1 := Nat.cast_injective hcast
  rcases Nat.pow_eq_one.mp h2 with h3 | h3
  · exact absurd h3 (by norm_num)
  · exact absurd h3 hn.ne'

/-- ★★★★★**#10 は零写像を排除する** —— `arithDivOfElt ≡ 0` を採ると
`units_eq_roots_of_unity` は偽になる。

★これが `Example 6.3` の 10 件のうち唯一の錨である。#10 を落とすと
残り 9 件は(#6 も指示関数で)すべて内容ゼロで埋まり、形式化はまるごと空になる。 -/
theorem zeroArithDiv_not_units_eq_roots_of_unity (L : Type u) [Field L] [NumberField L]
    (h : ∀ f : Lˣ, (zeroArithDiv L).div f = 0 → ∃ n : ℕ, 0 < n ∧ (f : L) ^ n = 1) : False := by
  obtain ⟨n, hn, hfn⟩ := h (Units.mk0 (2 : L) (by norm_num)) rfl
  exact two_pow_ne_one L hn (by simpa using hfn)

/-! ## ★5. ★★新種——「素朴に退化を直すと偽になる」

#2 を素直に強めた形、すなわち
`{x : L // x ≠ 0 ∧ v x ≤ 1}` の(単元による)商が `ℝ≥0` と同型、は偽である。
可算性で示す。 -/

/-- ★数体は可算である(`ℚ` 上有限次元だから)。 -/
theorem countable_of_numberField (L : Type u) [Field L] [NumberField L] : Countable L := by
  have b := Module.Free.chooseBasis ℚ L
  exact (b.equivFun.toEquiv).countable_iff.mpr inferInstance

/-- ★`ℝ≥0` は非可算(`exp : ℝ → ℝ≥0` が単射だから)。 -/
theorem not_countable_nnreal : ¬ Countable NNReal := by
  intro h
  have hinj : Function.Injective (fun x : ℝ => Real.toNNReal (Real.exp x)) := by
    intro a b hab
    refine Real.exp_injective ?_
    have h' := congrArg NNReal.toReal hab
    simpa [Real.coe_toNNReal _ (Real.exp_pos a).le,
      Real.coe_toNNReal _ (Real.exp_pos b).le] using h'
  haveI : Countable ℝ := @Function.Injective.countable ℝ NNReal h _ hinj
  exact Cardinal.not_countable_real Set.countable_univ

/-- ★★**`{x : L // x ≠ 0 ∧ v x ≤ 1}` から `ℝ≥0` への全射は存在しない**。

左辺は可算(`L` が可算)、右辺は非可算。 -/
theorem ordMonArch_not_surjective (L : Type u) [Field L] [NumberField L] (v : InfinitePlace L)
    (g : {x : L // x ≠ 0 ∧ v x ≤ 1} → NNReal) : ¬ Function.Surjective g := by
  intro hg
  haveI : Countable L := countable_of_numberField L
  exact not_countable_nnreal hg.countable

/-- ★★★★★★**#2 を「素直に」強めると偽になる**(新種)。

`π` を任意の商写像(とくに単元による商 `O_v^▷ → O_v^▷/O_v^×`)とするとき、
その商が `ℝ≥0` と同型になることはありえない。

★理由:`L` は可算なので `{|x|_v : x ∈ L^×}` は `ℝ_{>0}` の可算部分群にすぎない。
原文の `O_v^▷` は完備化 `F_v` の中の対象であって `L` の部分集合ではない。

★★これは「退化した statement を強めれば直る」という直観が間違っていることを示す。
強める方向を間違えると、`sorry` を消せない(消せば `False` が出る)`Skeleton` ができる。
正しい直し方は `Found/Divisor/ArithOrd.lean` のように
`ordArch (v : InfinitePlace L) (x : v.Completion) : ℝ` と完備化の側で書くことである。 -/
theorem ordMon_arch_equiv_strengthened_false (L : Type u) [Field L] [NumberField L]
    (v : InfinitePlace L) {Q : Type*} (π : {x : L // x ≠ 0 ∧ v x ≤ 1} → Q)
    (hπ : Function.Surjective π) (e : Q ≃ NNReal) : False :=
  ordMonArch_not_surjective L v (fun x => e (π x)) (e.surjective.comp hπ)

/-- ★商を取らない場合(`π = id`)—— `{x : L // x ≠ 0 ∧ v x ≤ 1} ≃ ℝ≥0` も偽。 -/
theorem ordMon_arch_equiv_nnreal_false (L : Type u) [Field L] [NumberField L]
    (v : InfinitePlace L) (e : {x : L // x ≠ 0 ∧ v x ≤ 1} ≃ NNReal) : False :=
  ordMon_arch_equiv_strengthened_false L v id Function.surjective_id e

/-- **`O_v^▷` を `L` の中で素朴に読んだもの** —— 乗法モノイド
`{x ∈ L | x ≠ 0, |x|_v ≤ 1}`。★原文の `O_v^▷` は完備化 `F_v` の中の対象なので、
これは読み違いの側である。 -/
def OvTriOld (L : Type u) [Field L] [NumberField L] (v : InfinitePlace L) : Submonoid L where
  carrier := {x : L | x ≠ 0 ∧ v x ≤ 1}
  one_mem' := ⟨one_ne_zero, by simp⟩
  mul_mem' := by
    rintro a b ⟨ha0, ha1⟩ ⟨hb0, hb1⟩
    exact ⟨mul_ne_zero ha0 hb0, by rw [map_mul]; exact mul_le_one₀ ha1 (apply_nonneg v b) hb1⟩

/-- ★★★**モノイド版** —— `O_v^▷`(素朴読み)のどんなモノイド商も `ℝ≥0` と
`≃*` にはならない。単元による商 `ord(O_v^▷) = O_v^▷/O_v^×` はその特別な場合である。 -/
theorem ordMon_arch_mulEquiv_false (L : Type u) [Field L] [NumberField L]
    (v : InfinitePlace L) {Q : Type*} [Monoid Q] (π : OvTriOld L v →* Q)
    (hπ : Function.Surjective π) (e : Q ≃* NNReal) : False :=
  ordMon_arch_equiv_strengthened_false L v (fun x => π x) hπ e.toEquiv

end ABC3.Check.FrdI
