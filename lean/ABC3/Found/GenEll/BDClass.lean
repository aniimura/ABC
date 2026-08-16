import Mathlib.Order.Filter.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import ABC3.Meta.Claim

/-!
# [GenEll] Definition 1.2, (ii) —— BD-class(bounded discrepancy class)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.5。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.5):
> (ii) Fix a subset F[scr] ⊆ X(Q[bb][bar]). If α, β : F[scr] → R[bb] are functions, then we shall write

原文 (GenEll p.5):
> if there exists a [“constant”] C ∈ R[bb] such that β(x) − α(x) ≤ C (respectively,

## ★なぜここだけ先に実装できるか

`Definition 1.2` は (i)(ii) の 2 条から成り、**(i) は算術直線束**(`Definition 1.1`)を要する。
一方 **(ii) は関数 `F → ℝ` 上の同値関係でしかなく、Arakelov 理論を一切要しない。**
そして **abc の結論そのものがこの語彙で書かれている**——
[IUTchIV] Corollary 2.3 は「an **inequality of "bounded discrepancy classes"**」と述べる。

★したがって本ファイルは `Definition 1.2` の **(ii) だけ**の実装である。
条なしの `.src`(命題全体の完全実装の主張)は付けない。

## ★★原典の記号の向きが食い違っている(2026-08-16、3 箇所を PDF 目視して確定)

| 出典 | 目視した内容 |
|---|---|
| [GenEll] p.5, Def 1.2, (ii) | `α ≲_F β` ⟺ ∃C, **`β(x) − α(x) ≤ C`**(= `β − α` が**上**に有界) |
| [GenEll] p.11, Thm 2.1, (i) | 表題が **“Effective Mordell/ABC/Vojta Conjecture”** で、内容は `ht ≲ (1+ε)(log-diff + log-cond)` |
| [IUTchIV] p.54, Cor 2.3 | 同じ式に対し「i.e., the function `(1+ε)(log-diff + log-cond) − ht` **is bounded below** by a constant」<br>しかも `[cf. [GenEll], Definition 1.2, (ii)]` と明示的に引いている |

★**ここまでが観測である。以下は推論であり、区別して書く。**

**観測**: (a) の印字どおりに `≲` を読むと、(b) の式は
「`(1+ε)(log-diff + log-cond) ≤ ht + C`」になる。
一方 (b) の**表題**は “Effective Mordell/ABC/Vojta Conjecture” であり、
(c) は同じ形の式を「`(1+ε)(…) − ht` が**下**に有界」と言い換えている。
**abc/Vojta の内容は `ht ≤ (1+ε)(…) + C`**(= `(1+ε)(…) − ht` が下に有界)である。

**推論**: (b) の表題と (c) の言い換えは互いに整合し、abc の内容とも整合する。
(a) の印字だけがそれと逆になる。**最も無理のない読みは、(a) で `α` と `β` が
入れ替わっている(誤植)というものである。** ★ただしこれは**推論であって観測ではない**。

★**この推論に依存しない形で実装する**——下の `BDle` は (a) の印字を逐語で写し、
abc の向きは `BDge` として別に持つ。どちらが「正しい」かを Lean は判定しない。

**トリアージ**(`PLAN.md` §5-2):

- **① 我々のモデル化の誤り** —— **読字については排除した**。3 箇所とも 260 dpi で目視し、
  `pdftotext` ではなく PDF の版面そのものを読んでいる(事実2 の罠を踏んでいない)。
  ★ただし「abc の内容は `ht ≤ (1+ε)(…) + C` である」という**我々の数学的了解**が
  誤っている可能性までは排除していない。
- **② 未構築の数学** —— 該当しない。
- **③ 原典側** —— 上の推論が正しければ該当する。ただしこれは `PLAN.md` §5 が想定した
  「**飛躍**」ではなく**記号定義の誤植**であり、`Gap` の機構(追加仮説として型に出す)は当たらない。

★**独立確認: 取れた**(2026-08-16)。文脈を持たない子に、**結論を伏せたまま**同じ 3 箇所を
読ませた(260 dpi 全頁 ＋ 該当箇所 500 dpi 拡大)。子の判定:

> **逆** —— 同じ差 `(1+ε)(log-diff + log-cond) − ht` を、
> (a) は「≤ C(上に有界)」と定義し、(c) は「bounded below by a constant(下に有界)」と
> 言い換えているから。

★これは `GAP` 型(否定的)の判定であり、`memory/challenger-audit-without-context.md` の
基準では `OK` より強い証拠である。**読字については独立に裏が取れた。**
(それでも「abc の内容についての我々の了解」は依然として共有前提であり、そこは
第三者の確認を経ていない。)

★**ゆえの扱い**: **どちらも型に出し、両者が互いに逆であることを定理にする。**
逐語どおりの定義を `BDle`(原典の印字)、下流が要する向きを `BDge` として書き、
`bdle_iff_bdge_swap` で関係を明示し、`bdle_ne_bdge` で
**両者が本当に別物である**ことを反例で示す(記号の読み替えを「無害な流儀の違い」と
片づけないため)。

★**abc の主張を書くときは `BDge`(= 原文の印字の `≳`)を使う**。
その選択は本 docstring と `bdle_ne_bdge` によって機械可読になっている。

## ★何を証明したか

原文は「The relation “≈” **clearly** defines an equivalence relation」と書く。
その「clearly」を実際に証明した(`bdeq_equivalence`)。
また原文の「Thus, α ≈ β if and only if α ≲ β and α ≳ β」も証明した(`bdeq_iff`)。

## ★★Challenger 監査で 1 件落ちた(2026-08-16)——記号を BD-class に適用していなかった

文脈を持たない子に、原文 p.5 と本ファイルを突き合わせさせた結果:

> **GAP: 4** —— 原文は「Finally, we observe that it makes sense to apply the notation
> “≳”, “≲”, “≈” to BD-classes.」と 3 記号すべてを **BD-class に適用する**と述べるが、
> Lean にあるのは関数レベルの両立性 `bdle_congr` / `bdge_congr` の 2 本だけで、
> `BDClass F` 上の関係は 1 つも無く、**記号は BD-class に適用されていない**。

★**指摘は正しかった。** 両立性(= well-defined 性の前提)を示しただけで、
**商の上に関係を定義していなかった**。`BDClass.le` / `BDClass.ge` /
`BDClass.mk_eq_mk` / `BDClass.eq_of_le_of_ge` を足して埋めた。

★**これが「完成宣言を文脈のない子に監査させる」規律の実例である**
(`memory/challenger-audit-without-context.md`)。親は「両立性を示したから適用できる」と
思っていたが、原文が言っているのは**適用そのもの**だった。
-/

namespace ABC3.Found.GenEll

variable {F : Type*}

/-- 原文 `α ≲_F β` の**逐語どおり**の定義。

原文 (GenEll p.5):
> if there exists a [“constant”] C ∈ R[bb] such that β(x) − α(x) ≤ C (respectively,

★**この向きは印字どおりである**。abc の向きではない(上の docstring を参照)。 -/
def BDle (α β : F → ℝ) : Prop := ∃ C : ℝ, ∀ x, β x - α x ≤ C

/-- 原文 `α ≳_F β`。逐語どおり `α(x) − β(x) ≤ C`。

★**abc の不等式はこちらの向き**である——`ht ≤ (1+ε)(…) + C`。 -/
def BDge (α β : F → ℝ) : Prop := ∃ C : ℝ, ∀ x, α x - β x ≤ C

/-- 原文 `α ≈_F β`。逐語どおり `|α(x) − β(x)| ≤ C`。 -/
def BDeq (α β : F → ℝ) : Prop := ∃ C : ℝ, ∀ x, |α x - β x| ≤ C

/-! ### ★2 つの向きの関係 -/

theorem bdle_iff_bdge_swap (α β : F → ℝ) : BDle α β ↔ BDge β α := Iff.rfl

/-- ★**`≲` と `≳` は同じではない。** 反例で示す。

`α = 0`、`β = −(n : ℝ)` とすると `β − α = −n ≤ 0` なので `BDle α β` は成り立つが、
`α − β = n` は非有界なので `BDge α β` は成り立たない。

★これを置く理由: 原典の向きの食い違い(docstring)を「無害な流儀の違い」と
片づけないため。**向きを取り違えると主張が別物になる**ことの機械的な証拠である。 -/
theorem bdle_ne_bdge :
    ∃ α β : ℕ → ℝ, BDle α β ∧ ¬ BDge α β := by
  refine ⟨fun _ => 0, fun n => -(n : ℝ), ⟨0, fun n => by simp⟩, ?_⟩
  rintro ⟨C, hC⟩
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have := hC n
  simp only [sub_neg_eq_add, zero_add] at this
  linarith

/-! ### ★原文「Thus, α ≈ β if and only if α ≲ β and α ≳ β」 -/

theorem bdeq_iff (α β : F → ℝ) : BDeq α β ↔ BDle α β ∧ BDge α β := by
  constructor
  · rintro ⟨C, hC⟩
    exact ⟨⟨C, fun x => by have := abs_le.1 (hC x); linarith [this.1]⟩,
           ⟨C, fun x => by have := abs_le.1 (hC x); linarith [this.2]⟩⟩
  · rintro ⟨⟨C₁, h₁⟩, ⟨C₂, h₂⟩⟩
    refine ⟨max C₁ C₂, fun x => abs_le.2 ⟨?_, ?_⟩⟩
    · linarith [h₁ x, le_max_left C₁ C₂]
    · linarith [h₂ x, le_max_right C₁ C₂]

/-! ### ★原文「The relation “≈” clearly defines an equivalence relation」

★原文の「clearly」を実際に証明する。 -/

theorem bdeq_refl (α : F → ℝ) : BDeq α α := ⟨0, fun x => by simp⟩

theorem bdeq_symm {α β : F → ℝ} (h : BDeq α β) : BDeq β α := by
  obtain ⟨C, hC⟩ := h
  exact ⟨C, fun x => by rw [abs_sub_comm]; exact hC x⟩

theorem bdeq_trans {α β γ : F → ℝ} (h₁ : BDeq α β) (h₂ : BDeq β γ) : BDeq α γ := by
  obtain ⟨C₁, h₁⟩ := h₁
  obtain ⟨C₂, h₂⟩ := h₂
  refine ⟨C₁ + C₂, fun x => abs_le.2 ⟨?_, ?_⟩⟩
  · have a1 := (abs_le.1 (h₁ x)).1
    have a2 := (abs_le.1 (h₂ x)).1
    linarith
  · have b1 := (abs_le.1 (h₁ x)).2
    have b2 := (abs_le.1 (h₂ x)).2
    linarith

theorem bdeq_equivalence : Equivalence (BDeq (F := F)) :=
  ⟨bdeq_refl, bdeq_symm, bdeq_trans⟩

/-- BD-class の台となる setoid。原文の「equivalence class relative to this
equivalence relation」がこの商である。 -/
def bdSetoid (F : Type*) : Setoid (F → ℝ) := ⟨BDeq, bdeq_equivalence⟩

/-- 原文 `[α]_F`(BD-class)。 -/
abbrev BDClass (F : Type*) : Type _ := Quotient (bdSetoid F)

/-- 原文 `[α]`。 -/
def BDClass.mk (α : F → ℝ) : BDClass F := Quotient.mk (bdSetoid F) α

/-! ### ★原文「it makes sense to apply the notation “≳”, “≲”, “≈” to BD-classes」

★原文は「makes sense」と書くだけで、**well-defined 性を示していない**。
実際に示す——`≲` は `≈` と両立する。 -/

theorem bdle_congr {α α' β β' : F → ℝ} (hα : BDeq α α') (hβ : BDeq β β')
    (h : BDle α β) : BDle α' β' := by
  obtain ⟨Ca, ha⟩ := hα
  obtain ⟨Cb, hb⟩ := hβ
  obtain ⟨C, hC⟩ := h
  refine ⟨C + Ca + Cb, fun x => ?_⟩
  have h1 : α x - α' x ≤ Ca := (abs_le.1 (ha x)).2
  have h2 : β' x - β x ≤ Cb := by
    have := (abs_le.1 (hb x)).1
    linarith
  linarith [hC x]

theorem bdge_congr {α α' β β' : F → ℝ} (hα : BDeq α α') (hβ : BDeq β β')
    (h : BDge α β) : BDge α' β' :=
  bdle_congr hβ hα h

/-! ### ★★記号を BD-class に**実際に**適用する

★**この段は Challenger 監査(2026-08-16)の指摘で足した。**
上の `bdle_congr` / `bdge_congr` は**関数レベルの両立性**であって、
**商の上に関係を定義してはいなかった**。原文は

> Finally, we observe that it makes sense to apply the notation “≳”, “≲”, “≈” to BD-classes.

と、3 つの記号を **BD-class に適用する**と述べている。「意味を成す」で済ませず降ろす。 -/

/-- 原文の `≲` を BD-class の上へ降ろしたもの。 -/
def BDClass.le : BDClass F → BDClass F → Prop :=
  Quotient.lift₂ BDle fun _ _ _ _ ha hb =>
    propext ⟨bdle_congr ha hb, bdle_congr (bdeq_symm ha) (bdeq_symm hb)⟩

/-- 原文の `≳` を BD-class の上へ降ろしたもの。 -/
def BDClass.ge : BDClass F → BDClass F → Prop :=
  Quotient.lift₂ BDge fun _ _ _ _ ha hb =>
    propext ⟨bdge_congr ha hb, bdge_congr (bdeq_symm ha) (bdeq_symm hb)⟩

@[simp] theorem BDClass.le_mk (α β : F → ℝ) :
    BDClass.le (BDClass.mk α) (BDClass.mk β) ↔ BDle α β := Iff.rfl

@[simp] theorem BDClass.ge_mk (α β : F → ℝ) :
    BDClass.ge (BDClass.mk α) (BDClass.mk β) ↔ BDge α β := Iff.rfl

/-- ★原文の `≈` を BD-class に適用したものは、**類の等号**である。 -/
theorem BDClass.mk_eq_mk (α β : F → ℝ) :
    BDClass.mk α = BDClass.mk β ↔ BDeq α β := Quotient.eq

/-- ★`bdeq_iff`(原文の「α ≈ β iff α ≲ β and α ≳ β」)の**類の側での形**。
`≲` と `≳` の両方が成り立つ 2 つの BD-class は**同じ類**である。 -/
theorem BDClass.eq_of_le_of_ge {a b : BDClass F}
    (hle : BDClass.le a b) (hge : BDClass.ge a b) : a = b := by
  induction a using Quotient.inductionOn with
  | _ α =>
    induction b using Quotient.inductionOn with
    | _ β => exact (BDClass.mk_eq_mk α β).2 ((bdeq_iff α β).2 ⟨hle, hge⟩)

theorem BDClass.le_refl (a : BDClass F) : BDClass.le a a := by
  induction a using Quotient.inductionOn with
  | _ α => exact ⟨0, fun x => by simp⟩

theorem BDClass.le_trans {a b c : BDClass F}
    (hab : BDClass.le a b) (hbc : BDClass.le b c) : BDClass.le a c := by
  induction a using Quotient.inductionOn with
  | _ α =>
    induction b using Quotient.inductionOn with
    | _ β =>
      induction c using Quotient.inductionOn with
      | _ γ =>
        obtain ⟨C₁, h₁⟩ := hab
        obtain ⟨C₂, h₂⟩ := hbc
        exact ⟨C₁ + C₂, fun x => by linarith [h₁ x, h₂ x]⟩

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2` は (i) が算術直線束(`Definition 1.1`)を要するため未実装。 -/

def BDle.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2, (ii)",
    sectionId := "genell-def-1-2-ii" }

end ABC3.Found.GenEll
