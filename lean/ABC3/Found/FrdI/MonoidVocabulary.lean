import Mathlib.GroupTheory.MonoidLocalization.Basic
import Mathlib.GroupTheory.Subgroup.Saturated
import Mathlib.Data.ENat.Basic
import Mathlib.Algebra.Order.Group.Nat
import Mathlib.Tactic.Abel
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Algebra.Group.Submonoid.Membership
import Mathlib.Data.ZMod.Basic
import Mathlib.GroupTheory.Congruence.Hom
import ABC3.Meta.Claim

/-!
# [FrdI] §0 のモノイド語彙 —— `sharp` / `integral` / `saturated`

原典: S. Mochizuki, *The Geometry of Frobenioids I: The General Theory* [FrdI]、
物理 p.11(全 126 ページ。**400 dpi 目視確認 2026-08-15**)。

★**原文は加法で書く**:

原文 (FrdI p.11):
> of M will be written additively. We shall denote by

原文 (FrdI p.11):
> the submonoid [which, in fact, forms a group] of invertible elements of M, by

原文 (FrdI p.11):
> the natural homomorphism from M to its groupification M gp. Thus, M gp is the

そして3語の定義は1文にまとまっている:

原文 (FrdI p.11):
> that M is torsion-free if M has no torsion elements; we shall say that M is sharp if

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;

原文 (FrdI p.11):
> we shall say that M is saturated if every a ∈M gp for which n · a lies in the image

原文 (FrdI p.11):
> of M for some n ∈N≥1 lies in the image of M.

★`M ±` は原文では **`M^±`**(上付き)であり、`pdftotext` が上付きを空白区切りの別トークンに
するために `M ±` と出る。`papers.json` の `FrdI.notationNotes` に実測を記録した。

## ★なぜここを作るか(層0の仕分けの結論)

2026-08-15 の層0(129節点)の仕分けで、**いま着手できるもの**として6件を名指しした。
その筆頭がこの3語である。理由は2つ:

1. **mathlib に土台がある。** `IsAffineMonoid = IsCancelMul + Monoid.FG + IsMulTorsionFree`
   (`Mathlib/Algebra/AffineMonoid/Basic.lean`)。
2. ★**これは Frobenioid の入口である。** `Frobenioid` は
   `[IUTchIII] Definition 3.8` の解消条件(中層待ち2件のうち1件)であり、
   [FrdI] は機械概算で最も密(到達45項目・辺184本・平均出次数 4.1)。
   **層0の最も安い部分が、中層の最大の壁の入口になっている。**

## ★mathlib の何を使ったか(車輪の再発明をしない)

- `M^gp`(groupification): mathlib の **`AddLocalization (⊤ : AddSubmonoid M)`**。
  ★原文の同値関係「`a₁ + b₂ + c = b₁ + a₂ + c` for some `c ∈ M`」は
  `AddLocalization.r` と**同じ形**である(原文 p.11 目視)。
  mathlib には `GrothendieckAddGroup M` という `abbrev` もあるが、
  2026-08-15 時点の本ビルドでは `import` しても識別子が見えなかった(原因未特定)。
  `abbrev` の右辺そのものを使っているので、実質は同じものである。
- `saturated`: **`AddSubmonoid.NSMulSaturated`**(`Mathlib/GroupTheory/Subgroup/Saturated.lean`)。
  下の `frdI_saturated_iff` で、原文の定義がこれと一致することを示す。
- `sharp` に相当する名前は mathlib に無い(2026-08-15 実測: `\bsharp\b` は 17 件あるが
  **全部「sharp bound」**で別物)。1行なので自前で書く。
-/

namespace ABC3.Found.FrdI

open AddLocalization

variable (M : Type*) [AddCommMonoid M]

/-! ### `M^gp` —— groupification -/

/-- **`M^gp`** —— 原文の groupification。mathlib の `AddLocalization` そのもの。

原文 (FrdI p.11):
> the natural homomorphism from M to its groupification M gp. Thus, M gp is the
-/
abbrev Gp : Type _ := AddLocalization (⊤ : AddSubmonoid M)

/-- 正準写像 `M → M^gp`。 -/
def toGp (m : M) : Gp M := AddLocalization.mk m 0

variable {M}

/-- `toGp` が等しいことの判定 —— 原文の同値関係そのもの
(`a₁ + b₂ + c = b₁ + a₂ + c` for some `c ∈ M`)。 -/
theorem toGp_eq_iff {m m' : M} : toGp M m = toGp M m' ↔ ∃ c : M, m + c = m' + c := by
  constructor
  · intro h
    rw [toGp, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at h
    obtain ⟨c, hc⟩ := h
    refine ⟨(c : M), ?_⟩
    simpa [add_comm] using hc
  · rintro ⟨c, hc⟩
    rw [toGp, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    refine ⟨⟨c, trivial⟩, ?_⟩
    simpa [add_comm] using hc

variable (M)

/-! ### 3語の定義 -/

/-- **[FrdI] §0 `sharp`** —— `M^± = 0`、すなわち可逆元が `0` だけ。

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;
-/
def IsSharp : Prop := ∀ a : M, IsAddUnit a → a = 0

/-- **[FrdI] §0 `integral`** —— `M → M^gp` が単射。

原文 (FrdI p.11):
> M ± = 0; we shall say that M is integral if the natural map M →M gp is injective;
-/
def IsIntegralMonoid : Prop := Function.Injective (toGp M)

/-- **[FrdI] §0 `saturated`** —— `n · a` が `M` の像に入るような `a ∈ M^gp`
(`n ≥ 1`)は、それ自身 `M` の像に入る。

原文 (FrdI p.11):
> we shall say that M is saturated if every a ∈M gp for which n · a lies in the image

原文 (FrdI p.11):
> of M for some n ∈N≥1 lies in the image of M.
-/
def IsSaturatedMonoid : Prop :=
  ∀ (a : Gp M) (n : ℕ), 0 < n → n • a ∈ Set.range (toGp M) → a ∈ Set.range (toGp M)

/-! ### mathlib との接続 -/

/-- ★**原文の `saturated` の条件は、mathlib の `AddSubmonoid.NSMulSaturated` そのものである。**

原文は「`n ≥ 1` かつ `n·a` が入るなら `a` も入る」と述べ、mathlib は
「`n • g ∈ H` なら `n = 0` または `g ∈ H`」と述べる。**同じ主張**である。

★これを任意の部分モノイドについて示しておくのは、`M^gp` の中の `M` の像に
そのまま当てられるようにするためである。 -/
theorem frdI_saturated_iff {G : Type*} [AddMonoid G] (H : AddSubmonoid G) :
    (∀ (a : G) (n : ℕ), 0 < n → n • a ∈ H → a ∈ H) ↔ H.NSMulSaturated := by
  constructor
  · intro h n g hng
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · exact Or.inl rfl
    · exact Or.inr (h g n hn hng)
  · intro h a n hn hmem
    rcases h hmem with h0 | h1
    · omega
    · exact h1

/-- **cancellative なら integral**(mathlib の局所化の単射性が要求するのと同じ条件)。 -/
theorem isIntegralMonoid_of_isCancelAdd [IsCancelAdd M] : IsIntegralMonoid M := by
  intro m m' h
  obtain ⟨c, hc⟩ := toGp_eq_iff.mp h
  exact add_right_cancel hc

/-- ★**逆も成り立つ —— integral なら cancellative**(2026-08-15 追加)。

★これは `[FrdI] Proposition 1.5` の証明の第1段を写すために測った問いである。
原文は `𝔽_Φ` が totally epimorphic であることを
「`𝒟` が totally epimorphic」「**pre-divisorial monoid が integral**」
「`Definition 1.1, (ii), (a)` の単射性」の3つから出す。
そこで **integral が cancellative を与えるか**が関門になる。**与える。**

`a + c = b + c` は `toGp_eq_iff` によりそのまま `toGp a = toGp b` を意味するので、
`toGp` の単射性が `a = b` を与える。 -/
theorem isCancelAdd_of_isIntegralMonoid (h : IsIntegralMonoid M) : IsCancelAdd M where
  add_left_cancel a b c hbc := h (toGp_eq_iff.mpr ⟨a, by
    rw [add_comm b a, add_comm c a]; exact hbc⟩)
  add_right_cancel a b c hac := h (toGp_eq_iff.mpr ⟨a, hac⟩)

/-- ★**integral ⟺ cancellative**。 -/
theorem isIntegralMonoid_iff_isCancelAdd : IsIntegralMonoid M ↔ Nonempty (IsCancelAdd M) :=
  ⟨fun h => ⟨isCancelAdd_of_isIntegralMonoid M h⟩,
   fun ⟨h⟩ => letI := h; isIntegralMonoid_of_isCancelAdd M⟩

/-! ### ★非退化 —— 3語それぞれに「満たす例」と「満たさない例」

★これが無いと「型は付くが中身が空」を排除できない。 -/

/-- `ℕ` は **sharp** —— 可逆元は `0` だけ。 -/
theorem isSharp_nat : IsSharp ℕ := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h := u.val_neg
  omega

/-- ★`ℤ` は **sharp でない** —— `1` が可逆で `1 ≠ 0`。 -/
theorem not_isSharp_int : ¬ IsSharp ℤ := by
  intro h
  have h1 : IsAddUnit (1 : ℤ) := ⟨⟨1, -1, by omega, by omega⟩, rfl⟩
  have h2 := h 1 h1
  omega

/-- `ℕ` は **integral** —— `ℕ` は cancellative。 -/
theorem isIntegralMonoid_nat : IsIntegralMonoid ℕ :=
  isIntegralMonoid_of_isCancelAdd ℕ

/-- ★`ℕ∞` は **integral でない** —— `1 + ⊤ = 2 + ⊤` なので `toGp 1 = toGp 2`、
しかし `1 ≠ 2`。 -/
theorem not_isIntegralMonoid_enat : ¬ IsIntegralMonoid ℕ∞ := by
  intro h
  have heq : toGp ℕ∞ 1 = toGp ℕ∞ 2 := toGp_eq_iff.mpr ⟨⊤, by simp⟩
  have := h heq
  simp at this

/-- ★`ℕ` の中の**偶数**のなす部分モノイドは **saturated でない** ——
`2 • 1 = 2` は入るが `1` は入らない。

`frdI_saturated_iff` により、これは原文の saturated の条件が破れている実例である。 -/
theorem not_saturated_evens :
    ¬ (∀ (a : ℕ) (n : ℕ), 0 < n → n • a ∈ AddSubmonoid.closure ({2} : Set ℕ) →
        a ∈ AddSubmonoid.closure ({2} : Set ℕ)) := by
  intro h
  have h2 : (2 : ℕ) ∈ AddSubmonoid.closure ({2} : Set ℕ) :=
    AddSubmonoid.subset_closure rfl
  have h1 : (1 : ℕ) ∈ AddSubmonoid.closure ({2} : Set ℕ) := by
    refine h 1 2 (by omega) ?_
    simpa using h2
  rw [AddSubmonoid.mem_closure_singleton] at h1
  obtain ⟨n, hn⟩ := h1
  simp only [smul_eq_mul] at hn
  omega

/-- `⊤` は **saturated**(自明に)。「全部落とす」定義ではないことの対照。 -/
theorem saturated_top :
    ∀ (a : ℕ) (n : ℕ), 0 < n → n • a ∈ (⊤ : AddSubmonoid ℕ) → a ∈ (⊤ : AddSubmonoid ℕ) :=
  fun _ _ _ _ => trivial

/-! ### ★独立性 —— どの1語も他の2語から従わない -/

/-- `M → M^gp` が全射なら saturated(自明に)。 -/
theorem isSaturatedMonoid_of_range_eq_univ (h : Set.range (toGp M) = Set.univ) :
    IsSaturatedMonoid M := fun a _ _ _ => h ▸ Set.mem_univ a

/-- `ℤ` では `M → M^gp` は全射。 -/
theorem range_toGp_int : Set.range (toGp ℤ) = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  induction x using AddLocalization.ind with
  | _ p =>
    obtain ⟨a, b⟩ := p
    refine ⟨a - (b : ℤ), ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simp⟩

/-- `ℕ∞` では `M → M^gp` は全射(`⊤` がすべてを潰すので `M^gp` は1点)。 -/
theorem range_toGp_enat : Set.range (toGp ℕ∞) = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  induction x using AddLocalization.ind with
  | _ p =>
    obtain ⟨a, b⟩ := p
    refine ⟨0, ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    refine ⟨⟨⊤, trivial⟩, ?_⟩
    simp

/-- ★**独立性 (1)**: `ℤ` は **integral かつ saturated だが sharp でない**。 -/
theorem independence_sharp :
    IsIntegralMonoid ℤ ∧ IsSaturatedMonoid ℤ ∧ ¬ IsSharp ℤ :=
  ⟨isIntegralMonoid_of_isCancelAdd ℤ,
   isSaturatedMonoid_of_range_eq_univ ℤ range_toGp_int,
   not_isSharp_int⟩

/-- `ℕ∞` は sharp —— 正準順序つき加法モノイドでは `a + b = 0 → a = 0`。 -/
theorem isSharp_enat : IsSharp ℕ∞ := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h := u.val_neg
  exact (add_eq_zero.mp h).1

/-- ★**独立性 (2)**: `ℕ∞` は **sharp かつ saturated だが integral でない**。 -/
theorem independence_integral :
    IsSharp ℕ∞ ∧ IsSaturatedMonoid ℕ∞ ∧ ¬ IsIntegralMonoid ℕ∞ :=
  ⟨isSharp_enat,
   isSaturatedMonoid_of_range_eq_univ ℕ∞ range_toGp_enat,
   not_isIntegralMonoid_enat⟩

/-! ### 独立性 (3) —— `saturated` が他の2語から従わないこと -/

/-- `⟨2,3⟩ ⊆ ℕ` が生成する部分モノイド。`1` を含まないことが効く。 -/
abbrev TwoThree : AddSubmonoid ℕ := AddSubmonoid.closure ({2, 3} : Set ℕ)

theorem two_mem_twoThree : (2 : ℕ) ∈ TwoThree :=
  AddSubmonoid.subset_closure (by simp)

theorem three_mem_twoThree : (3 : ℕ) ∈ TwoThree :=
  AddSubmonoid.subset_closure (by simp)

theorem one_not_mem_twoThree : (1 : ℕ) ∉ TwoThree := by
  intro h
  rw [AddSubmonoid.mem_closure_pair] at h
  obtain ⟨m, n, hmn⟩ := h
  simp only [smul_eq_mul] at hmn
  omega

/-- ★**独立性 (3)**: `⟨2,3⟩` は **sharp かつ integral だが saturated でない**。

`a := mk 3 2`(= `1`)は `2 • a = mk 6 4`(= `2`)が像に入るのに、
それ自身は像に入らない —— `1 ∉ ⟨2,3⟩` だから。 -/
theorem independence_saturated :
    IsSharp (TwoThree) ∧ IsIntegralMonoid (TwoThree) ∧ ¬ IsSaturatedMonoid (TwoThree) := by
  refine ⟨?_, isIntegralMonoid_of_isCancelAdd _, ?_⟩
  · intro a ha
    obtain ⟨u, rfl⟩ := ha
    have h := u.val_neg
    have hsum : ((u : TwoThree) : ℕ) + ((u.neg : TwoThree) : ℕ) = 0 := congrArg Subtype.val h
    refine Subtype.ext ?_
    show ((u : TwoThree) : ℕ) = 0
    omega
  · intro hsat
    set a : Gp TwoThree :=
      AddLocalization.mk ⟨3, three_mem_twoThree⟩ ⟨⟨2, two_mem_twoThree⟩, trivial⟩ with ha
    have hmem : (2 : ℕ) • a ∈ Set.range (toGp TwoThree) := by
      refine ⟨⟨2, two_mem_twoThree⟩, ?_⟩
      rw [ha, AddLocalization.mk_nsmul, toGp, AddLocalization.mk_eq_mk_iff,
        AddLocalization.r_iff_exists]
      refine ⟨0, ?_⟩
      ext
      simp
    obtain ⟨m, hm⟩ := hsat a 2 (by omega) hmem
    rw [ha, toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at hm
    obtain ⟨c, hc⟩ := hm
    have : ((m : TwoThree) : ℕ) + 2 = 3 := by
      have := congrArg Subtype.val hc
      simp at this
      omega
    refine one_not_mem_twoThree ?_
    have hm1 : ((m : TwoThree) : ℕ) = 1 := by omega
    simp [← hm1]

/-! ### `torsion-free` と、原文が主張する含意

原文 (FrdI p.11):
> inductive system I∗determines a natural morphism

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if

★**「which」の正体**: 直前の `M → M^pf`(**perfection** への自然な射)である。
`M^gp` ではない。原文は `M^pf` を「乗法 `n·` による帰納極限」として定める:

原文 (FrdI p.11):
> the perfection of M, that is to say, the inductive limit of the inductive system I∗of

★★**原文は「torsion element」を定義していない**(2026-08-15 実測: [FrdI] 全体で
`torsion` は 3 箇所——上の 2 箇所と、無関係な `torsion subgroup` の言及 1 箇所のみ)。
したがって `torsion-free` の意味は**原文からは決まらない**。標準的な読みを採り、
その旨をここに明記する。

★mathlib の `IsAddTorsionFree` は **「`n • ·` が単射」**であり、素朴な読み
(「`n • a = 0` なら `a = 0`」)とは**別の条件**である。群では一致するが、
モノイドでは一般に別物なので、両方を置いて関係を証明する。 -/

/-- **[FrdI] §0 `torsion-free`**(素朴な読み) —— torsion 元を持たない。

★原文はこの語を定義していない。ここは標準的な読み
「`n ≥ 1` かつ `n • a = 0` なら `a = 0`」を採る。 -/
def IsTorsionFreeNaive : Prop := ∀ (a : M) (n : ℕ), 0 < n → n • a = 0 → a = 0

variable {M}

/-- ★★**原文の「hence, in particular」の中身** —— **sharp なら torsion-free**。

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if

原文 (FrdI p.11):
> M is sharp, integral, and saturated. We shall say that M is perfect if multiplication

★証明は1行で、**`sharp` だけで足りる**(`integral` も `saturated` も要らない):
`n • a = 0` かつ `n ≥ 1` なら `a + (n-1) • a = 0` なので `a` は可逆、
`sharp` より `a = 0`。 -/
theorem isTorsionFreeNaive_of_isSharp (h : IsSharp M) : IsTorsionFreeNaive M := by
  intro a n hn hna
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  refine h a ⟨⟨a, m • a, ?_, ?_⟩, rfl⟩
  · calc a + m • a = m • a + a := add_comm _ _
      _ = (m + 1) • a := (succ_nsmul a m).symm
      _ = 0 := hna
  · calc m • a + a = (m + 1) • a := (succ_nsmul a m).symm
      _ = 0 := hna

/-- mathlib の `IsAddTorsionFree`(`n • ·` が単射)なら、素朴な読みも成り立つ。

★逆は一般には成り立たない。**同じ語で別の条件**である。 -/
theorem isTorsionFreeNaive_of_isAddTorsionFree [IsAddTorsionFree M] :
    IsTorsionFreeNaive M := by
  intro a n hn hna
  have h0 : n • a = n • (0 : M) := by simpa using hna
  exact IsAddTorsionFree.nsmul_right_injective (by omega) h0

/-! ### `M^±` と `of characteristic type`

原文 (FrdI p.11):
> the submonoid [which, in fact, forms a group] of invertible elements of M, by

原文 (FrdI p.11):
> the quotient monoid of M by M ±, which we shall refer to as the characteristic of

原文 (FrdI p.11):
> characteristic type if the fibers of the natural map M →M char are torsors over M ±.

★**`M^char` そのもの(商モノイド)は作っていない。** 作るにはモノイドの合同関係
(`AddCon`)が要る。ここで要るのは「**同じファイバーに属する**」という関係だけなので、
それを直接書く——`M^± ` で生成される合同関係は
「`∃ u v ∈ M^±, a + u = b + v`」である。**この選択は原文の言い換えであって、
原文が書いていないことを足してはいない。** -/

variable (M)

/-- **`M^±`** —— 可逆元のなす部分モノイド。 -/
def unitsSubmonoid : AddSubmonoid M where
  carrier := {a : M | IsAddUnit a}
  zero_mem' := isAddUnit_zero
  add_mem' := IsAddUnit.add

variable {M}

/-- `M → M^char` の**同じファイバー**に属する関係
(= `M^±` が生成する合同関係)。 -/
def CharRel (a b : M) : Prop := ∃ u v : M, IsAddUnit u ∧ IsAddUnit v ∧ a + u = b + v

theorem charRel_refl (a : M) : CharRel a a := ⟨0, 0, isAddUnit_zero, isAddUnit_zero, rfl⟩

variable (M)

/-- **[FrdI] §0 `of characteristic type`** —— `M → M^char` のファイバーが
`M^±` 上の**トーサー**(単純推移的)であること。

原文 (FrdI p.11):
> characteristic type if the fibers of the natural map M →M char are torsors over M ±.

★「トーサー」を「推移的かつ自由」に開いて書く。 -/
def IsOfCharacteristicType : Prop :=
  ∀ a b : M, CharRel a b → ∃! u : M, IsAddUnit u ∧ b = a + u

/-! ### 非退化 -/

/-- `ℕ` は **of characteristic type** —— 可逆元が `0` だけなのでファイバーは1点。 -/
theorem isOfCharacteristicType_nat : IsOfCharacteristicType ℕ := by
  intro a b hab
  obtain ⟨u, v, hu, hv, huv⟩ := hab
  obtain ⟨w, rfl⟩ := hu
  obtain ⟨w', rfl⟩ := hv
  have hw : ((w : ℕ)) = 0 := by have := w.val_neg; omega
  have hw' : ((w' : ℕ)) = 0 := by have := w'.val_neg; omega
  refine ⟨0, ⟨isAddUnit_zero, by omega⟩, ?_⟩
  rintro y ⟨hy, hy2⟩
  obtain ⟨z, rfl⟩ := hy
  have := z.val_neg
  omega

/-- ★`WithTop (ZMod 2)` は **of characteristic type でない** ——
`⊤` のファイバー上で `M^±` の作用が**自由でない**(`⊤ + 0 = ⊤ + 1 = ⊤`)。 -/
theorem not_isOfCharacteristicType_withTop_zmod2 :
    ¬ IsOfCharacteristicType (WithTop (ZMod 2)) := by
  intro h
  have h1 : IsAddUnit ((1 : ZMod 2) : WithTop (ZMod 2)) :=
    ⟨⟨((1 : ZMod 2) : WithTop (ZMod 2)), ((1 : ZMod 2) : WithTop (ZMod 2)), by decide, by decide⟩, rfl⟩
  obtain ⟨u, hu, huniq⟩ := h ⊤ ⊤ (charRel_refl ⊤)
  have e0 : (0 : WithTop (ZMod 2)) = u := huniq 0 ⟨isAddUnit_zero, by simp⟩
  have e1 : ((1 : ZMod 2) : WithTop (ZMod 2)) = u := huniq _ ⟨h1, by simp⟩
  rw [← e0] at e1
  exact absurd e1 (by decide)

/-! ### ★[FrdI] Definition 1.1, (i) —— `pre-divisorial` / `divisorial`

原典: [FrdI] 物理 p.19(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.19):
> (i) We shall say that M ∈Ob(Mon) is pre-divisorial if it is integral [cf. §0],

原文 (FrdI p.19):
> saturated [cf. §0], and of characteristic type [cf. §0]. Suppose that M is pre-

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is

★**Definition 1.1 (i) が §0 から引く語は `integral` / `saturated` /
`of characteristic type` の3つで、いま全部揃った。**

★ただし同じ (i) が定める `group-like`(「`M^char` が零」)は**書けていない**——
`M^char` を商モノイドとして作っていないため。 -/

/-- **[FrdI] Definition 1.1, (i)** —— `pre-divisorial`。

原文 (FrdI p.19):
> (i) We shall say that M ∈Ob(Mon) is pre-divisorial if it is integral [cf. §0],
-/
def IsPreDivisorial : Prop :=
  IsIntegralMonoid M ∧ IsSaturatedMonoid M ∧ IsOfCharacteristicType M

/-- **[FrdI] Definition 1.1, (i)** —— `divisorial` = pre-divisorial かつ sharp。

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is
-/
def IsDivisorial : Prop := IsPreDivisorial M ∧ IsSharp M

variable {M}

/-- 群は **of characteristic type** —— `M^± = M` なのでファイバーは全体、
`u = b - a` が一意。 -/
theorem isOfCharacteristicType_int : IsOfCharacteristicType ℤ := by
  intro a b _
  refine ⟨b - a, ⟨⟨⟨b - a, a - b, by omega, by omega⟩, rfl⟩, by omega⟩, ?_⟩
  rintro y ⟨-, hy⟩
  omega

/-- ★`ℤ` は **pre-divisorial**(Definition 1.1 (i) を満たす実例)。 -/
theorem isPreDivisorial_int : IsPreDivisorial ℤ :=
  ⟨isIntegralMonoid_of_isCancelAdd ℤ,
   isSaturatedMonoid_of_range_eq_univ ℤ range_toGp_int,
   isOfCharacteristicType_int⟩

/-- ★`ℤ` は **divisorial ではない** —— sharp でないから。
`pre-divisorial` と `divisorial` が別物であることの確認。 -/
theorem not_isDivisorial_int : ¬ IsDivisorial ℤ := fun h => not_isSharp_int h.2

/-! ### ★出典の紐付け(`.src`)

`Skeleton/` の宣言は `.src` で原典項目に紐づいている。`Found/` にはその規約が
及んでいなかったため、**実装を原典に対して数えられなかった**(器具が 0 件と印字していた)。
ここで主要な定義に `.src` を付ける。★`Skeleton/` と**同じ検査**(論文の登記・
物理ページの範囲・`1_Structured` の `sectionId` の実在)が走る。 -/

def IsSharp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — sharp",
    sectionId := "frdi-s0-sharp-integral" }

def IsIntegralMonoid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — integral",
    sectionId := "frdi-s0-sharp-integral" }

def IsSaturatedMonoid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — saturated",
    sectionId := "frdi-s0-saturated" }

def IsOfCharacteristicType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — of characteristic type",
    sectionId := "frdi-s0-char-type" }

def MChar.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — M^char = M/M^±",
    sectionId := "frdi-s0-char-type" }

def IsCharacteristicallyInjective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — characteristically injective",
    sectionId := "frdi-s0-char-inj" }

def IsPrimaryElt.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12, item := "§0 Monoids — primary",
    sectionId := "frdi-s0-primary" }

def IsPreDivisorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (i) — pre-divisorial",
    sectionId := "frdi-def-1-1-i" }

def IsDivisorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (i) — divisorial",
    sectionId := "frdi-def-1-1-i" }

def IsGroupLike.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (i) — group-like",
    sectionId := "frdi-def-1-1-i" }

def IsPerfectMonoid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — perfect",
    sectionId := "frdi-s0-perfect" }

/-! ### ★`M^char = M/M^±` を**商モノイドとして**作る(2026-08-15 追加)

前回はここを `CharRel`(同じファイバーに属する関係)で止めた。今回 `Definition 1.1` を
閉じるために、**商モノイドそのもの**を作る。これが要るのは2箇所である:

* `Definition 1.1, (i)` の `group-like`(「`M^char` が零」)
* `Definition 1.1, (ii)` (a) の `characteristically injective`(§0)

原文 (FrdI p.11):
> the quotient monoid of M by M ±, which we shall refer to as the characteristic of
-/

variable (M)

/-- `CharRel` は**加法合同関係**である。

推移律: `a + u = b + v`、`b + u' = c + v'` から
`a + (u+u') = (a+u) + u' = (b+v) + u' = (b+u') + v = (c+v') + v = c + (v'+v)`。
-/
def charCon : AddCon M where
  r := CharRel
  iseqv :=
    { refl := charRel_refl
      symm := fun ⟨u, v, hu, hv, h⟩ => ⟨v, u, hv, hu, h.symm⟩
      trans := fun ⟨u, v, hu, hv, h⟩ ⟨u', v', hu', hv', h'⟩ =>
        ⟨u + u', v' + v, hu.add hu', hv'.add hv, by
          rw [← add_assoc, h, add_right_comm, h', add_assoc]⟩ }
  add' := fun ⟨u, v, hu, hv, h⟩ ⟨u', v', hu', hv', h'⟩ =>
    ⟨u + u', v + v', hu.add hu', hv.add hv', by
      rw [add_add_add_comm, h, h', add_add_add_comm]⟩

/-- **[FrdI] §0** —— **`M^char = M/M^±`**、商モノイドとして。 -/
abbrev MChar : Type _ := (charCon M).Quotient

variable {M}

/-- 自然な全射 `M ↠ M^char`。 -/
def toChar : M →+ MChar M := (charCon M).mk'

@[simp] theorem toChar_eq_iff {a b : M} : toChar a = toChar b ↔ CharRel a b :=
  AddCon.eq _

variable (M)

/-- **[FrdI] Definition 1.1, (i)** —— `group-like` = 「`M^char` が零」。

原文 (FrdI p.19):
> divisorial. Then we shall say that M is group-like if M char is zero; we shall say that
-/
def IsGroupLike : Prop := ∀ x : MChar M, x = 0

/-- ★`group-like` ⟺ **M が群である**(全元が可逆)。

原文の「`M^char` が零」を、`M` の側の条件に翻訳したもの。★これが
`M^char` を**商モノイドとして**作った見返りである——`CharRel` のままでは
「零である」という条件そのものが書けなかった。 -/
theorem isGroupLike_iff : IsGroupLike M ↔ ∀ a : M, IsAddUnit a := by
  constructor
  · intro h a
    have : (toChar a : MChar M) = toChar 0 := by rw [h (toChar a), map_zero]
    obtain ⟨u, v, hu, hv, huv⟩ := toChar_eq_iff.mp this
    obtain ⟨V, rfl⟩ := hv
    refine ⟨⟨a, u + V.neg, ?_, ?_⟩, rfl⟩
    · rw [← add_assoc, huv, zero_add, V.val_neg]
    · rw [add_comm, ← add_assoc, huv, zero_add, V.val_neg]
  · intro h x
    induction x using AddCon.induction_on with
    | H a =>
      obtain ⟨A, rfl⟩ := h a
      show (toChar (A : M) : MChar M) = toChar 0
      exact toChar_eq_iff.mpr ⟨A.neg, 0, ⟨-A, rfl⟩, isAddUnit_zero, by simp⟩

/-! ### 非退化(`group-like`) -/

/-- `ℤ` は **group-like**。 -/
theorem isGroupLike_int : IsGroupLike ℤ :=
  (isGroupLike_iff ℤ).mpr fun a => ⟨⟨a, -a, by omega, by omega⟩, rfl⟩

/-- ★`ℕ` は **group-like でない** —— `1` が可逆でない。 -/
theorem not_isGroupLike_nat : ¬ IsGroupLike ℕ := by
  intro h
  obtain ⟨A, hA⟩ := (isGroupLike_iff ℕ).mp h 1
  have := A.val_neg
  omega

/-! ### `ℕ≥1` の乗法モノイドの構造

★原文が `Proposition 1.7` の証明で「the well-known structure of the multiplicative
monoid `N≥1`」と呼ぶもの。使うのは**単元が 1 だけ**という一点である。 -/

/-- `ℕ+` で `a * b = 1` なら `a = 1`。 -/
theorem pnat_left_eq_one {a b : ℕ+} (h : a * b = 1) : a = 1 := by
  have h' : ((a : ℕ+) : ℕ) * ((b : ℕ+) : ℕ) = 1 := by
    exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) h
  exact PNat.coe_eq_one_iff.mp (Nat.dvd_one.mp ⟨_, h'.symm⟩)

/-- `ℕ+` で `a * b = 1` なら `b = 1`。 -/
theorem pnat_right_eq_one {a b : ℕ+} (h : a * b = 1) : b = 1 :=
  pnat_left_eq_one (by rw [mul_comm]; exact h)

/-! ### `perfect` モノイド

原文 (FrdI p.11):
> M is sharp, integral, and saturated. We shall say that M is perfect if multiplication

原文 (FrdI p.11):
> by any element of N≥1 on M is bijective. Thus, M pf is always perfect; M is perfect
-/

/-- **[FrdI] §1 (p.11)** —— `perfect` モノイド。

`ℕ≥1` の任意の元による**乗法(= `n` 倍)が全単射**。 -/
def IsPerfectMonoid : Prop := ∀ n : ℕ+, Function.Bijective (fun a : M => (n : ℕ) • a)

/-! ### ★★`M^pf` —— **perfection**(2026-08-16 追加)

原文 (FrdI p.11):
> the perfection of M, that is to say, the inductive limit of the inductive system I∗of

★★**帰納系の実体**: 添字は `(ℕ≥1, |)`、`I_a = M`、`a | b` の射は `(b/a)·`。

★★**濾過的帰納極限なので、対 `(m, a)`(＝「`m/a`」)の同値類として書ける**:
```
(m, a) ~ (m', a')  ⟺  ∃ k : ℕ+, (k*a') • m = (k*a) • m'
```
★**`k` を自由に取れるのが濾過性**であり、**推移律はそこで効く**
(中間の分母を `k` に取り込む)。★**`k` を落として `a' • m = a • m'` にすると
推移律が消去律を要求してしまう** —— 一般のモノイドでは成り立たない。

★**これが無いために `Φ^pf`(Definition 1.1)と `Frobenius-compact`
(Definition 1.2, (iv))の両方が書けなかった**(2026-08-16 の監査)。 -/

/-- perfection を定める同値関係。★`k` が濾過性を担う。 -/
def PfRel : (M × ℕ+) → (M × ℕ+) → Prop :=
  fun x y => ∃ k : ℕ+, ((k : ℕ) * (y.2 : ℕ)) • x.1 = ((k : ℕ) * (x.2 : ℕ)) • y.1

theorem pfRel_refl (x : M × ℕ+) : PfRel M x x := ⟨1, rfl⟩

theorem pfRel_symm {x y : M × ℕ+} (h : PfRel M x y) : PfRel M y x :=
  ⟨h.choose, h.choose_spec.symm⟩

theorem pfRel_trans {x y z : M × ℕ+} (h1 : PfRel M x y) (h2 : PfRel M y z) :
    PfRel M x z := by
  obtain ⟨k₁, e₁⟩ := h1
  obtain ⟨k₂, e₂⟩ := h2
  refine ⟨k₁ * k₂ * y.2, ?_⟩
  have s₁ : (((k₂ : ℕ) * (z.2 : ℕ)) * ((k₁ : ℕ) * (y.2 : ℕ))) • x.1
      = (((k₂ : ℕ) * (z.2 : ℕ)) * ((k₁ : ℕ) * (x.2 : ℕ))) • y.1 := by
    have e₁' := e₁
    simp only [mul_smul] at e₁' ⊢
    rw [e₁']
  have s₂ : (((k₁ : ℕ) * (x.2 : ℕ)) * ((k₂ : ℕ) * (z.2 : ℕ))) • y.1
      = (((k₁ : ℕ) * (x.2 : ℕ)) * ((k₂ : ℕ) * (y.2 : ℕ))) • z.1 := by
    have e₂' := e₂
    simp only [mul_smul] at e₂' ⊢
    rw [e₂']
  push_cast
  calc ((k₁ : ℕ) * (k₂ : ℕ) * (y.2 : ℕ) * (z.2 : ℕ)) • x.1
      = (((k₂ : ℕ) * (z.2 : ℕ)) * ((k₁ : ℕ) * (y.2 : ℕ))) • x.1 := by ring_nf
    _ = (((k₂ : ℕ) * (z.2 : ℕ)) * ((k₁ : ℕ) * (x.2 : ℕ))) • y.1 := s₁
    _ = (((k₁ : ℕ) * (x.2 : ℕ)) * ((k₂ : ℕ) * (z.2 : ℕ))) • y.1 := by ring_nf
    _ = (((k₁ : ℕ) * (x.2 : ℕ)) * ((k₂ : ℕ) * (y.2 : ℕ))) • z.1 := s₂
    _ = ((k₁ : ℕ) * (k₂ : ℕ) * (y.2 : ℕ) * (x.2 : ℕ)) • z.1 := by ring_nf

/-- perfection の setoid。 -/
def pfSetoid : Setoid (M × ℕ+) where
  r := PfRel M
  iseqv := ⟨pfRel_refl M, pfRel_symm M, pfRel_trans M⟩

/-- **[FrdI] §0** `M^pf` —— `M` の **perfection**。 -/
def Pf : Type _ := Quotient (pfSetoid M)

namespace Pf

variable {M}

/-- `M^pf` の元 `m/a`。 -/
def mk (m : M) (a : ℕ+) : Pf M := Quotient.mk (pfSetoid M) (m, a)

theorem sound {m m' : M} {a a' : ℕ+} (k : ℕ+)
    (h : ((k : ℕ) * (a' : ℕ)) • m = ((k : ℕ) * (a : ℕ)) • m') : mk m a = mk m' a' :=
  Quotient.sound ⟨k, h⟩

@[elab_as_elim]
theorem inductionOn {p : Pf M → Prop} (x : Pf M) (h : ∀ (m : M) (a : ℕ+), p (mk m a)) : p x :=
  Quotient.inductionOn x fun y => h y.1 y.2

/-- `m/a + m'/a' = (a'·m + a·m')/(a·a')`。 -/
instance : Add (Pf M) :=
  ⟨Quotient.map₂ (fun x y : M × ℕ+ => ((y.2 : ℕ) • x.1 + (x.2 : ℕ) • y.1, x.2 * y.2))
    (by
      rintro ⟨m₁, a₁⟩ ⟨n₁, b₁⟩ ⟨k₁, e₁⟩ ⟨m₂, a₂⟩ ⟨n₂, b₂⟩ ⟨k₂, e₂⟩
      refine ⟨k₁ * k₂, ?_⟩
      simp only [PNat.mul_coe, smul_add, smul_smul]
      have h₁ : ((k₂ : ℕ) * b₂ * a₂) • (((k₁ : ℕ) * b₁) • m₁)
          = ((k₂ : ℕ) * b₂ * a₂) • (((k₁ : ℕ) * a₁) • n₁) := by rw [e₁]
      have h₂ : ((k₁ : ℕ) * b₁ * a₁) • (((k₂ : ℕ) * b₂) • m₂)
          = ((k₁ : ℕ) * b₁ * a₁) • (((k₂ : ℕ) * a₂) • n₂) := by rw [e₂]
      simp only [smul_smul] at h₁ h₂
      calc ((k₁ : ℕ) * (k₂ : ℕ) * ((b₁ : ℕ) * (b₂ : ℕ)) * (a₂ : ℕ)) • m₁
            + ((k₁ : ℕ) * (k₂ : ℕ) * ((b₁ : ℕ) * (b₂ : ℕ)) * (a₁ : ℕ)) • m₂
          = ((k₂ : ℕ) * b₂ * a₂ * ((k₁ : ℕ) * b₁)) • m₁
            + ((k₁ : ℕ) * b₁ * a₁ * ((k₂ : ℕ) * b₂)) • m₂ := by ring_nf
        _ = ((k₂ : ℕ) * b₂ * a₂ * ((k₁ : ℕ) * a₁)) • n₁
            + ((k₁ : ℕ) * b₁ * a₁ * ((k₂ : ℕ) * a₂)) • n₂ := by rw [h₁, h₂]
        _ = ((k₁ : ℕ) * (k₂ : ℕ) * ((a₁ : ℕ) * (a₂ : ℕ)) * (b₂ : ℕ)) • n₁
            + ((k₁ : ℕ) * (k₂ : ℕ) * ((a₁ : ℕ) * (a₂ : ℕ)) * (b₁ : ℕ)) • n₂ := by ring_nf)⟩

@[simp] theorem mk_add_mk (m m' : M) (a a' : ℕ+) :
    mk m a + mk m' a' = mk ((a' : ℕ) • m + (a : ℕ) • m') (a * a') := rfl

instance : Zero (Pf M) := ⟨mk 0 1⟩

@[simp] theorem zero_def : (0 : Pf M) = mk 0 1 := rfl

instance : AddCommMonoid (Pf M) where
  add_assoc x y z := by
    induction x using Pf.inductionOn with | _ m₁ a₁ =>
    induction y using Pf.inductionOn with | _ m₂ a₂ =>
    induction z using Pf.inductionOn with | _ m₃ a₃ =>
    refine (Pf.sound 1 ?_).trans (congrArg (mk _) (mul_assoc a₁ a₂ a₃))
    push_cast
    simp only [smul_add, smul_smul, one_mul]
    ring_nf <;> abel
  zero_add x := by
    induction x using Pf.inductionOn with | _ m a =>
    refine (Pf.sound 1 ?_).trans (congrArg (mk _) (one_mul a))
    push_cast
    simp only [smul_add, smul_smul, one_mul, smul_zero, zero_add]
    ring_nf <;> abel
  add_zero x := by
    induction x using Pf.inductionOn with | _ m a =>
    refine (Pf.sound 1 ?_).trans (congrArg (mk _) (mul_one a))
    push_cast
    simp only [smul_add, smul_smul, one_mul, smul_zero, add_zero]
    ring_nf <;> abel
  add_comm x y := by
    induction x using Pf.inductionOn with | _ m₁ a₁ =>
    induction y using Pf.inductionOn with | _ m₂ a₂ =>
    refine (Pf.sound 1 ?_).trans (congrArg (mk _) (mul_comm a₁ a₂))
    push_cast
    simp only [smul_add, smul_smul, one_mul]
    ring_nf <;> abel
  nsmul n x := nsmulRec n x

@[simp] theorem nsmul_mk (n : ℕ) (m : M) (a : ℕ+) : n • mk m a = mk (n • m) a := by
  induction n with
  | zero => exact (Pf.sound 1 (by simp)).symm
  | succ n ih =>
    rw [succ_nsmul, ih, mk_add_mk]
    refine Pf.sound 1 ?_
    push_cast
    simp only [smul_add, smul_smul, succ_nsmul, one_mul]
    ring_nf <;> abel

/-- ★**自然な射 `M → M^pf`**。

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if
-/
def of : M →+ Pf M where
  toFun m := mk m 1
  map_zero' := rfl
  map_add' m m' := (Pf.sound 1 (by simp)).symm

@[simp] theorem of_apply (m : M) : (of m : Pf M) = mk m 1 := rfl

/-- ★★**`M^pf` は常に perfect**。

原文 (FrdI p.11):
> by any element of N≥1 on M is bijective. Thus, M pf is always perfect; M is perfect

★**全射**は分母に `n` を掛ければよい(`m/a = (m/(n·a)) を n 倍したもの`)。
★**単射**は `k` を `k·n` に取り替えるだけ —— ★**濾過性がここで効く。** -/
theorem isPerfectMonoid_pf : IsPerfectMonoid (Pf M) := by
  intro n
  constructor
  · intro x y h
    induction x using Pf.inductionOn with | _ m a =>
    induction y using Pf.inductionOn with | _ m' a' =>
    simp only [nsmul_mk] at h
    obtain ⟨k, e⟩ := Quotient.exact h
    refine Pf.sound (k * n) ?_
    push_cast at e ⊢
    simp only [smul_smul] at e
    calc ((k : ℕ) * (n : ℕ) * (a' : ℕ)) • m
        = ((k : ℕ) * (a' : ℕ) * (n : ℕ)) • m := by ring_nf
      _ = ((k : ℕ) * (a : ℕ) * (n : ℕ)) • m' := e
      _ = ((k : ℕ) * (n : ℕ) * (a : ℕ)) • m' := by ring_nf
  · intro x
    induction x using Pf.inductionOn with | _ m a =>
    refine ⟨mk m (n * a), ?_⟩
    simp only [nsmul_mk]
    refine Pf.sound 1 ?_
    push_cast
    simp only [smul_smul]
    ring_nf

/-! ### ★★perfection の関手性(`Φ^pf` を `MonoidOn` にするための部品)

原文 (FrdI p.19):
> is a monoid on D, then Φ determines monoids “Φchar”, “Φgp”, Φpf” on D [i.e.,

★`Φ^pf` を `monoid on D` にするには、原文 (ii) の 2 条件
(a) characteristically injective の保存、(b) FSM-morphism で同型、が要る。
★**ここでは関手性と、単射性・全射性の保存までを作る。** -/

/-- ★**perfection は関手的** —— `m/a ↦ f(m)/a`。 -/
def map {N : Type*} [AddCommMonoid N] (f : M →+ N) : Pf M →+ Pf N where
  toFun := Quotient.map (fun x : M × ℕ+ => (f x.1, x.2))
    (by
      rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, e⟩
      exact ⟨k, by simpa only [← map_nsmul] using congrArg f e⟩)
  map_zero' := by
    show mk (f 0) 1 = mk 0 1
    rw [map_zero]
  map_add' x y := by
    induction x using Pf.inductionOn with | _ m a =>
    induction y using Pf.inductionOn with | _ m' a' =>
    show mk (f ((a' : ℕ) • m + (a : ℕ) • m')) (a * a')
      = mk ((a' : ℕ) • f m + (a : ℕ) • f m') (a * a')
    rw [map_add, map_nsmul, map_nsmul]

@[simp] theorem map_mk {N : Type*} [AddCommMonoid N] (f : M →+ N) (m : M) (a : ℕ+) :
    map f (mk m a) = mk (f m) a := rfl

@[simp] theorem map_id : map (AddMonoidHom.id M) = AddMonoidHom.id (Pf M) := by
  ext x
  induction x using Pf.inductionOn with | _ m a => rfl

theorem map_comp {N O : Type*} [AddCommMonoid N] [AddCommMonoid O]
    (f : M →+ N) (g : N →+ O) : map (g.comp f) = (map g).comp (map f) := by
  ext x
  induction x using Pf.inductionOn with | _ m a => rfl

/-- ★**単射性は perfection で保たれる**。

★`f((k a')•m) = f((k a)•m')` から単射性で `(k a')•m = (k a)•m'` が出る ——
★**`k` はそのまま使える**(濾過性を壊さない)。 -/
theorem map_injective {N : Type*} [AddCommMonoid N] {f : M →+ N}
    (hf : Function.Injective f) : Function.Injective (map f) := by
  intro x y h
  induction x using Pf.inductionOn with | _ m a =>
  induction y using Pf.inductionOn with | _ m' a' =>
  obtain ⟨k, e⟩ := Quotient.exact h
  refine Pf.sound k ?_
  apply hf
  simpa only [map_nsmul] using e

/-- ★**全射性も保たれる** —— 分母はそのままでよい。 -/
theorem map_surjective {N : Type*} [AddCommMonoid N] {f : M →+ N}
    (hf : Function.Surjective f) : Function.Surjective (map f) := by
  intro y
  induction y using Pf.inductionOn with | _ n a =>
  obtain ⟨m, rfl⟩ := hf n
  exact ⟨mk m a, rfl⟩

/-! ### ★★`ℚ>0` 倍作用を**構成せずに書く**（2026-08-16）

原文 (FrdI p.23) の `Frobenius-compact` の第3節は
「every element of Aut_C(C) that acts on O×(C)pf via multiplication by an element
∈Q>0 in fact acts trivially」である。

★★**ここで `ℚ>0` を型として構成する必要はない。**
「`q = c/d` 倍として作用する」は、`ℚ` を一切使わずに

    ∀ x, (d : ℕ) • σ x = (c : ℕ) • x       （`c d : ℕ+`）

と書ける。★**`Pf M` は perfect なので `d •` は全単射であり、
この式は `σ x` を一意に決める** —— したがって両者は同値である。

★**監査は「`ℚ>0` 作用が未実装で、これが最大の未着手部分」と見積もったが、
★★実際には「作用を書く」のには不要だった。** 作るのは式であって型ではない。 -/

/-- ★**`Pf M` では `n •` は単射** —— `isPerfectMonoid_pf` の半分を取り出した形。

★「`d • y = c • x` が `y` を一意に決める」の中身である。 -/
theorem nsmul_injective (n : ℕ+) :
    Function.Injective (fun x : Pf M => ((n : ℕ+) : ℕ) • x) :=
  (isPerfectMonoid_pf (M := M) n).1

/-- ★★**「`c/d` 倍として作用する」は写像を一意に決める**。

★`σ` と `τ` がどちらも `c/d` 倍として作用するなら、両者は一致する。
★**これが「`ℚ>0` を構成せずに済む」ことの根拠である。** -/
theorem eq_of_qact {c d : ℕ+} {σ τ : Pf M → Pf M}
    (hσ : ∀ x, ((d : ℕ+) : ℕ) • σ x = ((c : ℕ+) : ℕ) • x)
    (hτ : ∀ x, ((d : ℕ+) : ℕ) • τ x = ((c : ℕ+) : ℕ) • x) : σ = τ := by
  funext x
  exact nsmul_injective d ((hσ x).trans (hτ x).symm)

/-- ★**`c = d` なら作用は自明** —— 原文の「acts trivially」の形。 -/
theorem eq_id_of_qact_one {c : ℕ+} {σ : Pf M → Pf M}
    (hσ : ∀ x, ((c : ℕ+) : ℕ) • σ x = ((c : ℕ+) : ℕ) • x) : σ = id := by
  funext x
  exact nsmul_injective c (hσ x)

end Pf

def Pf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11, item := "§0 Monoids — perfection M^pf",
    sectionId := "frdi-s0-perfect" }

/-! ### 非退化(`perfect`) -/

/-- `ℚ` は **perfect** —— `n` 倍は `1/n` 倍を逆に持つ。 -/
theorem isPerfectMonoid_rat : IsPerfectMonoid ℚ := by
  intro n
  have hn : (n : ℚ) ≠ 0 := by
    exact_mod_cast (n : ℕ+).ne_zero
  constructor
  · intro a b hab
    simp only [nsmul_eq_mul] at hab
    exact mul_left_cancel₀ hn hab
  · intro b
    refine ⟨b / (n : ℚ), ?_⟩
    show (n : ℕ) • (b / (n : ℚ)) = b
    rw [nsmul_eq_mul, mul_comm, div_mul_cancel₀ _ hn]

/-- ★`ℕ` は **perfect でない** —— `2` 倍は全射でない(`1` が像にない)。 -/
theorem not_isPerfectMonoid_nat : ¬ IsPerfectMonoid ℕ := by
  intro h
  obtain ⟨a, ha⟩ := (h 2).2 1
  simp only [PNat.val_ofNat, smul_eq_mul] at ha
  omega

/-- ★**group-like でも perfect とは限らない** —— `ℤ` は群だが `2` 倍が全射でない。

★負の対照を2つ並べる理由: `perfect` は `group-like` の帰結ではない。 -/
theorem not_isPerfectMonoid_int : ¬ IsPerfectMonoid ℤ := by
  intro h
  obtain ⟨a, ha⟩ := (h 2).2 1
  simp only [PNat.val_ofNat, nsmul_eq_mul] at ha
  push_cast at ha
  omega

/-! ### ★`M^char` の性質 —— `Proposition 1.5` の `𝔽_Φ → 𝔽_{Φ^char}` に要る

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-

★`Proposition 1.5` は `Φ` が **pre-divisorial**(sharp とは限らない)のとき、
pre-Frobenioid 構造を **`Φ^char` 経由**で入れる。`PreFrobenioid` は `Φ` が
**divisorial**(= pre-divisorial + sharp)であることを要求するので、
`Φ^char` が sharp であることが鍵になる。 -/

/-- ★**`M^char` は sharp** —— `Definition 1.1, (i)` の但し書き
「Thus, if M is pre-divisorial, then M char is divisorial」の sharp の部分。

`toChar a` が可逆なら `a` 自身が可逆で、そのとき `CharRel a 0` が成り立つ。 -/
theorem toChar_surjective : Function.Surjective (toChar : M → MChar M) := fun x => by
  induction x using AddCon.induction_on with
  | H a => exact ⟨a, rfl⟩

theorem isSharp_mChar : IsSharp (MChar M) := by
  intro x hx
  obtain ⟨a, rfl⟩ := toChar_surjective M x
  obtain ⟨X, hX⟩ := hx
  obtain ⟨b, hb⟩ := toChar_surjective M (X.neg : MChar M)
  have h0 : (toChar a : MChar M) + toChar b = 0 := by
    rw [hb, ← hX]; exact X.val_neg
  have hab : (toChar (a + b) : MChar M) = toChar 0 := by
    rw [map_add, map_zero]; exact h0
  obtain ⟨u, v, hu, hv, huv⟩ := toChar_eq_iff.mp hab
  obtain ⟨V, rfl⟩ := hv
  -- `a + (b + u + V.neg) = 0` なので `a` は可逆
  have ha : IsAddUnit a := by
    refine ⟨⟨a, b + u + V.neg, ?_, ?_⟩, rfl⟩
    · have e : a + (b + u + V.neg) = (a + b + u) + V.neg := by simp only [add_assoc]
      rw [e, huv, zero_add, V.val_neg]
    · have e : b + u + V.neg + a = (a + b + u) + V.neg := by
        simp [add_comm, add_left_comm, add_assoc]
      rw [e, huv, zero_add, V.val_neg]
  obtain ⟨A, rfl⟩ := ha
  show (toChar (A : M) : MChar M) = toChar 0
  exact toChar_eq_iff.mpr ⟨A.neg, 0, ⟨-A, rfl⟩, isAddUnit_zero, by simp⟩

/-- ★**sharp なら `CharRel` は等号** —— したがって sharp なモノイドでは
`M^char` は `M` と同じである。 -/
theorem charRel_iff_eq (h : IsSharp M) {a b : M} : CharRel a b ↔ a = b := by
  constructor
  · rintro ⟨u, v, hu, hv, huv⟩
    rw [h u hu, h v hv, add_zero, add_zero] at huv
    exact huv
  · rintro rfl
    exact charRel_refl a

/-! ### ★★`M` が pre-divisorial なら `M^char` は divisorial

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is

★原文は「Thus」と書くだけで証明を置かない。**これは原文の主張であって我々の証明ではない**ので、
1条ずつ測る。`sharp` は上の `isSharp_mChar` で済んだ。残りは
**integral / saturated / of characteristic type** の3つである。 -/

/-- ★**簡約性は `M^char` に遺伝する**。

`[a+c] = [b+c]` は `∃u v 可逆, a+c+u = b+c+v` を意味し、`M` の簡約性で `c` を消せば
`a+u = b+v`、すなわち `[a] = [b]`。 -/
theorem isCancelAdd_mChar (h : IsCancelAdd M) : IsCancelAdd (MChar M) where
  add_left_cancel x y z hyz := by
    obtain ⟨a, rfl⟩ := toChar_surjective M x
    obtain ⟨b, rfl⟩ := toChar_surjective M y
    obtain ⟨c, rfl⟩ := toChar_surjective M z
    have hyz2 : (toChar (a + b) : MChar M) = toChar (a + c) := by rw [map_add, map_add]; exact hyz
    obtain ⟨u, v, hu, hv, huv⟩ := toChar_eq_iff.mp hyz2
    refine toChar_eq_iff.mpr ⟨u, v, hu, hv, ?_⟩
    have e1 : a + b + u = (b + u) + a := by simp [add_comm, add_left_comm, add_assoc]
    have e2 : a + c + v = (c + v) + a := by simp [add_comm, add_left_comm, add_assoc]
    rw [e1, e2] at huv
    letI := h
    exact add_right_cancel huv
  add_right_cancel x y z hyz := by
    obtain ⟨a, rfl⟩ := toChar_surjective M x
    obtain ⟨b, rfl⟩ := toChar_surjective M y
    obtain ⟨c, rfl⟩ := toChar_surjective M z
    have hyz2 : (toChar (b + a) : MChar M) = toChar (c + a) := by rw [map_add, map_add]; exact hyz
    obtain ⟨u, v, hu, hv, huv⟩ := toChar_eq_iff.mp hyz2
    refine toChar_eq_iff.mpr ⟨u, v, hu, hv, ?_⟩
    have e1 : b + a + u = (b + u) + a := by simp [add_comm, add_left_comm, add_assoc]
    have e2 : c + a + v = (c + v) + a := by simp [add_comm, add_left_comm, add_assoc]
    rw [e1, e2] at huv
    letI := h
    exact add_right_cancel huv

/-- ★**1. `M` が integral なら `M^char` も integral**。 -/
theorem isIntegralMonoid_mChar (h : IsIntegralMonoid M) : IsIntegralMonoid (MChar M) :=
  letI := isCancelAdd_mChar M (isCancelAdd_of_isIntegralMonoid M h)
  isIntegralMonoid_of_isCancelAdd (MChar M)

/-- ★**3. sharp なら of characteristic type**(一般の事実)。

sharp では可逆元が `0` だけなので `CharRel` は等号(`charRel_iff_eq`)、
したがってファイバーは1点で、トーサー条件は `u = 0` の一意性に帰着する。 -/
theorem isOfCharacteristicType_of_isSharp (h : IsSharp M) : IsOfCharacteristicType M := by
  intro a b hab
  rw [charRel_iff_eq M h] at hab
  subst hab
  refine ⟨0, ⟨isAddUnit_zero, (add_zero a).symm⟩, ?_⟩
  rintro u ⟨hu, -⟩
  exact h u hu

/-- ★**3'. `M^char` は of characteristic type**(`isSharp_mChar` の系)。 -/
theorem isOfCharacteristicType_mChar : IsOfCharacteristicType (MChar M) :=
  isOfCharacteristicType_of_isSharp (MChar M) (isSharp_mChar M)

/-! ### ★§0 `characteristically injective`

`Definition 1.1, (ii)` (a) が引く語。**`M^char` を作ったので初めて書ける。**

原文 (FrdI p.11):
> φ : M1 →M2 is a morphism of Mon, then we shall say that φ is characteristically
原文 (FrdI p.11):
> injective if φ is injective, and, moreover, the morphism M char
-/

variable {M} {N : Type*} [AddCommMonoid N]

/-- `φ : M →+ N` が誘導する `M^char →+ N^char`。 -/
def charMap (φ : M →+ N) : MChar M →+ MChar N :=
  AddCon.lift _ (toChar.comp φ) <| by
    rintro a b ⟨u, v, hu, hv, huv⟩
    exact toChar_eq_iff.mpr ⟨φ u, φ v, hu.map φ, hv.map φ, by
      simpa using congrArg φ huv⟩

@[simp] theorem charMap_toChar (φ : M →+ N) (a : M) :
    charMap φ (toChar a) = toChar (φ a) := rfl

/-- **[FrdI] §0** —— `characteristically injective`。

★原文は「`φ` が単射で、**かつ** `M₁^char → M₂^char` も単射」と**2つ**要求する。
片方だけでは足りないことは下の非退化で示す。 -/
def IsCharacteristicallyInjective (φ : M →+ N) : Prop :=
  Function.Injective φ ∧ Function.Injective (charMap φ)

/-! ### 非退化(`characteristically injective`) -/

/-- `AddMonoidHom.id` は characteristically injective。 -/
theorem isCharacteristicallyInjective_id :
    IsCharacteristicallyInjective (AddMonoidHom.id M) := by
  refine ⟨fun a b h => h, ?_⟩
  intro x y h
  induction x using AddCon.induction_on with
  | H a =>
    induction y using AddCon.induction_on with
    | H b => exact toChar_eq_iff.mpr (toChar_eq_iff.mp h)

/-- ★**単射だが characteristically injective でない射がある** ——
`ℕ →+ ℤ`(包含)。`ℕ^char = ℕ`(可逆元が 0 のみ)だが `ℤ^char = 0` なので、
誘導される `ℕ^char → ℤ^char` は単射でない。

★これが「2つ要求する」ことの意味である。 -/
theorem not_isCharacteristicallyInjective_natCast :
    ¬ IsCharacteristicallyInjective ((Nat.castAddMonoidHom ℤ)) := by
  rintro ⟨-, h2⟩
  have h0 : charMap (Nat.castAddMonoidHom ℤ) (toChar (0 : ℕ)) = charMap (Nat.castAddMonoidHom ℤ) (toChar (1 : ℕ)) := by
    rw [charMap_toChar, charMap_toChar]
    exact (isGroupLike_int _).trans (isGroupLike_int _).symm
  have : (0 : ℕ) = 1 := by
    have := toChar_eq_iff.mp (h2 h0)
    obtain ⟨u, v, hu, hv, huv⟩ := this
    obtain ⟨U, rfl⟩ := hu
    obtain ⟨V, rfl⟩ := hv
    have := U.val_neg; have := V.val_neg
    omega
  exact absurd this (by decide)

/-! ### ★`ℕ` は divisorial である

`Definition 1.2` の非退化 witness に**具体的な pre-Frobenioid** が要る。
その divisor monoid として `ℕ` を使うので、`Definition 1.1, (i)` の条件を
全部確かめておく。残っていたのは `saturated` である。 -/

/-- `Gp ℕ` の元 `mk p q` が `ℕ` の像に入る ⟺ `q ≤ p`。 -/
theorem mem_range_toGp_nat (p : ℕ) (q : (⊤ : AddSubmonoid ℕ)) :
    AddLocalization.mk p q ∈ Set.range (toGp ℕ) ↔ (q : ℕ) ≤ p := by
  constructor
  · rintro ⟨m, hm⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at hm
    obtain ⟨c, hc⟩ := hm
    simp only [AddSubmonoid.coe_zero] at hc
    omega
  · intro h
    refine ⟨p - (q : ℕ), ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simp only [AddSubmonoid.coe_zero]; omega⟩

/-- ★`ℕ` は **saturated**。

`n • mk p q = mk (n·p) (n·q)` が `ℕ` の像に入るなら `n·q ≤ n·p`、
`n ≥ 1` なので `q ≤ p`、よって `mk p q` 自身も像に入る。 -/
theorem isSaturatedMonoid_nat : IsSaturatedMonoid ℕ := by
  intro a n hn
  induction a using AddLocalization.ind with
  | _ y =>
    obtain ⟨p, q⟩ := y
    intro hmem
    rw [AddLocalization.mk_nsmul] at hmem
    have h1 := (mem_range_toGp_nat _ _).mp hmem
    simp only [AddSubmonoidClass.coe_nsmul, smul_eq_mul] at h1
    exact (mem_range_toGp_nat p q).mpr (Nat.le_of_mul_le_mul_left h1 hn)

/-- ★`ℕ` は **pre-divisorial**。 -/
theorem isPreDivisorial_nat : IsPreDivisorial ℕ :=
  ⟨isIntegralMonoid_nat, isSaturatedMonoid_nat, isOfCharacteristicType_nat⟩

/-- ★`ℕ` は **divisorial** —— `Definition 1.1, (iv)` が `Φ` に課す条件を満たす
具体的なモノイド。 -/
theorem isDivisorial_nat : IsDivisorial ℕ := ⟨isPreDivisorial_nat, isSharp_nat⟩

/-! ### ★§0 の順序と `primary`(物理 p.12、**400 dpi 目視確認 2026-08-15**)

`[FrdI] Definition 1.2, (iii)` の `primary pre-step` が引く語。

原文 (FrdI p.12):
> if ∃c ∈M such that a + c = b and

原文 (FrdI p.12):
> if ∃n ∈N≥1 such that a ≤n · b. If a subset S ⊆M satisfies the property that there

原文 (FrdI p.12):
> is irreducible if any equation a = b + c in M, where b, c ∈M, implies that b = 0

原文 (FrdI p.12):
> Denote by Primary(M) the set of primary

★`primary` の定義文そのもの(「for any `M ∋ b ⪯ a`, where `b ≠ 0`, it holds that
`a ⪯ b`」)は逐語照合に掛けていない。`≠` が pdftotext では合成文字
(スラッシュ + 等号)になり、照合器が拾えないためである。**書き換えずに、
照合できない事実として記す。**
-/

/-- **[FrdI] §0** —— `a ≤ b`。 -/
def MLe (a b : M) : Prop := ∃ c : M, a + c = b

/-- **[FrdI] §0** —— `a ⪯ b`。 -/
def MPrec (a b : M) : Prop := ∃ n : ℕ, 0 < n ∧ MLe a (n • b)

/-- **[FrdI] §0** —— `irreducible`。 -/
def IsIrreducibleElt (a : M) : Prop := a ≠ 0 ∧ ∀ b c : M, a = b + c → b = 0 ∨ c = 0

/-- **[FrdI] §0** —— `primary`。 -/
def IsPrimaryElt (a : M) : Prop := a ≠ 0 ∧ ∀ b : M, b ≠ 0 → MPrec b a → MPrec a b

/-! #### 非退化 -/

/-- `1 : ℕ` は irreducible。 -/
theorem isIrreducibleElt_one_nat : IsIrreducibleElt (1 : ℕ) := ⟨one_ne_zero, fun b c h => by omega⟩

/-- ★`2 : ℕ` は **irreducible でない**(`2 = 1 + 1`)。 -/
theorem not_isIrreducibleElt_two_nat : ¬ IsIrreducibleElt (2 : ℕ) := by
  rintro ⟨-, h⟩
  rcases h 1 1 rfl with h1 | h1 <;> omega

/-- `1 : ℕ` は primary —— `ℕ` では `0` でない元は互いに `⪯` で結ばれる。 -/
theorem isPrimaryElt_one_nat : IsPrimaryElt (1 : ℕ) := by
  refine ⟨one_ne_zero, fun b hb _ => ⟨1, one_pos, b - 1, ?_⟩⟩
  simp only [one_smul]
  omega

/-- ★`(1,1) : ℕ × ℕ` は **primary でない**。

`(1,0) ⪯ (1,1)` は成り立つが、`(1,1) ⪯ (1,0)` は第2成分で破れる
(`1 + c = 0` は `ℕ` で解けない)。★`primary` が自明な条件でないことの実例。 -/
theorem not_isPrimaryElt_nat_prod : ¬ IsPrimaryElt ((1, 1) : ℕ × ℕ) := by
  rintro ⟨-, h⟩
  obtain ⟨n, -, c, hc⟩ := h (1, 0) (by simp) ⟨1, one_pos, (0, 1), by simp⟩
  have := congrArg Prod.snd hc
  simp at this

/-! ### ★★`non-dilating`(2026-08-16 追加)

原文 (FrdI p.19):
> we shall say that α is non-dilating if the endomorphism αchar of M char induced by

★**引用を選び直した記録(事故 #3 の 6 度目)**: 続く行
「α is the identity endomorphism of M char whenever αchar(a) ≼a for all primary [cf.」は
★**`≼` が抽出で落ちるため引用できない**(49/68 文字で停止)。
★**`′`・`▷` に続く 3 つ目の「目では見えるが抽出に無い文字」である。

★★**「whenever」は含意である** —— 「すべての primary `a` について
`α^char(a) ≼ a`」**ならば**「`α^char = id`」。★条件そのものではなく**含意**が
定義であることに注意(原文の語順が逆なので読み違えやすい)。

★**監査(2026-08-16)で `Definition 1.1, (i)` の未実装 1 件として挙がった。** -/

/-- **[FrdI] Definition 1.1, (i)** —— 自己準同型が `non-dilating`。 -/
def IsNonDilating (α : M →+ M) : Prop :=
  (∀ a : MChar M, IsPrimaryElt a → MPrec (charMap α a) a) →
    charMap α = AddMonoidHom.id (MChar M)

/-! #### ★非退化を両側から(我々の作法) -/

/-- ★**非退化(下)**: 恒等射は non-dilating。 -/
theorem isNonDilating_id : IsNonDilating (AddMonoidHom.id M) := by
  intro _
  ext a
  rfl

variable (M)

/-! #### 2. saturated —— ★簡約的なモノイドでは初等的な形に落ちる -/

/-- 簡約的なモノイドでは `mk x y` が `M` の像に入る ⟺ `y ≤ x`。 -/
theorem mem_range_toGp_iff (hc : IsCancelAdd M) (x : M) (y : (⊤ : AddSubmonoid M)) :
    AddLocalization.mk x y ∈ Set.range (toGp M) ↔ MLe (y : M) x := by
  letI := hc
  constructor
  · rintro ⟨m, hm⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at hm
    obtain ⟨c, hc'⟩ := hm
    simp only [AddSubmonoid.coe_zero, zero_add] at hc'
    refine ⟨m, ?_⟩
    exact add_left_cancel hc'
  · rintro ⟨w, hw⟩
    refine ⟨w, ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    refine ⟨0, ?_⟩
    simp only [AddSubmonoid.coe_zero, zero_add]
    rw [add_comm ((y : M)) w] at hw ⊢
    exact hw

/-- ★**簡約的なら `saturated` は初等的な形と同値**。

`∀ a b n>0, n•b ≤ n•a → b ≤ a`。★`Gp` を経由しないので扱いやすい。 -/
theorem isSaturatedMonoid_iff_mle (hc : IsCancelAdd M) :
    IsSaturatedMonoid M ↔ ∀ (a b : M) (n : ℕ), 0 < n → MLe (n • b) (n • a) → MLe b a := by
  constructor
  · intro hsat a b n hn hle
    have hmem : n • (AddLocalization.mk a (⟨b, trivial⟩ : (⊤ : AddSubmonoid M)))
        ∈ Set.range (toGp M) := by
      rw [AddLocalization.mk_nsmul]
      refine (mem_range_toGp_iff M hc _ _).mpr ?_
      simpa [AddSubmonoidClass.coe_nsmul] using hle
    exact (mem_range_toGp_iff M hc _ _).mp (hsat _ n hn hmem)
  · intro h X n hn
    induction X using AddLocalization.ind with
    | _ y =>
      obtain ⟨p, q⟩ := y
      intro hmem
      rw [AddLocalization.mk_nsmul] at hmem
      have h1 := (mem_range_toGp_iff M hc _ _).mp hmem
      rw [AddSubmonoidClass.coe_nsmul] at h1
      exact (mem_range_toGp_iff M hc _ _).mpr (h p (q : M) n hn h1)

/-- ★**2. `M` が integral かつ saturated なら `M^char` も saturated**。

初等形に落としてから移す: `n•[b] ≤ n•[a]` は `∃u v 可逆, n•b+z+u = n•a+v` を意味し、
`v` が可逆なので `n•b ≤ n•a` が `M` の中で成り立つ。`M` の saturated から `b ≤ a`、
それを `toChar` で送れば `[b] ≤ [a]`。 -/
theorem isSaturatedMonoid_mChar (hint : IsIntegralMonoid M) (hsat : IsSaturatedMonoid M) :
    IsSaturatedMonoid (MChar M) := by
  letI hcM := isCancelAdd_of_isIntegralMonoid M hint
  letI hcC := isCancelAdd_mChar M hcM
  rw [isSaturatedMonoid_iff_mle (MChar M) hcC]
  have hM := (isSaturatedMonoid_iff_mle M hcM).mp hsat
  intro X Y n hn hle
  obtain ⟨a, rfl⟩ := toChar_surjective M X
  obtain ⟨b, rfl⟩ := toChar_surjective M Y
  obtain ⟨Z, hZ⟩ := hle
  obtain ⟨z, rfl⟩ := toChar_surjective M Z
  -- `[n•b + z] = [n•a]` なので `CharRel (n•b + z) (n•a)`
  have hcr : CharRel (n • b + z) (n • a) := by
    refine toChar_eq_iff.mp ?_
    rw [map_add, map_nsmul, map_nsmul]
    exact hZ
  obtain ⟨u, v, hu, hv, huv⟩ := hcr
  obtain ⟨V, rfl⟩ := hv
  -- `v` が可逆なので `n•b ≤ n•a`
  have hle' : MLe (n • b) (n • a) := by
    refine ⟨z + u + V.neg, ?_⟩
    have e : n • b + (z + u + V.neg) = (n • b + z + u) + V.neg := by
      simp only [add_assoc]
    rw [e, huv, add_assoc, V.val_neg, add_zero]
  obtain ⟨w, hw⟩ := hM a b n hn hle'
  exact ⟨toChar w, by rw [← map_add, hw]⟩

/-- ★★**原文の「Thus」を証明した** —— `M` が pre-divisorial なら `M^char` は divisorial。

原文 (FrdI p.19):
> M is divisorial if M is sharp [cf. §0]. [Thus, if M is pre-divisorial, then M char is

★これが `Proposition 1.5` で `𝔽_Φ → 𝔽_{Φ^char}` を使う理由である ——
`PreFrobenioid` は `Φ` が **divisorial**(= pre-divisorial + sharp)を要求するが、
`Φ` は pre-divisorial なだけかもしれない。`Φ^char` なら**無条件に**要求を満たす。 -/
theorem isDivisorial_mChar (h : IsPreDivisorial M) : IsDivisorial (MChar M) :=
  ⟨⟨isIntegralMonoid_mChar M h.1, isSaturatedMonoid_mChar M h.1 h.2.1,
    isOfCharacteristicType_mChar M⟩, isSharp_mChar M⟩


/-! #### ★負の対照 —— `M^char` の integral / saturated は仮定が効いている

★`isSharp_mChar` と `isOfCharacteristicType_mChar` は**仮定なし**で成り立つが、
`isIntegralMonoid_mChar` と `isSaturatedMonoid_mChar` は `M` の側の仮定を使う。
**その仮定が効いていること**を実物で示す。 -/

variable {M}

/-- sharp なモノイドでは `toChar` は単射(`charRel_iff_eq` の言い換え)。 -/
theorem toChar_injective_of_isSharp (h : IsSharp M) :
    Function.Injective (toChar : M → MChar M) :=
  fun _ _ hab => (charRel_iff_eq M h).mp (toChar_eq_iff.mp hab)

variable (M)

/-- ★**`M^char` が integral であるには `M` が integral である必要がある。**

`ℕ∞` は sharp なので `toChar` は単射、しかも全射なので `ℕ∞^char` は `ℕ∞` と
「同じ」である。`ℕ∞` は integral でない(`1 + ⊤ = 2 + ⊤`)ので、
`ℕ∞^char` も integral でない。

★したがって `isIntegralMonoid_mChar` の仮定は**落とせない**。 -/
private theorem isCancelAdd_enat_of_mChar (h : IsCancelAdd (MChar ℕ∞)) :
    IsCancelAdd ℕ∞ where
  add_left_cancel a b c hbc := by
    refine toChar_injective_of_isSharp isSharp_enat ?_
    have e : (toChar (a + b) : MChar ℕ∞) = toChar (a + c) := congrArg (⇑toChar) hbc
    rw [map_add, map_add] at e
    letI := h
    exact add_left_cancel e
  add_right_cancel a b c hac := by
    refine toChar_injective_of_isSharp isSharp_enat ?_
    have e : (toChar (b + a) : MChar ℕ∞) = toChar (c + a) := congrArg (⇑toChar) hac
    rw [map_add, map_add] at e
    letI := h
    exact add_right_cancel e

theorem not_isIntegralMonoid_mChar_enat : ¬ IsIntegralMonoid (MChar ℕ∞) := fun h =>
  letI := isCancelAdd_enat_of_mChar (isCancelAdd_of_isIntegralMonoid (MChar ℕ∞) h)
  not_isIntegralMonoid_enat (isIntegralMonoid_of_isCancelAdd ℕ∞)

/-- ★**まとめ** —— `M^char` の4条件のうち、
**2つは無条件、2つは `M` の仮定が効いている**。 -/
theorem mChar_divisorial_hypotheses_are_load_bearing :
    (IsSharp (MChar ℕ∞) ∧ IsOfCharacteristicType (MChar ℕ∞)) ∧
      ¬ IsIntegralMonoid (MChar ℕ∞) :=
  ⟨⟨isSharp_mChar ℕ∞, isOfCharacteristicType_mChar ℕ∞⟩, not_isIntegralMonoid_mChar_enat⟩

/-- ★**`toChar a = 0` ⟺ `a` が可逆**。

`𝔽_Φ → 𝔽_{Φ^char}` の下では `Div φ = toChar φ.div` なので、
**isometric であることは「零因子が可逆」を意味する**——`0` であることではない。
★この差が `Proposition 1.5` と我々の模型 `wP` を分ける。 -/
theorem toChar_eq_zero_iff {A : Type*} [AddCommMonoid A] {a : A} :
    (toChar a : MChar A) = 0 ↔ IsAddUnit a := by
  constructor
  · intro h
    have h0 : (toChar a : MChar A) = toChar 0 := by rw [h, map_zero]
    obtain ⟨u, v, hu, hv, huv⟩ := toChar_eq_iff.mp h0
    obtain ⟨V, rfl⟩ := hv
    refine ⟨⟨a, u + V.neg, ?_, ?_⟩, rfl⟩
    · rw [← add_assoc, huv, zero_add, V.val_neg]
    · rw [add_comm, ← add_assoc, huv, zero_add, V.val_neg]
  · rintro ⟨A0, rfl⟩
    show (toChar (A0 : A) : MChar A) = toChar 0
    exact toChar_eq_iff.mpr ⟨A0.neg, 0, ⟨-A0, rfl⟩, isAddUnit_zero, by simp⟩


/-- ★**非退化(上)**: `ℕ` の零自己準同型は non-dilating **でない**。

★`ℕ` は sharp なので `toChar` は単射。primary 元 `1` について
`0 ≼ 1` は成り立つので仮定は満たされるが、`α^char = 0 ≠ id` である。 -/
theorem not_isNonDilating_zero_nat : ¬ IsNonDilating (0 : ℕ →+ ℕ) := by
  intro h
  have hprem : ∀ a : MChar ℕ, IsPrimaryElt a → MPrec (charMap (0 : ℕ →+ ℕ) a) a := by
    intro a _
    refine ⟨1, one_pos, ?_⟩
    induction a using Quotient.inductionOn with | _ x =>
    refine ⟨toChar x, ?_⟩
    show (charMap (0 : ℕ →+ ℕ)) (toChar x) + toChar x = (1 : ℕ) • toChar x
    rw [charMap_toChar]
    simp
  have hid := h hprem
  have h1 : charMap (0 : ℕ →+ ℕ) (toChar (1 : ℕ)) = toChar (1 : ℕ) := by
    rw [hid]; rfl
  rw [charMap_toChar] at h1
  exact one_ne_zero (toChar_injective_of_isSharp isSharp_nat h1.symm)

/-! #### ★`charMap` の単射性・全射性 —— `Φ^char` を `MonoidOn` にするための部品 -/

/-- `charMap` は全射性を保つ。 -/
theorem charMap_surjective {A B : Type*} [AddCommMonoid A] [AddCommMonoid B]
    {g : A →+ B} (hg : Function.Surjective g) : Function.Surjective (charMap g) := by
  intro y
  obtain ⟨b, rfl⟩ := toChar_surjective B y
  obtain ⟨a, rfl⟩ := hg b
  exact ⟨toChar a, charMap_toChar g a⟩

/-- ★**行き先が sharp なら、`charMap` は単射性を保つ**。

`M^char` はつねに sharp(`isSharp_mChar`)なので、これは
`(M^char)^char` の段でそのまま使える。 -/
theorem charMap_injective_of_sharp {A B : Type*} [AddCommMonoid A] [AddCommMonoid B]
    (hB : IsSharp B) {g : A →+ B} (hg : Function.Injective g) :
    Function.Injective (charMap g) := by
  intro x y hxy
  obtain ⟨a, rfl⟩ := toChar_surjective A x
  obtain ⟨b, rfl⟩ := toChar_surjective A y
  rw [charMap_toChar, charMap_toChar] at hxy
  exact congrArg (⇑toChar) (hg (toChar_injective_of_isSharp hB hxy))

/-- ★**`M^char` の間の射は、単射なら characteristically injective**。

`Definition 1.1, (ii), (a)` を `Φ^char` について確かめるのに使う。 -/
theorem isCharacteristicallyInjective_of_injective_mChar
    {A B : Type*} [AddCommMonoid A] [AddCommMonoid B]
    {g : MChar A →+ MChar B} (hg : Function.Injective g) :
    IsCharacteristicallyInjective g :=
  ⟨hg, charMap_injective_of_sharp (isSharp_mChar B) hg⟩

/-! ### ★`≤` の反対称性 —— `Order(M)` が半順序になる条件

原文 (FrdI p.12):
> observe that the relation “≤” on elements of M determines a category

★原文は `Order(M)` を**前順序**の圏として作るだけで、半順序とは言わない。
実際に反対称性が要るかは使う側の問題である。ここでは
**divisorial(= integral かつ sharp)なら反対称的**であることと、
**どちらの仮定も落とせない**ことを示す。
-/

/-- ★★**`M` が integral かつ sharp なら `≤` は反対称的**。

★**2つの仮定がそれぞれ別の役割を持つ**:
`a + c = b`, `b + d = a` から `a + (c + d) = a + 0` を作り、
**integral(= 消約律)**で `c + d = 0`、すなわち `c` は可逆、
**sharp** で `c = 0`。 -/
theorem mle_antisymm {M : Type*} [AddCommMonoid M] (hint : IsIntegralMonoid M)
    (hsh : IsSharp M) {a b : M} (h1 : MLe a b) (h2 : MLe b a) : a = b := by
  obtain ⟨c, hc⟩ := h1
  obtain ⟨d, hd⟩ := h2
  letI := isCancelAdd_of_isIntegralMonoid M hint
  have hstab : a + (c + d) = a + 0 := by
    rw [← add_assoc, hc, hd, add_zero]
  have hcd : c + d = 0 := add_left_cancel hstab
  have hc0 : c = 0 := hsh c ⟨⟨c, d, hcd, by rw [add_comm]; exact hcd⟩, rfl⟩
  rw [← hc, hc0, add_zero]

/-- ★**負の対照 (1): `sharp` は落とせない** —— `ℤ` は integral だが sharp でなく、
`0 ≤ 1` かつ `1 ≤ 0` なのに `0 ≠ 1`。 -/
theorem not_mle_antisymm_int : MLe (0 : ℤ) 1 ∧ MLe (1 : ℤ) 0 ∧ (0 : ℤ) ≠ 1 :=
  ⟨⟨1, by ring⟩, ⟨-1, by ring⟩, by decide⟩

/-! #### ★負の対照 (2) 用の模型 —— `ℕ/(n ∼ n+2, n ≥ 1)`

3元 `{0, o, e}`(`o` = 奇数の類、`e` = 正の偶数の類)。
`o + o = e`、`o + e = o`、`e + e = e`。 -/

/-- `ℕ` を「0 / 正の奇数 / 正の偶数」で潰したモノイド。 -/
inductive Par : Type
  | zero : Par
  | odd : Par
  | even : Par
  deriving DecidableEq, Repr

namespace Par

instance : Fintype Par where
  elems := {Par.zero, Par.odd, Par.even}
  complete := by intro x; cases x <;> simp

instance : Zero Par := ⟨Par.zero⟩

/-- `Par` の加法。 -/
def add : Par → Par → Par
  | zero, y => y
  | x, zero => x
  | odd, odd => even
  | odd, even => odd
  | even, odd => odd
  | even, even => even

instance : Add Par := ⟨Par.add⟩

instance : AddCommMonoid Par where
  add_assoc := by decide
  zero_add := by decide
  add_zero := by decide
  add_comm := by decide
  nsmul := nsmulRec

/-- `Par` では `v + w = 0` なら `v = 0`。 -/
theorem add_eq_zero_left : ∀ v w : Par, v + w = 0 → v = 0 := by decide

end Par

/-- ★`Par` は **sharp** —— `0` 以外に可逆元は無い。 -/
theorem isSharp_par : IsSharp Par := by
  rintro a ⟨u, rfl⟩
  exact Par.add_eq_zero_left _ _ u.val_neg

/-- ★`Par` は **integral でない** —— `o + o = e + e` だが `o ≠ e`。 -/
theorem not_isIntegralMonoid_par : ¬ IsIntegralMonoid Par := by
  intro h
  letI := isCancelAdd_of_isIntegralMonoid Par h
  have : Par.zero = Par.even :=
    add_left_cancel (a := Par.odd) (show Par.odd + Par.zero = Par.odd + Par.even from by decide)
  exact absurd this (by decide)

/-- ★★**負の対照 (2): `integral` も落とせない** —— `Par` は sharp だが integral でなく、
`o ≤ e`(`o + o = e`)かつ `e ≤ o`(`e + o = o`)なのに `o ≠ e`。

★これで `mle_antisymm` の2つの仮定が**どちらも本質的**であることが確定する。 -/
theorem not_mle_antisymm_par :
    IsSharp Par ∧ MLe Par.odd Par.even ∧ MLe Par.even Par.odd ∧ Par.odd ≠ Par.even :=
  ⟨isSharp_par, ⟨Par.odd, by decide⟩, ⟨Par.odd, by decide⟩, by decide⟩

end ABC3.Found.FrdI
