import Mathlib.GroupTheory.MonoidLocalization.Basic
import Mathlib.GroupTheory.Subgroup.Saturated
import Mathlib.Data.ENat.Basic
import Mathlib.Algebra.Order.Group.Nat
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

end ABC3.Found.FrdI
