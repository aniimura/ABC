import Mathlib.CategoryTheory.Category.Preorder
import Mathlib.CategoryTheory.IsConnected
import Mathlib.CategoryTheory.Types.Basic
import Mathlib.CategoryTheory.Endomorphism
import Mathlib.CategoryTheory.MorphismProperty.Basic
import Mathlib.CategoryTheory.Comma.Basic
import Mathlib.CategoryTheory.ObjectProperty.FullSubcategory
import ABC3.Meta.Claim

/-!
# [FrdI] §0 の圏の語彙 —— `totally epimorphic` / `fiberwise-surjective` / `FSM-morphism`

原典: S. Mochizuki, *The Geometry of Frobenioids I: The General Theory* [FrdI]、
物理 p.14 と p.15(**400 dpi 目視確認 2026-08-15**)。

## なぜここを作るか

`[FrdI] Definition 1.1, (ii)/(iv)` は圏の側の語彙を §0 から引く:

* `(ii)` (b): 「if α is an **FSM-morphism** [cf. §0] of 𝒟, then α* is an isomorphism」
* `(iv)`: 「𝒞, 𝒟 are connected, **totally epimorphic** categories [cf. §0]」

★前回作った `MonoidVocabulary.lean` は `Definition 1.1, (i)`(モノイドの側)だけを
閉じたものである。**`Definition 1.2` が依拠するのは (ii)(iii)(iv)、すなわち圏の側**であり、
そこはまだ何も無かった。

## ★mathlib の実測(2026-08-15)

| 原文の語 | mathlib | 判定 |
|---|---|---|
| `totally epimorphic` | **0 件** | 無い。自前で書く(1行) |
| `fiberwise-surjective` | **0 件** | 無い。自前で書く |
| `FSM-morphism` | **0 件** | 無い |
| epimorphism / monomorphism | `CategoryTheory.Epi` / `.Mono` | ★**使う** |
| 前順序を圏と見る | `Preorder.smallCategory`、`homOfLE` / `leOfHom` | ★**使う** |
| `Type` の epi ⟺ 全射 | `CategoryTheory.epi_iff_surjective` | ★**使う** |

★語で 0 件だったので、S2(概念)でも探した。`fiberwise-surjective` は
「引き戻しが空でない」という形の条件だが、mathlib の `Limits.HasPullback` 系は
**極限の存在**を要求するので別物である(原文の条件は極限を要求しない)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u

variable (C : Type u) [Category.{v} C]

/-! ### 定義 -/

/-- **[FrdI] §0 `totally epimorphic`** —— 圏のすべての射が epimorphism。

原文 (FrdI p.15):
> nected objects are always quasi-connected. We shall say that a category C is totally

原文 (FrdI p.15):
> (respectively, almost totally) epimorphic if every morphism in C whose domain is

原文 (FrdI p.15):
> arbitrary (respectively, nonempty) and whose codomain is arbitrary (respectively,

原文 (FrdI p.15):
> connected) is an epimorphism.

★原文は括弧で `almost totally` と二重に述べるが、`totally` の側は
「domain も codomain も**任意**」なので、条件は「すべての射が epi」に尽きる。
-/
def IsTotallyEpimorphic : Prop := ∀ (A B : C) (f : A ⟶ B), Epi f

variable {C}

/-- **[FrdI] §0 `fiberwise-surjective`**。

原文 (FrdI p.14):
> We shall say that an arrow β : B →A of a category C is fiberwise-surjective

原文 (FrdI p.14):
> if, for every arrow γ : C →A of C, there exist arrows δB : D →B, δC : D →C

原文 (FrdI p.14):
> such that β ◦δB = γ ◦δC. An arrow of a category which is a fiberwise-surjective

★原文の `β ◦ δB` は Lean の `≫` では `δB ≫ β` である(合成の向きが逆)。
-/
def IsFiberwiseSurjective {B A : C} (β : B ⟶ A) : Prop :=
  ∀ ⦃Z : C⦄ (γ : Z ⟶ A), ∃ (D : C) (δB : D ⟶ B) (δZ : D ⟶ Z), δB ≫ β = δZ ≫ γ

/-- **[FrdI] §0 `FSM-morphism`** —— fiberwise-surjective な monomorphism。

原文 (FrdI p.14):
> monomorphism will be referred to as an FSM-morphism. One verifies immediately
-/
def IsFSMMorphism {B A : C} (β : B ⟶ A) : Prop := IsFiberwiseSurjective β ∧ Mono β

variable (C)

/-- **[FrdI] §0 `category of FSM-type`** —— FSM-morphism がすべて同型。

原文 (FrdI p.14):
> which satisfies the property that every FSM-morphism of C is, in fact, an isomor-

原文 (FrdI p.14):
> phism will be referred to as a category of FSM-type.
-/
def IsOfFSMType : Prop := ∀ (B A : C) (β : B ⟶ A), IsFSMMorphism β → IsIso β

/-! ### ★原文が主張していることを証明する

原文 (FrdI p.14):
> that every composite of FSM-morphisms is again an FSM-morphism. A category C

★原文は「One verifies immediately that …」と書くだけで証明を置かない。
**それを実際に示す**のがここの目的である。
-/

variable {C}

/-- 恒等射は fiberwise-surjective。 -/
theorem isFiberwiseSurjective_id (A : C) : IsFiberwiseSurjective (𝟙 A) :=
  fun _ γ => ⟨_, γ, 𝟙 _, by simp⟩

/-- fiberwise-surjective は合成で閉じる。

`γ : Z ⟶ A'` に対し、まず `β'` の性質で `D₁, δA, δZ` を取り、
その `δA : D₁ ⟶ A` に `β` の性質を当てて `D₂, δB, δD₁` を取る。すると
`δB ≫ (β ≫ β') = (δD₁ ≫ δZ) ≫ γ`。 -/
theorem IsFiberwiseSurjective.comp {B A A' : C} {β : B ⟶ A} {β' : A ⟶ A'}
    (hβ : IsFiberwiseSurjective β) (hβ' : IsFiberwiseSurjective β') :
    IsFiberwiseSurjective (β ≫ β') := by
  intro Z γ
  obtain ⟨D₁, δA, δZ, h₁⟩ := hβ' γ
  obtain ⟨D₂, δB, δD₁, h₂⟩ := hβ δA
  refine ⟨D₂, δB, δD₁ ≫ δZ, ?_⟩
  calc δB ≫ β ≫ β' = (δB ≫ β) ≫ β' := by simp
    _ = (δD₁ ≫ δA) ≫ β' := by rw [h₂]
    _ = δD₁ ≫ (δA ≫ β') := by simp
    _ = δD₁ ≫ (δZ ≫ γ) := by rw [h₁]
    _ = (δD₁ ≫ δZ) ≫ γ := by simp

/-- ★**原文の主張** —— FSM-morphism の合成は FSM-morphism。

原文 (FrdI p.14):
> that every composite of FSM-morphisms is again an FSM-morphism. A category C
-/
theorem IsFSMMorphism.comp {B A A' : C} {β : B ⟶ A} {β' : A ⟶ A'}
    (hβ : IsFSMMorphism β) (hβ' : IsFSMMorphism β') : IsFSMMorphism (β ≫ β') :=
  ⟨hβ.1.comp hβ'.1, @mono_comp _ _ _ _ _ _ hβ.2 _ hβ'.2⟩

/-! ### ★非退化 —— 各語に「満たす例」と「満たさない例」 -/

/-- hom 集合が subsingleton な圏(= 前順序の圏)は totally epimorphic。 -/
theorem isTotallyEpimorphic_of_subsingleton_hom
    (hs : ∀ X Y : C, Subsingleton (X ⟶ Y)) : IsTotallyEpimorphic C :=
  fun _ _ _ => ⟨fun g h _ => (hs _ _).elim g h⟩

/-- `Type` の中の非全射な射(`PUnit → Bool`、値は `true` のみ)。 -/
private def punitToBool : PUnit ⟶ Bool := TypeCat.ofHom fun _ => true

/-- ★`Type` は **totally epimorphic でない** —— `PUnit → Bool` は全射でないので epi でない。 -/
theorem not_isTotallyEpimorphic_type : ¬ IsTotallyEpimorphic (Type) := by
  intro h
  have he : Epi punitToBool := h _ _ _
  rw [CategoryTheory.epi_iff_surjective] at he
  obtain ⟨_, hx⟩ := he false
  exact Bool.noConfusion hx

/-- ★★**始対象を持つ圏では「すべての射」が fiberwise-surjective である。**

始対象を `D` に取ると、条件 `δB ≫ β = δZ ≫ γ` が**始対象からの射の一意性**だけで
満たされてしまうからである。

★これは原文の定義についての測定結果である: `fiberwise-surjective` は
**始対象を持つ圏では退化する**。だから原文は `Definition 1.1, (iv)` で
「**connected**, totally epimorphic categories」と、圏の側に条件を課している。

★ここでは `Limits` を import せず、始対象の性質(全対象への射があり、それが一意)を
仮定としてそのまま書く。 -/
theorem isFiberwiseSurjective_of_initial (I : C) (arr : ∀ X : C, I ⟶ X)
    (huniq : ∀ (X : C) (f g : I ⟶ X), f = g) {B A : C} (β : B ⟶ A) :
    IsFiberwiseSurjective β :=
  fun _ _ => ⟨I, arr _, arr _, huniq _ _ _⟩

/-- ★`Type` では `PEmpty` が始対象なので、**すべての射**が fiberwise-surjective。 -/
theorem isFiberwiseSurjective_of_type {X Y : Type} (f : X ⟶ Y) :
    IsFiberwiseSurjective f :=
  isFiberwiseSurjective_of_initial PEmpty (fun _ => TypeCat.ofHom fun e => e.elim)
    (fun _ _ _ => by ext e; exact e.elim) f

/-! #### ★V 字の半順序 —— 下限を持たない圏

`left` と `right` はともに `top` の下にあるが、互いに比較不能で、
**共通の下界を持たない**。これが `fiberwise-surjective` を破る最小の形である。 -/

/-- 3点の半順序 `left ≤ top`、`right ≤ top`、`left` と `right` は比較不能。 -/
inductive Vee : Type
  | left : Vee
  | right : Vee
  | top : Vee
  deriving DecidableEq

instance : Preorder Vee where
  le a b := a = b ∨ b = Vee.top
  le_refl _ := Or.inl rfl
  le_trans a b c hab hbc := by
    rcases hab with rfl | rfl
    · exact hbc
    · rcases hbc with rfl | h
      · exact Or.inr rfl
      · exact Or.inr h

/-- `Vee` は totally epimorphic(前順序の圏だから)。 -/
theorem isTotallyEpimorphic_vee : IsTotallyEpimorphic Vee :=
  isTotallyEpimorphic_of_subsingleton_hom fun X Y => Preorder.subsingleton_hom X Y

/-- ★**`Vee` は connected**(2026-08-16 に追加)。

★`left → top ← right` の V 字なので、`top` から両側へ辺がある。
★原文 `Definition 1.1, (iv)` の `connected` を witness が満たすことに要る。 -/
instance isConnected_vee : IsConnected Vee := by
  refine IsConnected.of_induct (j₀ := Vee.top) ?_
  intro p hp0 hstep j
  have hl : (Vee.left : Vee) ∈ p :=
    (hstep (homOfLE (show Vee.left ≤ Vee.top from Or.inr rfl))).mpr hp0
  have hr : (Vee.right : Vee) ∈ p :=
    (hstep (homOfLE (show Vee.right ≤ Vee.top from Or.inr rfl))).mpr hp0
  cases j
  · exact hl
  · exact hr
  · exact hp0

/-- `Vee` で `top` へ向かう射のうち、始点が `top` でないものは
**fiberwise-surjective でない**。

`X` と別の `Z`(ともに `top` でない)を取ると、`D ≤ X` かつ `D ≤ Z` は
`D = X` かつ `D = Z` を強いるので矛盾する。 -/
theorem not_isFiberwiseSurjective_vee (X Z : Vee) (hXZ : X ≠ Z)
    (hX : X ≠ Vee.top) (hZ : Z ≠ Vee.top) (β : X ⟶ Vee.top) :
    ¬ IsFiberwiseSurjective β := by
  intro h
  obtain ⟨D, δB, δZ, -⟩ := h (homOfLE (show Z ≤ Vee.top from Or.inr rfl))
  rcases leOfHom δB with h1 | h1
  · rcases leOfHom δZ with h2 | h2
    · exact hXZ (h1.symm.trans h2)
    · exact hZ h2
  · exact hX h1

/-- ★`Vee` の `left ⟶ top` は **fiberwise-surjective でない**(具体例)。 -/
theorem not_isFiberwiseSurjective_vee_left :
    ¬ IsFiberwiseSurjective (homOfLE (show Vee.left ≤ Vee.top from Or.inr rfl)) :=
  not_isFiberwiseSurjective_vee Vee.left Vee.right (by decide) (by decide) (by decide) _

/-- `Vee` の自己射はすべて同型(hom が subsingleton だから)。 -/
theorem isIso_vee_self {X : Vee} (f : X ⟶ X) : IsIso f :=
  ⟨⟨𝟙 X, Subsingleton.elim _ _, Subsingleton.elim _ _⟩⟩

/-- ★`Vee` は **of FSM-type** —— FSM-morphism は恒等射しかない。 -/
theorem isOfFSMType_vee : IsOfFSMType Vee := by
  intro X Y β hβ
  rcases leOfHom β with rfl | rfl
  · exact isIso_vee_self β
  · obtain (rfl | rfl | rfl) : X = Vee.left ∨ X = Vee.right ∨ X = Vee.top := by
      cases X
      · exact Or.inl rfl
      · exact Or.inr (Or.inl rfl)
      · exact Or.inr (Or.inr rfl)
    · exact absurd hβ.1
        (not_isFiberwiseSurjective_vee Vee.left Vee.right (by decide) (by decide) (by decide) β)
    · exact absurd hβ.1
        (not_isFiberwiseSurjective_vee Vee.right Vee.left (by decide) (by decide) (by decide) β)
    · exact isIso_vee_self β

/-! ### ★★退化の射程 —— 「始対象があると `fiberwise-surjective` が潰れる」の帰結

`isFiberwiseSurjective_of_initial` は「始対象があれば全射が fiberwise-surjective」
と言う。すると `FSM-morphism` は単なる monomorphism に潰れ、
`of FSM-type` は「すべての mono が同型」という非常に強い条件になる。

**では原文の設定でそれが起きるか。** 測定した結果は次のとおりである。

★**§0 は圏についての `connected` を定義している**(物理 p.16、400 dpi 目視確認)。

原文 (FrdI p.16):
> as a connected component of C. In particular, we shall say that C is connected if

そして**それはグラフの連結性であって、始対象を排除しない**——始対象を持つ圏は
つねにグラフとして連結である。したがって `Definition 1.1, (iv)` の
「connected, totally epimorphic categories」の `connected` では退化を防げない。

★★**防いでいるのは `totally epimorphic` の側である。** 下の
`subsingleton_hom_of_isTotallyEpimorphic_of_initial` が示すとおり、
**totally epimorphic な圏が始対象を持てば、その圏は前順序になる**
(hom がすべて subsingleton)。したがって

* `𝒟` が前順序でない(= 実際の設定)なら、`𝒟` に始対象は無く、退化は起きない
* `𝒟` が始対象を持つなら、`𝒟` は前順序であり、Frobenioid 論は自明化する

★**原文はこの条件を明示していない**(「𝒟 は始対象を持たない」とは書いていない)。
読者が `totally epimorphic` から導く必要がある。ここではその導出を定理として置く。
-/

/-- ★★**totally epimorphic な圏が始対象を持てば、hom はすべて subsingleton**
(= その圏は前順序である)。

始対象からの一意な射 `i : I ⟶ X` は epi である。`f g : X ⟶ Y` に対し
`i ≫ f` と `i ≫ g` はどちらも `I ⟶ Y` なので始対象の一意性から等しく、
`i` が epi なので `f = g`。 -/
theorem subsingleton_hom_of_isTotallyEpimorphic_of_initial
    (h : IsTotallyEpimorphic C) (I : C) (arr : ∀ X : C, I ⟶ X)
    (huniq : ∀ (X : C) (f g : I ⟶ X), f = g) (X Y : C) : Subsingleton (X ⟶ Y) := by
  refine ⟨fun f g => ?_⟩
  haveI : Epi (arr X) := h _ _ _
  exact (cancel_epi (arr X)).mp (huniq Y _ _)

/-- ★**原文が「in passing」と述べる主張** —— totally epimorphic な圏では、
`α ∘ β` が同型なら `α` と `β` も同型。

原文 (FrdI p.16):
> We observe in passing that if C is a totally epimorphic category, and α ◦β

原文 (FrdI p.16):
> [where α, β ∈Arr(C)] is an isomorphism, then α, β are isomorphisms.

★原文は証明を置かない。**それを実際に示す。**
`β ≫ (α ≫ (β ≫ α)⁻¹) = 𝟙` なので `β` は split mono、さらに `β` が epi なので
逆向きの合成も `𝟙` になり `β` は同型。`α` はその後 `β⁻¹ ≫ (β ≫ α)` として同型。 -/
theorem isIso_of_isIso_comp (h : IsTotallyEpimorphic C) {A B E : C}
    (β : A ⟶ B) (α : B ⟶ E) [IsIso (β ≫ α)] : IsIso β ∧ IsIso α := by
  haveI : Epi β := h _ _ _
  have hβ : IsIso β := by
    refine ⟨α ≫ inv (β ≫ α), ?_, ?_⟩
    · rw [← Category.assoc, IsIso.hom_inv_id]
    · refine (cancel_epi β).mp ?_
      rw [← Category.assoc, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp,
        Category.comp_id]
  refine ⟨hβ, ?_⟩
  haveI := hβ
  have : α = inv β ≫ (β ≫ α) := by rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  rw [this]
  infer_instance

/-- ★`Type` は **of FSM-type でない** —— `PUnit → Bool` は
単射(mono)で、かつ `Type` では自動的に fiberwise-surjective なので FSM-morphism だが、
全単射でないので同型でない。

★`isFiberwiseSurjective_of_type` の退化がここに効いている。 -/
theorem not_isOfFSMType_type : ¬ IsOfFSMType (Type) := by
  intro h
  have hm : Mono punitToBool := by
    rw [CategoryTheory.mono_iff_injective]
    intro a b _
    exact Subsingleton.elim a b
  have := h _ _ _ ⟨isFiberwiseSurjective_of_type _, hm⟩
  rw [CategoryTheory.isIso_iff_bijective] at this
  obtain ⟨_, hx⟩ := this.2 false
  exact Bool.noConfusion hx

/-! ### §0 —— `sub-automorphism`(原文 p.14)

`Definition 1.2, (iv)` の `Aut^sub-ample` が土台に要求する語である。

原文 (FrdI p.14):
> shall say that an endomorphism α ∈EndC(A) of C is a sub-automorphism if there

原文 (FrdI p.14):
> exists an arrow φ : B →A of C and an automorphism β ∈AutC(B) such that

原文 (FrdI p.14):
> φ ◦β = α ◦φ; write
-/

/-- **§0** `sub-automorphism` —— 自己射 `α ∈ End_𝒞(A)` であって、
ある射 `φ : B ⟶ A` と自己同型 `β ∈ Aut_𝒞(B)` が `φ ∘ β = α ∘ φ` を満たすもの。

★原文の `φ ◦ β = α ◦ φ` は「先に `β`」なので、Lean では `β ≫ φ = φ ≫ α`。 -/
def IsSubAutomorphism {A : C} (α : End A) : Prop :=
  ∃ (B : C) (φ : B ⟶ A) (β : End B), IsIso β ∧ (β ≫ φ : B ⟶ A) = φ ≫ (α : A ⟶ A)

/-- **§0** `Aut^sub_𝒞(A)` —— sub-automorphism の集合。

★原文は `(Aut_𝒞(A) ⊆) Aut^sub_𝒞(A) ⊆ End_𝒞(A)` と書き、**部分集合**としか言わない
(部分モノイドとは主張していない)ので、`Set` として写す。

原文 (FrdI p.14):
> for the subset of EndC(A) determined by the sub-automorphisms of A. We shall say
-/
def SubAut (A : C) : Set (End A) := {α | IsSubAutomorphism α}

@[simp] theorem mem_subAut {A : C} {α : End A} :
    α ∈ SubAut A ↔ IsSubAutomorphism α := Iff.rfl

/-- ★原文の括弧書き `(Aut_𝒞(A) ⊆)` —— **自己同型は sub-automorphism**。

`φ := 𝟙_A`、`β := α` に取ればよい。 -/
theorem isSubAutomorphism_of_isIso {A : C} (α : End A) (h : IsIso α) :
    IsSubAutomorphism α :=
  ⟨A, 𝟙 A, α, h, by simp⟩

/-- ★★**始対象があると `Aut^sub` は `End` 全体に潰れる**。

`φ : B ⟶ A` を始対象からの一意な射、`β := 𝟙_B` に取ると、
`β ≫ φ` と `φ ≫ α` はともに始対象から出る射なので**必ず一致する**。

★これは原文の定義の**退化条件**である。§0 の `totally epimorphic` が
始対象を持つと前順序に潰れる(`subsingleton_hom_of_isTotallyEpimorphic_of_initial`)
のと**同じ形の退化**であり、`Aut^sub-ample` を非自明に使うには
始対象を持たない圏が要る。

★`isFiberwiseSurjective_of_initial` と同じく、`Limits` を import せず
始対象の性質を仮定としてそのまま書く。 -/
theorem isSubAutomorphism_of_initial {A : C} (I : C) (arr : ∀ X : C, I ⟶ X)
    (huniq : ∀ (X : C) (f g : I ⟶ X), f = g) (α : End A) :
    IsSubAutomorphism α :=
  ⟨I, arr A, 𝟙 I, inferInstance, huniq _ _ _⟩

/-! ### 非退化(`sub-automorphism`)

★**満たす例**は `isSubAutomorphism_of_isIso`(自己同型)。
★**自己同型でない sub-automorphism** は下の `isSubAutomorphism_of_initial` を
`Type` に当てたもの —— `PEmpty` が始対象なので、**定数写像すら** sub-automorphism になる。
★**満たさない例**は `ElementaryFrobenioid.lean` の
`not_isSubAutomorphism_of_deg_ne_one` —— `Vee` のような前順序圏では
`End A = {𝟙}` で反例が作れず、`𝔽_Φ` の Frobenius 次数 `≥ 2` の自己射を要した。 -/

/-- ★**自己同型でない sub-automorphism が実在する** —— `Type` の定数写像
`Bool → Bool`。`PEmpty` が始対象であることが効いている。 -/
def constTrue : (Bool : Type) ⟶ Bool := TypeCat.ofHom fun _ => true

theorem isSubAutomorphism_constTrue : IsSubAutomorphism constTrue :=
  isSubAutomorphism_of_initial PEmpty (fun _ => TypeCat.ofHom fun e => e.elim)
    (fun _ _ _ => by ext e; exact e.elim) _

/-- ★その定数写像は**自己同型ではない** —— 全射でない。 -/
theorem not_isIso_constTrue : ¬ IsIso constTrue := by
  intro h
  rw [CategoryTheory.isIso_iff_bijective] at h
  obtain ⟨_, hx⟩ := h.2 false
  exact Bool.noConfusion hx

/-! ### §0 —— `minimal-adjoint` / `minimal-coadjoint` / `mid-adjoint`(原文 p.17)

`Proposition 1.7, (ii)(iii)(iv)` がこの3語で述べられている。

原文 (FrdI p.17):
> Let C be a category; S a collection of arrows in C; φ ∈Arr(C). Then we shall say

原文 (FrdI p.17):
> that φ is minimal-adjoint to S (respectively, minimal-coadjoint to S; mid-adjoint

原文 (FrdI p.17):
> to S) if every factorization φ = α ◦β (respectively, φ = β ◦α; φ = α ◦β ◦γ) of φ in

原文 (FrdI p.17):
> C such that β lies in S satisfies the property that β is, in fact, an isomorphism. If φ

★原文の「a collection of arrows in `C`」は mathlib の `MorphismProperty C` である。
★合成の向き: `α ◦ β` は「先に `β`」なので `φ = β ≫ α`。
したがって **minimal-adjoint では `S` の元が最初の因子**、
**minimal-coadjoint では最後の因子**、**mid-adjoint では真ん中**である。 -/

/-- **§0** `minimal-adjoint to S` —— `φ = β ≫ α`(`β ∈ S`)なる分解では `β` が同型。 -/
def IsMinimalAdjoint (S : MorphismProperty C) {A B : C} (φ : A ⟶ B) : Prop :=
  ∀ (X : C) (β : A ⟶ X) (α : X ⟶ B), φ = β ≫ α → S β → IsIso β

/-- **§0** `minimal-coadjoint to S` —— `φ = α ≫ β`(`β ∈ S`)なる分解では `β` が同型。 -/
def IsMinimalCoadjoint (S : MorphismProperty C) {A B : C} (φ : A ⟶ B) : Prop :=
  ∀ (X : C) (α : A ⟶ X) (β : X ⟶ B), φ = α ≫ β → S β → IsIso β

/-- **§0** `mid-adjoint to S` —— `φ = γ ≫ β ≫ α`(`β ∈ S`)なる分解では `β` が同型。 -/
def IsMidAdjoint (S : MorphismProperty C) {A B : C} (φ : A ⟶ B) : Prop :=
  ∀ (X Y : C) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B), φ = γ ≫ β ≫ α → S β → IsIso β

/-! ### §0 —— `subordinate` / `irreducible` / `FSMI` / `FSMFF-type`(原文 p.17–18)

原文 (FrdI p.17):
> admits a factorization φ = α ◦β ◦γ in C, then we shall say that β is subordinate to

原文 (FrdI p.17):
> that either α or β is an isomorphism, then we shall say that φ is irreducible. We

原文 (FrdI p.17):
> shall refer to an FSM-morphism which is irreducible as an FSMI-morphism. Thus,

原文 (FrdI p.17):
> not an isomorphism factors as a composite of finitely many FSMI-morphisms; (b)

原文 (FrdI p.18):
> that n ≤N. Thus, if C is of FSM-type, then it is of FSMFF-type. Also, we observe

★★**この 4 語は同じ段落で、順に積み上がっている**(`irreducible` → `FSMI` →
`FSMFF-type`、`subordinate` は `mid-adjoint` と同じ分解の形)。
★**`Proposition 1.14` はこの 4 語すべてを使う。**

★★**`irreducible` は当初 `Prop110.lean` に置いていたが、ここへ移した**(2026-08-16)。
★**§0 の語彙は §0 の場所に置く** —— `Proposition 1.14` が同じ段落の残り 3 語を
要求したことで、分散していたことが問題になった。
-/

/-- **§0** `subordinate to φ` —— `φ = γ ≫ β ≫ α` と分解できるとき `β` は `φ` に従属。

★`mid-adjoint` とちょうど同じ分解の形である ——
★**`φ` が `S` に mid-adjoint とは「`S` に属する従属射はすべて同型」**ということ。 -/
def IsSubordinateTo {A B X Y : C} (β : X ⟶ Y) (φ : A ⟶ B) : Prop :=
  ∃ (γ : A ⟶ X) (α : Y ⟶ B), φ = γ ≫ β ≫ α

/-- **§0** `irreducible`(射)。

原文の `φ = α ◦ β` は **「先に `β`、次に `α`」**であり、Lean では `φ = β ≫ α` と書く。
★**逆に写すと別の定義になる** —— `α` と `β` の役割は非対称だからである
(minimal-adjoint と minimal-coadjoint の違いがまさにそれ)。

★**「分解できない」ではなく「自明にしか分解できない」。** -/
def IsIrreducibleMor {A B : C} (φ : A ⟶ B) : Prop :=
  ¬ IsIso φ ∧ ∀ (X : C) (β : A ⟶ X) (α : X ⟶ B), φ = β ≫ α → IsIso β ∨ IsIso α

/-- ★**非退化(下から)**: 同型は irreducible でない。定義の第1条件そのものだが、
**「irreducible が空でない」を主張する前に「全体でもない」を押さえる**ために書く。 -/
theorem not_isIrreducibleMor_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] :
    ¬ IsIrreducibleMor φ := fun h => h.1 inferInstance

/-- ★**非退化(上から)**: `φ` が irreducible なら、`φ = β ≫ α` の分解で
**両方が非同型になることはない**。定義の言い換えだが、
★**使うときはこの向き**(分解を与えて片方の同型性を得る)である。 -/
theorem IsIrreducibleMor.isIso_left_or_right {A B X : C} {φ : A ⟶ B}
    (h : IsIrreducibleMor φ) (β : A ⟶ X) (α : X ⟶ B) (hf : φ = β ≫ α) :
    IsIso β ∨ IsIso α := h.2 X β α hf

/-- ★**irreducible は同型との合成で保たれる(右から)**。

★`Proposition 1.14, (ii)` の十分性で、底の分解を `𝒞` へ持ち上げるときに要る ——
持ち上げは**同型を除いてしか**底を再現しないからである。 -/
theorem IsIrreducibleMor.comp_isIso {A B E : C} {f : A ⟶ B}
    (h : IsIrreducibleMor f) (g : B ⟶ E) [IsIso g] : IsIrreducibleMor (f ≫ g) := by
  refine ⟨fun hiso => h.1 ?_, fun X u v huv => ?_⟩
  · haveI := hiso
    have hf : f = (f ≫ g) ≫ inv g := by simp
    rw [hf]
    infer_instance
  · have hf : f = u ≫ (v ≫ inv g) := by rw [← Category.assoc, ← huv]; simp
    rcases h.2 X u (v ≫ inv g) hf with h1 | h2
    · exact Or.inl h1
    · refine Or.inr ?_
      haveI := h2
      have hv : v = (v ≫ inv g) ≫ g := by simp
      rw [hv]
      infer_instance

/-- **§0** `FSMI-morphism` —— irreducible な FSM 射。 -/
def IsFSMI {A B : C} (φ : A ⟶ B) : Prop := IsFSMMorphism φ ∧ IsIrreducibleMor φ

/-- ★**同型は FSM 射**。 -/
theorem isFSMMorphism_of_isIso {A B : C} (f : A ⟶ B) [IsIso f] : IsFSMMorphism f :=
  ⟨fun Z γ => ⟨Z, γ ≫ inv f, 𝟙 Z, by
    rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id, Category.id_comp]⟩,
   inferInstance⟩

/-- ★**irreducible は同型との合成で保たれる(左から)**。 -/
theorem IsIrreducibleMor.isIso_comp {A A' B : C} (e : A ⟶ A') [IsIso e]
    {f : A' ⟶ B} (h : IsIrreducibleMor f) : IsIrreducibleMor (e ≫ f) := by
  refine ⟨fun hiso => h.1 ?_, fun X u v huv => ?_⟩
  · haveI := hiso
    have hf : f = inv e ≫ (e ≫ f) := by simp
    rw [hf]
    infer_instance
  · have hf : f = (inv e ≫ u) ≫ v := by rw [Category.assoc, ← huv, ← Category.assoc]; simp
    rcases h.2 X (inv e ≫ u) v hf with h1 | h2
    · refine Or.inl ?_
      haveI := h1
      have hu : u = e ≫ (inv e ≫ u) := by simp
      rw [hu]
      infer_instance
    · exact Or.inr h2

/-- ★**FSMI は同型との合成で保たれる(左から)**。 -/
theorem IsFSMI.isIso_comp {A A' B : C} (e : A ⟶ A') [IsIso e]
    {f : A' ⟶ B} (h : IsFSMI f) : IsFSMI (e ≫ f) :=
  ⟨IsFSMMorphism.comp (isFSMMorphism_of_isIso e) h.1, h.2.isIso_comp e⟩

/-- ★**原文の「Thus, a category of FSM-type does not contain any FSMI-morphisms.」** ——
FSM 型の圏には FSMI 射が無い。★**原文が `Thus` で片づける一歩**だが、
中身は「FSM 射は同型 ⟹ irreducible の第1条件に反する」だけである。 -/
theorem not_isFSMI_of_isOfFSMType (h : IsOfFSMType C) {A B : C} (φ : A ⟶ B) :
    ¬ IsFSMI φ := fun hf => hf.2.1 (h A B φ hf.1)

/-- **§0** `FSMI 射の合成鎖` —— `φ : A ⟶ B` が長さ `n` の FSMI 射の合成であること。

★原文の `φ_n ◦ … ◦ φ_1` を、**長さを添字として持つ帰納的述語**で書く。
★**`FSMFF-type` の条件 (a)(b) はどちらもこの述語で述べられる** ——
(a) は「ある `n` がある」、(b) は「`n` に上界がある」。 -/
inductive IsFSMIChain : ℕ → ∀ {A B : C}, (A ⟶ B) → Prop
  | nil {A : C} : IsFSMIChain 0 (𝟙 A)
  | cons {n : ℕ} {A B E : C} {φ : A ⟶ B} {ψ : B ⟶ E} :
      IsFSMI φ → IsFSMIChain n ψ → IsFSMIChain (n + 1) (φ ≫ ψ)

/-- ★**長さ 1 以上の鎖は、最初の因子として FSMI 射を切り出せる**。

★★**`FSMFF-type` の条件 (a) を「従属する FSMI 射がある」に変換する道具**である
(原文の「Base(α) admits a subordinate FSMI-morphism」)。 -/
theorem IsFSMIChain.exists_first {n : ℕ} (hn : 0 < n) {A B : C} {f : A ⟶ B}
    (h : IsFSMIChain n f) :
    ∃ (X : C) (u : A ⟶ X) (v : X ⟶ B), f = u ≫ v ∧ IsFSMI u := by
  cases h with
  | nil => omega
  | cons hu _ => exact ⟨_, _, _, rfl, hu⟩

/-- ★**FSMI 鎖の合成は FSM 射**。 -/
theorem IsFSMIChain.isFSMMorphism {n : ℕ} {A B : C} {f : A ⟶ B}
    (h : IsFSMIChain n f) : IsFSMMorphism f := by
  induction h with
  | nil => exact isFSMMorphism_of_isIso _
  | cons hu _ ih => exact IsFSMMorphism.comp hu.1 ih

/-- ★★**`φ ≫ ψ` が FSM で `ψ` が mono なら `φ` も FSM**。

★★**原文の「since ψ◦φ and ψ are FSM-morphisms, it thus follows formally that
φ is also an FSM-morphism」がこれである。**
★fiberwise-surjectivity は `γ ≫ ψ` に当てて `ψ` の mono で消す。 -/
theorem isFSMMorphism_of_comp {A B E : C} (φ : A ⟶ B) (ψ : B ⟶ E)
    (h : IsFSMMorphism (φ ≫ ψ)) (hψ : Mono ψ) : IsFSMMorphism φ := by
  haveI := hψ
  haveI := h.2
  refine ⟨fun Z γ => ?_, mono_of_mono φ ψ⟩
  obtain ⟨Dd, δA, δZ, hsq⟩ := h.1 (γ ≫ ψ)
  refine ⟨Dd, δA, δZ, (cancel_mono ψ).mp ?_⟩
  rw [Category.assoc, Category.assoc]
  exact hsq

/-- ★**鎖に同型を前置しても長さは変わらない**(長さ 1 以上のとき)。

★★**`Proposition 1.14, (iii)` の「(a)(b) 型は底が同型なので隣に吸収できる」**が
これである。★**長さ 0 のときは成り立たない**(`𝟙` に同型を前置しても `𝟙` にならない)ので、
そこは別に扱う。 -/
theorem IsFSMIChain.isIso_comp {n : ℕ} {A A' B : C} (e : A ⟶ A') [IsIso e]
    {f : A' ⟶ B} (h : IsFSMIChain n f) (hn : 0 < n) : IsFSMIChain n (e ≫ f) := by
  cases h with
  | nil => omega
  | cons hu hv =>
      rw [← Category.assoc]
      exact IsFSMIChain.cons (hu.isIso_comp e) hv

variable (C) in
/-- **§0** `category of FSMFF-type`(FSM-finitely factorizable type)。

原文の 2 条件:
- **(a)** 同型でない FSM 射は、**有限個の FSMI 射の合成**に分解する
- **(b)** 各対象 `A` について、`A` を域とする FSMI 射の合成鎖の長さに**一様な上界**がある

★★**(b) の量化子の順序が要点**である —— 上界 `N` は **`A` ごと**に取れればよく、
圏全体で一様である必要はない。★**「for every A ∈ Ob(C), there exists a natural
number N」の順序をそのまま写す。** -/
def IsOfFSMFFType : Prop :=
  (∀ (A B : C) (φ : A ⟶ B), IsFSMMorphism φ → ¬ IsIso φ →
      ∃ n : ℕ, 0 < n ∧ IsFSMIChain n φ) ∧
  (∀ A : C, ∃ N : ℕ, ∀ (n : ℕ) (B : C) (φ : A ⟶ B), IsFSMIChain n φ → n ≤ N)

/-- ★**原文の「Thus, if C is of FSM-type, then it is of FSMFF-type.」** ——
FSM 型なら FSMFF 型。

★★**原文は `Thus` で済ませるが、2 条件それぞれに理由が要る**:
- **(a)**: FSM 型では同型でない FSM 射が**存在しない**ので、含意は空虚に真
- **(b)**: FSMI 射が 1 本も無いので、鎖の長さは必ず `0`。したがって `N = 0` で足りる

★**(b) のほうは「空虚に真」ではない** —— 鎖 `IsFSMIChain n φ` から
`n = 0` を**導く**必要がある(`cons` の場合に FSMI 射の非存在で矛盾を出す)。 -/
theorem isOfFSMFFType_of_isOfFSMType (h : IsOfFSMType C) : IsOfFSMFFType C := by
  constructor
  · exact fun A B φ hφ hni => absurd (h A B φ hφ) hni
  · refine fun A => ⟨0, ?_⟩
    intro n B φ hchain
    cases hchain with
    | nil => exact Nat.le_refl 0
    | cons hfsmi _ => exact absurd hfsmi (not_isFSMI_of_isOfFSMType h _)

/-- ★**原文の「no endomorphism of an object of a category of FSMFF-type is an
FSMI-morphism」** —— FSMFF 型の圏では自己射は FSMI でない。

★★**原文は「by condition (b)」と一言で片づけるが、中身は
「FSMI な自己射があれば、それを何度でも繰り返して任意の長さの鎖が作れる」**である。
★`N + 1` 回繰り返した鎖が `N` を超えるので矛盾する。 -/
theorem not_isFSMI_endo_of_isOfFSMFFType (h : IsOfFSMFFType C) {A : C} (φ : End A) :
    ¬ IsFSMI (φ : A ⟶ A) := by
  intro hφ
  obtain ⟨N, hN⟩ := h.2 A
  have hchain : ∀ n : ℕ, IsFSMIChain n ((φ : A ⟶ A) ^ n : End A) := by
    intro n
    induction n with
    | zero => exact IsFSMIChain.nil
    | succ m ih =>
        have hpow : ((φ : A ⟶ A) ^ (m + 1) : End A) = (φ : A ⟶ A) ≫ ((φ : A ⟶ A) ^ m : End A) := by
          rw [pow_succ]; rfl
        rw [hpow]
        exact IsFSMIChain.cons hφ ih
  exact absurd (hN (N + 1) A _ (hchain (N + 1))) (by omega)

/-! ### 非退化(3語)

★**満たす例**と**満たさない例**を、`S = ⊤`(すべての射)について両方与える。
★この2つを合わせると **totally epimorphic な圏では
「`⊤` に minimal-adjoint」⟺「同型」**が出る —— 定義が空でも全体でもないことの証拠。 -/

/-- ★**同型でない射は `⊤` に minimal-adjoint でない** —— `φ = φ ≫ 𝟙` を取ればよい。

★coadjoint / mid-adjoint も同じ分解で落ちる。 -/
theorem not_isMinimalAdjoint_top_of_not_isIso {A B : C} (φ : A ⟶ B) (h : ¬ IsIso φ) :
    ¬ IsMinimalAdjoint (⊤ : MorphismProperty C) φ :=
  fun hm => h (hm B φ (𝟙 B) (by simp) trivial)

theorem not_isMinimalCoadjoint_top_of_not_isIso {A B : C} (φ : A ⟶ B) (h : ¬ IsIso φ) :
    ¬ IsMinimalCoadjoint (⊤ : MorphismProperty C) φ :=
  fun hm => h (hm A (𝟙 A) φ (by simp) trivial)

theorem not_isMidAdjoint_top_of_not_isIso {A B : C} (φ : A ⟶ B) (h : ¬ IsIso φ) :
    ¬ IsMidAdjoint (⊤ : MorphismProperty C) φ :=
  fun hm => h (hm A B (𝟙 A) φ (𝟙 B) (by simp) trivial)

/-- ★**totally epimorphic な圏では、同型は `⊤` に minimal-adjoint** ——
`isIso_of_isIso_comp`(原文 p.16 の「in passing」)の直接の帰結。 -/
theorem isMinimalAdjoint_top_of_isIso (h : IsTotallyEpimorphic C) {A B : C} (φ : A ⟶ B)
    [IsIso φ] : IsMinimalAdjoint (⊤ : MorphismProperty C) φ := by
  rintro X β α rfl -
  exact (isIso_of_isIso_comp h β α).1

/-- ★同上(coadjoint 側)。 -/
theorem isMinimalCoadjoint_top_of_isIso (h : IsTotallyEpimorphic C) {A B : C} (φ : A ⟶ B)
    [IsIso φ] : IsMinimalCoadjoint (⊤ : MorphismProperty C) φ := by
  rintro X α β rfl -
  exact (isIso_of_isIso_comp h α β).2

/-- ★★**したがって、totally epimorphic な圏では
「`⊤` に minimal-adjoint」は「同型」とちょうど一致する。**

★定義が退化していない(空でも全体でもない)ことの証拠であり、
同時に **`Proposition 1.7, (ii)(iii)` が「`S` を小さく取ると条件が弱まる」形の主張**
であることを示している。 -/
theorem isMinimalAdjoint_top_iff (h : IsTotallyEpimorphic C) {A B : C} (φ : A ⟶ B) :
    IsMinimalAdjoint (⊤ : MorphismProperty C) φ ↔ IsIso φ := by
  refine ⟨fun hm => ?_, fun hi => by haveI := hi; exact isMinimalAdjoint_top_of_isIso h φ⟩
  by_contra hc
  exact not_isMinimalAdjoint_top_of_not_isIso φ hc hm

/-- ★任意の射は `isomorphisms` に minimal-adjoint(定義から)。
`S` を小さく取ると条件が空になる、という**もう一方の端**。 -/
theorem isMinimalAdjoint_isomorphisms {A B : C} (φ : A ⟶ B) :
    IsMinimalAdjoint (MorphismProperty.isomorphisms C) φ :=
  fun _ _ _ _ hs => hs

/-! ### §0 —— `of Aut-type` / `End-equivalence` / `abstractly equivalent`

`Proposition 1.8, (i)(ii)` がこの3語で述べられている。

原文 (FrdI p.14):
> that A is Aut-saturated (respectively, Autsub-saturated; of Aut-type) if AutC(A) =

原文 (FrdI p.14):
> object of C is Aut-saturated (respectively, Autsub-saturated; of Aut-type), then we

原文 (FrdI p.14):
> shall say that an arrow A →B of C is an End-equivalence if there exists an arrow

原文 (FrdI p.14):
> B →A in C.
-/

/-- **§0** 対象 `A` が `Aut-saturated` —— `Aut_𝒞(A) = Aut^sub_𝒞(A)`。

★`Aut ⊆ Aut^sub` はつねに成り立つ(`isSubAutomorphism_of_isIso`)ので、
条件は「sub-automorphism がすべて同型」に尽きる。 -/
def IsAutSaturatedObj (A : C) : Prop := ∀ α : End A, IsSubAutomorphism α → IsIso α

/-- **§0** 対象 `A` が `Aut^sub-saturated` —— `Aut^sub_𝒞(A) = End_𝒞(A)`。 -/
def IsAutSubSaturatedObj (A : C) : Prop := ∀ α : End A, IsSubAutomorphism α

/-- **§0** 対象 `A` が `of Aut-type` —— `Aut_𝒞(A) = End_𝒞(A)`、
すなわち**自己射がすべて同型**。 -/
def IsOfAutTypeObj (A : C) : Prop := ∀ α : End A, IsIso α

variable (C) in
/-- **§0** 圏 `𝒞` が `of Aut-type` —— すべての対象が `of Aut-type`。 -/
def IsOfAutType : Prop := ∀ A : C, IsOfAutTypeObj A

/-- **§0** `End-equivalence`。

★★**測定**: 原文は「an arrow `A → B` … is an End-equivalence if there exists an
arrow `B → A`」と**射について**述べるが、条件に現れるのは**始域と終域だけ**で、
その射自身は一切現れない。**射の性質ではなく対象の対の性質である。**
写す側でこれを勝手に直さず、原文どおり射を引数に取る形で写す。 -/
def IsEndEquivalence {A B : C} (_φ : A ⟶ B) : Prop := Nonempty (B ⟶ A)

/-! ### §0 —— `abstractly equivalent`(原文 p.17)

原文 (FrdI p.17):
> — where the horizontal arrows are isomorphisms in C — as an abstract equivalence

原文 (FrdI p.17):
> from f1 to f2. If there exists an abstract equivalence from f1 to f2, then we shall

原文 (FrdI p.17):
> say that f1, f2 are abstractly equivalent.
-/

/-- **§0** `abstractly equivalent` —— 横向きが同型である可換四角形で結ばれること。 -/
def AbstractlyEquivalent {A₁ B₁ A₂ B₂ : C} (f₁ : A₁ ⟶ B₁) (f₂ : A₂ ⟶ B₂) : Prop :=
  ∃ (a : A₁ ≅ A₂) (b : B₁ ≅ B₂), f₁ ≫ b.hom = a.hom ≫ f₂

theorem AbstractlyEquivalent.refl {A B : C} (f : A ⟶ B) : AbstractlyEquivalent f f :=
  ⟨Iso.refl A, Iso.refl B, by simp⟩

theorem AbstractlyEquivalent.symm {A₁ B₁ A₂ B₂ : C} {f₁ : A₁ ⟶ B₁} {f₂ : A₂ ⟶ B₂}
    (h : AbstractlyEquivalent f₁ f₂) : AbstractlyEquivalent f₂ f₁ := by
  obtain ⟨a, b, hab⟩ := h
  refine ⟨a.symm, b.symm, ?_⟩
  show f₂ ≫ b.inv = a.inv ≫ f₁
  rw [Iso.eq_inv_comp, ← Category.assoc, ← hab]
  simp

/-- ★`f ≫ g` と `g ≫ f` がともに同型なら `f` は同型。

`f ≫ g` から `f` が split mono、`g ≫ f` から `f` が split epi になり、
mathlib の `isIso_of_mono_of_isSplitEpi` で閉じる。

★`Proposition 1.8, (i)` の逆向き(`𝒟` が `of Aut-type` なら linear End-equivalence は
pre-step)で使う。 -/
theorem isIso_of_comp_isIso_both {X Y : C} (f : X ⟶ Y) (g : Y ⟶ X)
    [IsIso (f ≫ g)] [IsIso (g ≫ f)] : IsIso f := by
  haveI : IsSplitMono f := ⟨⟨g ≫ inv (f ≫ g), by rw [← Category.assoc, IsIso.hom_inv_id]⟩⟩
  haveI : IsSplitEpi f := ⟨⟨inv (g ≫ f) ≫ g, by rw [Category.assoc, IsIso.inv_hom_id]⟩⟩
  exact isIso_of_mono_of_isSplitEpi f

/-! ### 非退化(`of Aut-type` / `End-equivalence`) -/

/-- ★前順序の圏は `of Aut-type` —— 自己射は `𝟙` しかない。 -/
theorem isOfAutType_of_subsingleton_hom (hs : ∀ X Y : C, Subsingleton (X ⟶ Y)) :
    IsOfAutType C := by
  intro A α
  refine ⟨𝟙 A, ?_, ?_⟩ <;> exact (hs A A).elim _ _

/-- `Vee` は `of Aut-type`。 -/
theorem isOfAutType_vee : IsOfAutType Vee :=
  isOfAutType_of_subsingleton_hom fun X Y => Preorder.subsingleton_hom X Y

/-- ★`Type` は **`of Aut-type` でない** —— `Bool` の定数写像は同型でない。 -/
theorem not_isOfAutType_type : ¬ IsOfAutType (Type) := fun h =>
  not_isIso_constTrue (h Bool constTrue)

/-- ★**同型は End-equivalence**(逆射を取ればよい)。 -/
theorem isEndEquivalence_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] :
    IsEndEquivalence φ := ⟨inv φ⟩

/-- ★`Vee` の `left ⟶ top` は **End-equivalence でない** —— `top ⟶ left` が無い。

★これが「射の性質ではなく対象の対の性質」であることの実例でもある。 -/
theorem not_isEndEquivalence_vee (β : Vee.left ⟶ Vee.top) : ¬ IsEndEquivalence β := by
  rintro ⟨g⟩
  rcases leOfHom g with h | h <;> exact Vee.noConfusion h

/-! ### ★§0 `CFP`(categorical fiber product)

原文 (FrdI p.17):
> are functors, then we define the “CFP” — i.e., “categorical fiber product” —

原文 (FrdI p.17):
> where Ai ∈Ob(Ci) (for i = 1, 2); α is an isomorphism of D; and whose morphisms

★★**mathlib の実測(2026-08-15)**: `IsoComma` は **0 件**。しかし
`CategoryTheory.Comma` はあり、その対象は三つ組
`⟨left, right, hom : L.obj left ⟶ R.obj right⟩`、射は
`L.map left ≫ Y.hom = X.hom ≫ R.map right` を満たす対である。

★原文の条件 `β ∘ Φ₁(γ₁) = Φ₂(γ₂) ∘ α` は、合成の向きを直すと
`Φ₁.map γ₁ ≫ β = α ≫ Φ₂.map γ₂` であり、**`Comma` の条件と字義どおり一致する**。
したがって **CFP は `Comma` の「`hom` が同型」という充満部分圏**である
(p.17 を 400 dpi で目視して確認した)。

★**mathlib にあるものを書き直さない** —— `Comma` はそのまま使い、
`ObjectProperty.FullSubcategory` を被せるだけにする。
-/

section CFP

universe v₁ u₁ v₂ u₂ v₃ u₃

variable {C₁ : Type u₁} [Category.{v₁} C₁] {C₂ : Type u₂} [Category.{v₂} C₂]
  {Dd : Type u₃} [Category.{v₃} Dd]

/-- CFP を切り出す対象の性質 —— 三つ組の `α` が**同型**であること。 -/
def cfpProp (Φ₁ : C₁ ⥤ Dd) (Φ₂ : C₂ ⥤ Dd) : ObjectProperty (Comma Φ₁ Φ₂) :=
  fun X => IsIso X.hom

/-- **§0 `𝒞₁ ×_𝒟 𝒞₂`** —— categorical fiber product。 -/
abbrev CFP (Φ₁ : C₁ ⥤ Dd) (Φ₂ : C₂ ⥤ Dd) : Type _ := (cfpProp Φ₁ Φ₂).FullSubcategory

/-- 第1射影 `𝒞₁ ×_𝒟 𝒞₂ ⥤ 𝒞₁`。 -/
def cfpFst (Φ₁ : C₁ ⥤ Dd) (Φ₂ : C₂ ⥤ Dd) : CFP Φ₁ Φ₂ ⥤ C₁ :=
  (cfpProp Φ₁ Φ₂).ι ⋙ Comma.fst Φ₁ Φ₂

/-- 第2射影 `𝒞₁ ×_𝒟 𝒞₂ ⥤ 𝒞₂`。 -/
def cfpSnd (Φ₁ : C₁ ⥤ Dd) (Φ₂ : C₂ ⥤ Dd) : CFP Φ₁ Φ₂ ⥤ C₂ :=
  (cfpProp Φ₁ Φ₂).ι ⋙ Comma.snd Φ₁ Φ₂

/-- ★**§0 の主張** —— `Φ₂` が圏同値なら第1射影も圏同値。

原文 (FrdI p.17):
> Φ2(γ2)◦α. One verifies easily that if Φ2 is an equivalence, then the natural projection

★**中身は「`α` が同型であること」が第2成分を完全に決める**という一点である:
`Φ₂(γ₂) = α⁻¹ ∘ Φ₁(γ₁) ∘ β` なので、`Φ₂` の忠実性から一意、充満性から存在が出る。 -/
theorem cfpFst_isEquivalence (Φ₁ : C₁ ⥤ Dd) (Φ₂ : C₂ ⥤ Dd) [Φ₂.IsEquivalence] :
    (cfpFst Φ₁ Φ₂).IsEquivalence := by
  haveI hfaith : (cfpFst Φ₁ Φ₂).Faithful := by
    constructor
    intro X Y f g hfg
    haveI hX : IsIso X.obj.hom := X.property
    have hl : f.hom.left = g.hom.left := hfg
    have hr : Φ₂.map f.hom.right = Φ₂.map g.hom.right := by
      rw [← cancel_epi X.obj.hom, ← f.hom.w, ← g.hom.w, hl]
    exact InducedCategory.hom_ext (CommaMorphism.ext hl (Φ₂.map_injective hr))
  haveI hfull : (cfpFst Φ₁ Φ₂).Full := by
    constructor
    intro X Y u
    haveI hX : IsIso X.obj.hom := X.property
    obtain ⟨v, hv⟩ := Φ₂.map_surjective (inv X.obj.hom ≫ Φ₁.map u ≫ Y.obj.hom)
    refine ⟨InducedCategory.homMk ⟨u, v, ?_⟩, rfl⟩
    show Φ₁.map u ≫ Y.obj.hom = X.obj.hom ≫ Φ₂.map v
    rw [hv, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    rfl
  haveI hess : (cfpFst Φ₁ Φ₂).EssSurj := by
    constructor
    intro A
    obtain ⟨e, he⟩ : ∃ e : Φ₂.obj (Φ₂.objPreimage (Φ₁.obj A)) ≅ Φ₁.obj A,
        e = Φ₂.objObjPreimageIso (Φ₁.obj A) := ⟨_, rfl⟩
    exact ⟨⟨⟨A, Φ₂.objPreimage (Φ₁.obj A), e.inv⟩, inferInstanceAs (IsIso e.inv)⟩,
      ⟨Iso.refl A⟩⟩
  exact ⟨hfaith, hfull, hess⟩

end CFP

/-! ### ★出典の紐付け(`.src`) -/

def IsOfAutTypeObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — of Aut-type",
    sectionId := "frdi-s0-aut-type" }

def IsAutSaturatedObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — Aut-saturated",
    sectionId := "frdi-s0-aut-type" }

def IsEndEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — End-equivalence",
    sectionId := "frdi-s0-end-equivalence" }

def AbstractlyEquivalent.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — abstractly equivalent",
    sectionId := "frdi-s0-abstract-equiv" }

def IsMinimalAdjoint.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — minimal-adjoint",
    sectionId := "frdi-s0-minimal-adjoint" }

def IsMinimalCoadjoint.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — minimal-coadjoint",
    sectionId := "frdi-s0-minimal-adjoint" }

def IsMidAdjoint.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — mid-adjoint",
    sectionId := "frdi-s0-minimal-adjoint" }

def IsSubordinateTo.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — subordinate",
    sectionId := "frdi-s0-irreducible" }

def IsIrreducibleMor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — irreducible",
    sectionId := "frdi-s0-irreducible" }

def IsFSMI.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — FSMI-morphism",
    sectionId := "frdi-s0-irreducible" }

def IsOfFSMFFType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — FSMFF-type",
    sectionId := "frdi-s0-fsmff" }

def IsSubAutomorphism.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — sub-automorphism",
    sectionId := "frdi-s0-subaut" }

def SubAut.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — Aut^sub",
    sectionId := "frdi-s0-subaut" }

def IsTotallyEpimorphic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 15, item := "§0 Categories — totally epimorphic",
    sectionId := "frdi-s0-tot-epi" }

def IsFiberwiseSurjective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — fiberwise-surjective",
    sectionId := "frdi-s0-fsm" }

def IsFSMMorphism.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — FSM-morphism",
    sectionId := "frdi-s0-fsm" }

def IsOfFSMType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14, item := "§0 Categories — category of FSM-type",
    sectionId := "frdi-s0-fsm" }

def CFP.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — CFP (categorical fiber product)",
    sectionId := "frdi-s0-cfp" }

end ABC3.Found.FrdI
