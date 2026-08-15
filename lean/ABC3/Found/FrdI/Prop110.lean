import ABC3.Found.FrdI.Prop17

/-!
# [FrdI] Proposition 1.10 —— Morphisms of Frobenius Type

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.34–35(**400 dpi 目視確認 2026-08-15**、親が p.34 / p.35 / p.17 を描画して照合)。

原文 (FrdI p.34):
> (Morphisms of Frobenius Type) Let Φ be a divisorial

原文 (FrdI p.34):
> monoid on a connected, totally epimorphic category D; C →FΦ a Frobenioid.

## ★この命題の規模(着手前の測定)

**6 条 (i)–(vi)、主張は 21**:

| 条 | 主張 | 内訳 |
|---|---|---|
| (i) | **10** | 一意存在 1 / `degFr` 1 / `Div` 1 / **7 タイプの保存** |
| (ii) | 3 | 分解 1 / `degFr` 1 / `Div` 1 |
| (iii) | 3 | 像のモノイド 1 / `𝒪^▷` 1 / `𝒪^×` 1 |
| (iv) | 2 | prime-Frobenius ⟺ irreducible 1 / 無限個の同型類 1 |
| (v) | 1 | Frobenius 型 ⟺ prime-Frobenius の合成 |
| (vi) | 2 | `𝒞^istr` が sub-quasi-Frobenius-trivial 1 / group-like ⟹ Frobenius-trivial 1 |

## ★★`Proposition 1.6` との対照

`1.6` は外から来る関手 `G : 𝒟' ⥤ 𝒟` に **FSM 保存しか仮定されていない**ことが原因で
未解決が 4 件出た。★**`1.10` は単一の Frobenioid の中の話**であり、
外から入ってくる未定義の対象がない。★**着手前の見立て: 未解決は出ない。**

## ★不足していた語彙

`Proposition 1.10` が要求する §0 の語彙のうち、**irreducible** が未実装だった。
ここで足す(物理 p.17、400 dpi 目視確認)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★§0 の語彙 —— `irreducible`

原文 (FrdI p.17):
> φ. If φ is not an isomorphism, but, for every factorization φ = α ◦β in C, it holds

原文 (FrdI p.17):
> that either α or β is an isomorphism, then we shall say that φ is irreducible. We
-/

section Vocab

/-- **[FrdI] §0** `irreducible`(射)。

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

end Vocab

def IsIrreducibleMor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 17, item := "§0 Categories — irreducible",
    sectionId := "frdi-s0-irreducible" }

/-! ## ★(i) —— 一意性

原文 (FrdI p.34):
> pre-step; pull-back morphism; co-angular morphism; base-isomorphism;

★**引用を選び直した記録**: 最初は「Then there exists a **unique** morphism φ′ : A′ →B′ …」を
引いたが、★**ゲートが落とした** ——`′` が layout モードで消えるため
(31/65 文字で停止)。**Lean の docstring には `data-txt` の逃げ道が無い**ので、
★**照合できる行に差し替えるしかない。**
★**これは事故 #3 の 3 度目**であり、**同じページの同じ `′` でも捕まる箇所と
捕まらない箇所がある**という前段の測定を、独立に再現したことになる。

★**証明(原文 p.35、目視)**:
> Observe that uniqueness follows from the fact that C is totally epimorphic.

★**我々の `IsTotallyEpimorphic C` は「すべての射が epi」**なので、
`α : A ⟶ A′` が epi であることから直ちに出る。★**Frobenius 型であることは
一意性には使わない** —— 使うのは存在のほうである。
-/

/-- **(i) の一意性** —— 四角形を可換にする `φ′` は高々1つ。

★原文は「totally epimorphic から従う」の一言。実際に使うのは
**`α` が epi であること**だけで、★**`α` が Frobenius 型であることは使わない。**
(これは原文が言っていないことである: 仮定の一部が一意性には効いていない。) -/
theorem prop_1_10_i_uniq {A B A' B' : C} (hC : IsTotallyEpimorphic C)
    (φ : A ⟶ B) (α : A ⟶ A') (β : B ⟶ B') (φ₁' φ₂' : A' ⟶ B')
    (h₁ : φ ≫ β = α ≫ φ₁') (h₂ : φ ≫ β = α ≫ φ₂') : φ₁' = φ₂' := by
  haveI : Epi α := hC _ _ _
  exact (cancel_epi α).mp (h₁ ▸ h₂)

include P in
/-- ★**上の一意性を Frobenioid の言葉で**(実際に使う形)。
`P.totEpiC` が仮定を供給する。

★**`include P in` が要る**(分類表 #7): 文に `P` が現れないのに本体で使うから。
★**これは「文に出ない仮定を本体が使っている」ことの機械的な検出**でもある ——
上の `prop_1_10_i_uniq` は仮定を文に出したので `include` が要らない。 -/
theorem prop_1_10_i_uniq' {A B A' B' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (β : B ⟶ B') (φ₁' φ₂' : A' ⟶ B')
    (h₁ : φ ≫ β = α ≫ φ₁') (h₂ : φ ≫ β = α ≫ φ₂') : φ₁' = φ₂' :=
  prop_1_10_i_uniq P.totEpiC φ α β φ₁' φ₂' h₁ h₂

/-! ## ★(i) —— 存在(Frobenius 型の場合)

★**原文が名指す材料**(p.35、目視):
> since morphisms of Frobenius type are closed under composition, with multiplying
> Frobenius degrees [cf. Proposition 1.7, (i); Remark 1.1.1], the existence of a
> morphism of Frobenius type φ′ as desired follows immediately from the existence
> and (essential) uniqueness of morphisms of Frobenius type of a given Frobenius
> degree [cf. Definition 1.3, (ii)].

★**我々の持ち物との対応**:
- 「closed under composition」→ `IsFrobeniusType.comp`(`Prop17.lean`)
- 「existence … of a given Frobenius degree」→ `F.frobDegSurj`
- 「(essential) uniqueness」→ `F.frobDegUniq`
-/

/-- **(i) の存在の骨格** —— `A'` からも `B` からも、指定された次数の
Frobenius 型射が取れる。★**原文が `Definition 1.3, (ii)` として名指す材料そのもの。**

★これ自体は `F.frobDegSurj` の言い換えだが、**(i) の存在証明で 2 回使う**ので
先に取り出しておく(★どちらの対象から取るかで役割が違う: `A'` からは `φ′`、
`B` からは `β` を作る)。 -/
theorem prop_1_10_i_frobType_available (F : FrobenioidCore P) (A : C) (n : ℕ+) :
    ∃ (B : C) (φ : A ⟶ B), IsFrobeniusType P φ ∧ P.degFr φ = n :=
  F.frobDegSurj A n

/-- **(i) の存在で使う可換性の型** —— 「`φ` に対し四角形を可換にする `φ′` が
Frobenius 型として存在する」を、**次数を明示した形**で述べたもの。

★**この宣言は「述べられること」の確認である。** 原文は 3 場合(Frobenius 型 /
pre-step / pull-back)に分けて存在を示すので、★**各場合で埋めるべき穴の形が
これで確定する。** -/
def Prop110iGoal {A B A' : C}
    (φ : A ⟶ B) (α : A ⟶ A') : Prop :=
  ∃ (B' : C) (β : B ⟶ B') (φ' : A' ⟶ B'),
    IsFrobeniusType P β ∧ P.degFr β = P.degFr α ∧ φ ≫ β = α ≫ φ'

end ABC3.Found.FrdI
