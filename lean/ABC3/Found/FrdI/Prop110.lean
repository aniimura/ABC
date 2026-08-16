import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import ABC3.Found.FrdI.Prop19

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

/-! ### ★(i) の存在 —— **第1の場合(`φ` が Frobenius 型)**

原文 (FrdI p.35):
> Frobenius type are closed under composition, with multiplying Frobenius degrees [cf.

★**原文が名指す 3 つの材料と、我々の持ち物の対応**:

| 原文 | 我々 |
|---|---|
| closed under composition | `IsFrobeniusType.comp`(`Prop17.lean`) |
| multiplying Frobenius degrees [Remark 1.1.1] | `PreFrobenioid.degFr_comp` |
| existence … of a given Frobenius degree [Def 1.3, (ii)] | `F.frobDegSurj` |
| (essential) uniqueness [Def 1.3, (ii)] | `F.frobDegUniq` |

★**原文は「follows immediately」と書くが、構成は 3 手である**:
1. `B` から次数 `degFr α` の Frobenius 型射 `β` を取る
2. `A'` から次数 `degFr φ` の Frobenius 型射 `ψ` を取る
3. `α ≫ ψ` と `φ ≫ β` は **`A` から出る同次数の Frobenius 型射**なので
   `frobDegUniq` が同型 `γ` を与える。★**`φ′ := ψ ≫ γ`。**

★★**次数が一致する理由が `mul_comm` である** —— `degFr ψ * degFr α` と
`degFr β * degFr φ` がそれぞれ `degFr φ * degFr α`、`degFr α * degFr φ` になる。
★**原文の「multiplying Frobenius degrees」の一言が、ここでは可換性として効く。**
-/

/-- **(i) の存在(第1の場合: `φ` が Frobenius 型)**。

★`φ′` も Frobenius 型として得られる —— これは (i) の後半「7 タイプの保存」の
**Frobenius 型の行**にあたる。 -/
theorem prop_1_10_i_exists_frobType (F : FrobenioidCore P) {A B A' : C}
    (φ : A ⟶ B) (hφ : IsFrobeniusType P φ)
    (α : A ⟶ A') (hα : IsFrobeniusType P α) :
    ∃ (B' : C) (β : B ⟶ B') (φ' : A' ⟶ B'),
      IsFrobeniusType P β ∧ P.degFr β = P.degFr α ∧
      IsFrobeniusType P φ' ∧ φ ≫ β = α ≫ φ' := by
  obtain ⟨B', β, hβ, hdβ⟩ := F.frobDegSurj B (P.degFr α)
  obtain ⟨Z, ψ, hψ, hdψ⟩ := F.frobDegSurj A' (P.degFr φ)
  have h1 : IsFrobeniusType P (α ≫ ψ) := IsFrobeniusType.comp P F hα hψ
  have h2 : IsFrobeniusType P (φ ≫ β) := IsFrobeniusType.comp P F hφ hβ
  have hdeg : P.degFr (α ≫ ψ) = P.degFr (φ ≫ β) := by
    rw [P.degFr_comp, P.degFr_comp, hdβ, hdψ, mul_comm]
  obtain ⟨γ, hγiso, hγ⟩ := F.frobDegUniq A Z B' (α ≫ ψ) (φ ≫ β) h1 h2 hdeg
  haveI : IsIso γ := hγiso
  refine ⟨B', β, ψ ≫ γ, hβ, hdβ,
    IsFrobeniusType.comp P F hψ (isFrobeniusType_of_isIso P γ), ?_⟩
  rw [← Category.assoc]
  exact hγ.symm

/-! ### ★(i) の `degFr` の等式

原文 (FrdI p.34):
> isometry; LB-invertible morphism), then the same is true of φ

★原文は `deg_Fr(φ) = deg_Fr(φ′)` と述べ、証明では
> The portion of assertion (i) concerning “degFr(−)”, “Div(−)” then follows immediately from
> Remark 1.1.1.
とだけ書く(p.35、目視)。

★**我々の側では、上の構成で `φ′ = ψ ≫ γ`(`γ` は同型)なので**
`degFr φ′ = degFr γ * degFr ψ = 1 * degFr φ = degFr φ` となる。
★**「同型の Frobenius 次数は 1」が要る** —— それは `isLinear_of_isIso`
(`IsLinear φ := degFr φ = 1`)である。
-/

/-- ★**同型の Frobenius 次数は 1**(`isLinear_of_isIso` の言い換え)。
★**名前が `IsLinear` なので `degFr` を探しても見つからない** ——
前段の教訓「推測せず `grep`」がここでも効いた。 -/
theorem degFr_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : P.degFr φ = 1 :=
  isLinear_of_isIso P φ

/-! ### ★★(i) の「7 タイプの保存」—— **4 つは独立な主張ではない**

原文 (FrdI p.34):
> pre-step; pull-back morphism; co-angular morphism; base-isomorphism;

★原文は 7 タイプを並べる: Frobenius type / pre-step / pull-back morphism /
co-angular morphism / base-isomorphism / isometry / LB-invertible morphism。

★★**ところが第1の場合(`φ` が Frobenius 型)では、後ろの 4 つは
「`φ′` が Frobenius 型である」ことから自動で従う。** 定義の入れ子:
```
IsFrobeniusType φ := IsLBInvertible φ ∧ IsBaseIsomorphism φ
IsLBInvertible  φ := IsCoAngular φ ∧ IsIsometric φ
```
★**つまり Frobenius 型 ⟹ LB-invertible ∧ base-isomorphism ∧ co-angular ∧ isometry。**

★**測定**: 原文の「7 タイプ」は**場合によって独立な主張の数が違う**。
第1の場合では **7 のうち 4 は自動**で、実質は 3 である。
★**列挙の長さは主張の数ではない** —— `Proposition 1.7` で我々が学んだことの再現。
-/

/-- **(i) の 4 タイプ** —— Frobenius 型から自動で従うぶん。

★**この宣言の値は「短いこと」にある。** 原文が「follows immediately from the
definitions」と書いた箇所が、本当に定義の射影だけであることを示している。 -/
theorem prop_1_10_i_four_types {A B : C} {φ : A ⟶ B} (h : IsFrobeniusType P φ) :
    IsLBInvertible P φ ∧ IsBaseIsomorphism P φ ∧ IsCoAngular P φ ∧ IsIsometric P φ :=
  ⟨h.1, h.2, h.1.1, h.1.2⟩

/-! ### ★★(i) の `Div` の公式 —— 原文が言っていない材料が要る

原文 (FrdI p.34):
> isomorphism α]. Finally, if φ is a morphism of Frobenius type (respectively,

★**引用を選び直した記録(事故 #3 の 4 度目)**: 公式そのものの行
(「In this situation, degFr(φ) = degFr(φ′); Div(φ′) = d · α∗(Div(φ)) …」)を引こうとしたが、
★**ゲートが落とした**(`′` が layout モードで消え、32/68 文字で停止)。
★**同じページの `′` でも、`α]. Finally,` の行には無いので通る。**
★**公式そのものは引用できない** —— だから**読みとして日本語で書き、
機械照合には隣接する行を使う**。★これは我々の器具の制約であって、
**「引用できない ＝ 原文に無い」ではない**(400 dpi 目視では読める)。

★原文の証明(p.35、目視)は
> The portion of assertion (i) concerning “degFr(−)”, “Div(−)” then follows
> immediately from Remark 1.1.1.
とだけ書く。

★★**実際に計算すると、`Remark 1.1.1`(＝ `Div_comp`)だけでは足りない。**
`Div_comp` を四角形の両辺に当てると
```
Div(α ≫ φ′) = Φ.map (Base α) (Div φ′) + (degFr φ′) • Div α
Div(φ ≫ β)  = Φ.map (Base φ) (Div β)  + (degFr β)  • Div φ
```
★**ここで `Div α = 0` と `Div β = 0` が要る** —— つまり
**`α` と `β` が isometry であること**。それは「Frobenius 型 ⟹ LB-invertible ⟹ isometry」
から出るが、★**原文はこの一歩を書いていない。**「Remark 1.1.1 から」だけでは
両辺の第2項・第1項が消える理由が出ない。

★結果は
```
Φ.map (Base α) (Div φ′) = (degFr β) • Div φ
```
であり、`Φ.map (Base α) : Φ(A′) → Φ(A)` の逆が原文の `α∗ : Φ(A) ⥲ Φ(A′)` なので、
これは `Div(φ′) = d · α∗(Div φ)` と同じことである
(★**我々は `α∗` を作らずに、`Φ.map (Base α)` を当てた形で述べる** ——
逆写像の構成は base-isomorphism の同型性から別に出るもので、
**公式そのものはこの形で完結している**)。
-/

/-! ### ★(i) の存在 —— **第2の場合(`φ` が pre-step)**

原文 (FrdI p.35):
> In the case where φ is a pre-step, the existence of a pre-step φ

★**原文が名指す材料**(p.35、目視): `Definition 1.3, (iv), (a)` の分解
[cf. also `Proposition 1.4, (iv)`]、`Definition 1.3, (ii)` の本質的一意性、
`Proposition 1.7, (i)` の co-angular の合成閉性。

★★**構成(4 手)**:
1. `B` から次数 `degFr α` の Frobenius 型射 `β` を取る
2. `φ ≫ β` に **`F.arbFactor`** を当てて `γ ≫ β₀ ≫ α₀` を得る
   (`γ` Frobenius 型 / `β₀` pre-step / `α₀` pull-back)
3. ★**`α₀` は同型である** —— `φ ≫ β` は base-isomorphism(`φ` も `β` も base-iso)で、
   `γ`・`β₀` も base-iso だから `Base α₀` は同型。
   **底が同型な pull-back は同型**(`isIso_of_isPullBack_of_isBaseIso`、`Prop19.lean`)。
4. 次数を数えると `degFr γ = degFr α` になるので `F.frobDegUniq` が同型 `δ` を与え、
   ★**`φ′ := inv δ ≫ β₀ ≫ α₀`。**

★★**原文が書いていないもの**: 手順 3 の「`α₀` が同型である」。
原文は分解を持ち出すだけで、**pull-back の因子が消える理由**を書いていない。
★我々の側では `Proposition 1.9` で作った `isIso_of_isPullBack_of_isBaseIso` が要る
——★**`Proposition 1.10` が `Proposition 1.9` に依存する**、という依存が 1 本増えた
(原文の引用リストには `1.9` は無い)。
-/

include P in
/-- **(i) の存在(第2の場合: `φ` が pre-step)**。

★`φ′` も pre-step として得られる。 -/
theorem prop_1_10_i_exists_preStep (F : FrobenioidCore P) {A B A' : C}
    (φ : A ⟶ B) (hφ : IsPreStep P φ)
    (α : A ⟶ A') (hα : IsFrobeniusType P α) :
    ∃ (B' : C) (β : B ⟶ B') (φ' : A' ⟶ B'),
      IsFrobeniusType P β ∧ P.degFr β = P.degFr α ∧
      IsPreStep P φ' ∧ φ ≫ β = α ≫ φ' := by
  obtain ⟨B', β, hβ, hdβ⟩ := F.frobDegSurj B (P.degFr α)
  obtain ⟨X, Y, γ, β₀, α₀, hfac, hγ, hβ₀, hα₀⟩ := F.arbFactor (φ ≫ β)
  -- `φ ≫ β` は base-isomorphism
  have hcomp_bi : IsBaseIsomorphism P (φ ≫ β) := isBaseIsomorphism_comp P hφ.2 hβ.2
  -- したがって `Base α₀` は同型 ⟹ `α₀` は同型
  have hα₀bi : IsBaseIsomorphism P α₀ := by
    -- ★`𝒟` が totally epimorphic であることを使って、底で左から簡約する
    have h1 : IsIso (P.Base (γ ≫ β₀ ≫ α₀)) := hfac ▸ hcomp_bi
    rw [P.Base_comp, P.Base_comp] at h1
    haveI := h1
    have h2 := (isIso_of_isIso_comp P.totEpiD (P.Base γ) (P.Base β₀ ≫ P.Base α₀)).2
    haveI := h2
    exact (isIso_of_isIso_comp P.totEpiD (P.Base β₀) (P.Base α₀)).2
  haveI : IsIso α₀ := isIso_of_isPullBack_of_isBaseIso P F α₀ hα₀ hα₀bi
  -- 次数を数える(★`α₀` が同型なので次数 1、`β₀` は pre-step なので次数 1)
  have hdγ : P.degFr γ = P.degFr α := by
    have h : P.degFr (φ ≫ β) = P.degFr (γ ≫ β₀ ≫ α₀) := by rw [hfac]
    rw [P.degFr_comp, P.degFr_comp, P.degFr_comp, hφ.1, hβ₀.1,
      degFr_of_isIso P α₀, hdβ] at h
    simpa using h.symm
  obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A X A' γ α hγ hα hdγ
  haveI : IsIso δ := hδiso
  refine ⟨B', β, inv δ ≫ β₀ ≫ α₀, hβ, hdβ, ?_, ?_⟩
  · exact IsPreStep.comp P (isPreStep_of_isIso P (inv δ))
      (IsPreStep.comp P hβ₀ (isPreStep_of_isIso P α₀))
  · rw [← hδ, Category.assoc, ← Category.assoc δ, IsIso.hom_inv_id,
      Category.id_comp]
    exact hfac

/-! ### ★(i) の存在 —— **第3の場合(`φ` が pull-back)**

原文 (FrdI p.35):
> sition [cf. Proposition 1.7, (i)], the existence of a pull-back morphism φ

★**原文が名指す材料**(p.35、目視): pull-back は LB-invertible [`Prop 1.4, (ii)`]、
LB-invertible は合成で閉じる [`Prop 1.7, (i)`]、**`Proposition 1.4, (v)` の分解**、
`Definition 1.3, (ii)` の本質的一意性。

★★**第2の場合と違い、ここでは分解が 2 項**(`Proposition 1.4, (v)`:
LB-invertible ⟹ `Frobenius 型 ≫ pull-back`)。だから **pull-back の因子が
残ってよい** —— `φ′` 自身が pull-back であってほしいのだから当然である。
★**第2の場合で `α₀` を消す必要があったのは、`φ′` を pre-step にしたかったから。**
★**同じ「分解して比べる」でも、目標のタイプによって因子の扱いが変わる。**
-/

include P in
/-- **(i) の存在(第3の場合: `φ` が pull-back)**。

★`φ′` も pull-back として得られる。 -/
theorem prop_1_10_i_exists_pullBack (F : FrobenioidCore P) {A B A' : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ)
    (α : A ⟶ A') (hα : IsFrobeniusType P α) :
    ∃ (B' : C) (β : B ⟶ B') (φ' : A' ⟶ B'),
      IsFrobeniusType P β ∧ P.degFr β = P.degFr α ∧
      IsPullBack P φ' ∧ φ ≫ β = α ≫ φ' := by
  obtain ⟨B', β, hβ, hdβ⟩ := F.frobDegSurj B (P.degFr α)
  obtain ⟨hφlb, hφlin⟩ := F.pullBackLB φ hφ
  -- `φ ≫ β` は LB-invertible
  have hlb : IsLBInvertible P (φ ≫ β) := IsLBInvertible.comp P F hφlb hβ.1
  -- `Proposition 1.4, (v)` の分解(2 項)
  obtain ⟨X, γ, α₀, hfac, hγ, hα₀⟩ := prop_1_4_v_mpr P F (φ ≫ β) hlb
  obtain ⟨_, hα₀lin⟩ := F.pullBackLB α₀ hα₀
  -- 次数を数える
  have hdγ : P.degFr γ = P.degFr α := by
    have h : P.degFr (φ ≫ β) = P.degFr (γ ≫ α₀) := by rw [hfac]
    rw [P.degFr_comp, P.degFr_comp, hφlin, hα₀lin, hdβ] at h
    simpa using h.symm
  obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A X A' γ α hγ hα hdγ
  haveI : IsIso δ := hδiso
  refine ⟨B', β, inv δ ≫ α₀, hβ, hdβ, ?_, ?_⟩
  · exact IsPullBack.comp P (isPullBack_of_isIso P (inv δ)) hα₀
  · rw [← hδ, Category.assoc, ← Category.assoc δ, IsIso.hom_inv_id,
      Category.id_comp]
    exact hfac

/-! ### ★★(i) の存在 —— **任意の `φ`**(3 場合を繋ぐ)

原文 (FrdI p.35):
> desired, first in the case where φ is a morphism of Frobenius type, then in the case

★**原文の構造**: 「it suffices to prove the existence of `φ′` as desired,
first in the case where `φ` is a morphism of Frobenius type, then in the case
where `φ` is a pre-step, and finally in the case where `φ` is a pull-back morphism
[cf. **the factorization of Definition 1.3, (iv), (a)**]」。

★★**「it suffices」の中身**: 任意の `φ` を `arbFactor` で `γ ≫ β₀ ≫ α₀` に分解し、
**3 つの場合を順に当てて、出てきた `β` を次の入力にする**。
```
γ  と α  → 第1の場合 → β₁, γ′    (γ  ≫ β₁ = α  ≫ γ′)
β₀ と β₁ → 第2の場合 → β₂, β₀′   (β₀ ≫ β₂ = β₁ ≫ β₀′)
α₀ と β₂ → 第3の場合 → β₃, α₀′   (α₀ ≫ β₃ = β₂ ≫ α₀′)
```
★**`φ′ := γ′ ≫ β₀′ ≫ α₀′`、`β := β₃`。**

★★**次数が全部 `degFr α` に揃う**のが要点である ——
第1の場合が `degFr β₁ = degFr α` を返し、それが第2の場合の入力になり……と
**バケツリレーで `d` が運ばれる**。★原文の「of Frobenius degree `d ∈ ℕ≥1`」という
一言が、ここでは**3 段の受け渡し**として現れる。
-/

include P in
/-- **(i) の存在(任意の `φ`)** —— ★**`Proposition 1.10, (i)` の存在部分の完成形**。

★3 つの場合を `F.arbFactor` の分解に沿って順に当てる。 -/
theorem prop_1_10_i_exists (F : FrobenioidCore P) {A B A' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (hα : IsFrobeniusType P α) :
    ∃ (B' : C) (β : B ⟶ B') (φ' : A' ⟶ B'),
      IsFrobeniusType P β ∧ P.degFr β = P.degFr α ∧ φ ≫ β = α ≫ φ' := by
  obtain ⟨X, Y, γ, β₀, α₀, hfac, hγ, hβ₀, hα₀⟩ := F.arbFactor φ
  obtain ⟨X', β₁, γ', hβ₁, hd₁, _, hsq₁⟩ :=
    prop_1_10_i_exists_frobType P F γ hγ α hα
  obtain ⟨Y', β₂, β₀', hβ₂, hd₂, _, hsq₂⟩ :=
    prop_1_10_i_exists_preStep P F β₀ hβ₀ β₁ hβ₁
  obtain ⟨B', β₃, α₀', hβ₃, hd₃, _, hsq₃⟩ :=
    prop_1_10_i_exists_pullBack P F α₀ hα₀ β₂ hβ₂
  refine ⟨B', β₃, γ' ≫ β₀' ≫ α₀', hβ₃, by rw [hd₃, hd₂, hd₁], ?_⟩
  calc φ ≫ β₃ = γ ≫ β₀ ≫ α₀ ≫ β₃ := by rw [hfac]; simp
    _ = γ ≫ β₀ ≫ β₂ ≫ α₀' := by rw [hsq₃]
    _ = γ ≫ (β₀ ≫ β₂) ≫ α₀' := by simp
    _ = γ ≫ (β₁ ≫ β₀') ≫ α₀' := by rw [hsq₂]
    _ = (γ ≫ β₁) ≫ β₀' ≫ α₀' := by simp
    _ = (α ≫ γ') ≫ β₀' ≫ α₀' := by rw [hsq₁]
    _ = α ≫ γ' ≫ β₀' ≫ α₀' := by simp

include P in
/-- ★★★**`Proposition 1.10, (i)` の原文の形** —— `β` を**与えられたもの**として受ける。

原文 (FrdI p.34):
> pre-step; pull-back morphism; co-angular morphism; base-isomorphism;

★★**監査で「量化子が逆」と指摘された項目**である。
原文は「Suppose that α, β **are** morphisms of Frobenius type」と `β` を与える（∀）が、
`prop_1_10_i_exists` は `β` を**自分で作る**（∃）。

★★**橋渡しは `Definition 1.3, (ii)` の本質的一意性（`frobDegUniq`）**である ——
同じ次数の Frobenius 型射は同型を除いて一意なので、
自分で作った `β₀` を与えられた `β` に合わせられる。
★一意性は `prop_1_10_i_uniq`（既存）。 -/
theorem prop_1_10_i_exists_given (F : FrobenioidCore P) {A B A' B' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (hα : IsFrobeniusType P α)
    (β : B ⟶ B') (hβ : IsFrobeniusType P β) (hd : P.degFr α = P.degFr β) :
    ∃! φ' : A' ⟶ B', φ ≫ β = α ≫ φ' := by
  obtain ⟨B₀, β₀, φ₀, hβ₀, hd₀, hsq₀⟩ := prop_1_10_i_exists P F φ α hα
  obtain ⟨θ, hθiso, hθ⟩ := F.frobDegUniq B B₀ B' β₀ β hβ₀ hβ (by rw [hd₀, hd])
  have hmain : φ ≫ β = α ≫ (φ₀ ≫ θ) := by
    rw [← hθ, ← Category.assoc, hsq₀, Category.assoc]
  refine ⟨φ₀ ≫ θ, hmain, ?_⟩
  intro ψ hψ
  exact prop_1_10_i_uniq P.totEpiC φ α β ψ (φ₀ ≫ θ) hψ hmain

/-! ### ★★監査で欠落が判明した主張(2026-08-16)

原文 (FrdI p.34):
> pre-step; pull-back morphism; co-angular morphism; base-isomorphism;

★**引用を選び直した記録(事故 #3 の 5 度目)**: 主張の行
「In this situation, degFr(φ) = degFr(φ′); Div(φ′) = d · α∗(Div(φ)) …」は
★**`′` が抽出で落ちるため引用できない**(32/68 文字で停止)。
★**同じ段落の、`′` を含まない行を引く。**

★★**以下は「与えられた `α`, `β`, `φ′` が可換方形をなす」という原文の状況そのもの**を扱う。
★`prop_1_10_i_exists` が `β` を**自分で作る**のに対し、
★★**こちらは `β` を与えられたものとして受け取る**(原文の量化子に合わせる)。

★**原文が「the same is true of φ′」と述べる 7 タイプのうち、
base-isomorphism と isometry をここで実装する。** -/

include P in
/-- ★★**原文「In this situation, degFr(φ) = degFr(φ′)」** ——
★**監査以前、この等式はファイルのどこにも存在しなかった。**

★`degFr` の乗法性と `ℕ+` の消去だけで出る。 -/
theorem prop_1_10_i_degFr_phi_eq {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hd : P.degFr α = P.degFr β) (hsq : φ ≫ β = α ≫ φ') :
    P.degFr φ = P.degFr φ' := by
  have h := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, hd] at h
  exact mul_left_cancel (a := P.degFr β) (by rw [h, mul_comm])

include P in
/-- ★**`φ` が base-isomorphism なら `φ′` もそう**(7 タイプのうちの 1 つ)。

★`Base` は関手的なので、方形の 3 辺が同型なら 4 辺目も同型。 -/
theorem prop_1_10_i_baseIso_of {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsBaseIsomorphism P α) (hβ : IsBaseIsomorphism P β)
    (hsq : φ ≫ β = α ≫ φ') (hφ : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism P φ' := by
  haveI : IsIso (P.Base α) := hα
  haveI : IsIso (P.Base β) := hβ
  haveI : IsIso (P.Base φ) := hφ
  have h := congrArg P.Base hsq
  rw [P.Base_comp, P.Base_comp] at h
  have hrw : P.Base φ' = inv (P.Base α) ≫ P.Base φ ≫ P.Base β := by
    rw [h, IsIso.inv_hom_id_assoc]
  show IsIso (P.Base φ')
  rw [hrw]
  infer_instance

include P in
/-- ★**`φ` が isometry なら `φ′` もそう**(7 タイプのうちの 1 つ)。

★`α`, `β` も isometry(Frobenius 型は LB-invertible、したがって isometric)なので
`Div` の公式の項が落ち、★**`Φ.map (Base α)` の単射性**(`charInj`)で閉じる。 -/
theorem prop_1_10_i_isometric_of {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsIsometric P α) (hβ : IsIsometric P β)
    (hsq : φ ≫ β = α ≫ φ') (hφ : IsIsometric P φ) :
    IsIsometric P φ' := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp] at h
  rw [show P.Div β = 0 from hβ, show P.Div φ = 0 from hφ, show P.Div α = 0 from hα] at h
  simp only [map_zero, smul_zero, add_zero] at h
  refine Φ.map_injective (P.Base α) ?_
  rw [map_zero]
  first
  | exact h
  | exact h.symm

include P in
/-- ★★**`φ` が co-angular なら `φ′` もそう**（7 タイプのうちの 1 つ）。

★★**鍵は「`φ ≫ β` が co-angular である」こと**である ——
`β` は Frobenius 型なので LB-invertible、したがって co-angularであり、
co-angular は合成で閉じる（`Definition 1.3, (iii), (a)`）。

★`φ′` の分解 `φ′ = γ ≫ β₀ ≫ α₀` に `α` を前合成すれば
`φ ≫ β` の分解になる。★**側面条件のうち base-isomorphism の選択肢だけが
手を要す** —— `γ` が base-isomorphism なら `α ≫ γ` もそうである
（`α` は Frobenius 型なので base-isomorphism）。 -/
theorem prop_1_10_i_coAngular_of (F : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsFrobeniusType P α) (hβ : IsFrobeniusType P β)
    (hsq : φ ≫ β = α ≫ φ') (hφ : IsCoAngular P φ) :
    IsCoAngular P φ' := by
  intro X Y γ β₀ α₀ hfac hlin hiso hstep hor
  have hcomp : IsCoAngular P (φ ≫ β) := F.coAngularComp φ β hφ hβ.1.1
  refine hcomp X Y (α ≫ γ) β₀ α₀ ?_ hlin hiso hstep ?_
  · rw [hsq, hfac]; simp
  · rcases hor with h | h
    · exact Or.inl h
    · exact Or.inr (isBaseIsomorphism_comp P hα.2 h)

include P in
/-- ★**`φ` が LB-invertible なら `φ′` もそう**（7 タイプのうちの 1 つ）。

★LB-invertible = co-angular ∧ isometric なので、上の 2 本の合わせ技。 -/
theorem prop_1_10_i_lbInvertible_of (F : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsFrobeniusType P α) (hβ : IsFrobeniusType P β)
    (hsq : φ ≫ β = α ≫ φ') (hφ : IsLBInvertible P φ) :
    IsLBInvertible P φ' :=
  ⟨prop_1_10_i_coAngular_of P F hα hβ hsq hφ.1,
   prop_1_10_i_isometric_of P hα.1.2 hβ.1.2 hsq hφ.2⟩

include P in
/-- ★**`φ` が linear なら `φ′` もそう**。★`degFr` の一致から直ちに出る。 -/
theorem prop_1_10_i_linear_of {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hd : P.degFr α = P.degFr β) (hsq : φ ≫ β = α ≫ φ') (hφ : IsLinear P φ) :
    IsLinear P φ' := by
  show P.degFr φ' = 1
  rw [← prop_1_10_i_degFr_phi_eq P hd hsq]
  exact hφ

include P in
/-- ★**`φ` が pre-step なら `φ′` もそう**（7 タイプのうちの 1 つ）。 -/
theorem prop_1_10_i_preStep_of {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsFrobeniusType P α) (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) (hsq : φ ≫ β = α ≫ φ') (hφ : IsPreStep P φ) :
    IsPreStep P φ' :=
  ⟨prop_1_10_i_linear_of P hd hsq hφ.1,
   prop_1_10_i_baseIso_of P hα.2 hβ.2 hsq hφ.2⟩

include P in
/-- ★**`φ` が Frobenius 型なら `φ′` もそう**（7 タイプのうちの 1 つ）。 -/
theorem prop_1_10_i_frobType_of (F : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsFrobeniusType P α) (hβ : IsFrobeniusType P β)
    (hsq : φ ≫ β = α ≫ φ') (hφ : IsFrobeniusType P φ) :
    IsFrobeniusType P φ' :=
  ⟨prop_1_10_i_lbInvertible_of P F hα hβ hsq hφ.1,
   prop_1_10_i_baseIso_of P hα.2 hβ.2 hsq hφ.2⟩

include P in
/-- ★★**`φ` が pull-back なら `φ′` もそう**（7 タイプの最後の 1 つ）。

★★**`Proposition 1.4, (ii)` を両向きに使う**:
`IsPullBack ⇔ IsLBInvertible ∧ IsLinear`。
★普遍性（hom 集合の全単射）を直接扱わずに済む。 -/
theorem prop_1_10_i_pullBack_of (F : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {α : A ⟶ A'} {β : B ⟶ B'} {φ' : A' ⟶ B'}
    (hα : IsFrobeniusType P α) (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) (hsq : φ ≫ β = α ≫ φ') (hφ : IsPullBack P φ) :
    IsPullBack P φ' := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F φ).mp hφ
  exact (prop_1_4_ii P F φ').mpr
    ⟨prop_1_10_i_lbInvertible_of P F hα hβ hsq hlb,
     prop_1_10_i_linear_of P hd hsq hlin⟩

/-! ## ★★(ii) —— **pre-step と Frobenius 型は順序を入れ替えられる**

原文 (FrdI p.34):
> (ii) Any composite morphism β ◦α of C, where α is a pre-step, and β is of

原文 (FrdI p.34):
> Frobenius type, may be written as a composite

★**3 本目の行(「where α′ is a pre-step, and β′ is of Frobenius type such that:」)は
引用できない**(事故 #3 の 5 度目。`′` が layout モードで消え 6/50 文字で停止)。
400 dpi 目視では読める。読みは下に日本語で書く。

★**合成の向き**: 原文の `β ◦ α` は「先に `α`、次に `β`」＝ Lean の `α ≫ β`。
`α′ ◦ β′` は「先に `β′`、次に `α′`」＝ `β′ ≫ α′`。
★**つまり「pre-step のあと Frobenius 型」を「Frobenius 型のあと pre-step」に組み替える。**

★原文の証明(p.35、目視)は
> assertion (ii) follows **formally** from assertion (i).
(直前に「in light of the existence of morphisms of Frobenius type of a given
Frobenius degree — cf. Definition 1.3, (ii)」と断ってある)。

★★**「formally」の中身は 3 手**:
1. `A` から次数 `degFr β` の Frobenius 型射 `β′` を取る(`frobDegSurj`)
2. **(i) の第2の場合**を `φ := α`(pre-step)と縦 `β′` に当てて、
   `α ≫ β'' = β′ ≫ α''`(`β''` は `X` から出る次数 `degFr β` の Frobenius 型)を得る
3. `β''` と `β` は **`X` から出る同次数の Frobenius 型射**なので `frobDegUniq` が同型 `ε` を与え、
   ★**`α′ := α'' ≫ ε`。**

★**原文が「formally」と言うのは正しい** —— 新しい材料は 1 つも要らず、
**(i) と `Definition 1.3, (ii)` だけで閉じる。**
★**ただし手数は 3 手ある。**「formally」は「短い」ではなく「新しい材料が要らない」の意味である。
-/

include P in
/-- **(ii)** —— `pre-step ≫ Frobenius 型` を `Frobenius 型 ≫ pre-step` に組み替える。

★次数は保たれる(`degFr β′ = degFr β`)。 -/
theorem prop_1_10_ii (F : FrobenioidCore P) {A X B : C}
    (α : A ⟶ X) (hα : IsPreStep P α) (β : X ⟶ B) (hβ : IsFrobeniusType P β) :
    ∃ (Y : C) (β' : A ⟶ Y) (α' : Y ⟶ B),
      IsFrobeniusType P β' ∧ P.degFr β' = P.degFr β ∧
      IsPreStep P α' ∧ β' ≫ α' = α ≫ β := by
  obtain ⟨Y, β', hβ', hdβ'⟩ := F.frobDegSurj A (P.degFr β)
  obtain ⟨X'', β'', α'', hβ'', hd'', hpre'', hsq⟩ :=
    prop_1_10_i_exists_preStep P F α hα β' hβ'
  -- `β''` と `β` は `X` から出る同次数の Frobenius 型射
  have hdeg : P.degFr β'' = P.degFr β := by rw [hd'', hdβ']
  obtain ⟨ε, hεiso, hε⟩ := F.frobDegUniq X X'' B β'' β hβ'' hβ hdeg
  haveI : IsIso ε := hεiso
  refine ⟨Y, β', α'' ≫ ε, hβ', hdβ',
    IsPreStep.comp P hpre'' (isPreStep_of_isIso P ε), ?_⟩
  rw [← Category.assoc, ← hsq, Category.assoc, hε]

/-! ### ★(ii) の `Div` の公式 —— **(i) のものがそのまま使える**

原文 (FrdI p.34):
> ∗for the bijection induced by applying the functor Φ to the base-

★**公式そのものの行(「Div(α′) = degFr(β) · β′∗(Div(α))」)は引用できない**
(事故 #3 の 6 度目、5/19 文字で停止)。隣接する `β′∗` の説明行を引く。

★★**(i) で作った `prop_1_10_i_Div_formula` が、四角形の向きを合わせるだけで
そのまま (ii) に効く。** 対応は
```
(i):  φ ≫ β_上 = α_左 ≫ φ′
(ii): α ≫ β    = β′   ≫ α′
```
で、`φ := α`(pre-step)、`β_上 := β`、`α_左 := β′`、`φ′ := α′`。

★**要る仮定は「`β′` と `β` が isometry」**で、どちらも Frobenius 型だから出る。
★**(i) で「原文が書いていない材料」として明示した isometry が、
(ii) でも同じ位置に要る** —— つまりあれは (i) 固有の抜けではなく、
**この形の公式に共通して要るもの**だった。

★**実体は下(`prop_1_10_ii_Div_formula`)に置く** —— (i) の公式を使うので、
**宣言の順序でそちらが先に要る**。★これは Lean の制約であって数学の順序ではない。
-/

include P in
/-- **(i) の `Div` の公式**(`Φ.map (Base α)` を当てた形)。

★**`α` と `β` が isometry であることを仮定に明示している** ——
原文は「Remark 1.1.1 から」と書くが、**実際にはこれが要る**。
★仮定を文に出したので、**何が効いているかが宣言の型から読める。** -/
theorem prop_1_10_i_Div_formula {A B A' B' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (β : B ⟶ B') (φ' : A' ⟶ B')
    (hα : IsIsometric P α) (hβ : IsIsometric P β)
    (hsq : φ ≫ β = α ≫ φ') :
    Φ.map (P.Base α) (P.Div φ') = (P.degFr β : ℕ) • P.Div φ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hα, hβ] at h
  simpa using h.symm

/-- **(i) の `degFr` の等式**(第1の場合の構成に対して)。

★上の `prop_1_10_i_exists_frobType` が返す `φ′` は `ψ ≫ γ` の形で、
`degFr ψ = degFr φ`、`γ` は同型。★**だから次数は保たれる。**

★**`FrobenioidCore` は要らない** —— `degFr` の関手性と同型の次数だけで出る。
★**「原文が (i) の材料として `Definition 1.3, (ii)` を挙げている」のは
存在のためであって、次数の等式のためではない。** -/
theorem prop_1_10_i_degFr_eq {A' Z B' : C}
    (ψ : A' ⟶ Z) (γ : Z ⟶ B') [IsIso γ] {n : ℕ+} (hψ : P.degFr ψ = n) :
    P.degFr (ψ ≫ γ) = n := by
  rw [P.degFr_comp, degFr_of_isIso, one_mul, hψ]

include P in
/-- **(ii) の `Div` の公式**。★(i) のものを向きを合わせて使うだけ。

★★**(i) で「原文が書いていない材料」として明示した isometry が、
(ii) でも同じ位置に要る** —— つまりあれは (i) 固有の抜けではなく、
**この形の公式に共通して要るもの**だった。 -/
theorem prop_1_10_ii_Div_formula {A X B Y : C}
    (α : A ⟶ X) (β : X ⟶ B) (β' : A ⟶ Y) (α' : Y ⟶ B)
    (hβ' : IsFrobeniusType P β') (hβ : IsFrobeniusType P β)
    (hsq : β' ≫ α' = α ≫ β) :
    Φ.map (P.Base β') (P.Div α') = (P.degFr β : ℕ) • P.Div α :=
  prop_1_10_i_Div_formula P α β' β α' hβ'.1.2 hβ.1.2 hsq.symm

include P in
/-- ★★★**原文の形の `Div` の公式**（2026-08-16）。

原文 (FrdI p.34):
> ∗for the bijection induced by applying the functor Φ to the base-

★**引用を選び直した記録（事故 #3 の 9 度目）**: 主張の行
「where α′ is a pre-step, and β′ is of Frobenius type such that:」は
★**`′` が抽出で落ちる**ため引用できない（6/50 文字で停止）。
★`′` は本ファイルで何度も落ちている既知の文字である。

★★**監査で「`β′` の base-isomorphism 性を仮定せず、原文の `β′∗`（全単射）の
形になっていない」と指摘されたもの**である。

★`β′` は Frobenius 型なので base-isomorphism、したがって
`Φ.map (Base β′)` は全単射であり、その逆を使って原文の向きに直せる。 -/
theorem prop_1_10_ii_Div_formula' {A X B Y : C}
    (α : A ⟶ X) (β : X ⟶ B) (β' : A ⟶ Y) (α' : Y ⟶ B)
    (hβ' : IsFrobeniusType P β') (hβ : IsFrobeniusType P β)
    (hsq : β' ≫ α' = α ≫ β) :
    haveI : IsIso (P.Base β') := hβ'.2
    P.Div α' = (P.degFr β : ℕ) • Φ.map (inv (P.Base β')) (P.Div α) := by
  haveI : IsIso (P.Base β') := hβ'.2
  have h := prop_1_10_ii_Div_formula P α β β' α' hβ' hβ hsq
  have h2 := congrArg (Φ.map (inv (P.Base β'))) h
  rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at h2
  exact h2


/-! ## ★(iii) —— `of perfect type` と像のモノイドの perfect 性

原文 (FrdI p.34):
> (iii) Suppose that C is of perfect type. Then the monoids in the image of

原文 (FrdI p.34):
> Φ are perfect. If, moreover, C is of isotropic and Frobenius-normalized type,

★**原文の証明**(p.35、目視):
> In light of the existence of morphisms of Frobenius type of a given Frobenius degree
> [cf. Definition 1.3, (ii)] and the equivalences of categories of Definition 1.3, (iii), (d),
> the fact that Φ(A) is perfect follows immediately [cf. Remark 1.1.1] from the fact
> that A is perfect [cf. Definition 1.2, (iv)].

★**不足していた語彙 `of perfect type` をここで足す** ——
`IsPerfectObj` は対象について定義済みで、「型」は**全対象がそれである**という形。
-/

/-- **[FrdI] Definition 1.3 系** `𝒞` が **of perfect type** —— 全対象が perfect。

★**「型」の定義はどれもこの形**(全対象がその性質を持つ)なので、
`IsPerfectObj` から機械的に作れる。★**原文が別に定義を置いていないのは
そのためだが、我々は名前が要る。** -/
def IsOfPerfectType : Prop := ∀ A : C, IsPerfectObj P A

/-! ### ★★(iii) の核心 —— `n` 倍が全射になる仕組み

★**原文は「follows immediately [cf. Remark 1.1.1]」と書く。**
その `Remark 1.1.1`(＝ `Div_comp`)を四角形
`φ₁ ≫ ψ' = ψ ≫ φ₂`(`φ₁`・`φ₂` は次数 `n` の Frobenius 型、`ψ`・`ψ'` は pre-step)
に当てると:
```
Div(φ₁ ≫ ψ') = Φ.map (Base φ₁) (Div ψ') + (degFr ψ') • Div φ₁
Div(ψ ≫ φ₂)  = Φ.map (Base ψ)  (Div φ₂) + (degFr φ₂)  • Div ψ
```
★**`Div φ₁ = Div φ₂ = 0`**(Frobenius 型 ⟹ isometry)なので、両辺は
```
Φ.map (Base φ₁) (Div ψ') = n • Div ψ
```
★★**これが「`n` 倍が全射」の正体である** ——
`IsPerfectObj` の条件 (b) が `ψ` の存在を与え、その `Div ψ` が
`Φ.map (Base φ₁) (Div ψ')` の `n` 分の 1 になる。

★**(i) の `Div` の公式とまったく同じ形**である(`prop_1_10_i_Div_formula`)。
★**同じ補題が (i)・(ii)・(iii) の 3 箇所で効いている** ——
原文が 3 箇所で別々に「Remark 1.1.1 から」と書いているものは、**1 本の等式**だった。
-/

include P in
/-- **(iii) の核心** —— 次数 `n` の Frobenius 型射に沿った四角形は、
`Div` を **`n` 倍の関係**に翻訳する。

★`prop_1_10_i_Div_formula` の直接の帰結。★**新しい証明は 1 行も要らない。** -/
theorem prop_1_10_iii_div_nsmul {B₁ B₁' B₂ B₂' : C}
    (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂') (ψ : B₁ ⟶ B₂) (ψ' : B₁' ⟶ B₂')
    (hφ₁ : IsFrobeniusType P φ₁) (hφ₂ : IsFrobeniusType P φ₂)
    {n : ℕ+} (hn : P.degFr φ₂ = n)
    (hsq : φ₁ ≫ ψ' = ψ ≫ φ₂) :
    Φ.map (P.Base ψ) (P.Div φ₂) + (n : ℕ) • P.Div ψ
      = Φ.map (P.Base φ₁) (P.Div ψ') + (P.degFr ψ' : ℕ) • P.Div φ₁ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hn] at h
  exact h.symm

include P in
/-- ★**上を、両辺の消える項を落とした形で**。★**これが原文の言う
「Remark 1.1.1 から従う」の実体**である。 -/
theorem prop_1_10_iii_div_nsmul' {B₁ B₁' B₂ B₂' : C}
    (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂') (ψ : B₁ ⟶ B₂) (ψ' : B₁' ⟶ B₂')
    (hφ₁ : IsIsometric P φ₁) (hφ₂ : IsIsometric P φ₂)
    {n : ℕ+} (hn : P.degFr φ₂ = n)
    (hsq : φ₁ ≫ ψ' = ψ ≫ φ₂) :
    Φ.map (P.Base φ₁) (P.Div ψ') = (n : ℕ) • P.Div ψ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hn, hφ₁, hφ₂] at h
  simpa using h

/-! ## ★(v) —— Frobenius 型 ⟺ prime-Frobenius の合成

原文 (FrdI p.35):
> (v) A morphism of C is a morphism of Frobenius type if and only if it is a

原文 (FrdI p.35):
> composite of prime-Frobenius morphisms.

★★**「合成」を帰納的に定義するとき、原文の曖昧さが 2 段階で出る。**

**第1段**: 原文は「a composite of prime-Frobenius morphisms」と書くだけで、
**空の合成を含むかを言っていない**。含めないと、次数 1 の Frobenius 型射が
合成として書けない。

**第2段(こちらが本質)**: ★**空の合成を「恒等射」に取るだけでは足りない。**
`Proposition 1.4, (iii)` により **次数 1 の Frobenius 型射は同型**であり、
同型は一般に **`𝟙 A` そのものではない**。(v) が iff である以上、
★**基底の場合は「同型」でなければならない。**

★★**これが原文の 3 つ目の曖昧さである。**
「composite of prime-Frobenius morphisms」は **同型を法として読む**しかない。
★**形式化しなければ「恒等射」と「同型」の差は見えない** ——
自然言語では「合成が無い場合」で済んでしまう。
-/

/-- **prime-Frobenius 射の合成**(帰納的)。

★**基底は同型**(恒等射ではない)。理由は上の docstring を見よ ——
次数 1 の Frobenius 型射は同型であって、`𝟙` とは限らない。 -/
inductive IsPrimeFrobComposite : ∀ {A B : C}, (A ⟶ B) → Prop
  | iso {A B : C} (φ : A ⟶ B) : IsIso φ → IsPrimeFrobComposite φ
  | cons {A B E : C} {φ : A ⟶ B} {ψ : B ⟶ E} :
      IsPrimeFrobenius P φ → IsPrimeFrobComposite ψ → IsPrimeFrobComposite (φ ≫ ψ)

include P in
/-- **(v) の片向き** —— prime-Frobenius の合成は Frobenius 型。

★**`Proposition 1.7, (i)` の閉性だけで出る。** 原文が (v) の証明で
何を使うかは書いていないが、★**この向きは合成閉性の帰納法である。** -/
theorem isFrobeniusType_of_isPrimeFrobComposite (F : FrobenioidCore P)
    {A B : C} {φ : A ⟶ B} (h : IsPrimeFrobComposite P φ) : IsFrobeniusType P φ := by
  induction h with
  | iso ψ hψ => exact @isFrobeniusType_of_isIso _ _ _ _ _ P _ _ ψ hψ
  | cons hφ _ ih => exact IsFrobeniusType.comp P F hφ.1 ih

include P in
/-- ★**非退化**: prime-Frobenius 自身は(恒等射との合成として)合成である。
★**「合成の集合が空でない」を押さえる。** -/
theorem isPrimeFrobComposite_of_isPrimeFrobenius {A B : C} {φ : A ⟶ B}
    (h : IsPrimeFrobenius P φ) : IsPrimeFrobComposite P φ := by
  have := IsPrimeFrobComposite.cons h
    (IsPrimeFrobComposite.iso (P := P) (𝟙 B) inferInstance)
  simpa using this

/-! ### ★(v) の逆向き —— Frobenius 型 ⟹ prime-Frobenius の合成

★**原文は証明を書いていない**((v) の証明は本文に無い)。
★**構成は次数についての強帰納法**である:
- `degFr φ = 1` のとき: ★**`Proposition 1.4, (iii)` により `φ` は同型**。基底の場合。
- `degFr φ = n > 1` のとき: `n` の素因数 `p` を取り `n = m * p` と書く。
  次数 `p` の Frobenius 型射 `φ₁` と次数 `m` の `φ₂` を `frobDegSurj` で取ると
  `φ₁ ≫ φ₂` の次数は `m * p = n` なので、`frobDegUniq` が同型 `δ` を与え
  ★**`φ = φ₁ ≫ (φ₂ ≫ δ)`**。`φ₂ ≫ δ` の次数は `m < n` なので帰納法が回る。

★★**「素因数を 1 つ剥がす」が帰納のステップである。** 原文の (v) は
「Frobenius 型 ⟺ prime-Frobenius の合成」と書くだけだが、
★**その中身は「次数の素因数分解」そのもの**である。
-/

include P in
/-- **(v) の逆向き** —— Frobenius 型射は prime-Frobenius 射の合成。

★次数についての強帰納法。基底は **同型**(次数 1)。 -/
theorem isPrimeFrobComposite_of_isFrobeniusType (F : FrobenioidCore P) :
    ∀ (n : ℕ) {A B : C} (φ : A ⟶ B), IsFrobeniusType P φ →
      ((P.degFr φ : ℕ+) : ℕ) = n → IsPrimeFrobComposite P φ := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro A B φ hφ hn
    rcases eq_or_ne n 1 with h1 | h1
    · -- 次数 1 ⟹ 同型
      have hlin : IsLinear P φ := by
        show P.degFr φ = 1
        exact PNat.coe_injective (by rw [hn, h1]; rfl)
      exact IsPrimeFrobComposite.iso φ (prop_1_4_iii P F φ hφ.1 ⟨hlin, hφ.2⟩)
    · -- 素因数を 1 つ剥がす
      obtain ⟨p, hp, k, hk⟩ := Nat.exists_prime_and_dvd h1
      have hnpos : 0 < n := by rw [← hn]; exact (P.degFr φ).pos
      have hppos : 0 < p := hp.pos
      have hkpos : 0 < k := Nat.pos_of_ne_zero (by rintro rfl; omega)
      set pp : ℕ+ := ⟨p, hppos⟩ with hpp
      set kk : ℕ+ := ⟨k, hkpos⟩ with hkk
      obtain ⟨Z, φ₁, hφ₁, hd₁⟩ := F.frobDegSurj A pp
      obtain ⟨W, φ₂, hφ₂, hd₂⟩ := F.frobDegSurj Z kk
      have hcomp : IsFrobeniusType P (φ₁ ≫ φ₂) := IsFrobeniusType.comp P F hφ₁ hφ₂
      have hdc : P.degFr (φ₁ ≫ φ₂) = P.degFr φ := by
        rw [P.degFr_comp, hd₁, hd₂]
        refine PNat.coe_injective ?_
        simp only [PNat.mul_coe, hpp, hkk, hn, hk]
        exact Nat.mul_comm k p
      obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A W B (φ₁ ≫ φ₂) φ hcomp hφ hdc
      haveI : IsIso δ := hδiso
      have hdk : ((P.degFr (φ₂ ≫ δ) : ℕ+) : ℕ) = k := by
        rw [P.degFr_comp, degFr_of_isIso P δ, one_mul, hd₂, hkk]
        rfl
      have hklt : k < n := by
        rw [hk]
        calc k = 1 * k := (Nat.one_mul k).symm
          _ < p * k := Nat.mul_lt_mul_of_lt_of_le hp.one_lt (le_refl k) hkpos
      have hrec := ih k hklt (φ₂ ≫ δ)
        (IsFrobeniusType.comp P F hφ₂ (isFrobeniusType_of_isIso P δ)) hdk
      have hprime : IsPrimeFrobenius P φ₁ := ⟨hφ₁, by rw [hd₁]; simpa [hpp] using hp⟩
      have : φ = φ₁ ≫ φ₂ ≫ δ := by rw [← hδ, Category.assoc]
      rw [this]
      exact IsPrimeFrobComposite.cons hprime hrec

include P in
/-- ★★**`Proposition 1.10, (v)` の完成形** —— Frobenius 型 ⟺ prime-Frobenius の合成。

★**原文はこの iff を証明なしで述べる。** 我々の側では
- `⟸` は `Proposition 1.7, (i)` の合成閉性の帰納法
- `⟹` は **次数の素因数分解による強帰納法**(基底は `Proposition 1.4, (iii)` で同型)
であり、★**「合成」の基底を同型に取らないと iff にならない**ことが分かった。 -/
theorem prop_1_10_v (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsFrobeniusType P φ ↔ IsPrimeFrobComposite P φ :=
  ⟨fun h => isPrimeFrobComposite_of_isFrobeniusType P F _ φ h rfl,
   isFrobeniusType_of_isPrimeFrobComposite P F⟩

/-! ### ★★(v) の基底 —— 「意図的な修復」を**反例**に格上げする(2026-08-16)

★★**取り下げ表の (v) の項目**: 基底を「任意の同型」に取ったことを
**意図的な修復**として記録していた。★しかし「修復」と「反例がある」は違う ——
★**素直な読み(基底は恒等射)が実際に偽であることを示せば、
同型に取るのは選択ではなく強制になる。**

★**ここでそれを示す。** 素直な読み `IsPrimeFrobCompositeId` を別に定義し、
★**恒等射でない同型は決してその形に書けない**ことを証明する。
`Proposition 1.7, (i)`(同型は Frobenius 型)により、そのような射は
**(v) の左辺を満たす**から、素直な読みでは iff が破れる。

★★**理由は次数だけである** —— prime-Frobenius 射の次数は素数(≥2)であり、
次数は合成で掛かる(`Remark 1.1.1`)。★空でない合成の次数は 2 以上だが、
同型の次数は 1 である。★**だから空の合成しか残らず、それは `𝟙` である。**

★★**これは「原文の言葉が形式化で 2 通りに分かれ、一方が偽になる」例である。**
自然言語の「composite of prime-Frobenius morphisms」は、
★**同型を法として読む以外にない。**
-/

/-- ★**(v) の素直な読み** —— 基底を**恒等射**に取った「prime-Frobenius の合成」。

★これが原文の字面どおりの読みである。★下でこれが**偽**であることを示す。 -/
inductive IsPrimeFrobCompositeId : ∀ {A B : C}, (A ⟶ B) → Prop
  | nil {A : C} : IsPrimeFrobCompositeId (𝟙 A)
  | cons {A B E : C} {φ : A ⟶ B} {ψ : B ⟶ E} :
      IsPrimeFrobenius P φ → IsPrimeFrobCompositeId ψ → IsPrimeFrobCompositeId (φ ≫ ψ)

include P in
/-- ★**次数 1 の「素直な合成」は空の合成しかない**。

★prime-Frobenius 射の次数は素数(≥2)で、次数は合成で掛かるので、
★**`cons` の場合は次数が 2 以上になり、1 にはなり得ない。** -/
theorem eq_id_of_isPrimeFrobCompositeId_of_degFr_one :
    ∀ {A B : C} {φ : A ⟶ B}, IsPrimeFrobCompositeId P φ →
      ((P.degFr φ : ℕ+) : ℕ) = 1 → ∃ h : A = B, φ = eqToHom h := by
  intro A B φ h
  induction h with
  | nil => exact fun _ => ⟨rfl, by simp⟩
  | @cons A' B' E' φ' ψ' hp _ _ =>
      intro hd
      exfalso
      rw [P.degFr_comp, PNat.mul_coe] at hd
      have h2 : 2 ≤ ((P.degFr φ' : ℕ+) : ℕ) := hp.2.two_le
      have h1 : ((P.degFr φ' : ℕ+) : ℕ) = 1 :=
        Nat.dvd_one.mp (Dvd.intro_left _ hd)
      omega

include P in
/-- ★★★**(v) の素直な読みは偽である** —— 恒等射でない自己同型は、
prime-Frobenius 射の(恒等射を基底とする)合成には決して書けない。

★★**しかし `Proposition 1.7, (i)` によりそれは Frobenius 型である。**
したがって素直な読みでは (v) の iff が破れる ——
★**基底を同型に取ったのは修復ではなく、強制であった。** -/
theorem not_isPrimeFrobCompositeId_of_isIso_of_ne_id {A : C} (u : A ⟶ A) [IsIso u]
    (hne : u ≠ 𝟙 A) : ¬ IsPrimeFrobCompositeId P u := by
  intro h
  obtain ⟨hAA, hu⟩ := eq_id_of_isPrimeFrobCompositeId_of_degFr_one P h
    (by rw [degFr_of_isIso P u]; rfl)
  exact hne (by simpa using hu)

/-! ## ★(vi) —— `𝒞^istr` と group-like 対象

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

原文 (FrdI p.35):
> every group-like object A ∈Ob(Cistr) is Frobenius-trivial.

★**後半(group-like ⟹ Frobenius-trivial)の証明**(p.36、目視):
> If, moreover, A is group-like, then [since Cistr is a Frobenioid — cf. Proposition 1.9, (v)]
> it follows from Definition 1.3, (i), (a), (b), that there exist [co-angular — cf.
> Proposition 1.4, (i)] pre-steps A′ →A, A′ →A′′, where A′′ is Frobenius-trivial.
> But by Proposition 1.4, (iii), these pre-steps are isomorphisms, so A is Frobenius-trivial.

★★**「これらの pre-step が同型である」に `Proposition 1.4, (iii)` を使うには
LB-invertible が要る。** 原文は `[co-angular — cf. Proposition 1.4, (i)]` で
**co-angular のほうだけ**括弧に入れており、★**isometric のほうを書いていない。**

★**isometric はどこから来るか**: `A` が group-like で、`Φ` が **divisorial(＝ sharp を含む)**
なら、★**`Φ(A)` のすべての元が 0** である。よって `A` から出る任意の射の `Div` は 0、
すなわち **isometric**。

★★**これは (i) の `Div` の公式・(ii)・(iii) で見つけた「原文が isometry を書かない」
というパターンの 4 例目**である。★**同じ省略が 4 箇所で起きている。**
-/

/-- ★**group-like ＋ sharp ⟹ モノイドは零**。

★`IsGroupLike M ↔ ∀ a, IsAddUnit a`(`isGroupLike_iff`)と
`IsSharp M := ∀ a, IsAddUnit a → a = 0` を合わせるだけ。
★**しかしこの一行が (vi) の後半の要点である。** -/
theorem eq_zero_of_isGroupLike_of_isSharp {M : Type w} [AddCommMonoid M]
    (hg : IsGroupLike M) (hs : IsSharp M) (a : M) : a = 0 :=
  hs a ((isGroupLike_iff M).mp hg a)

include P in
/-- ★★**group-like 対象から出る射はすべて isometric**。

★**原文が書いていない一歩**(上の docstring を見よ)。
`Φ` が divisorial であることの **sharp の部分**がここで効く ——
★**divisorial の定義に sharp が入っている理由の 1 つがこれである。** -/
theorem isIsometric_of_isGroupLikeObj {A B : C} (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) : IsIsometric P φ :=
  eq_zero_of_isGroupLike_of_isSharp hA (P.divisorial _).2 (P.Div φ)

include P in
/-- ★**(vi) の後半の要点** —— isotropic かつ group-like な対象から出る
pre-step は同型。

★原文の「by Proposition 1.4, (iii), these pre-steps are isomorphisms」の中身。
★**co-angular は isotropic から(`Proposition 1.4, (i)`)、
isometric は group-like から(上の補題)。原文は後者を書いていない。** -/
theorem isIso_of_preStep_of_isGroupLikeObj (F : FrobenioidCore P) {A B : C}
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X) (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) (hφ : IsPreStep P φ) : IsIso φ :=
  prop_1_4_iii P F φ ⟨prop_1_4_i P φ hiso, isIsometric_of_isGroupLikeObj P hA φ⟩ hφ

include P in
/-- ★**`𝒞` が isotropic 型のときの形**(`𝒞^istr` で使うのはこちら)。

★`Proposition 1.9` の `𝒞^istr` は全対象が isotropic なので、仮定はそこから供給される。 -/
theorem isIso_of_preStep_of_isGroupLikeObj' (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) (hφ : IsPreStep P φ) : IsIso φ :=
  isIso_of_preStep_of_isGroupLikeObj P F (fun X _ => hiso X) hA φ hφ

/-! ## ★(iv) —— prime-Frobenius ⟺ irreducible

原文 (FrdI p.34):
> nius morphism if and only if it is irreducible [cf. §0]. In particular, if A ∈Ob(C)

★**原文の証明**(p.36、目視):
> assertion (iv) follows immediately from Proposition 1.7, (v), and
> **the well-known structure of the multiplicative monoid N≥1**
> [cf. also Definition 1.3, (ii)]

★★**「the well-known structure of the multiplicative monoid ℕ≥1」とは何か。**
原文はこの一句で (iv) を片づける。★**その中身は「乗法モノイド `(ℕ≥1, ×)` において
既約な元とは素数のことである」**である。

★**mathlib には `PNat.Prime`(＝ `(p : ℕ).Prime`)はあるが、
「乗法モノイド `ℕ+` の既約元」という形の補題は無い**(`Mathlib/Data/PNat/Prime.lean` を確認)。
★**だからここで書く。** ——★**原文が「well-known」と呼ぶものが、
我々の道具立てには存在しなかった、という測定でもある。**

★**そしてこれが (iv) の骨格である**: 射 `φ` の既約性は
「`φ = β ≫ α` の分解で片方が同型」であり、Frobenius 型射に限れば
**次数の分解 `degFr φ = degFr α * degFr β`** に写る。同型は次数 1 なので、
★**射の既約性は次数の既約性に、そのまま翻訳される。**
-/

/-- ★★**乗法モノイド `ℕ≥1` の既約元は素数**(原文の「well-known structure」)。

★`n : ℕ+` について「1 でなく、どう 2 つに分けても片方が 1」⟺「素数」。
★**mathlib に無いので書いた。** -/
theorem pnat_irreducible_iff_prime (n : ℕ+) :
    (n ≠ 1 ∧ ∀ a b : ℕ+, n = a * b → a = 1 ∨ b = 1) ↔ Nat.Prime (n : ℕ) := by
  have hn0 : (n : ℕ) ≠ 0 := n.pos.ne'
  constructor
  · rintro ⟨hne, hsplit⟩
    rw [← Nat.irreducible_iff_nat_prime]
    refine ⟨fun h => hne (PNat.coe_injective ?_), fun a b hab => ?_⟩
    · rw [Nat.isUnit_iff] at h; rw [h]; rfl
    · have ha0 : a ≠ 0 := by rintro rfl; rw [Nat.zero_mul] at hab; exact hn0 hab
      have hb0 : b ≠ 0 := by rintro rfl; rw [Nat.mul_zero] at hab; exact hn0 hab
      have hsp := hsplit ⟨a, Nat.pos_of_ne_zero ha0⟩ ⟨b, Nat.pos_of_ne_zero hb0⟩
        (PNat.coe_injective (by simpa using hab))
      rcases hsp with h | h
      · exact Or.inl (Nat.isUnit_iff.mpr (congrArg PNat.val h))
      · exact Or.inr (Nat.isUnit_iff.mpr (congrArg PNat.val h))
  · intro hp
    have hirr : Irreducible (n : ℕ) := (Nat.irreducible_iff_nat_prime _).mpr hp
    refine ⟨fun h => ?_, fun a b hab => ?_⟩
    · rw [h] at hp; exact Nat.not_prime_one (by simpa using hp)
    · have hab' : (n : ℕ) = (a : ℕ) * (b : ℕ) := by rw [hab]; rfl
      rcases hirr.isUnit_or_isUnit hab' with h | h
      · exact Or.inl (PNat.coe_injective (by simpa using Nat.isUnit_iff.mp h))
      · exact Or.inr (PNat.coe_injective (by simpa using Nat.isUnit_iff.mp h))

/-! ### ★(iv) の片向き —— irreducible ⟹ prime-Frobenius

★**射の分解が次数の分解に写る**のが要点:
`φ = φ₁ ≫ (φ₂ ≫ δ)` の形を `frobDegSurj` ＋ `frobDegUniq` で作り、
既約性から片方が同型と分かると、★**その次数が 1 になる**。
つまり **`degFr φ = a * b` ⟹ `a = 1` または `b = 1`**。
そこに `pnat_irreducible_iff_prime` を当てて素数性を得る。

★**`isotropic` の仮定が効く場所**: 「次数 1 ⟹ 同型」に `Proposition 1.4, (iii)` を
使うには LB-invertible が要り、その co-angular の部分を
**`Proposition 1.4, (i)`(isotropic なら全射が co-angular)**が供給する。
★**原文が `[cf. §0]` としか書かない (iv) の中で、isotropic はここに効いている。**
-/

include P in
/-- **(iv) の片向き** —— isotropic な域を持つ Frobenius 型射が irreducible なら
prime-Frobenius。 -/
theorem prop_1_10_iv_mpr (F : FrobenioidCore P) {A B : C}
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X)
    (φ : A ⟶ B) (hφ : IsFrobeniusType P φ) (hirr : IsIrreducibleMor φ) :
    IsPrimeFrobenius P φ := by
  refine ⟨hφ, ?_⟩
  rw [← pnat_irreducible_iff_prime]
  refine ⟨fun h => hirr.1 (prop_1_4_iii P F φ hφ.1 ⟨h, hφ.2⟩), fun a b hab => ?_⟩
  obtain ⟨Z, φ₁, hφ₁, hd₁⟩ := F.frobDegSurj A a
  obtain ⟨W, φ₂, hφ₂, hd₂⟩ := F.frobDegSurj Z b
  have hcomp : IsFrobeniusType P (φ₁ ≫ φ₂) := IsFrobeniusType.comp P F hφ₁ hφ₂
  have hdc : P.degFr (φ₁ ≫ φ₂) = P.degFr φ := by
    rw [P.degFr_comp, hd₁, hd₂, hab, mul_comm]
  obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A W B (φ₁ ≫ φ₂) φ hcomp hφ hdc
  haveI : IsIso δ := hδiso
  have hfac : φ = φ₁ ≫ (φ₂ ≫ δ) := by rw [← hδ, Category.assoc]
  rcases hirr.2 Z φ₁ (φ₂ ≫ δ) hfac with h | h
  · haveI := h
    exact Or.inl (by rw [← hd₁]; exact (degFr_of_isIso P φ₁).symm ▸ rfl)
  · haveI := h
    refine Or.inr ?_
    have hd := degFr_of_isIso P (φ₂ ≫ δ)
    rw [P.degFr_comp, degFr_of_isIso P δ, one_mul, hd₂] at hd
    exact hd

/-! ### ★sharp モノイドの 2 つの補題 —— (iv) の逆向きに要る

★**どちらも「和が 0 なら可逆、可逆なら 0」という sharp の一行の帰結**だが、
**使う形が違う**ので別々に置く。

★**原文はこれらを一度も書かない。** `Definition 1.1, (i)` で `divisorial` に
`sharp` を入れておいて、以後は暗黙に使う。★**「sharp を仮定に入れる」という
一回の判断が、下流の何箇所で効いているかは、形式化するまで見えない。**
-/

/-- ★**sharp なら、和が 0 の両項は 0**。 -/
theorem eq_zero_of_add_eq_zero_of_isSharp {M : Type w} [AddCommMonoid M]
    (hs : IsSharp M) {a b : M} (h : a + b = 0) : a = 0 ∧ b = 0 :=
  ⟨hs a ⟨⟨a, b, h, by rw [add_comm]; exact h⟩, rfl⟩,
   hs b ⟨⟨b, a, by rw [add_comm]; exact h, h⟩, rfl⟩⟩

/-- ★**sharp なら、`n • a = 0`(`n ≥ 1`)から `a = 0`**。

★`n • a = a + (n-1) • a` なので `a` は可逆。 -/
theorem eq_zero_of_nsmul_eq_zero_of_isSharp {M : Type w} [AddCommMonoid M]
    (hs : IsSharp M) {n : ℕ} (hn : 0 < n) {a : M} (h : n • a = 0) : a = 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [succ_nsmul] at h
  exact (eq_zero_of_add_eq_zero_of_isSharp hs h).2

/-! ### ★(iv) の逆向き —— prime-Frobenius ⟹ irreducible

★**任意の分解 `φ = β ≫ α` に対して片方が同型**を示す。分解は**任意の射**なので、
`β`・`α` が Frobenius 型であることは**与えられていない**。そこで 3 つを別々に出す:

1. **base-isomorphism**: `Base φ = Base β ≫ Base α` が同型で `𝒟` は totally epimorphic
   なので、**両方が同型**(`isIso_of_isIso_comp`)
2. **isometric**: `Div φ = 0`(Frobenius 型)を `Div_comp` で開くと
   `Φ.map (Base β) (Div α) + (degFr α) • Div β = 0`。
   ★**sharp なので両項が 0**、さらに `Φ.map` は**常に単射**なので `Div α = 0`、
   `n • x = 0` から `Div β = 0`。
3. **co-angular**: `𝒞` が isotropic 型なので `Proposition 1.4, (i)` から自動

★あとは次数。`degFr α * degFr β = degFr φ` が素数なので
`pnat_irreducible_iff_prime` の `⟸` 向きから **片方が 1**。
次数 1 ＋ LB-invertible ＋ base-iso なら `Proposition 1.4, (iii)` で同型。

★★**ここで sharp が 2 回効いている**(和の分解と `n` 倍の消去)。
★**原文は (iv) を「Proposition 1.7, (v) と ℕ≥1 の well-known structure から
immediately」と書くだけで、この 3 段は書いていない。**
-/

include P in
/-- **(iv) の逆向き** —— prime-Frobenius なら irreducible。

★`𝒞` が isotropic 型であることを使う(原文の `Proposition 1.4, (i)` 経由)。 -/
theorem prop_1_10_iv_mp (F : FrobenioidCore P) {A B : C} (hA : IsIsotropic P A)
    (φ : A ⟶ B) (hp : IsPrimeFrobenius P φ) : IsIrreducibleMor φ := by
  refine ⟨fun h => ?_, fun X β α hfac => ?_⟩
  · haveI := h
    have h2 := hp.2
    rw [degFr_of_isIso P φ] at h2
    exact Nat.not_prime_one (by simpa using h2)
  · -- base-isomorphism を両方に降ろす
    have hb : IsIso (P.Base β ≫ P.Base α) := by
      rw [← P.Base_comp, ← hfac]; exact hp.1.2
    haveI := hb
    obtain ⟨hbβ, hbα⟩ := isIso_of_isIso_comp P.totEpiD (P.Base β) (P.Base α)
    -- isometric を両方に降ろす
    have hsum : Φ.map (P.Base β) (P.Div α) + (P.degFr α : ℕ) • P.Div β = 0 := by
      rw [← P.Div_comp, ← hfac]; exact hp.1.1.2
    obtain ⟨e1, e2⟩ := eq_zero_of_add_eq_zero_of_isSharp (P.divisorial _).2 hsum
    have hdβ : P.Div β = 0 :=
      eq_zero_of_nsmul_eq_zero_of_isSharp (P.divisorial _).2 (P.degFr α).pos e2
    have hdα : P.Div α = 0 :=
      Φ.map_injective (P.Base β) (by rw [e1, map_zero])
    -- 次数の分解
    have hd : P.degFr φ = P.degFr α * P.degFr β := by rw [hfac, P.degFr_comp]
    rcases ((pnat_irreducible_iff_prime (P.degFr φ)).mpr hp.2).2 _ _ hd with h | h
    · -- degFr α = 1 ⟹ α は同型
      -- ★`α : X ⟶ B` の域 `X` は `β : A ⟶ X` を通じて `A` から到達するので isotropic
      exact Or.inr (prop_1_4_iii P F α
        ⟨prop_1_4_i P α (fun _ f => F.isotropicClosed f (F.isotropicClosed β hA)), hdα⟩
        ⟨h, hbα⟩)
    · -- degFr β = 1 ⟹ β は同型
      exact Or.inl (prop_1_4_iii P F β
        ⟨prop_1_4_i P β (fun _ f => F.isotropicClosed f hA), hdβ⟩ ⟨h, hbβ⟩)

include P in
/-- ★★**`Proposition 1.10, (iv)` の前半の完成形**(iff)。

★★**仮定は「域 `A` が isotropic」だけ**である。

原文 (FrdI p.34):
> (iv) A morphism of Frobenius type with isotropic domain is a prime-Frobe-

★**監査以前は `∀ X : C, IsIsotropic P X`(圏全体)を仮定していた** ——
原文が要求しない仮定を足していた。`F.isotropicClosed` で域から伝播させれば足りる。 -/
theorem prop_1_10_iv (F : FrobenioidCore P) {A B : C} (hA : IsIsotropic P A)
    (φ : A ⟶ B) (hφ : IsFrobeniusType P φ) :
    IsPrimeFrobenius P φ ↔ IsIrreducibleMor φ :=
  ⟨prop_1_10_iv_mp P F hA φ,
   prop_1_10_iv_mpr P F (fun _ f => F.isotropicClosed f hA) φ hφ⟩

/-! ### ★(vi) の前半の核心 —— base-identity 自己射を pre-step に沿って降ろす

原文 (FrdI p.36、目視):
> since φC is a base-identity endomorphism of Frobenius degree d, and γ is a pre-
> step, it follows [cf. Remark 1.1.1] that φB is also a base-identity endomorphism of
> Frobenius degree d.

★**この一段だけを取り出す**。四角形 `γ ≫ φ_B = φ_C ≫ γ` があり、
`φ_C` が base-identity で次数 `d`、`γ` が pre-step のとき、
★**`φ_B` も base-identity で次数 `d`**。

★**証明の中身**:
- **次数**: `degFr (γ ≫ φ_B) = degFr φ_B * degFr γ = degFr φ_B`(`γ` は pre-step ⟹ 次数 1)、
  `degFr (φ_C ≫ γ) = degFr γ * degFr φ_C = d`。よって `degFr φ_B = d`。
- **base-identity**: `Base γ ≫ Base φ_B = Base φ_C ≫ Base γ = Base γ`(`φ_C` が base-identity)。
  ★**`𝒟` が totally epimorphic なので `Base γ` は epi**、よって `Base φ_B = 𝟙`。

★★**原文の「[cf. Remark 1.1.1]」は次数の部分だけを指している。**
base-identity の部分は **`𝒟` の totally epimorphic 性**から出るもので、
`Remark 1.1.1`(合成則)ではない。★**1 つの括弧が 2 つの理由を覆っていた。**
-/

include P in
/-- **(vi) の前半の核心** —— base-identity 自己射を pre-step に沿って降ろす。

★次数は `Remark 1.1.1`(合成則)から、base-identity は **`𝒟` の totally epimorphic 性**から。
★**原文の括弧 `[cf. Remark 1.1.1]` は前者だけを指している。** -/
theorem prop_1_10_vi_descend {B Cc : C} (γ : B ⟶ Cc) (hγ : IsPreStep P γ)
    (φC : Cc ⟶ Cc) (φB : B ⟶ B) (hbC : IsBaseIdentity P φC)
    {d : ℕ+} (hdC : P.degFr φC = d)
    (hsq : φB ≫ γ = γ ≫ φC) :
    IsBaseIdentity P φB ∧ P.degFr φB = d := by
  haveI : IsIso (P.Base γ) := hγ.2
  constructor
  · -- base-identity: ★`γ` が pre-step ⟹ `Base γ` は**同型**なので右から消せる
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base φC = P.Base (𝟙 Cc) from hbC,
      P.Base_id, Category.comp_id] at h
    show P.Base φB = P.Base (𝟙 B)
    rw [P.Base_id]
    exact (cancel_mono (P.Base γ)).mp (by rw [h, Category.id_comp])
  · -- 次数: 合成則から
    have h := congrArg P.degFr hsq
    rw [P.degFr_comp, P.degFr_comp, hγ.1, hdC, one_mul, mul_one] at h
    exact h

/-! ## ★★`Definition 1.3, (iii), (d)` の圏同値を**消費する**

★**測定(2026-08-16)**: `coaPreUnderEquiv` / `coaPreOverEquiv` は
`Prop15`(`𝔽_Φ` の場合)と `Prop19`(`𝒞^istr` への移送)で**構成**されているが、
★**下流で一度も消費されていない**(`grep` で確認)。

★**`Proposition 1.10` の残り 2 件((iii) の後半、(vi) の前半)は
どちらもこの圏同値を待っている。** だから**橋を 1 本架ける**。

★**橋の形**: `(𝒞^coa-pre)_{Cc} ⥤ Order(Φ(Cc))^opp` が**充満**であることから、
`Div` の順序関係 `γ の不変量 ≤ γ′ の不変量` を、
★**射の因子分解 `β ≫ γ = γ′`** に変換する。

★**なぜ `op` が要るか**: 関手は反変(`Order(...)^opp` へ行く)。
`F(Z_{γ′}) ⟶ F(Z_γ)` は `op γ′_inv ⟶ op γ_inv`、すなわち
`γ_inv ⟶ γ′_inv` の `op`、すなわち **`γ_inv ≤ γ′_inv`**。
★**「大きいほうから小さいほうへ射がある」の向きが、`op` で 1 回ひっくり返る。**
-/

include P in
/-- ★★**`Definition 1.3, (iii), (d)` の第2の圏同値の初の消費**。

co-angular pre-step の `Div` 不変量に順序関係があれば、**因子分解が存在する**。

★原文(p.36)の
> the portion of assertion (ii) concerning the relationship between Div(γ), Div(γ′)
> implies, in light of the second equivalence of categories of Definition 1.3, (iii), (d),
> that γ′ factors through γ
の「in light of the second equivalence」の実体である。

★**`[IsIso (P.Base γ)]` を instance binder で取っている理由**: `inv (P.Base γ)` が
**文の中**に現れるので、`hγs.2` を証明の中で `haveI` しても遅い。
★**「仮定に含まれている事実を、インスタンスとしても渡す」必要がある** ——
`IsPreStep` は構造体の連言であってインスタンスではないから。
★これは分類表 #5(インスタンス合成)の**設計側の現れ**である。 -/
theorem coaPre_factor_of_mle (G : Frobenioid P) {Cc B B' : C}
    (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (γ' : B' ⟶ Cc) (hγ'c : IsCoAngular P γ') (hγ's : IsPreStep P γ')
    [IsIso (P.Base γ)] [IsIso (P.Base γ')]
    (hle : MLe (Φ.map (inv (P.Base γ)) (P.Div γ) : Φ.val (P.toElem.obj Cc).base)
      (Φ.map (inv (P.Base γ')) (P.Div γ'))) :
    ∃ β : B' ⟶ B, IsCoAngular P β ∧ IsPreStep P β ∧ β ≫ γ = γ' := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv Cc
  set Zγ : Over (⟨Cc⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨Cc⟩ from ⟨γ, hγc, hγs⟩) with hZγ
  set Zγ' : Over (⟨Cc⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨Cc⟩ from ⟨γ', hγ'c, hγ's⟩) with hZγ'
  have hmor : (coaPreOverFunctor P Cc).obj Zγ' ⟶ (coaPreOverFunctor P Cc).obj Zγ :=
    (homOfLE (show (toOrderCat (Φ.map (inv (P.Base γ)) (P.Div γ)))
      ≤ toOrderCat (Φ.map (inv (P.Base γ')) (P.Div γ')) from hle)).op
  let f := (coaPreOverFunctor P Cc).preimage hmor
  refine ⟨f.left.hom, f.left.property.1, f.left.property.2, ?_⟩
  have hw := Over.w f
  exact congrArg (fun t => t.hom) hw

/-! ### ★(vi) の前半の組み立て

原文(p.36)の段取りは 7 段:
1. `Definition 1.3, (i), (a), (b)` で co-angular pre-step `α : B ⟶ A`、`γ : B ⟶ C`
   (`C` は Frobenius-trivial)を得る
2. `C` が Frobenius-trivial なので、各 `d` に base-identity 自己射 `φ_C`(次数 `d`)がある
3. **(ii)** で `γ ≫ φ_C = ψ ≫ γ′` と組み替える(`ψ` は Frobenius 型、`γ′` は co-angular pre-step)
4. **`Div` の関係 ＋ 圏同値**で `γ′` が `γ` を経由する: `β ≫ γ = γ′`
5. `φ_B := ψ ≫ β` と置くと `φ_B ≫ γ = γ ≫ φ_C`
6. `prop_1_10_vi_descend` で `φ_B` は base-identity・次数 `d`
7. よって `B` は quasi-Frobenius-trivial、`α` により `A` は sub-quasi-Frobenius-trivial

★★**ここでは 6・7 を実装し、3–5 が作るべきものを `hstep` として仮定に出す。**
★**こうすると「残りは何か」が宣言の型として確定する** ——
`Prop110iGoal` で使ったのと同じやり方である。
★**未完を散文で書くのではなく、型で書く。**
-/

include P in
/-- **(vi) 前半の 6・7 段** —— `γ` に沿って base-identity 自己射が降りるなら
`B` は quasi-Frobenius-trivial。

★仮定 `hstep` が原文の 3–5 段(**(ii) と圏同値**)が作るべきものである。
★**`prop_1_10_vi_descend` がそれを受けて 6 段を実行する。** -/
theorem prop_1_10_vi_quasi {B Cc : C} (γ : B ⟶ Cc) (hγs : IsPreStep P γ)
    (hC : IsFrobeniusTrivial P Cc)
    (hstep : ∀ (n : ℕ+) (φC : Cc ⟶ Cc), IsBaseIdentity P φC → IsFrobeniusType P φC →
      P.degFr φC = n → ∃ φB : B ⟶ B, φB ≫ γ = γ ≫ φC) :
    IsQuasiFrobeniusTrivial P B := by
  intro n
  obtain ⟨ζ, hdeg, hprop⟩ := hC
  obtain ⟨φB, hsq⟩ := hstep n (ζ n) (hprop n).1 (hprop n).2 (hdeg n)
  obtain ⟨hb, hd⟩ := prop_1_10_vi_descend P γ hγs (ζ n) φB (hprop n).1 (hdeg n) hsq
  exact ⟨φB, hb, hd⟩

include P in
/-- **(vi) 前半の 7 段** —— quasi-Frobenius-trivial な対象へ co-angular pre-step が
あれば sub-quasi-Frobenius-trivial。

★**定義そのものだが、原文の「hence that A is sub-quasi-Frobenius-trivial」の
一言に対応する**ので明示する。 -/
theorem prop_1_10_vi_subQuasi {A B : C} (α : B ⟶ A)
    (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (hq : IsQuasiFrobeniusTrivial P B) : IsSubQuasiFrobeniusTrivial P A :=
  ⟨B, α, hαc, hαs, hq⟩

/-! ### ★MLe の 2 つの補助 —— 穴を塞ぐのに要る

★どちらも `MLe a b := ∃ c, a + c = b` の定義から直に出るが、
**使う形が違う**ので別々に置く。
-/

/-- ★**`x ≤ n • x`**(`n ≥ 1`)。★`n • x = x + (n-1) • x` から。 -/
theorem mle_nsmul_self {M : Type w} [AddCommMonoid M] {n : ℕ} (hn : 0 < n) (x : M) :
    MLe x (n • x) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  exact ⟨m • x, by rw [succ_nsmul, add_comm]⟩

/-- ★**加法準同型は `MLe` を保つ**。 -/
theorem MLe.map {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N)
    {a b : M} (h : MLe a b) : MLe (f a) (f b) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨f c, by rw [← map_add, hc]⟩

/-! ### ★★(vi) の最後の穴を塞ぐ —— 3–5 段

★**道具は全部そろっている**: (ii)(`prop_1_10_ii`)、`Div` の公式
(`prop_1_10_ii_Div_formula`)、橋(`coaPre_factor_of_mle`)。
★**残っていたのは `Φ.map` の記帳だけ**だった。

★**記帳の中身**:
- 四角形 `β′ ≫ α′ = γ ≫ φ_C` から `Base γ = Base β′ ≫ Base α′`(`φ_C` は base-identity)
- `Div` の公式から `Φ.map (Base β′) (Div α′) = d • Div γ`
- 両辺に `Φ.map (inv (Base β′))` を当てて ★**`Div α′ = d • Φ.map (inv (Base β′)) (Div γ)`**
- `x ≤ d • x` なので `MLe (Φ.map (inv (Base β′)) (Div γ)) (Div α′)`
- `Φ.map (inv (Base α′))` で押し出すと、★**橋が要求する形**になる
-/

include P in
/-- ★★**(vi) 前半の 3–5 段** —— `prop_1_10_vi_quasi` が要求する `hstep` を作る。

★これで `Proposition 1.10, (vi)` の前半の穴が塞がる。 -/
theorem prop_1_10_vi_step (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {B Cc : C} (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (n : ℕ+) (φC : Cc ⟶ Cc) (hbC : IsBaseIdentity P φC) (hftC : IsFrobeniusType P φC)
    (hdC : P.degFr φC = n) :
    ∃ φB : B ⟶ B, φB ≫ γ = γ ≫ φC := by
  -- 3 段: (ii) で組み替える
  obtain ⟨Y, β', α', hβ', hdβ', hα's, hsq⟩ := prop_1_10_ii P G.core γ hγs φC hftC
  haveI hbβ' : IsIso (P.Base β') := hβ'.2
  haveI hbα' : IsIso (P.Base α') := hα's.2
  haveI hbγ : IsIso (P.Base γ) := hγs.2
  -- `Div` の公式
  have hdiv : Φ.map (P.Base β') (P.Div α') = (P.degFr φC : ℕ) • P.Div γ :=
    prop_1_10_ii_Div_formula P γ φC β' α' hβ' hftC hsq
  -- `Div α′ = d • Φ.map (inv (Base β′)) (Div γ)`
  have hdiv' : P.Div α' = (P.degFr φC : ℕ) • Φ.map (inv (P.Base β')) (P.Div γ) := by
    have := congrArg (Φ.map (inv (P.Base β'))) hdiv
    rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at this
    exact this
  -- `MLe` を作る
  have hmle0 : MLe (Φ.map (inv (P.Base β')) (P.Div γ)) (P.Div α') := by
    rw [hdiv']
    exact mle_nsmul_self (P.degFr φC).pos _
  -- 橋が要求する形へ押し出す
  -- ★`inv` の下を書き換えると motive エラー(分類表 #2)になるので、
  --   `inv_eq_of_hom_inv_id` で**外から**特徴づける
  have hb : P.Base β' ≫ P.Base α' = P.Base γ := by
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base φC = P.Base (𝟙 Cc) from hbC,
      P.Base_id, Category.comp_id] at h
    exact h
  have hbase : inv (P.Base γ) = inv (P.Base α') ≫ inv (P.Base β') := by
    refine CategoryTheory.IsIso.inv_eq_of_hom_inv_id ?_
    rw [← hb, Category.assoc, ← Category.assoc (P.Base α'), IsIso.hom_inv_id,
      Category.id_comp, IsIso.hom_inv_id]
  have hmle : MLe (Φ.map (inv (P.Base γ)) (P.Div γ))
      (Φ.map (inv (P.Base α')) (P.Div α')) := by
    rw [hbase, Φ.map_comp]
    exact MLe.map _ hmle0
  -- 4 段: 橋で因子分解を得る
  obtain ⟨β, _, _, hβγ⟩ :=
    coaPre_factor_of_mle P G γ hγc hγs α' (prop_1_4_i P α' (fun Y' _ => hiso Y')) hα's hmle
  -- 5 段: `φ_B := β′ ≫ β`
  exact ⟨β' ≫ β, by rw [Category.assoc, hβγ, hsq]⟩

include P in
/-- ★★**`Proposition 1.10, (vi)` 前半の完成形** ——
isotropic 型の Frobenioid において、Frobenius-trivial な対象へ co-angular pre-step で
つながる対象は quasi-Frobenius-trivial であり、そこへ co-angular pre-step で
つながる対象は sub-quasi-Frobenius-trivial。

★原文(p.36)の 7 段をすべて実装した形である。
★**`Definition 1.3, (i), (a), (b)` から `α`・`γ`・`C` を得る 1・2 段は仮定に出している**
——それは `𝒞^istr` が Frobenioid であること(`Proposition 1.9, (v)`)から供給されるもので、
**この定理の内容ではない**。 -/
theorem prop_1_10_vi_first (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B Cc : C} (α : B ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (γ : B ⟶ Cc) (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (hC : IsFrobeniusTrivial P Cc) :
    IsSubQuasiFrobeniusTrivial P A :=
  prop_1_10_vi_subQuasi P α hαc hαs
    (prop_1_10_vi_quasi P γ hγs hC
      (fun n φC hbC hftC hdC =>
        prop_1_10_vi_step P G hiso γ hγc hγs n φC hbC hftC hdC))

include P in
/-- ★★★**`α`・`γ`・`C` を `Definition 1.3, (i), (a), (b)` から導く**（2026-08-16）。

★★**監査で「原文は導いているのに仮定に逃がしている」と指摘された段**である。

★**導き方**（3 段）:
1. `baseSurj` を `Base A` に当て、**Frobenius-trivial な `A₀` と底の同型**を得る
2. `preStepSpan` をその同型に当て、**pre-step の span** `A ← X → A₀` を得る
3. ★`𝒞` が isotropic 型なので、`Proposition 1.4, (i)` により
   **すべての射が co-angular** —— span の両辺が co-angular pre-step になる

★★**これで `prop_1_10_vi_first` の仮定がすべて供給される**。 -/
theorem prop_1_10_vi_data (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X) (A : C) :
    ∃ (B Cc : C) (α : B ⟶ A) (γ : B ⟶ Cc),
      IsCoAngular P α ∧ IsPreStep P α ∧ IsCoAngular P γ ∧ IsPreStep P γ ∧
      IsFrobeniusTrivial P Cc := by
  obtain ⟨A₀, hA₀ft, ⟨e⟩⟩ := G.core.baseSurj (P.toElem.obj A).base
  haveI : IsIso e.inv := inferInstance
  obtain ⟨X, φ, ψ, hφ, hψ, -⟩ := G.core.preStepSpan A₀ A e.hom e.isIso_hom
  exact ⟨X, A₀, ψ, φ,
    prop_1_4_i P ψ (fun Y _ => hiso Y), hψ,
    prop_1_4_i P φ (fun Y _ => hiso Y), hφ, hA₀ft⟩

include P in
/-- ★★★**`Proposition 1.10, (vi)` 前半の「型」の形**（2026-08-16）。

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

★★**監査で「「type」の ∀ が無い」と指摘された形**である。
`prop_1_10_vi_data`（`Definition 1.3, (i), (a), (b)` からの導出）と
`prop_1_10_vi_first`（本体）を合成するだけで出る。

★**isotropic 型の Frobenioid は sub-quasi-Frobenius-trivial 型である** ——
★`𝒞^istr` は isotropic 型なので、これを `istrPre` に当てれば原文の形になる。 -/
theorem prop_1_10_vi_ofType (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X) :
    IsOfSubQuasiFrobeniusTrivialType P := by
  intro A
  obtain ⟨B, Cc, α, γ, hαc, hαs, hγc, hγs, hC⟩ := prop_1_10_vi_data P G hiso A
  exact prop_1_10_vi_first P G hiso α hαc hαs γ hγc hγs hC


/-! ## ★★(iii) —— 第1の圏同値も消費する

★**前段では第2(スライス)の圏同値を消費した**(`coaPre_factor_of_mle`)。
(iii) では **第1(コスライス)の圏同値**を使う。
関手は `_A(𝒞^coa-pre) ⥤ Order(Φ(A))`、`Z ↦ Div (Z.hom)`。

★**本質的全射性**は「任意の `c ∈ Φ(A)` は、ある co-angular pre-step の `Div` に**同型**」
を与える。★**`Order` は前順序圏なので、同型 = 両向きの `MLe`。**
★**そこに `mle_antisymm`(integral ＋ sharp)を当てると等号になる。**
★★**`divisorial` の 2 つの条件(integral と sharp)が、ここで一緒に効く。**
-/

include P in
/-- ★★**第1の圏同値の初の消費** —— `Φ(A)` の任意の元は
co-angular pre-step の `Div` として実現できる。

★**前順序圏の同型は両向きの `MLe`** なので、`mle_antisymm` で等号に直す。 -/
theorem coaPre_realize (G : Frobenioid P) (A : C)
    (c : Φ.val (P.toElem.obj A).base) :
    ∃ (X : C) (ψ : A ⟶ X), IsCoAngular P ψ ∧ IsPreStep P ψ ∧ P.Div ψ = c := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := coaPreUnderFunctor P A)
    (toOrderCat c)
  refine ⟨Z.right.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
  refine mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 ?_ ?_
  · exact leOfHom e.hom
  · exact leOfHom e.inv

/-! ### ★(iii) の全射性 —— `n` 倍が全射であること

★**組み立て**:
1. `c ∈ Φ(A)` を `coaPre_realize` で co-angular pre-step `ψ : A ⟶ X`(`Div ψ = c`)に実現
2. `IsPerfectObj A` の (a) で、`A` と `X` の上に**次数 `n` の Frobenius 型射**
   `φ₁ : B₀ ⟶ A`、`φ₂ : B₂ ⟶ X` を取る(`ψ` が pre-step なので `X` は `A` と base-isomorphic)
3. `IsPerfectObj A` の (b) で `ψ` を降ろし、pre-step `ψ₀ : B₀ ⟶ B₂` と
   四角形 `φ₁ ≫ ψ = ψ₀ ≫ φ₂` を得る
4. ★**核心の等式**(`prop_1_10_iii_div_nsmul'`)から `Φ.map (Base φ₁) (Div ψ) = n • Div ψ₀`
5. `Φ.map (inv (Base φ₁))` を当てて ★**`c = n • Φ.map (inv (Base φ₁)) (Div ψ₀)`**

★**「A は perfect」から「Φ(A) は perfect」への橋が、これで渡り切る。**
-/

include P in
/-- ★★**(iii) の全射性** —— `𝒞` の対象が perfect なら、`Φ(A)` で `n` 倍は全射。 -/
theorem prop_1_10_iii_nsmul_surjective (G : Frobenioid P) {A : C}
    (hperf : IsPerfectObj P A) (n : ℕ+) :
    Function.Surjective (fun a : Φ.val (P.toElem.obj A).base => (n : ℕ) • a) := by
  intro c
  obtain ⟨X, ψ, hψc, hψs, hψd⟩ := coaPre_realize P G A c
  haveI hbψ : IsIso (P.Base ψ) := hψs.2
  have hAX : BaseIsomorphic P A X := ⟨asIso (P.Base ψ)⟩
  have hAA : BaseIsomorphic P A A := ⟨Iso.refl _⟩
  obtain ⟨B₀, φ₁, hφ₁, hd₁⟩ := (hperf n).1 A hAA
  obtain ⟨B₂, φ₂, hφ₂, hd₂⟩ := (hperf n).1 X hAX
  haveI hbφ₁ : IsIso (P.Base φ₁) := hφ₁.2
  -- `B₀`・`B₂` は `A` と base-isomorphic
  have hB₀A : BaseIsomorphic P B₀ A := ⟨asIso (P.Base φ₁)⟩
  have hB₂A : BaseIsomorphic P B₂ A := by
    haveI : IsIso (P.Base φ₂) := hφ₂.2
    exact ⟨asIso (P.Base φ₂) ≪≫ (asIso (P.Base ψ)).symm⟩
  obtain ⟨ψ₀, ⟨hψ₀s, hsq⟩, _⟩ :=
    (hperf n).2 B₀ A B₂ X φ₁ φ₂ hφ₁ hd₁ hφ₂ hd₂ hB₀A hB₂A ψ hψs
  -- 核心の等式
  have hcore : Φ.map (P.Base φ₁) (P.Div ψ) = (n : ℕ) • P.Div ψ₀ :=
    prop_1_10_iii_div_nsmul' P φ₁ φ₂ ψ₀ ψ hφ₁.1.2 hφ₂.1.2 hd₂ hsq
  refine ⟨Φ.map (inv (P.Base φ₁)) (P.Div ψ₀), ?_⟩
  have h := congrArg (Φ.map (inv (P.Base φ₁))) hcore
  rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at h
  show (n : ℕ) • Φ.map (inv (P.Base φ₁)) (P.Div ψ₀) = c
  rw [← h, hψd]

/-! ### ★★(iii) の単射性 —— **親の前段の見立ては誤っていた。そして次の壁が見えた**

前段で親はこう書いた:
> 単射性は `IsPerfectObj` の `∃!` から出るはずで、`divisorial` からは出ない
> (`Gp M` の捻れを排除する条件が無い)。

★★**誤りである。紙の上では `divisorial` から出る。** `Gp M` は**捻れなし**である:

`x ∈ Gp M` が `n • x = 0`(`n ≥ 1`)を満たすとする。
- `n • x = 0 = toGp 0` は `M` の像に入るので、★**saturated** から `x = toGp m` と書ける
- 同様に `n • (-x) = 0` も像に入るので `-x = toGp m′`
- `toGp (m + m′) = x + (-x) = 0 = toGp 0` で、★**integral**(`toGp` が単射)から `m + m′ = 0`
- よって `m` は可逆、★**sharp** から `m = 0`、したがって `x = 0`

★★**`divisorial` の 3 条件(saturated / integral / sharp)が、この 1 つの主張で
一度に効いている。**

## ★★しかし Lean では書けなかった —— 次の壁

★**`Gp M := AddLocalization (⊤ : AddSubmonoid M)` に群構造が無い。**
上の議論は `x - y` と `-x` を使うが、我々の `Gp M` は
**`AddCommMonoid` としてしか実装されていない**(`HSub`・`Neg` のインスタンスが無い)。
★**`saturated` の定義自体が `Gp M` の元を量化しているのに、
`Gp M` を群として使う準備ができていない。**

★**したがって単射性はここで止まる。** `sorry` は置かない(`Found/` の規律)。
★**次に当たるときの仕事は明確**: `Gp M` に `AddCommGroup` インスタンスを付ける
(`neg (mk a s) := mk s a`、`S = ⊤` なので常に定義できる)。

## ★★記録: 否定的な予測は検証が甘くなる

親は前段で「`divisorial` からは出ない」と書いた。★**1 本も試さずに書いた。**
実際に試したら**紙の上では出た**(そして Lean では別の理由で止まった)。

★★**否定的な予測(「〜からは出ない」)は、肯定的な予測より検証が甘くなる。**
「出る」は 1 本示せば済むが、「出ない」は**全ての道を塞いだこと**を要する。
★**我々は「証明できなかった」と「反例がある」を区別してきたが、
「試していない」と「試して駄目だった」の区別も要る。**
★今回の親の発言は **3 つ目のカテゴリ(試していない)**だった。
-/

/-! ### ★★壁は import 1 行だった —— `Gp M` の群構造

前段で親は「`Gp M` に群構造が無いので単射性は書けない」と記録した。
★**確かめたら、mathlib に**あった**。**

`Mathlib/GroupTheory/MonoidLocalization/GrothendieckGroup.lean`:
```
abbrev AddGrothendieckGroup : Type _ := AddLocalization (⊤ : AddSubmonoid M)
instance instAddCommGroup : AddCommGroup (AddGrothendieckGroup M)
```
★★**我々の `Gp M := AddLocalization (⊤ : AddSubmonoid M)` と同じもの**である。
**import していなかっただけ。**

★**これは `comp_div` と同じ形の壁である**(第24段): 症状は「難しい」に見えたが、
原因は★**「1 行足りない」**だった。
★**前回は `@[simp]` 補題が 1 本、今回は `import` が 1 行。**
★★**「壁の高さ」と「原因の大きさ」は関係がない。**

★**そして 2 回とも、原因は「揃っているはずのものを並べて見る」で見つかった** ——
前回は `@[simp]` 補題の一覧、今回は mathlib のファイル一覧
(`ls` して `GrothendieckGroup.lean` が目に入った)。
-/

/-- ★★**divisorial なモノイドでは `n` 倍が単射**(`n ≥ 1`)。

★`Gp M` が**捻れなし**であることの帰結。saturated / integral / sharp の 3 つが要る。 -/
theorem nsmul_injective_of_isDivisorial {M : Type w} [AddCommMonoid M]
    (hd : IsDivisorial M) {n : ℕ} (hn : 0 < n) :
    Function.Injective (fun a : M => n • a) := by
  obtain ⟨⟨hint, hsat, _⟩, hsh⟩ := hd
  -- ★`AddLocalization.mk_add`(mathlib、`mk_mul` の `to_additive` 版)で足し算を開く
  have htoGp_add : ∀ x y : M, toGp M (x + y) = toGp M x + toGp M y := by
    intro x y
    rw [toGp, toGp, toGp, AddLocalization.mk_add]
    congr 1
    exact (add_zero (0 : (⊤ : AddSubmonoid M))).symm
  have htoGp_nsmul : ∀ (k : ℕ) (x : M), toGp M (k • x) = k • toGp M x := by
    intro k x
    induction k with
    | zero => simp [toGp, AddLocalization.mk_zero]
    | succ m ih => rw [succ_nsmul, succ_nsmul, htoGp_add, ih]
  intro a b hab
  set x : Gp M := toGp M a - toGp M b with hx
  have hnx : n • x = 0 := by
    rw [hx, smul_sub, ← htoGp_nsmul, ← htoGp_nsmul, sub_eq_zero]
    exact congrArg (toGp M) hab
  have hmem : ∀ y : Gp M, n • y = 0 → y ∈ Set.range (toGp M) := by
    intro y hy
    refine hsat y n hn ?_
    rw [hy]
    exact ⟨0, by simp [toGp, AddLocalization.mk_zero]⟩
  obtain ⟨m, hm⟩ := hmem x hnx
  obtain ⟨m', hm'⟩ := hmem (-x) (by rw [smul_neg, hnx, neg_zero])
  have hsum : m + m' = 0 := by
    refine hint ?_
    rw [htoGp_add, hm, hm', add_neg_cancel]
    simp [toGp, AddLocalization.mk_zero]
  have hm0 : m = 0 := hsh m ⟨⟨m, m', hsum, by rw [add_comm]; exact hsum⟩, rfl⟩
  have hx0 : x = 0 := by rw [← hm, hm0]; simp [toGp, AddLocalization.mk_zero]
  refine hint ?_
  rwa [hx, sub_eq_zero] at hx0

include P in
/-- ★★**`Proposition 1.10, (iii)` の第1主張 完成** ——
`𝒞` の対象 `A` が perfect なら、`Φ(A)` は perfect なモノイド。

★**全射性**は `IsPerfectObj`(原文が名指す `Definition 1.2, (iv)`)から、
★**単射性**は `divisorial`(saturated ＋ integral ＋ sharp)から。

★★**原文は「the fact that Φ(A) is perfect follows immediately … from the fact that
A is perfect」と書くが、実際には 2 つの別々の源から来ている** ——
全射性だけが `A` の perfect 性から来て、**単射性は `Φ` が divisorial であることから来る**。
★**原文の「from the fact that A is perfect」は半分しか説明していない。** -/
theorem prop_1_10_iii_perfect (G : Frobenioid P) {A : C} (hperf : IsPerfectObj P A) :
    IsPerfectMonoid (Φ.val (P.toElem.obj A).base) := fun n =>
  ⟨nsmul_injective_of_isDivisorial (P.divisorial _) n.pos,
   prop_1_10_iii_nsmul_surjective P G hperf n⟩

include P in
/-- ★**(iii) 第1主張の「型」版** —— `𝒞` が of perfect type なら
`Φ` の像のモノイドはすべて perfect。★**原文の言い回しそのもの。** -/
theorem prop_1_10_iii_image_perfect (G : Frobenioid P) (hpt : IsOfPerfectType P) (A : C) :
    IsPerfectMonoid (Φ.val (P.toElem.obj A).base) :=
  prop_1_10_iii_perfect P G (hpt A)

include P in
/-- ★★**(iii) 第1主張 ——「the monoids in the image of Φ」の本当の範囲**。

★★**取り下げの原因の1つを埋める**(2026-08-16)。
`prop_1_10_iii_image_perfect` は `𝒞` の対象の底、すなわち
**`Base` の像**についてしか述べていなかった。★原文の「the image of Φ」は
`Φ : 𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫` の像であり、★**`Ob(𝒟)` 全体に渡る。**

★**渡し方**: `Definition 1.3, (i), (a)`(`baseSurj`)が、任意の `Y : 𝒟` に対し
`Base(A) ≅ Y` なる(しかも Frobenius-trivial な)`A : 𝒞` を与える。
`Φ.map` はその同型に沿って全単射(`MonoidOn.map_bijective_of_iso`)であり、
perfect 性は加法的全単射で移る(`isPerfectMonoid_of_bijective`)。

★★**ここで初めて `𝒟` 側と `𝒞` 側の隔たりが埋まった** ——
`Proposition 1.10` の他の条はすべて `𝒞` の中の話だが、
★**(iii) の第1主張だけは `𝒟` の全対象について述べている。**
その差を原文は書いていない(「the image of Φ」の一語で済ませている)。 -/
theorem prop_1_10_iii_image_perfect_of_base (F : FrobenioidCore P) (G : Frobenioid P)
    (hpt : IsOfPerfectType P) (Y : D) : IsPerfectMonoid (Φ.val Y) := by
  obtain ⟨A, -, ⟨e⟩⟩ := F.baseSurj Y
  exact isPerfectMonoid_of_bijective (Φ.map e.symm.hom) (Φ.map_bijective_of_iso e.symm)
    (prop_1_10_iii_image_perfect P G hpt A)

/-! ## ★(iii) の「moreover」—— `𝒪^▷(A)` と `𝒪^×(A)` が perfect

原文 (FrdI p.35):
> immediately [cf. Remark 1.1.1] from the fact that A is perfect [cf. Definition 1.2,

★**主張そのものの行(「then the monoids O▷(A) and O×(A) are perfect.」)は
引用できない**(事故 #2 の再来。`▷` が layout モードで消え 15/38 文字で停止)。
★**`▷` は事故 #2 で既に記録した文字である** —— 同じ文字が別の場所で再び落ちた。
証明側の隣接行を引く。

★**まず語彙の問題**: `𝒪^▷(A)` は `Submonoid (End A)` であり**乗法**モノイドである。
我々の `IsPerfectMonoid` は**加法**モノイド用(`∀ n, Bijective (n • ·)`)。
★**乗法版が要る。** 原文は両方を同じ「perfect」で呼ぶが、
★**`Φ(A)` は加法、`𝒪^▷(A)` は乗法**で、**同じ語が 2 つの構造を指している。**

★**原文の証明**(p.35、目視):
> the fact that the monoids O▷(A) and O×(A) are perfect follows immediately from the
> fact that A is perfect [cf. Definition 1.2, (iv), **applied to the base-identity
> endomorphisms of Frobenius type of the Frobenius-trivial object A**] and Frobenius-normalized.
-/

/-- ★**乗法モノイドの perfect** —— `n` 乗が全単射。

★`IsPerfectMonoid`(加法)の乗法版。★**原文が同じ語で呼ぶ 2 つを、我々は別の述語で書く。** -/
def IsPerfectMulMonoid (M : Type*) [Monoid M] : Prop :=
  ∀ n : ℕ+, Function.Bijective (fun a : M => a ^ (n : ℕ))

/-- ★**非退化(下)**: 自明モノイドは perfect。 -/
theorem isPerfectMulMonoid_punit : IsPerfectMulMonoid PUnit :=
  fun _ => ⟨fun _ _ _ => rfl, fun _ => ⟨PUnit.unit, rfl⟩⟩

/-- ★**非退化(上)**: `Multiplicative ℕ` は perfect でない —— 2 乗写像が全射でない
(`ofAdd 1` は平方の像に入らない)。 -/
theorem not_isPerfectMulMonoid_nat : ¬ IsPerfectMulMonoid (Multiplicative ℕ) := by
  intro h
  obtain ⟨a, ha⟩ := (h 2).2 (Multiplicative.ofAdd 1)
  have h1 : Multiplicative.toAdd (a ^ (2 : ℕ)) = 1 := congrArg Multiplicative.toAdd ha
  rw [show (2 : ℕ) = 1 + 1 from rfl, pow_add, pow_one] at h1
  have h2 : Multiplicative.toAdd a + Multiplicative.toAdd a = 1 := h1
  omega

/-! ### ★★(iii) の moreover —— 残る穴を型で書く

★**構成**(原文の 3 段、目視で確認):
1. `A` が Frobenius-trivial なので、各 `n` に base-identity・Frobenius 型の
   自己射 `ζ n`(次数 `n`)がある
2. `IsPerfectObj` の (b) を `φ₁ = φ₂ = ζ n`、`ψ′ := α ∈ 𝒪^▷(A)` に当てると、
   ★**一意の pre-step `β` が `ζ n ≫ α = β ≫ ζ n` を満たす**
3. **Frobenius-normalized** を `φ := ζ n`、`α := β` に当てると
   `ζ n ≫ β^n = β ≫ ζ n = ζ n ≫ α`

★★**残る穴はただ 1 つ**: 上から `β^n = α` を結論するには
★**`ζ n` で左から消せること(mono)**が要る。
★**`Definition 1.3` は pre-step が mono であること(`preStepMono`)しか言っていない。**
`ζ n` は次数 `n` なので pre-step ではない。

★★**したがってこれが「原文の `follows immediately` が使っている、
書かれていない事実」である。** 3 通りのどれかである:
  (a) Frobenioid では Frobenius 型射は常に mono —— **原文にも我々の公理にも無い**
  (b) isotropic 型ではそうなる —— **未確認**
  (c) 別の議論で mono を回避できる —— **未発見**

★**`sorry` は置かない。** 3 通りのどれかを確かめるのが次の仕事である。
★**「試していない」と明記する**(前段で立てた 3 分類の第 3 category)。

★**語彙は先に置いた**(`IsPerfectMulMonoid`)ので、
**残りは上の 1 点だけ**である。
-/

/-! ### ★★穴が 1 段浅くなった —— mono ではなく `faithfulUpToUnits` だった

前段で親は「残る穴は `ζ n` が mono かどうか」と書いた。★**それは違った。**

★**`Definition 1.3, (vi)`(`faithfulUpToUnits`)が使える**:
`α` と `β^n` はどちらも `𝒪^▷(A)` の元(base-identity ＋ linear)であり、
- **base-equivalent**: どちらも `Base = 𝟙`
- **metrically equivalent**: `ζ n ≫ α = ζ n ≫ β^n` に `Div_comp` を当てると、
  ★**`ζ n` が Frobenius 型(`Div = 0`)かつ base-identity(`Φ.map 𝟙 = id`)**なので
  両辺が `Div α`・`Div (β^n)` に潰れる ⟹ **`Div α = Div (β^n)`**

★したがって `faithfulUpToUnits` から ★**`α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)**が出る。

★★**穴は「mono か」から「単元 `u` を吸収できるか」に変わった。**
★**これは前進である**: mono は**公理に無い**が、単元の吸収は
**`𝒪^×(A)` が perfect であること**(＝ (iii) の moreover のもう半分)から出るはずで、
★**同じ条の中で閉じる可能性がある。**

★**測定**: 「残る穴」を型で書いておいたおかげで、★**穴の位置が動いたことがすぐ分かった。**
散文で「あとは mono が要る」と書いていたら、`faithfulUpToUnits` を試すきっかけが無かった。
-/

include P in
/-- ★★**`α` と `β^n` は単元の違いしかない** —— `faithfulUpToUnits` の帰結。

★前段で「mono が要る」と書いた穴が、ここまで浅くなった。 -/
theorem prop_1_10_iii_otri_upToUnit (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X)
    (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (α β : End A) (hαb : IsBaseIdentity P α) (hαl : IsLinear P α)
    (hβb : IsBaseIdentity P β) (hβl : IsLinear P β)
    (heq : (ζn : A ⟶ A) ≫ (α : A ⟶ A) = (ζn : A ⟶ A) ≫ (β : A ⟶ A)) :
    ∃ u : End A, u ∈ OTimes P A ∧ (α : A ⟶ A) = (β : A ⟶ A) ≫ (u : A ⟶ A) := by
  -- `Div α = Div β`
  have hd : P.Div (α : A ⟶ A) = P.Div (β : A ⟶ A) := by
    have h := congrArg P.Div heq
    rw [P.Div_comp, P.Div_comp, show P.Div (ζn : A ⟶ A) = 0 from hζf.1.2,
      show P.Base (ζn : A ⟶ A) = P.Base (𝟙 A) from hζb, P.Base_id] at h
    simpa using h
  -- `α`・`β` は co-angular pre-step
  have hpre : ∀ γ : End A, IsBaseIdentity P γ → IsLinear P γ → IsPreStep P (γ : A ⟶ A) := by
    intro γ hb hl
    refine ⟨hl, ?_⟩
    show IsIso (P.Base (γ : A ⟶ A))
    rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
    infer_instance
  have hcoa : ∀ γ : End A, IsCoAngular P (γ : A ⟶ A) :=
    fun γ => prop_1_4_i P (γ : A ⟶ A) (fun X _ => hiso X)
  exact F.faithfulUpToUnits (α : A ⟶ A) (β : A ⟶ A)
    (show P.Base (α : A ⟶ A) = P.Base (β : A ⟶ A) by
      rw [show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hαb,
        show P.Base (β : A ⟶ A) = P.Base (𝟙 A) from hβb])
    hd (hcoa α) (hpre α hαb hαl) (hcoa β) (hpre β hβb hβl)

/-! ### ★★単元の吸収が難しい理由が構造から分かった

★**残る穴**: `α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)から `α = γ^n` を出したい。
`𝒪^×(A)` も perfect なら `u = v^n` と書けるが、
★**`β^n ≫ v^n = (β ≫ v)^n` には `β` と `v` の可換性が要る。**

★★**そして `𝒪^▷(A)` は一般に可換ではない。**

原典の序文(p.1、目視):
> kind of Frobenioid "FM" is the **non-commutative** monoid given by forming the
> "semi-direct product monoid" of a given commutative monoid M with the multi-
> plicative monoid N≥1

★**`𝔽_M` は明示的に非可換**である(`M` と `ℕ≥1` の半直積)。
ただし `𝒪^▷(A)` は **linear**(次数 1)な base-identity 自己射の集まりなので、
半直積の `M` 側にあたり、**素朴には可換に見える**。

★★**しかし一般の Frobenioid ではそうとは限らない。** 我々が持っている情報は
`faithfulUpToUnits`(`Definition 1.3, (vi)`)だけで、それが言うのは
★**「base と `Div` が一致する co-angular pre-step は単元の違いしかない」** ——
すなわち ★**`𝒪^▷(A)` は単元を法としてのみ `Φ(A)` に埋まる。**

★★**したがって「単元の吸収」は、この命題の中で最も構造に近い場所にある。**
`𝒪^×(A)` が何であるかは `Definition 1.3` が直接には規定していない。

★**次の仕事(3 通り、まだ試していない)**:
  (a) `𝒪^▷(A)` の可換性を Frobenioid の公理から導く —— **未確認**
  (b) `𝒪^×(A)` が中心に入ることを示す —— **未確認**
  (c) 可換性を使わない別の議論 —— **未発見**

★**「試していない」と明記する**(3 分類の第 3 カテゴリ)。
★**ただし前段と違い、今回は「なぜ難しいか」が構造から言えている** ——
★**原典が序文で「non-commutative」と宣言している当のものに触れている。**
-/

/-! ### ★★選択肢 (a) を試した —— `𝒪^▷(A)` は**単元を法として**可換

前段で 3 通りの「まだ試していない」を挙げた。★**(a) を試した。**

★**結果**: `𝒪^▷(A)` は**単元を法として可換**である。**完全な可換性は出ない。**

★**証明**: `α, β ∈ 𝒪^▷(A)` に対し `α ≫ β` と `β ≫ α` は
- **base-equivalent**: どちらも `Base = 𝟙`
- **metrically equivalent**: `Div_comp` で
  `Div (α ≫ β) = Φ.map 𝟙 (Div β) + 1 • Div α = Div β + Div α`、
  `Div (β ≫ α) = Div α + Div β`。★**`Φ(A)` が可換なので一致する。**

したがって `faithfulUpToUnits` から ★**`α ≫ β = (β ≫ α) ≫ u`(`u ∈ 𝒪^×(A)`)**。

★★**これは「試して、部分的に出た」である**(3 分類の第 2 カテゴリへ移った)。
★**単元の吸収に必要な完全な可換性は、ここからは出ない。**

★★**測定: 可換性の破れは `Φ(A)` ではなく `𝒪^×(A)` に局在している。**
`Div` は可換に振る舞い、ずれはすべて**単元**に吸い込まれる。
★**原典が序文で「non-commutative」と言うのは半直積の `ℕ≥1` 側の話だが、
`𝒪^▷(A)` の非可換性は `𝒪^×(A)` に由来する** —— **別の出どころ**である。
-/

include P in
/-- ★★**`𝒪^▷(A)` は単元を法として可換**。

★`Div` は可換に振る舞う(`Φ(A)` が可換だから)ので、
ずれはすべて `faithfulUpToUnits` の単元に吸い込まれる。 -/
theorem otri_comm_upToUnit (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X) (α β : OTri P A) :
    ∃ u : End A, u ∈ OTimes P A ∧
      ((α : End A) ≫ (β : End A)) = ((β : End A) ≫ (α : End A)) ≫ (u : A ⟶ A) := by
  have hpre : ∀ γ : OTri P A, IsPreStep P (γ : End A) := by
    intro γ
    refine ⟨γ.2.2, ?_⟩
    show IsIso (P.Base (γ : End A))
    rw [show P.Base (γ : End A) = P.Base (𝟙 A) from γ.2.1, P.Base_id]
    infer_instance
  have hcoa : ∀ f : A ⟶ A, IsCoAngular P f := fun f => prop_1_4_i P f (fun X _ => hiso X)
  have hbase : ∀ γ : OTri P A, P.Base (γ : End A) = 𝟙 _ := by
    intro γ
    rw [show P.Base (γ : End A) = P.Base (𝟙 A) from γ.2.1, P.Base_id]
  -- `Div` が一致する
  have hd : P.Div ((α : End A) ≫ (β : End A)) = P.Div ((β : End A) ≫ (α : End A)) := by
    rw [P.Div_comp, P.Div_comp, hbase α, hbase β, Φ.map_id, Φ.map_id,
      show P.degFr (β : End A) = 1 from β.2.2, show P.degFr (α : End A) = 1 from α.2.2]
    simp [add_comm]
  -- base が一致する
  have hb : P.Base ((α : End A) ≫ (β : End A)) = P.Base ((β : End A) ≫ (α : End A)) := by
    rw [P.Base_comp, P.Base_comp, hbase α, hbase β]
  exact F.faithfulUpToUnits _ _ hb hd (hcoa _)
    (IsPreStep.comp P (hpre α) (hpre β)) (hcoa _) (IsPreStep.comp P (hpre β) (hpre α))

/-! ### ★選択肢 (b) の一部 —— `β^n` が同型なら `β` も同型

★**`𝒪^×(A)` の側を先に片づけようとして出た副産物**である。

★**主張**: isotropic 型で、`β` が pre-step、`β^n` が同型なら `β` は同型。
★**証明**: `β^n` が同型 ⟹ `Div (β^n) = 0`。`β` は base-identity 型の pre-step なので
`Div (β^n) = n • Div β`(★`Div_comp` を `n` 回、`Base β = 𝟙` と `degFr β = 1` で潰す)。
★**sharp から `Div β = 0`** ⟹ `β` は isometric。isotropic から co-angular なので
LB-invertible な pre-step、`Proposition 1.4, (iii)` で**同型**。

★★**これは `𝒪^▷(A)` の中で `𝒪^×(A)` が「`n` 乗で閉じている」ことの逆向きである。**
★**単元の吸収そのものは片づかないが、「どの元が単元か」の判定は 1 つ増えた。**
-/

include P in
/-- ★`γ^k` の 3 つの不変量を一度に。★**`End A` の乗法は `x * y = y ≫ x`**
(mathlib の規約)なので、`pow_succ` を開くと `γ ≫ γ^m` になる。 -/
theorem otri_pow_invariants {A : C} (γ : End A) (hb : IsBaseIdentity P γ) (hl : IsLinear P γ) :
    ∀ k : ℕ, P.Base ((γ ^ k : End A) : A ⟶ A) = P.Base (𝟙 A) ∧
      P.degFr ((γ ^ k : End A) : A ⟶ A) = 1 ∧
      P.Div ((γ ^ k : End A) : A ⟶ A) = k • P.Div (γ : A ⟶ A) := by
  intro k
  induction k with
  | zero => exact ⟨rfl, by simpa using P.degFr_id A, by simpa using P.Div_id A⟩
  | succ m ih =>
      obtain ⟨ihb, ihd, ihv⟩ := ih
      have hmul : ((γ ^ (m + 1) : End A) : A ⟶ A) = (γ : End A) ≫ (γ ^ m : End A) := by
        rw [pow_succ]; rfl
      refine ⟨?_, ?_, ?_⟩
      · rw [hmul, P.Base_comp, ihb, P.Base_id, Category.comp_id]
        rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
      · rw [hmul, P.degFr_comp, ihd, hl, one_mul]
      · rw [hmul, P.Div_comp, ihv, ihd,
          show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id, Φ.map_id]
        simp [succ_nsmul, add_comm]

include P in
/-- ★★**`β^n` が同型なら `β` も同型**(isotropic 型、`β` は `𝒪^▷(A)` の元)。

★sharp が「`n • x = 0 ⟹ x = 0`」として効く。 -/
theorem isIso_of_pow_isIso (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X) (β : End A)
    (hb : IsBaseIdentity P β) (hl : IsLinear P β)
    {n : ℕ} (hn : 0 < n) (hpow : IsIso ((β ^ n : End A) : A ⟶ A)) :
    IsIso (β : A ⟶ A) := by
  haveI := hpow
  have hdiv0 : P.Div ((β ^ n : End A) : A ⟶ A) = 0 :=
    (isIsometric_of_isIso P ((β ^ n : End A) : A ⟶ A))
  have hdβ : P.Div (β : A ⟶ A) = 0 :=
    eq_zero_of_nsmul_eq_zero_of_isSharp (P.divisorial _).2 hn
      (by rw [← (otri_pow_invariants P β hb hl n).2.2]; exact hdiv0)
  have hbi : IsIso (P.Base (β : A ⟶ A)) := by
    rw [show P.Base (β : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
    infer_instance
  exact prop_1_4_iii P F (β : A ⟶ A)
    ⟨prop_1_4_i P (β : A ⟶ A) (fun X _ => hiso X), hdβ⟩ ⟨hl, hbi⟩

/-! ### ★★選択肢 (c) を試した —— **公理の非対称性が穴の正体だった**

★**(c)「可換性を使わない議論」を探した。** 見つからなかったが、
★★**なぜ見つからないかが公理の形から言える。**

我々に必要なのは:
> `ζ n ≫ α = ζ n ≫ β^n`(`α, β^n ∈ 𝒪^▷(A)`)から `α = β^n` を出す

すなわち ★**「Frobenius 型射 `ζ n` に沿った `𝒪^▷` の対応が単射である」**こと。

★**`Definition 1.3` が与えるものを並べる**:

| 公理 | 沿う射 | `∃!` の位置 |
|---|---|---|
| `otriFwd` (iii)(c) | ★**co-angular pre-step** | 終域側(`β ∈ 𝒪^▷(B)` が一意) |
| `otriBwd` (iii)(c) | ★**co-angular pre-step** | 始域側(`α ∈ 𝒪^▷(A)` が一意) |
| `IsPerfectObj` (b) | ★**Frobenius 型射** | ★**始域側のみ**(`ψ` が一意) |

★★**pre-step に沿っては両側の `∃!` があり(＝全単射)、
Frobenius 型射に沿っては片側の `∃!` しかない。**

★**我々が要るのは終域側の一意性**(`ζ n ≫ α = ζ n ≫ β^n ⟹ α = β^n`)であり、
★**それは Frobenius 型射については公理に無い。**

★★**これが穴の正体である。** 「原文が書いていない一歩」でも
「原文が挙げていない依存」でもなく、★**公理の非対称性**である。

★**3 通りの結論**:
- (a) 可換性 —— 試した。**単元を法としてのみ成立**
- (b) 単元が中心 —— 試した。**一部(`β^n` 同型 ⟹ `β` 同型)のみ**
- (c) 可換性なしの議論 —— ★**試した。公理の非対称性により、この形では出ない**

★★**したがって、原文の「follows immediately」は
`Definition 1.3` に無いものを使っている**か、あるいは
★**我々が `Definition 1.3` の写しで何かを落としている**かのどちらかである。
★**後者の可能性を否定していない** —— それが次に確かめることである。
-/

/-! ### ★★★原典自身が `𝒪^×(C)` の可換性を**仮定**として置いている

★前段で「我々が `Definition 1.3` の写しで何かを落としているかもしれない」と書いた。
★**原典 p.22–23 を目視で読み直した。結果:**

**(1) 我々の写しは正しい。**
- `perfect object`: 原文 `ψ′ ◦ φ₁ = φ₂ ◦ ψ` ＝ Lean の `φ₁ ≫ ψ′ = ψ ≫ φ₂` ✓
- `Frobenius-normalized`: 原文 `α^d ◦ φ = φ ◦ α` ＝ Lean の `φ ≫ α^d = α ≫ φ` ✓
- `group-like object`: 原文 `Φ(C) = 0` ✓(我々の `IsGroupLike` はその同値形)

**(2) ★★しかし別の定義に決定的なものがあった。**

原文 (FrdI p.23):
> A Frobenius-compact object of C is defined to be an object C

原文 (FrdI p.23):
> such that O×(C) is commutative,

★**続きは引用できない**(事故 #4 の再来。同じ行の `̸=` が layout モードで消える)。
400 dpi 目視では「`O×(C)pf ≠ 0`、および `Aut_𝒞(C)` の元で `𝒪^×(C)^pf` に
`ℚ>0` の元の乗法として作用するものは実際には自明に作用する」と続く。

★★★**原典は `Frobenius-compact` の定義の中で「`𝒪^×(C)` が可換であること」を
明示的な仮定として課している。**
★**すなわち `𝒪^×(C)` は一般には可換ではない** —— 可換なら仮定に書く必要がない。

## ★★これが我々の穴の位置を確定させる

(iii) の moreover を我々は
> `α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)から `α = γ^n` を出したい
まで縮めた。これには `u = v^n` と `β`・`v` の可換性が要る。

★★**原典自身が `𝒪^×` の可換性を別の場所で仮定として置いている以上、
それは `Definition 1.3` からは出ない。**

★**したがって 2 つのうちどちらかである**:
- **(α)** (iii) の moreover は、原文が (iii) で挙げていない仮定
  (`𝒪^×(A)` の可換性、あるいは Frobenius-compact 性)を使っている
- **(β)** 可換性を回避する議論があり、我々がまだ見つけていない

★**(β) を否定していない。** ただし ★**(α) の側に、原典自身が置いた証拠が 1 つある。**

★★**これは「原文が挙げる根拠が実際の根拠の一部」(4 種類目)の、
最も強い実例である** —— **原文は同じ論文の別の場所で、必要な仮定を明示している。**
-/

/-! ### ★★★選択肢 (β) が出た —— **探していた道具は最初から手元にあった**

★我々は「`ζ n` で**左から**消せること(mono)」を探し、`Definition 1.3` が
pre-step にしか mono を与えないことで詰まっていた。★**必要だったのは 2 つとも別物だった。**

★★**全射性**: `α ∈ 𝒪^▷(A)` に対し `IsPerfectObj` の (b) が pre-step `β` を与え
`ζ n ≫ α = β ≫ ζ n`。**Frobenius-normalized** から `ζ n ≫ β^n = β ≫ ζ n = ζ n ≫ α`。
★★**`𝒞` は totally epimorphic なので `ζ n` は epi** —— `Epi f` は
**`f ≫ x = f ≫ y ⟹ x = y`**、すなわち**左に置かれた `f` を消す**。
★**これがちょうど要る向きだった** ⟹ **`β^n = α`**。

★★**単射性**: `β₁^n = β₂^n = α` なら、Frobenius-normalized から
`ζ n ≫ α = βᵢ ≫ ζ n`(`i = 1, 2`)。
★★**(b) の `∃!` が「`ψ′ = α` に対する `ψ` は一意」と言っている** ⟹ **`β₁ = β₂`**。

★★★**我々は `faithfulUpToUnits` を経由して単元の吸収に苦しんでいたが、
`totEpiC` と (b) の `∃!` だけで足りた。**
★**「読み落とし」(6 種類目のギャップ)の 3 例目である。**
★★**3 回とも、必要なものは既に手元にあった。**
(1: §0 の torsion-free、2: `Gp M` の群構造の import、3: これ)

★**そして `𝒪^×(A)` の可換性は要らなかった** ——
原典が `Frobenius-compact` で可換性を仮定しているのは**別の目的**である。
★**我々が「穴」と呼んでいたものは、遠回りの経路にしか存在しなかった。**
-/

include P in
/-- ★★★**`Proposition 1.10, (iii)` の moreover(`𝒪^▷(A)` の側)完成** ——
`𝒪^▷(A)` では `n` 乗根が**一意に存在する**。

★**全射性**は `𝒞` の totally epimorphic 性(epi)から、
★**単射性**は `IsPerfectObj` の (b) の `∃!` から。
★**`𝒪^×(A)` の可換性は要らない。** -/
theorem prop_1_10_iii_otri_perfect {A : C}
    (hperf : IsPerfectObj P A) (hfn : IsFrobeniusNormalized P A)
    (n : ℕ+) (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (hζd : P.degFr ζn = n) (α : End A) (hα : α ∈ OTri P A) :
    ∃! β : End A, β ∈ OTri P A ∧ (β ^ (n : ℕ) : End A) = α := by
  haveI : Epi (ζn : A ⟶ A) := P.totEpiC _ _ _
  have hαs : IsPreStep P (α : A ⟶ A) := by
    refine ⟨hα.2, ?_⟩
    show IsIso (P.Base (α : A ⟶ A))
    rw [show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1, P.Base_id]
    infer_instance
  have hAA : BaseIsomorphic P A A := ⟨Iso.refl _⟩
  obtain ⟨β, ⟨hβs, hsq⟩, huniq⟩ :=
    (hperf n).2 A A A A ζn ζn hζf hζd hζf hζd hAA hAA (α : End A) hαs
  -- `β` は base-identity
  have hβb : IsBaseIdentity P β := by
    show P.Base β = P.Base (𝟙 A)
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1,
      P.Base_id, Category.comp_id, show P.Base ζn = P.Base (𝟙 A) from hζb,
      P.Base_id, Category.comp_id] at h
    rw [P.Base_id]
    simpa using h.symm
  have hβm : β ∈ OTri P A := ⟨hβb, hβs.1⟩
  -- Frobenius-normalized で `ζ n ≫ β^n = β ≫ ζ n`
  have hfnb := hfn ζn hζb β hβm
  rw [hζd] at hfnb
  refine ⟨β, ⟨hβm, ?_⟩, ?_⟩
  · -- 全射性: `ζ n` が epi
    refine (cancel_epi (ζn : A ⟶ A)).mp ?_
    rw [hfnb, hsq]
  · -- 単射性: (b) の `∃!`
    rintro γ ⟨hγm, hγp⟩
    refine huniq γ ⟨⟨hγm.2, ?_⟩, ?_⟩
    · show IsIso (P.Base (γ : A ⟶ A))
      rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hγm.1, P.Base_id]
      infer_instance
    · have hfng := hfn ζn hζb γ hγm
      rw [hζd, hγp] at hfng
      exact hfng

include P in
/-- ★★★**`Proposition 1.10, (iii)` の moreover 完成** —— `𝒪^×(A)` も perfect。

★`𝒪^▷(A)` の結果に、★**「`β^n` が同型なら `β` も同型」**(`isIso_of_pow_isIso`)を
足すだけ。★**単元の側は `𝒪^▷` の側から自動で従う。** -/
theorem prop_1_10_iii_otimes_perfect (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X)
    (hperf : IsPerfectObj P A) (hfn : IsFrobeniusNormalized P A)
    (n : ℕ+) (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (hζd : P.degFr ζn = n) (α : End A) (hα : α ∈ OTimes P A) :
    ∃! β : End A, β ∈ OTimes P A ∧ (β ^ (n : ℕ) : End A) = α := by
  obtain ⟨β, ⟨hβm, hβp⟩, huniq⟩ :=
    prop_1_10_iii_otri_perfect P hperf hfn n ζn hζb hζf hζd α hα.1
  -- `β^n = α` は単元なので `β` も単元
  haveI hpow : IsIso ((β ^ (n : ℕ) : End A) : A ⟶ A) := by
    rw [hβp]
    exact (CategoryTheory.isUnit_iff_isIso (α : End A)).mp hα.2
  have hβiso : IsIso ((β : End A) : A ⟶ A) :=
    isIso_of_pow_isIso P F hiso β hβm.1 hβm.2 n.pos hpow
  refine ⟨β, ⟨⟨hβm, (CategoryTheory.isUnit_iff_isIso (β : End A)).mpr hβiso⟩, hβp⟩, ?_⟩
  rintro γ ⟨hγm, hγp⟩
  exact huniq γ ⟨hγm.1, hγp⟩

/-! ### ★★(iii) の還元 —— 原文「we may assume that A is Frobenius-trivial」

原文 (FrdI p.35):
> of Frobenius-trivial objects [cf. Definition 1.3, (i), (a), (b); the isomorphism of Def-

原文 (FrdI p.35):
> inition 1.3, (iii), (c)], we may assume that A is Frobenius-trivial. Now the fact that

★★**取り下げの原因の1つを埋める**(2026-08-16)。
`prop_1_10_iii_otri_perfect` は「次数 `n` の base-identity Frobenius 型自己射 `ζ n` が
ある」を**仮定に置いていた** —— それは `A` が Frobenius-trivial だということである。
★原文はそこを**還元して**いる。その還元をここで実装する。

★**原文が名指す 3 つが、そのまま 3 段になる**:
1. `Definition 1.3, (i), (a)`(`baseSurj`) —— `Base(A)` を底に持つ
   **Frobenius-trivial な `A₀`** がある
2. `Definition 1.3, (i), (b)`(`preStepSpan`) —— その底の同型は
   **pre-step の span** `A₀ ← X → A` に持ち上がる
3. `Definition 1.3, (iii), (c)`(`otriFwd` / `otriBwd`) —— **co-angular** pre-step は
   `𝒪^▷` の全単射を与える

★★**`isotropic` の仮定がここで効く**。原文の角括弧
「so all morphisms of C are co-angular — cf. Proposition 1.4, (i)」がそれである ——
span の pre-step が co-angular だと分かって初めて 3 が使える。
★**原文は「isotropic」を moreover の仮定として置いているが、
その最初の用途はこの還元である**(perfect 性の議論そのものではない)。

★★**「一意存在が移る」ことだけが要る**。全単射そのもの(`MulEquiv`)は要らない ——
`otriFwd` / `otriBwd` の 2 つの `∃!` と、四角形が積と可換すること
(`φ ≫ (β₂ ≫ β₁) = (α₂ ≫ α₁) ≫ φ`)があれば、`n` 乗根の一意存在は渡る。
-/

section PowTransfer

variable {M N : Type*} [Monoid M] [Monoid N]

/-- ★**全単射の graph に沿って「`n` 乗根の一意存在」が移る**。

★`R` は `Definition 1.3, (iii), (c)` の四角形「`φ ≫ β = α ≫ φ`」であり、
`otriFwd` / `otriBwd` がちょうど**両向きの一意性**を与える。
★積と可換する(`hmul`)ことから `n` 乗とも可換し、両側の一意性で挟める。

★**`MulEquiv` を作らないのは、`R` が `∃!` の形でしか与えられないからである** ——
写像を取り出すには選択が要るが、主張は選択なしで書ける。 -/
theorem exists_unique_pow_transfer (R : M → N → Prop)
    (hfwd : ∀ a : M, ∃! b : N, R a b) (hbwd : ∀ b : N, ∃! a : M, R a b)
    (hone : R 1 1)
    (hmul : ∀ (a₁ a₂ : M) (b₁ b₂ : N), R a₁ b₁ → R a₂ b₂ → R (a₁ * a₂) (b₁ * b₂))
    (n : ℕ) (h : ∀ a : M, ∃! x : M, x ^ n = a) (b : N) : ∃! y : N, y ^ n = b := by
  have hpow : ∀ (a : M) (b : N), R a b → ∀ k : ℕ, R (a ^ k) (b ^ k) := by
    intro a b hab k
    induction k with
    | zero => rw [pow_zero, pow_zero]; exact hone
    | succ k ih => rw [pow_succ, pow_succ]; exact hmul _ _ _ _ ih hab
  obtain ⟨a, hab, hau⟩ := hbwd b
  obtain ⟨x, hxa, hxu⟩ := h a
  obtain ⟨y, hxy, hyu⟩ := hfwd x
  obtain ⟨z, -, hzu⟩ := hfwd (x ^ n)
  refine ⟨y, ?_, ?_⟩
  · show y ^ n = b
    rw [hzu _ (hpow x y hxy n), hzu b (by rw [hxa]; exact hab)]
  · intro y' hy'
    have hy'' : y' ^ n = b := hy'
    obtain ⟨a', ha'y, -⟩ := hbwd y'
    have ha'n : a' ^ n = a := hau _ (by rw [← hy'']; exact hpow a' y' ha'y n)
    exact hyu y' (by rw [← hxu a' ha'n]; exact ha'y)

end PowTransfer

/-- ★`𝒪^▷` の「`n` 乗根の一意存在」の 2 つの綴りが同値であること ——
`End A` の元＋所属で書いた形と、部分モノイドの元で書いた形。

★**前者は原文の言い回しに近く、後者は移送の補題に渡しやすい。** -/
theorem otri_existsUnique_pow_iff {A : C} (n : ℕ) :
    (∀ α ∈ OTri P A, ∃! β : End A, β ∈ OTri P A ∧ (β ^ n : End A) = α)
      ↔ (∀ a : OTri P A, ∃! x : OTri P A, x ^ n = a) := by
  constructor
  · intro h a
    obtain ⟨β, ⟨hβm, hβp⟩, hβu⟩ := h (a : End A) a.2
    refine ⟨⟨β, hβm⟩, Subtype.ext (by simpa using hβp), fun y hy => ?_⟩
    exact Subtype.ext (hβu (y : End A) ⟨y.2, by simpa using congrArg Subtype.val hy⟩)
  · intro h α hα
    obtain ⟨x, hx, hxu⟩ := h ⟨α, hα⟩
    refine ⟨(x : End A), ⟨x.2, by simpa using congrArg Subtype.val hx⟩, ?_⟩
    rintro β ⟨hβm, hβp⟩
    exact congrArg Subtype.val (hxu ⟨β, hβm⟩ (Subtype.ext (by simpa using hβp)))

/-- ★★**`Definition 1.3, (iii), (c)` を「perfect 性の移送」として使う**。

co-angular pre-step `φ : A ⟶ B` に沿って、`𝒪^▷` の `n` 乗根の一意存在は
★**両向きに**移る(`otriFwd` と `otriBwd` の両方があるため)。 -/
theorem otri_pow_transfer (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) (n : ℕ) :
    (∀ a : OTri P A, ∃! x : OTri P A, x ^ n = a) ↔
      (∀ b : OTri P B, ∃! y : OTri P B, y ^ n = b) := by
  set R : OTri P A → OTri P B → Prop :=
    fun a b => (φ ≫ (b : End B) : A ⟶ B) = (a : End A) ≫ φ with hR
  have hfwd : ∀ a : OTri P A, ∃! b : OTri P B, R a b := by
    intro a
    obtain ⟨b, ⟨hbm, hbr⟩, hbu⟩ := F.otriFwd φ hc hs (a : End A) a.2
    exact ⟨⟨b, hbm⟩, hbr, fun y hy => Subtype.ext (hbu (y : End B) ⟨y.2, hy⟩)⟩
  have hbwd : ∀ b : OTri P B, ∃! a : OTri P A, R a b := by
    intro b
    obtain ⟨a, ⟨ham, har⟩, hau⟩ := F.otriBwd φ hc hs (b : End B) b.2
    exact ⟨⟨a, ham⟩, har, fun y hy => Subtype.ext (hau (y : End A) ⟨y.2, hy⟩)⟩
  have hone : R 1 1 := by
    show (φ ≫ (𝟙 B : End B) : A ⟶ B) = (𝟙 A : End A) ≫ φ
    simp
  have hmul : ∀ (a₁ a₂ : OTri P A) (b₁ b₂ : OTri P B),
      R a₁ b₁ → R a₂ b₂ → R (a₁ * a₂) (b₁ * b₂) := by
    intro a₁ a₂ b₁ b₂ h₁ h₂
    show (φ ≫ ((b₂ : End B) ≫ (b₁ : End B)) : A ⟶ B)
      = ((a₂ : End A) ≫ (a₁ : End A)) ≫ φ
    rw [← Category.assoc φ, show (φ ≫ (b₂ : End B) : A ⟶ B) = (a₂ : End A) ≫ φ from h₂,
      Category.assoc, show (φ ≫ (b₁ : End B) : A ⟶ B) = (a₁ : End A) ≫ φ from h₁,
      ← Category.assoc]
  exact ⟨fun h => exists_unique_pow_transfer R hfwd hbwd hone hmul n h,
    fun h => exists_unique_pow_transfer (fun b a => R a b) hbwd hfwd hone
      (fun b₁ b₂ a₁ a₂ h₁ h₂ => hmul a₁ a₂ b₁ b₂ h₁ h₂) n h⟩

/-- ★★★**`Proposition 1.10, (iii)` の moreover —— 還元まで込めた完成形**(`𝒪^▷`)。

★原文の仮定そのもの(「`𝒞` が perfect・isotropic・Frobenius-normalized 型」)から、
★**任意の対象 `A`** について `𝒪^▷(A)` で `n` 乗根が一意に存在する。
`ζ n` の存在はもはや仮定ではなく、`Definition 1.3, (i), (a)` から**導かれる**。 -/
theorem prop_1_10_iii_otri_perfect_of_type (F : FrobenioidCore P)
    (hpt : IsOfPerfectType P) (hiso : ∀ X : C, IsIsotropic P X)
    (hfnt : ∀ X : C, IsFrobeniusNormalized P X) (A : C) (n : ℕ+)
    (α : End A) (hα : α ∈ OTri P A) :
    ∃! β : End A, β ∈ OTri P A ∧ (β ^ (n : ℕ) : End A) = α := by
  -- 1 段: `Definition 1.3, (i), (a)` —— 同じ底を持つ Frobenius-trivial な `A₀`
  obtain ⟨A₀, ⟨ζ, hζd, hζbf⟩, ⟨e⟩⟩ := F.baseSurj (P.toElem.obj A).base
  have hA₀ : ∀ a : OTri P A₀, ∃! x : OTri P A₀, x ^ (n : ℕ) = a :=
    (otri_existsUnique_pow_iff P (n : ℕ)).mp (fun β hβ =>
      prop_1_10_iii_otri_perfect P (hpt A₀) (hfnt A₀) n (ζ n)
        (hζbf n).1 (hζbf n).2 (hζd n) β hβ)
  -- 2 段: `Definition 1.3, (i), (b)` —— pre-step の span `A₀ ← X → A`
  obtain ⟨X, φ, ψ, hφs, hψs, -⟩ := F.preStepSpan A₀ A e.hom inferInstance
  -- 3 段: isotropic ⟹ co-angular(`Proposition 1.4, (i)`)、`(iii)(c)` で両向きに移送
  have hφc : IsCoAngular P φ := prop_1_4_i P φ (fun Y _ => hiso Y)
  have hψc : IsCoAngular P ψ := prop_1_4_i P ψ (fun Y _ => hiso Y)
  exact (otri_existsUnique_pow_iff P (n : ℕ)).mpr
    ((otri_pow_transfer P F ψ hψc hψs (n : ℕ)).mp
      ((otri_pow_transfer P F φ hφc hφs (n : ℕ)).mpr hA₀)) α hα

/-- ★★★**`Proposition 1.10, (iii)` の moreover —— 還元まで込めた完成形**(`𝒪^×`)。

★`𝒪^▷` の側に「`β^n` が同型なら `β` も同型」を足すだけ。 -/
theorem prop_1_10_iii_otimes_perfect_of_type (F : FrobenioidCore P)
    (hpt : IsOfPerfectType P) (hiso : ∀ X : C, IsIsotropic P X)
    (hfnt : ∀ X : C, IsFrobeniusNormalized P X) (A : C) (n : ℕ+)
    (α : End A) (hα : α ∈ OTimes P A) :
    ∃! β : End A, β ∈ OTimes P A ∧ (β ^ (n : ℕ) : End A) = α := by
  obtain ⟨β, ⟨hβm, hβp⟩, huniq⟩ :=
    prop_1_10_iii_otri_perfect_of_type P F hpt hiso hfnt A n α hα.1
  haveI hpow : IsIso ((β ^ (n : ℕ) : End A) : A ⟶ A) := by
    rw [hβp]
    exact (CategoryTheory.isUnit_iff_isIso (α : End A)).mp hα.2
  have hβiso : IsIso ((β : End A) : A ⟶ A) :=
    isIso_of_pow_isIso P F hiso β hβm.1 hβm.2 n.pos hpow
  refine ⟨β, ⟨⟨hβm, (CategoryTheory.isUnit_iff_isIso (β : End A)).mpr hβiso⟩, hβp⟩, ?_⟩
  rintro γ ⟨hγm, hγp⟩
  exact huniq γ ⟨hγm.1, hγp⟩

/-! ### ★(iv) の「In particular」—— 無限個の同型類

原文 (FrdI p.34):
> that arise from irreducible arrows with domain A.

★**構成**: 各素数 `p` に対し `frobDegSurj` が次数 `p` の Frobenius 型射 `φ_p : A ⟶ B_p` を与え、
(iv) によりそれは **irreducible**。
★**異なる素数は非同型な対象を与える**: `_A𝒞` の同型 `θ : B_p ⟶ B_q`(`φ_p ≫ θ = φ_q`)が
あれば `degFr φ_q = degFr θ * degFr φ_p = 1 * p = p`。★**`p ≠ q` に矛盾。**

★**素数が無限にあることは `Nat.exists_infinite_primes`(mathlib)から。**
★**「無限個の同型類」は「素数から同型類への単射がある」として書く** ——
原文の "infinitely many isomorphism classes" の忠実な形である。
-/

include P in
/-- ★**異なる次数の Frobenius 型射は `_A𝒞` で非同型**。

★これが「無限個の同型類」の核心。次数が同型類の不変量になっている。 -/
theorem frobType_not_iso_of_degFr_ne {A Bp Bq : C}
    (φp : A ⟶ Bp) (φq : A ⟶ Bq) (hne : P.degFr φp ≠ P.degFr φq)
    (θ : Bp ⟶ Bq) [IsIso θ] : φp ≫ θ ≠ φq := by
  intro h
  apply hne
  have : P.degFr (φp ≫ θ) = P.degFr φq := by rw [h]
  rwa [P.degFr_comp, degFr_of_isIso P θ, one_mul] at this

include P in
/-- ★★**(iv) の「In particular」** —— isotropic な `A` について、
**各素数ごとに irreducible な射があり、異なる素数のものは非同型**。

★**これが「無限個の同型類」の内容である**(素数は無限にあるので)。 -/
theorem prop_1_10_iv_infinitely_many (F : FrobenioidCore P)
    {A : C} (hA : IsIsotropic P A) :
    (∀ p : ℕ+, Nat.Prime (p : ℕ) →
       ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ P.degFr φ = p) ∧
    (∀ (Bp Bq : C) (φp : A ⟶ Bp) (φq : A ⟶ Bq), P.degFr φp ≠ P.degFr φq →
       ∀ θ : Bp ≅ Bq, φp ≫ θ.hom ≠ φq) := by
  constructor
  · intro p hp
    obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A p
    refine ⟨B, φ, ?_, hd⟩
    exact prop_1_10_iv_mp P F hA φ ⟨hφ, by rw [hd]; exact hp⟩
  · intro Bp Bq φp φq hne θ
    haveI : IsIso θ.hom := inferInstance
    exact frobType_not_iso_of_degFr_ne P φp φq hne θ.hom

include P in
/-- ★★**原文「there exist infinitely many isomorphism classes」そのもの**。

原文 (FrdI p.34):
> is isotropic, then there exist infinitely many isomorphism classes of objects of AC

★**監査以前は「無限個」を述べていなかった**(素数ごとの存在と非同型性まで)。
★**次数がいくらでも大きく取れる**という形で述べる ——
`prop_1_10_iv_infinitely_many` の第2主張(次数が違えば非同型)と合わせて、
★同型類が無限個あることを与える。 -/
theorem prop_1_10_iv_unbounded (F : FrobenioidCore P) {A : C} (hA : IsIsotropic P A) :
    ∀ n : ℕ, ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ n ≤ (P.degFr φ : ℕ) := by
  intro n
  obtain ⟨p, hpn, hp⟩ := Nat.exists_infinite_primes n
  have hppos : 0 < p := hp.pos
  obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A ⟨p, hppos⟩
  refine ⟨B, φ, prop_1_10_iv_mp P F hA φ ⟨hφ, ?_⟩, ?_⟩
  · rw [hd]; exact hp
  · rw [hd]; exact hpn

include P in
/-- ★★★**原文の「infinitely many isomorphism classes」を集合の無限性で述べる**。

原文 (FrdI p.34):
> is isotropic, then there exist infinitely many isomorphism classes of objects of AC

★★**監査で「「無限個」そのものを述べていない」と指摘されたもの**。

★**実現される次数の集合が無限**であることを述べる。
★`prop_1_10_iv_infinitely_many` の第2主張（次数が違えば非同型）と合わせれば、
★★**同型類が無限個ある**ことになる（次数が全単射の代わりをする）。 -/
theorem prop_1_10_iv_degrees_infinite (F : FrobenioidCore P) {A : C}
    (hA : IsIsotropic P A) :
    {n : ℕ | ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ ((P.degFr φ : ℕ+) : ℕ) = n}.Infinite := by
  refine Set.Infinite.mono ?_ Nat.infinite_setOf_prime
  intro p hp
  obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A ⟨p, hp.pos⟩
  exact ⟨B, φ, prop_1_10_iv_mp P F hA φ ⟨hφ, by rw [hd]; exact hp⟩, by rw [hd]; rfl⟩


/-! ### ★(vi) 後半の最終段 —— `IsFrobeniusTrivial` を同型に沿って移す

★**原文(p.36)の最後の 1 文**:
> But by Proposition 1.4, (iii), these pre-steps are isomorphisms, so A is Frobenius-trivial.

★**「so A is Frobenius-trivial」の中身**は「同型に沿った移送」である。
★**mathlib に `End` の共役の準同型が無い**ので自作する
(`grep` で確認: `Mathlib/CategoryTheory/Endomorphism.lean` に `conj` は無い)。

★**`End` の乗法が `x * y = y ≫ x`** なので、共役 `f ↦ θ⁻¹ ≫ f ≫ θ` が
準同型であることの計算も**その向きで**行う。
-/

/-- ★**同型による `End` の共役**(`End A →* End B`)。

★mathlib に無いので置いた。★`End` の乗法は `x * y = y ≫ x` である。 -/
@[simps] def endConj {A B : C} (θ : A ≅ B) : End A →* End B where
  toFun f := θ.inv ≫ f ≫ θ.hom
  map_one' := by simp
  map_mul' x y := by
    simp only [End.mul_def]
    simp

include P in
/-- ★★**共役は `𝒪^×` を保つ**（2026-08-16 追加）。

★`Definition 1.2, (iv)` の `Frobenius-compact` の第3節で
`Aut_𝒞(C)` の `𝒪^×(C)` への作用を書くのに要る。

★**3 成分とも `Base` / `degFr` の関手性と同型の linear 性から出る**。 -/
theorem endConj_mem_otimes {A B : C} (θ : A ≅ B) {f : End A} (hf : f ∈ OTimes P A) :
    endConj θ f ∈ OTimes P B := by
  refine ⟨⟨?_, ?_⟩, IsUnit.map (endConj θ) hf.2⟩
  · show P.Base (θ.inv ≫ (f : A ⟶ A) ≫ θ.hom) = P.Base (𝟙 B)
    rw [P.Base_comp, P.Base_comp,
      show P.Base (f : A ⟶ A) = P.Base (𝟙 A) from hf.1.1,
      P.Base_id, Category.id_comp, ← P.Base_comp, θ.inv_hom_id]
  · show P.degFr (θ.inv ≫ (f : A ⟶ A) ≫ θ.hom) = 1
    rw [P.degFr_comp, P.degFr_comp,
      show P.degFr (f : A ⟶ A) = 1 from hf.1.2,
      show P.degFr θ.hom = 1 from isLinear_of_isIso P θ.hom,
      show P.degFr θ.inv = 1 from isLinear_of_isIso P θ.inv]
    simp

/-! ### ★★`Frobenius-compact`（2026-08-16 追加）

原文 (FrdI p.23):
> C such that O×(C) is commutative, O×(C)pf

★**引用を切った記録（事故 #3 の 7 度目）**: 続く `̸= 0, and every element of AutC(C)` は
★**`̸=`（打ち消し付きの等号）が抽出で落ちる**ため引用できない（35/63 文字で停止）。
★`′`・`▷`・`≼` に続く 4 つ目の「目では見えるが抽出に無い文字」である。

★★**`Pf` を経由せずに書く。** 理由は 2 つある。

1. `𝒪^×(C)` は **乗法**の `Submonoid (End A)` であり、`Pf` は
   `[AddCommMonoid M]` を要求する。★**しかも可換性は①として
   仮定される側**であって既知ではないので、型を作るのに仮定が要る。
2. ★★**同値な式がある**。

★**②の翻訳**（`Pf.mk_eq_zero_iff` が根拠）:
`M^pf = 0` ⇔ 「すべての `m` に `∃ k, k • m = 0`」。
乗法で書き直すと `M^pf ≠ 0` ⇔ ★**「ある `u ∈ 𝒪^×(C)` が無限位数を持つ」**。

★**③の翻訳**:
「`q = c/d` 倍として作用する」は `ℚ` を使わずに
`∀ x, d • σ x = c • x` と書ける。`M → M^pf` の像の上で書き直すと
`∃ k, (σ u)^(d k) = u^(c k)`。「acts trivially」は `∃ k, (σ u)^k = u^k`。

★★**像への制限が情報を落とさない理由**:
任意の `x = mk m a` は `a • x = of m` を満たす。誘導写像は加法的なので
`a • (d • σ x) = d • σ (a • x) = d • σ (of m) = c • (of m) = a • (c • x)`。
★**`Pf` では `a •` が単射**（`Pf.nsmul_injective`）なので、
像の上の式から `M^pf` 全体の式が出る。

★一度ここを「`Pf.eq_of_qact` が根拠」と書いたが、それは「同じ `c/d` 倍として
作用する写像は一意」という**別の主張**であり、ここでは効かない。
★★**検証役の再監査で訂正した** —— 根拠として挙げた補題と、実際に効く補題が違っていた。

★★**原文の三節をそのまま三連言にする。** -/

include P in
/-- **[FrdI] Definition 1.2, (iv)** `Frobenius-compact object`。

★上の docstring の翻訳に従う。★**共役が `𝒪^×` を保つ**こと
（`endConj_mem_otimes`）が③を書くのに要る。 -/
def IsFrobeniusCompact (A : C) : Prop :=
  (∀ x y : End A, x ∈ OTimes P A → y ∈ OTimes P A → x * y = y * x) ∧
  (∃ u : End A, u ∈ OTimes P A ∧ ∀ k : ℕ+, (u ^ ((k : ℕ+) : ℕ) : End A) ≠ 1) ∧
  (∀ (θ : A ≅ A) (c d : ℕ+),
    (∀ u : End A, u ∈ OTimes P A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) →
    ∀ u : End A, u ∈ OTimes P A → ∃ k : ℕ+,
      ((endConj θ u) ^ ((k : ℕ+) : ℕ) : End A) = (u ^ ((k : ℕ+) : ℕ) : End A))

include P in
/-- **[FrdI] Definition 1.2, (v)** `of Frobenius-compact type`。 -/
def IsOfFrobeniusCompactType : Prop := ∀ A : C, IsFrobeniusCompact P A

def IsFrobeniusCompact.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact object",
    sectionId := "frdi-def-1-2-iv" }

include P in
/-- ★★**`IsFrobeniusTrivial` は同型に沿って移る**。

★原文の「so A is Frobenius-trivial」の中身。 -/
theorem isFrobeniusTrivial_of_iso (F : FrobenioidCore P) {A B : C} (θ : A ≅ B)
    (h : IsFrobeniusTrivial P A) : IsFrobeniusTrivial P B := by
  obtain ⟨ζ, hd, hp⟩ := h
  refine ⟨(endConj θ).comp ζ, ?_, ?_⟩
  · intro n
    show P.degFr (θ.inv ≫ (ζ n : A ⟶ A) ≫ θ.hom) = n
    rw [P.degFr_comp, P.degFr_comp, degFr_of_isIso P θ.hom, degFr_of_isIso P θ.inv,
      hd n, one_mul, mul_one]
  · intro n
    constructor
    · show P.Base (θ.inv ≫ (ζ n : A ⟶ A) ≫ θ.hom) = P.Base (𝟙 B)
      rw [P.Base_comp, P.Base_comp, show P.Base (ζ n : A ⟶ A) = P.Base (𝟙 A) from (hp n).1,
        P.Base_id, Category.id_comp, ← P.Base_comp, θ.inv_hom_id, P.Base_id]
    · exact IsFrobeniusType.comp P F (isFrobeniusType_of_isIso P θ.inv)
        (IsFrobeniusType.comp P F (hp n).2 (isFrobeniusType_of_isIso P θ.hom))

include P in
/-- ★**group-like は base-isomorphism に沿って移る**。

★`Φ.map (Base α)` は `Base α` が同型なので全単射(`Φ.map (inv (Base α))` が逆)。
`Φ(A)` の元がすべて 0 なら `Φ(A′)` の元もすべて 0。 -/
theorem isGroupLikeObj_of_baseIso {A A' : C} (α : A' ⟶ A) (hbi : IsBaseIsomorphism P α)
    (hA : IsGroupLikeObj P A) : IsGroupLikeObj P A' := by
  haveI : IsIso (P.Base α) := hbi
  show IsGroupLike (Φ.val (P.toElem.obj A').base)
  rw [isGroupLike_iff]
  intro y
  have hy : y = Φ.map (P.Base α) (Φ.map (inv (P.Base α)) y) := by
    rw [← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
  have h0 : Φ.map (inv (P.Base α)) y = 0 :=
    eq_zero_of_isGroupLike_of_isSharp hA (P.divisorial _).2 _
  rw [hy, h0, map_zero]
  exact isAddUnit_zero

include P in
/-- ★★★**`Proposition 1.10, (vi)` 後半 完成** ——
isotropic 型で group-like な対象は Frobenius-trivial。

★原文(p.36)の最後の 3 文:
「`Definition 1.3, (i), (a), (b)` から co-angular pre-step `A′ → A`、`A′ → A″`
(`A″` は Frobenius-trivial)がある。`Proposition 1.4, (iii)` によりこれらは同型。
よって `A` は Frobenius-trivial。」

★**`A′` が group-like であること**(`isGroupLikeObj_of_baseIso`)が
★**原文が書いていない一歩**である —— pre-step が同型であることに要る isometry は
`A′` の group-like 性から来るのであって、`A` のそれからではない。 -/
theorem prop_1_10_vi_groupLike (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A A' A'' : C} (hA : IsGroupLikeObj P A)
    (α : A' ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (γ : A' ⟶ A'') (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (hA'' : IsFrobeniusTrivial P A'') : IsFrobeniusTrivial P A := by
  -- `A′` も group-like
  have hA' : IsGroupLikeObj P A' := isGroupLikeObj_of_baseIso P α hαs.2 hA
  -- `α`・`γ` は同型
  haveI hαi : IsIso α := isIso_of_preStep_of_isGroupLikeObj' P F hiso hA' α hαs
  haveI hγi : IsIso γ := isIso_of_preStep_of_isGroupLikeObj' P F hiso hA' γ hγs
  -- `A″ ≅ A′ ≅ A`
  exact isFrobeniusTrivial_of_iso P F ((asIso γ).symm ≪≫ asIso α) hA''

/-! ### ★★(vi) の主語を置く —— `𝒞^istr` への具体化(2026-08-16)

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

原文 (FrdI p.35):
> every group-like object A ∈Ob(Cistr) is Frobenius-trivial.

★★**取り下げ表の (vi) の残り 2 件**。`prop_1_10_vi_ofType` /
`prop_1_10_vi_groupLike` はどちらも「isotropic 型の Frobenioid」について述べており、
★**原文の主語 `𝒞^istr` がコメントの外に一度も現れていなかった。**

★**具体化そのものは 1 手である** —— `Proposition 1.9` が
`istr_frobenioid`(`𝒞^istr` は Frobenioid)と `istr_isotropic`(全対象が isotropic)を
既に与えているので、それを当てるだけでよい。
★★**「主語が不在」という欠落は、材料が揃っていても消えない** ——
一般形を証明したことと、原文が名指した対象について述べたことは別である。
-/

/-- ★★★**`Proposition 1.10, (vi)` 前半 —— 原文の主語で**。
`𝒞^istr` は sub-quasi-Frobenius-trivial 型である。 -/
theorem prop_1_10_vi_istr (F : FrobenioidCore P) (G : Frobenioid P) :
    IsOfSubQuasiFrobeniusTrivialType (istrPre P F) :=
  prop_1_10_vi_ofType (istrPre P F) (istr_frobenioid P F G) (istr_isotropic P F)

/-- ★★★**`Proposition 1.10, (vi)` 後半 —— 原文の主語で**。
`𝒞^istr` の group-like な対象は Frobenius-trivial である。

★`α`・`γ`・`Cc` は `prop_1_10_vi_data` が `Definition 1.3, (i), (a), (b)` から
導くので、ここでも仮定に置かない。 -/
theorem prop_1_10_vi_istr_groupLike (F : FrobenioidCore P) (G : Frobenioid P)
    (A : Istr P) (hA : IsGroupLikeObj (istrPre P F) A) :
    IsFrobeniusTrivial (istrPre P F) A := by
  obtain ⟨_, _, α, γ, hαc, hαs, hγc, hγs, hCc⟩ :=
    prop_1_10_vi_data (istrPre P F) (istr_frobenioid P F G) (istr_isotropic P F) A
  exact prop_1_10_vi_groupLike (istrPre P F) (istr_frobenioidCore P F)
    (istr_isotropic P F) hA α hαc hαs γ hγc hγs hCc

/-! ## ★★★出典の紐付け(`.src`) —— `Proposition 1.10` は **21 主張すべて完成**

★**`.src` は「その原典項目を完全に実装した」という主張である**(2026-08-15 に明文化した規則)。
`Proposition 1.10` は 6 条 21 主張:

| 条 | 主張 | 実装 |
|---|---|---|
| (i) | 10 | 一意性 / 存在(3 場合＋任意) / `degFr` / `Div` / 7 タイプ |
| (ii) | 3 | 組み替え / `degFr` / `Div` |
| (iii) | 3 | 像のモノイド / `𝒪^▷(A)` / `𝒪^×(A)` |
| (iv) | 2 | iff / 無限個の同型類 |
| (v) | 1 | iff |
| (vi) | 2 | sub-quasi-Frobenius-trivial / group-like ⟹ Frobenius-trivial |

★★**この数え方が誤っていた**(下の取り下げ記録を見よ)。
-/

/-! ## ★★`.src` を**取り下げた**記録(2026-08-16)

★私の文脈を持たない検証役(Challenger)に監査を依頼し、
★★**6 個すべてが規則(「`.src` は完全実装の主張」)を満たしていない**ことが判明した。
★指摘はすべて**原文 (p.34–35) と Lean を私自身が直接読んで確認済み**である。

| 条 | 欠けているもの | 該当宣言 |
|---|---|---|
| (i) | ~~★`β` の量化子が逆~~ → ★**実装した**。`Definition 1.3, (ii)` の本質的一意性（`frobDegUniq`）で、自分で作った `β₀` を与えられた `β` に合わせる | `prop_1_10_i_exists_given` |
| (i) | ~~原文「In this situation, degFr(φ) = degFr(φ′)」が**ファイルに存在しない**~~ → ★**実装した** | `prop_1_10_i_degFr_phi_eq` |
| (i) | ★★原文「then the same is true of **φ′**」の 7 タイプ。`prop_1_10_i_four_types` は `φ` についての主張で `φ′` のものではない。★★**4 つすべて実装した**(`prop_1_10_i_baseIso_of` / `_isometric_of` / `_coAngular_of` / `_lbInvertible_of`)。★co-angular の鍵は「`φ ≫ β` が co-angular」であることだった —— 引き戻す必要はなく、**分解の側を前合成で延ばせばよかった**。★★**7 タイプすべて実装完了**（上の 4 本 ＋ `prop_1_10_i_linear_of` / `prop_1_10_i_preStep_of` / `prop_1_10_i_frobType_of` / `prop_1_10_i_pullBack_of`）。★pull-back は `Proposition 1.4, (ii)` を両向きに使って普遍性を避けた | `prop_1_10_i_four_types` |
| (ii) | ~~`Div` の式が `β′` の base-isomorphism 性を仮定せず、原文の `β′∗`(全単射)の形になっていない~~ → ★**実装した**。`β′` は Frobenius 型ゆえ base-isomorphism なので、`Φ.map (Base β′)` の逆で原文の向きに直せる | `prop_1_10_ii_Div_formula'` |
| (iii) | ~~★原文は証明中で「`A` を Frobenius-trivial としてよい」と**還元している**が、その還元が未実装で、仮定に逃がしている~~ → ★**実装した**。原文が名指す 3 つ(`Definition 1.3, (i), (a)` / `(i), (b)` / `(iii), (c)`)がそのまま 3 段になった。★`isotropic` の最初の用途が**span の pre-step を co-angular にすること**だったと分かった | `prop_1_10_iii_otri_perfect_of_type`, `prop_1_10_iii_otimes_perfect_of_type` |
| (iii) | ~~「the monoids in the image of Φ」は Ob(𝒟) 全体の像だが、実装は `Base` の像のみ~~ → ★**実装した**。`baseSurj` の同型を `MonoidOn.map_bijective_of_iso` で渡り、perfect 性を `isPerfectMonoid_of_bijective` で移した | `prop_1_10_iii_image_perfect_of_base` |
| (iv) | ~~★原文は「**域が** isotropic」だけを要求するが、実装は `∀ X : C, IsIsotropic P X`(**圏全体**)を仮定~~ → ★**実装した**。`isotropicClosed` で域の isotropic 性から下流へ伝播させれば足りた | `prop_1_10_iv`, `prop_1_10_iv_mp` |
| (iv) | ~~「infinitely many isomorphism classes」そのものを述べていない(素数ごとの存在と非同型性まで)~~ → ★**実装した**。★**次数の集合が無限**であることを述べ、非同型性と合わせて同型類の無限性にした | `prop_1_10_iv_degrees_infinite` |
| (v) | ~~`IsPrimeFrobComposite` の基底を「任意の同型」に取った——**意図的な修復**だが、原文項目そのものの実装ではない~~ → ★★**反例に格上げした**。素直な読み(基底＝恒等射)を `IsPrimeFrobCompositeId` として別に定義し、★**恒等射でない同型はその形に書けない**ことを示した。同型は Frobenius 型なので、素直な読みでは iff が破れる —— ★**基底を同型に取ったのは修復ではなく強制であった** | `not_isPrimeFrobCompositeId_of_isIso_of_ne_id` |
| (vi) | ~~★★**`𝒞^istr` がコメント外に一度も現れない**。(vi) の主語が不在~~ → ★**実装した**。`Proposition 1.9` の `istr_frobenioid` / `istr_isotropic` を当てるだけだった | `prop_1_10_vi_istr` |
| (vi) | ~~★`α`,`γ`,`Cc` の存在を**仮定に置いている**が、原文は `Definition 1.3, (i), (a), (b)` から**導いている**~~ → ★**実装した** | `prop_1_10_vi_data`, `prop_1_10_vi_ofType` |
| (vi) | ~~「Moreover, every group-like object … is Frobenius-trivial.」が別宣言で、そちらも `𝒞^istr` への具体化がない~~ → ★**実装した** | `prop_1_10_vi_istr_groupLike` |

★★**ここにある定理はどれも真であり、`sorry` も無い。**
★失われたのは「完全である」という**主張**だけである。
★**上の表がそのまま次の作業のチェックリストになった。**

★★**2026-08-16: 表の 12 行すべてが埋まった。**
★ただし ★**`.src` を付けるのは、私の文脈を持たない検証役が改めて全条を監査してからである** ——
過去 3 回、「埋めた」と数えた項目に**別種の**欠落が見つかっている。
★**表が埋まったことは「もう欠落が無い」ことを意味しない。**

★★**教訓**: 自分の文脈を継いだ検証は、自分の誤読も継ぐ。
「21/21」と数えたのは私であり、その数え方は「**実装したものを数えた**」ものであって、
★**「原文が述べたものを数えた」ものではなかった。**
-/

end ABC3.Found.FrdI
