/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import Mathlib.NumberTheory.PrimeCounting
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.MonoidPrime

/-!
# Prop110 —— §0 の語彙・(i) の一意性と存在

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
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

★★**2026-08-16 に `CategoryVocabulary.lean` へ移した。**
`Proposition 1.14` が同じ段落の `FSMI` / `FSMFF-type` / `subordinate` を要求し、
★**それらは `irreducible` の上に建つ**ので、§0 の語彙は 1 か所にまとめた。
名前空間は同じ(`ABC3.Found.FrdI`)なので、ここでの参照はそのまま通る。
-/

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

end ABC3.Found.FrdI
