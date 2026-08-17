import ABC3.Found.FrdI.CategoryAnchor
import ABC3.Found.FrdI.Prop110

/-!
# [FrdI] Definition 3.1, (i) —— quasi-isotropic 型 / standard 型 / Frobenius-slim

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.56–p.57。

原文 (FrdI p.56):
> (i) We shall say that C is of quasi-isotropic type if it holds that A

## ★この定義の規模(測定)

`Definition 3.1` は **4 条**あり、主張は **7**:

| 条 | # | 内容 | 状態 |
|---|---|---|---|
| (i) | 1 | `quasi-isotropic 型` | ★**ここで実装** |
| (i) | 2 | `standard 型`((a)–(e)) | ★**ここで実装** |
| (i) | 3 | `Frobenius-slim` な圏 | ★**ここで実装** |
| (ii) | 4 | `𝒞^Fr-tp` / `𝒞^bi-Fr` と `Hom^pf_𝒞(A,B)`(**帰納極限**) | ★**`Def31Pf.lean`** |
| (iii) | 5 | `𝒞^pf`(Frobenioid の**完備化**)と `𝒞 → 𝒞^pf` | ★**`Def31Pf.lean`** |
| (iv) | 6 | `unit-equivalent` と `Hom^un-tr` | ★**関係の定義はここ**、`Hom^un-tr` は `UnTr.lean` |
| (iv) | 7 | `𝒞^un-tr` | ★**`UnTr.lean`** |

★★**7 主張すべてが揃った(2026-08-17)。**
★条なしの `.src`(`definition_3_1.src`)は **`Def31Pf.lean` に 1 つ**置いてある ——
実装が 3 ファイルに分かれているので、対応表と一緒に 1 か所へ集約した。

## ★`Frobenius-slim` の `F` について

原文の `F` は **`Definition 1.1, (iii)` の standard Frobenioid `𝔽 = 𝔽_{ℤ≥0}`**
(`ElemFrob.Standard`)であり、`F ↠ N≥1` は次数を取る準同型
(`ElemFrob.degHom`)である。★**原文は `F` の定義を `Definition 1.1, (iii)` に
送っているだけなので、そこを引かないと読めない。**

★`E_A → E` はスライス `Over A` からの忘却関手であり、
`Aut(E_A → E)` はその**関手としての自己同型群**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★`Frobenius-slim`

原文 (FrdI p.56):
> Frobenius-slim if every homomorphism of monoids
-/

/-- ★★**[FrdI] Definition 3.1, (i)** —— 圏 `E` が **Frobenius-slim** であること。

★`𝔽 = 𝔽_{ℤ≥0}`(`ElemFrob.Standard`)からの任意のモノイド準同型
`𝔽 → Aut(E_A → E)` が、次数を取る全射 `𝔽 ↠ ℕ≥1` を経由する。 -/
def IsFrobeniusSlim (E : Type u2) [Category.{v2} E] : Prop :=
  ∀ (A : E) (f : ElemFrob.Standard →* Aut (Over.forget A)),
    ∃ g : ℕ+ →* Aut (Over.forget A), f = g.comp ElemFrob.degHom

/-- ★**非退化(上)** —— `Aut(E_A → E)` が自明な圏は Frobenius-slim。

★**`slim` な圏はすべて Frobenius-slim**(原文の「Thus, every slim category is
Frobenius-slim.」)の、最も退化した場合である。 -/
theorem isFrobeniusSlim_of_subsingleton (E : Type u2) [Category.{v2} E]
    (h : ∀ A : E, Subsingleton (Aut (Over.forget A))) : IsFrobeniusSlim E := by
  intro A f
  exact ⟨1, MonoidHom.ext fun x => (h A).elim _ _⟩

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★`quasi-isotropic 型`

原文 (FrdI p.56):
> is non-isotropic if and only if it is an iso-subanchor [cf.

★**「`A` が non-isotropic ⟺ `A` が iso-subanchor」**である。
`iso-subanchor` は §0 の語彙(`CategoryAnchor.lean`)。
-/

variable (C) in
/-- ★★★**[FrdI] Definition 3.1, (i)** —— `𝒞` が **quasi-isotropic 型**。 -/
def IsOfQuasiIsotropicType : Prop :=
  ∀ A : C, ¬ IsIsotropic P A ↔ IsIsoSubanchor C A

/-! ## ★`standard 型`

原文 (FrdI p.56):
> shall say that C is of standard type if the following conditions are satisfied: (a)

★**(a)–(e) の 5 条件**をそのまま構造体の 5 フィールドにする。
-/

variable (C D) in
/-- ★★★**[FrdI] Definition 3.1, (i)** —— `𝒞` が **standard 型**。

| 原文 | フィールド |
|---|---|
| (a) quasi-isotropic かつ Frobenius-isotropic 型 | `quasiIsotropic` / `frobIsotropic` |
| (b) group-like 型なら `𝒞^istr` が Frobenius-compact 対象を持つ | `groupLikeCompact` |
| (c) Frobenius-normalized 型 | `frobNormalized` |
| (d) `𝒟` が FSMFF-type | `baseFSMFF` |
| (e) `Φ` が non-dilating | `phiNonDilating` | -/
structure IsOfStandardType (F : FrobenioidCore P) : Prop where
  /-- **(a)** quasi-isotropic 型。 -/
  quasiIsotropic : IsOfQuasiIsotropicType C P
  /-- **(a)** Frobenius-isotropic 型。 -/
  frobIsotropic : IsOfFrobeniusIsotropicType P
  /-- **(b)** group-like 型なら `𝒞^istr` が Frobenius-compact 対象を持つ。 -/
  groupLikeCompact : IsOfGroupLikeType P → ∃ A : Istr P, IsFrobeniusCompact (istrPre P F) A
  /-- **(c)** Frobenius-normalized 型。 -/
  frobNormalized : IsOfFrobeniusNormalizedType P
  /-- **(d)** `𝒟` が FSMFF-type。 -/
  baseFSMFF : IsOfFSMFFType D
  /-- **(e)** `Φ` が non-dilating。 -/
  phiNonDilating : MonoidOn.IsNonDilatingOn Φ

/-! ## ★(iv) —— `unit-equivalent`

原文 (FrdI p.57):
> unit-equivalent if there exist morphisms

★**`α₁ = β ∘ γ`、`α₂ = β ∘ δ ∘ γ`** —— Lean の `≫` では
`α₁ = γ ≫ β`、`α₂ = γ ≫ δ ≫ β`。

★★**同値関係であることと合成で閉じることは `Proposition 3.3, (ii)` が与える**
(原文が明示的に前方参照している)。ここでは関係の定義だけを置く。
-/

/-- ★★★**[FrdI] Definition 3.1, (iv)** —— `𝒞^istr` の co-objective な 2 射が
**unit-equivalent** であること。

★中間対象 `Cc` と `γ : A ⟶ Cc`、`β : Cc ⟶ B`、`δ ∈ 𝒪^×(Cc)` があって
`α₁ = γ ≫ β`、`α₂ = γ ≫ δ ≫ β` となること。 -/
def IsUnitEquivalent {A B : C} (α₁ α₂ : A ⟶ B) : Prop :=
  ∃ (Cc : C) (γ : A ⟶ Cc) (β : Cc ⟶ B) (δ : End Cc), δ ∈ OTimes P Cc ∧
    α₁ = γ ≫ β ∧ α₂ = γ ≫ ((δ : Cc ⟶ Cc)) ≫ β

/-- ★**反射的** —— `δ = 1` を取ればよい。 -/
theorem isUnitEquivalent_refl {A B : C} (α : A ⟶ B) : IsUnitEquivalent P α α :=
  ⟨A, 𝟙 A, α, 1, (OTimes P A).one_mem, by simp, by simp⟩

/-- ★**単元は `𝒪^×` の元** —— `unit-equivalent` の `δ` は同型。 -/
theorem isIso_of_mem_otimes {Cc : C} {δ : End Cc} (h : δ ∈ OTimes P Cc) :
    IsIso ((δ : Cc ⟶ Cc)) := (CategoryTheory.isUnit_iff_isIso _).mp h.2

end ABC3.Found.FrdI
