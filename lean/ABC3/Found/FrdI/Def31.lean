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
| (ii) | 4 | `𝒞^Fr-tp` / `𝒞^bi-Fr` と `Hom^pf_𝒞(A,B)`(**帰納極限**) | 未 |
| (iii) | 5 | `𝒞^pf`(Frobenioid の**完備化**)と `𝒞 → 𝒞^pf` | 未 |
| (iv) | 6 | `unit-equivalent` と `Hom^un-tr` | 未(`Proposition 3.3, (ii)` を前方参照する) |
| (iv) | 7 | `𝒞^un-tr` | 未 |

★**(ii)(iii)(iv) が入るまで `.src` は付けない。**

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

end ABC3.Found.FrdI
