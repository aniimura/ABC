/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42

/-!
# [FrdI] Corollary 4.11 —— 初等 Frobenioid への関手の圏論性 I(仮定の型)

原文 (FrdI p.91):
> Corollary 4.11. (Category-theoreticity of the Functor to an Elementary

★★本ファイルは `Corollary 4.11` の**仮定を型にする**ところまでを担う。

## ★仮定はすべて在庫にあった(2026-08-19)

原文 (FrdI p.91):
> connected, totally epimorphic category Di which is Div-slim [with respect to

原文 (FrdI p.91):
> a Frobenioid of standard type;

| 原文 | 在庫 |
|---|---|
| perf-factorial | `MonoidOn.IsPerfFactorialOn`(`Thm42.lean` で置いた) |
| divisorial | `MonoidOn.IsDivisorialOn` |
| connected, totally epimorphic な `𝒟` | ★`PreFrobenioid` の**フィールド** |
| `Div-slim` | `IsDivSlim`(`Def45.lean`) |
| standard 型 | `IsOfStandardType`(`Def31.lean`) |

★★**`Div-slim` は在庫にあった** —— 節点の note が「`IsSlimCat` とは別の概念なので
先に検索せよ」と書いていたとおり別物で、`IsDivSlim Φ := ∀ A, Injective (overPhiAut Φ A)`。
★`isDivSlim_of_isSlim` で `slim ⟹ Div-slim` も既にある。

## ★★`Ψ` に依る条 —— group-like のときだけ増える

原文 (FrdI p.91):
> an equivalence of categories. If C1, C2 are of group-like type, then we also

★★この条は **`Ψ` を含む**ので `Hyp` には入らない。
`Theorem 3.4, (v)` の仮定 (b) と**同じ形**であり、そこでも個別の仮定として渡している。
★★★したがって `HypPair` として**両側 ＋ `Ψ`** をまとめた型を別に置く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3 w3 u4 v4

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

namespace Cor411

variable (D C) in
/-- ★★★**[FrdI] Corollary 4.11 の片側の仮定**。

原文 (FrdI p.91):
> a Frobenioid of standard type;
-/
structure Hyp (P : PreFrobenioid C Φ) (F : FrobenioidCore P) : Prop where
  /-- `Φ` は perf-factorial -/
  perfFactorial : Φ.IsPerfFactorialOn
  /-- `𝒟` は `Φ` に関して Div-slim -/
  divSlim : IsDivSlim Φ
  /-- `𝒞` は standard 型 -/
  standard : IsOfStandardType D C P F

/-- ★`divisorial` は `PreFrobenioid` から取れる。 -/
theorem divisorial (P : PreFrobenioid C Φ) : Φ.IsDivisorialOn := P.divisorial

/-- ★`𝒟` が slim なら `Div-slim` は自動。 -/
theorem divSlim_of_slim (hslim : IsSlimCat D) : IsDivSlim Φ :=
  isDivSlim_of_isSlim Φ hslim

section Pair

variable {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
  {Φ₂ : MonoidOn.{v3, u3, w3} D₂}

variable (D C D₂ C₂) in
/-- ★★★**両側 ＋ `Ψ` の仮定** —— group-like のときだけ増える条を含む。

原文 (FrdI p.91):
> an equivalence of categories. If C1, C2 are of group-like type, then we also

★★`Theorem 3.4, (v)` の仮定 (b) と同じ形である ——
「両方が group-like 型**なら**、`Ψ` とその quasi-inverse が base-isomorphism を保つ」。
★前件が偽(＝どちらかが group-like でない)なら**何も要求しない**。 -/
structure HypPair (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂) (Ψ : C ≌ C₂) : Prop where
  /-- 第 1 の Frobenioid の仮定 -/
  left : Hyp D C P F
  /-- 第 2 の Frobenioid の仮定 -/
  right : Hyp D₂ C₂ P₂ F₂
  /-- ★group-like 型のときだけ効く条 -/
  groupLikeBaseIso : IsOfGroupLikeType P → IsOfGroupLikeType P₂ →
    (∀ {A B : C} (φ : A ⟶ B), IsBaseIsomorphism P φ →
        IsBaseIsomorphism P₂ (Ψ.functor.map φ)) ∧
    (∀ {A B : C₂} (φ : A ⟶ B), IsBaseIsomorphism P₂ φ →
        IsBaseIsomorphism P (Ψ.inverse.map φ))

/-- ★前件が偽のときは `groupLikeBaseIso` は自動的に満たされる。 -/
theorem hypPair_groupLikeBaseIso_of_not
    {P : PreFrobenioid C Φ} {P₂ : PreFrobenioid C₂ Φ₂} {Ψ : C ≌ C₂}
    (h : ¬ IsOfGroupLikeType P) :
    IsOfGroupLikeType P → IsOfGroupLikeType P₂ →
      (∀ {A B : C} (φ : A ⟶ B), IsBaseIsomorphism P φ →
          IsBaseIsomorphism P₂ (Ψ.functor.map φ)) ∧
      (∀ {A B : C₂} (φ : A ⟶ B), IsBaseIsomorphism P₂ φ →
          IsBaseIsomorphism P (Ψ.inverse.map φ)) :=
  fun hg _ => absurd hg h

end Pair

def Hyp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91, item := "Corollary 4.11 — 仮定",
    sectionId := "frdi-cor-4-11" }

def HypPair.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91, item := "Corollary 4.11 — Ψ に依る条",
    sectionId := "frdi-cor-4-11" }

end Cor411

end ABC3.Found.FrdI
