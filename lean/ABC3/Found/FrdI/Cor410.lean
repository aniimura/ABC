/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def45

/-!
# [FrdI] Corollary 4.10 —— birationalization の圏論性(仮定の型)

原文 (FrdI p.90):
> Corollary 4.10. (Category-theoreticity of the Birationalization) For

★★本ファイルは `Corollary 4.10` の**仮定を型にする**ところまでを担う。

## ★仮定はすべて在庫にあった(2026-08-19)

原文 (FrdI p.90):
> be a divisorial monoid on a connected, totally epimorphic category

原文 (FrdI p.90):
> a Frobenioid of quasi-isotropic type;

| 原文 | 在庫 |
|---|---|
| divisorial | `MonoidOn.IsDivisorialOn` |
| connected, totally epimorphic な `𝒟` | ★`PreFrobenioid` の**フィールド** |
| `𝒟` が FSMFF 型 | `IsOfFSMFFType`(`CategoryVocabulary.lean`) |
| quasi-isotropic 型 | `IsOfQuasiIsotropicType`(`Def31.lean`) |

★★**新規の定義は 1 つも要らない** —— `Theorem 4.2` と違い perf-factorial が要らず、
`𝒟` 側の条件が `FSMFF` 型に替わるだけである。

## ★後段の仮定(rigidity のときだけ増える)

原文 (FrdI p.91):
> an equivalence of categories. Then there exists a 1-unique functor

★原文は最後に「`𝒟₁`, `𝒟₂` が slim で、`𝒞₁`, `𝒞₂` が birationally Frobenius-normalized
型なら合成関手は rigid」と足す。★その 2 条件は `HypRigid` に分けた ——
**本体(1-一意な `Ψ^birat` の存在)には要らない**からである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

namespace Cor410

variable (D C) in
/-- ★★★**[FrdI] Corollary 4.10 の仮定**。

原文 (FrdI p.90):
> a Frobenioid of quasi-isotropic type;

★`divisorial` / `connected` / `totally epimorphic` は `PreFrobenioid` が持っている。 -/
structure Hyp (P : PreFrobenioid C Φ) : Prop where
  /-- `𝒟` は FSMFF 型 -/
  fsmff : IsOfFSMFFType D
  /-- `𝒞` は quasi-isotropic 型 -/
  quasiIsotropic : IsOfQuasiIsotropicType C P

variable (D C) in
/-- ★★**rigidity のときだけ増える 2 条件**。

原文 (FrdI p.91):
> an equivalence of categories. Then there exists a 1-unique functor
-/
structure HypRigid (P : PreFrobenioid C Φ) (G : Frobenioid P) : Prop where
  /-- `𝒟` は slim -/
  slim : IsSlimCat D
  /-- `𝒞` は birationally Frobenius-normalized 型 -/
  biratFrobNormalized : IsOfBirationallyFrobeniusNormalizedType C P G

/-- ★`divisorial` は `PreFrobenioid` から取れる。 -/
theorem divisorial (P : PreFrobenioid C Φ) : Φ.IsDivisorialOn := P.divisorial

/-- ★FSMFF 型からは「自己射に FSMI なものが無い」が直ちに出る。 -/
theorem not_isFSMI_endo {P : PreFrobenioid C Φ} (h : Hyp D C P) {A : D} (φ : End A) :
    ¬ IsFSMI φ := not_isFSMI_endo_of_isOfFSMFFType h.fsmff φ

def Hyp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 90, item := "Corollary 4.10 — 仮定",
    sectionId := "frdi-cor-4-10" }

def HypRigid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91, item := "Corollary 4.10 — rigidity の追加仮定",
    sectionId := "frdi-cor-4-10" }

end Cor410

end ABC3.Found.FrdI
