/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def45

/-!
# [FrdI] Theorem 4.2 —— primary step の圏論性(仮定の型)

原文 (FrdI p.77):
> Theorem 4.2. (Category-theoreticity of Primary Steps) For i = 1, 2,

★★本ファイルは `Theorem 4.2` の**仮定を型にする**ところまでを担う。
主張 (i)(ii) の中身は別ファイルへ回す。

## ★仮定はすべて在庫にあった(2026-08-19)

原文 (FrdI p.77):
> be a perf-factorial divisorial monoid on a connected, totally epimorphic

原文 (FrdI p.77):
> a Frobenioid of standard and isotropic type, which is

| 原文 | 在庫 |
|---|---|
| perf-factorial | `IsPerfFactorial`(`Def24.lean`) |
| divisorial | `MonoidOn.IsDivisorialOn`(`ElementaryFrobenioid.lean`) |
| connected, totally epimorphic な `𝒟` | ★`PreFrobenioid` の**フィールド**(`connectedD` / `totEpiD`) |
| standard 型 | `IsOfStandardType`(`Def31.lean`) |
| isotropic 型 | `IsOfIsotropicType`(`MorphismTypes.lean`) |
| group-like 型でない | `¬ IsOfGroupLikeType`(`MorphismTypes.lean`) |

★★**新規は `IsPerfFactorialOn` の 1 本だけ**である ——
`IsDivisorialOn` と同じ形(`∀ A : 𝒟, …`)で `𝒟` の上へ持ち上げる。

★`connected` と `totally epimorphic` を `Hyp` に入れないのは、
**`PreFrobenioid` を持っている時点で既に持っているから**である
(`P.connectedD` / `P.totEpiD`)。★仮定を二重に書かない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★`perf-factorial` を `𝒟` の上へ -/

variable (Φ) in
/-- ★★**`Φ` が `𝒟` の上で perf-factorial** —— `IsDivisorialOn` と同じ形。

原文 (FrdI p.77):
> be a perf-factorial divisorial monoid on a connected, totally epimorphic
-/
def MonoidOn.IsPerfFactorialOn : Prop := ∀ A : D, IsPerfFactorial (Φ.val A)

def MonoidOn.IsPerfFactorialOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77, item := "Theorem 4.2 — perf-factorial divisorial monoid",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★`Theorem 4.2` の仮定 -/

namespace Thm42

variable (D C) in
/-- ★★★**[FrdI] Theorem 4.2 の仮定**。

原文 (FrdI p.77):
> a Frobenioid of standard and isotropic type, which is

★`connected` / `totally epimorphic` は `PreFrobenioid` が既に持っているので入れない。 -/
structure Hyp (P : PreFrobenioid C Φ) (F : FrobenioidCore P) : Prop where
  /-- `Φ` は perf-factorial -/
  perfFactorial : Φ.IsPerfFactorialOn
  /-- `𝒞` は standard 型 -/
  standard : IsOfStandardType D C P F
  /-- `𝒞` は isotropic 型 -/
  isotropic : IsOfIsotropicType P
  /-- `𝒞` は group-like 型**ではない** -/
  notGroupLike : ¬ IsOfGroupLikeType P

/-- ★`divisorial` は `PreFrobenioid` から取れる。 -/
theorem divisorial (P : PreFrobenioid C Φ) : Φ.IsDivisorialOn := P.divisorial

/-- ★`𝒞` の対象はすべて isotropic。 -/
theorem isIsotropic {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
    (h : Hyp D C P F) (A : C) : IsIsotropic P A := h.isotropic A

/-- ★group-like でない対象が**存在する**。 -/
theorem exists_not_isGroupLikeObj {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
    (h : Hyp D C P F) : ∃ A : C, ¬ IsGroupLikeObj P A := by
  by_contra hc
  exact h.notGroupLike fun A => not_not.mp fun hA => hc ⟨A, hA⟩

def Hyp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77, item := "Theorem 4.2 — 仮定",
    sectionId := "frdi-thm-4-2" }

end Thm42

end ABC3.Found.FrdI
