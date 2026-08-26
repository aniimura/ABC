/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Change

/-!
# [FrdI] Theorem 6.2, (i) の「formally」—— 引き戻しから `Ψ : 𝒞₁ ⥤ 𝒞₂` が出る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> [well-defined up to isomorphism] that is compatible with Frobenius degrees, the

## ★★原文が「formally」で畳んだ段

原文 `Theorem 6.2, (i)` は、支配射 `ψ : V₂ → V₁` から

* 底の関手 `𝒟₁ → 𝒟₂`(体の包含 `K₁ → K₂` が定めるもの)
* 因子の引き戻し `Φ₁ → Φ₂|𝒟₁`
* 有理函数の引き戻し `B₁ → B₂|𝒟₁`

を作ったあと、「**model Frobenioid の圏の定義から formally に**」`Ψ : 𝒞₁ → 𝒞₂` が
出ると言う。★★**その `formally` の中身が本ファイルである。**

★★★在庫の `ModelDataHomOver`(`Thm52Change.lean`)がまさにその入力の形であり、
`Ψ := F.functor` が出る。★残るのは原文が主張する**3 つの両立**であり、
**どれも `rfl` である**:

| 両立 | 宣言 |
|---|---|
| Frobenius 次数 | `thm62Psi_degFr` |
| 底の関手 `𝒟₁ → 𝒟₂` | `thm62Psi_comp_proj` |
| 零因子の引き戻し | `thm62Psi_Div` |
| 有理函数の引き戻し | `thm62Psi_u` |

★★**`rfl` で済むのは `ModelDataHomOver.hom` が 4 成分をそのまま運ぶからである** ——
原文の「formally」はこの事実そのものを指している。

## ★入力について(記録)

本ファイルは**入力(3 つの引き戻し)を仮定として受け取る**。
それを幾何から作るのは依存グラフの節点 `thm62-i-pull` であり、
そこは「正規スキームからの支配射は相対正規化を経由する」
(節点 `normalization-universal-normal`)を要する。
★★**(i) の一般形だけがそれを要求する** —— (ii)(iii)(iv) は `V₁ = V₂` なので
`normMap`(在庫)と同じ筋で出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

variable {D₁ : Type u} [Category.{v} D₁] {D₂ : Type u} [Category.{v} D₂]
  {M₁ : ModelData.{v, u, w} D₁} {M₂ : ModelData.{v, u, w} D₂} {ΨB : D₁ ⥤ D₂}

/-- ★★★★★**[FrdI] Theorem 6.2, (i)** —— 引き戻しの組から出る関手 `Ψ : 𝒞₁ ⥤ 𝒞₂`。

★入力は原文どおり「底の関手 `Ψ_𝒟`、`Φ₁ → Φ₂|𝒟₁`、`B₁ → B₂|𝒟₁`、
`Div_B` との両立」の 4 つ(＝ `ModelDataHomOver` の 4 フィールド)。 -/
noncomputable def thm62Psi (F : ModelDataHomOver ΨB M₁ M₂) :
    ModelData.Obj M₁ ⥤ ModelData.Obj M₂ := F.functor

/-! ## ★原文が主張する 3 つの両立 -/

/-- ★★**Frobenius 次数と両立する**。 -/
@[simp] theorem thm62Psi_degFr (F : ModelDataHomOver ΨB M₁ M₂)
    (h₁ : ModelData.Hyp M₁) (h₂ : ModelData.Hyp M₂)
    {A B : ModelData.Obj M₁} (φ : A ⟶ B) :
    (ModelData.modelPre h₂).degFr ((thm62Psi F).map φ)
      = (ModelData.modelPre h₁).degFr φ := rfl

/-- ★★**底の関手 `𝒟₁ ⥤ 𝒟₂` と両立する**。 -/
theorem thm62Psi_comp_proj (F : ModelDataHomOver ΨB M₁ M₂)
    (h₁ : ModelData.Hyp M₁) (h₂ : ModelData.Hyp M₂) :
    thm62Psi F ⋙ (ModelData.modelPre h₂).proj
      = (ModelData.modelPre h₁).proj ⋙ ΨB := rfl

/-- ★★**零因子の引き戻しと両立する**。 -/
@[simp] theorem thm62Psi_Div (F : ModelDataHomOver ΨB M₁ M₂)
    (h₁ : ModelData.Hyp M₁) (h₂ : ModelData.Hyp M₂)
    {A B : ModelData.Obj M₁} (φ : A ⟶ B) :
    (ModelData.modelPre h₂).Div ((thm62Psi F).map φ)
      = F.phiHom A.base ((ModelData.modelPre h₁).Div φ) := rfl

/-- ★★**有理函数の引き戻しと両立する**。 -/
@[simp] theorem thm62Psi_u (F : ModelDataHomOver ΨB M₁ M₂)
    {A B : ModelData.Obj M₁} (φ : A ⟶ B) :
    ((thm62Psi F).map φ).u = F.bmonHom A.base φ.u := rfl

/-- ★対象の底も `Ψ_𝒟` で送られる。 -/
@[simp] theorem thm62Psi_obj_base (F : ModelDataHomOver ΨB M₁ M₂) (A : ModelData.Obj M₁) :
    ((thm62Psi F).obj A).base = ΨB.obj A.base := rfl

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.2, (i)` の「formally」の段
(★**条つき**: 入力の 3 つの引き戻しは節点 `thm62-i-pull` が供給する)。 -/
def thm62Psi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (i) — 引き戻しから Ψ : 𝒞₁ ⥤ 𝒞₂ が formally に出る",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.FrdI
