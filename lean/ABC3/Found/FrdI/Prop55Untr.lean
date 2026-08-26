/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def31Pf
import ABC3.Found.FrdI.HomColimQuot
import ABC3.Found.FrdI.UnTr

/-!
# [FrdI] Proposition 5.5, (ii) —— 完全化と unit-trivial 化は交換する(射の集合の層)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> erality that C is of isotropic type. Then it follows immediately from the definition of

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★原文が「直ちに従う」と畳んだものの中身

原文は「`Definition 3.1, (iii)` の完全化の定義から、**2 対象の像のあいだの射の集合の
自然な全単射**を作れば足りる」と言う。★その全単射は 2 段でできている:

1. **帰納極限は商と交換する** —— `HomColimQuot.lean` の `quotEquiv`。
   `(colim_Z Hom(Z₁,Z₂))/∼` と `colim_Z (Hom(Z₁,Z₂)/∼)` が一致する。
2. **各段の商が `Hom^un-tr` である** —— `stageEquiv`。

★★本ファイルは **1 と 2** を閉じる。

## ★★★原文が畳んだ計算 —— 遷移写像は unit-equivalence を保つ

★1 の `compat`(遷移写像が商へ降りること)は自明ではない。
`Definition 3.1, (ii)` の遷移写像は `Proposition 1.10, (i)` の四角形

  `φ ≫ β = α ≫ φ′`

で定まる。`φ₁ ∼ φ₂`(すなわち `𝔽_Φ` で同じ像)から `φ₁′ ∼ φ₂′` を出すには、
**`𝔽_Φ` の中で `α` の像を左から消約する**必要がある。

★★**それができる**理由(`elemFrob_cancel_of_baseIso`):
`α` は Frobenius 型だから `Base(α)` は同型で、`𝔽_Φ` の射の 3 成分がそれぞれ消約できる。

| 成分 | 消約の根拠 |
|---|---|
| `deg` | `ℕ≥1` は簡約的 |
| `base` | `Base(α)` が同型 ⟹ エピ |
| `div` | `Φ(A)` が簡約的(`Φ` は divisorial ⟹ integral)＋ `Φ.map` が単射 |

★★★**`div` 成分に `Φ` の integral 性が効く** —— ここが原文の「immediately」の実体である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `elemFrob_cancel_of_baseIso` | `𝔽_Φ` で base が同型な射は**エピ** |
| `frobTransport_toElem_congr` | ★遷移写像は unit-equivalence を保つ |
| `unTrQuotPf` | `Hom^pf` の帰納系の上の同値関係の族 |
| `stageEquiv` | 各段の商は `Hom^un-tr` |
| `homPfUnTrEquiv` | ★★★**`Hom^pf(A,B)/∼ ≃ colim_Z Hom^un-tr(Z₁,Z₂)`** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

/-! ## ★1. `𝔽_Φ` の中の消約 -/

section ElemCancel

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

/-- ★★★**`𝔽_Φ` で base が同型な射は左から消約できる**。

★3 成分をそれぞれ消約する:
`deg` は `ℕ≥1` の簡約性、`base` は `Base(α)` が同型であること、
`div` は `Φ(A)` の簡約性と `Φ.map` の単射性。 -/
theorem elemFrob_cancel_of_baseIso (hcanc : ∀ A : D, IsCancelAdd (Φ.val A))
    {X Y Z : ElemFrobCat Φ} (a : X ⟶ Y) (ha : IsIso a.base) {f g : Y ⟶ Z}
    (h : a ≫ f = a ≫ g) : f = g := by
  haveI := ha
  have hdeg : f.deg = g.deg := by
    have hd : f.deg * a.deg = g.deg * a.deg := by
      simpa using congrArg ElemFrobCat.Hom.deg h
    exact mul_right_cancel hd
  refine ElemFrobCat.Hom.ext ?_ ?_ hdeg
  · have hb : a.base ≫ f.base = a.base ≫ g.base := by
      simpa using congrArg ElemFrobCat.Hom.base h
    exact (cancel_epi a.base).mp hb
  · have hdv : Φ.map a.base f.div + (f.deg : ℕ) • a.div
        = Φ.map a.base g.div + (g.deg : ℕ) • a.div := by
      simpa using congrArg ElemFrobCat.Hom.div h
    rw [hdeg] at hdv
    letI := hcanc X.base
    exact Φ.map_injective _ (add_right_cancel hdv)

end ElemCancel

/-! ## ★2. 遷移写像は unit-equivalence を保つ -/

section Transport

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-- ★`Φ` は divisorial なので各 `Φ(A)` は簡約的。 -/
theorem phi_isCancelAdd (P : PreFrobenioid C Φ) (A : D) : IsCancelAdd (Φ.val A) :=
  isCancelAdd_of_isIntegralMonoid _ (P.divisorial A).1.1

/-- ★★★**遷移写像は `𝔽_Φ` での像の等しさを保つ**。

★★原文の `Definition 3.1, (ii)` の遷移写像は四角形 `φ ≫ β = α ≫ φ′` で定まる。
`φ₁`, `φ₂` が `𝔽_Φ` で同じ像なら、四角形の左辺が一致し、
`α` の像が**エピ**(`elemFrob_cancel_of_baseIso`)なので右辺も一致する。 -/
theorem frobTransport_toElem_congr {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
    {A' B' A'' B'' : C}
    (α : A' ⟶ A'') (hα : IsFrobeniusType P α) (β : B' ⟶ B'') (hβ : IsFrobeniusType P β)
    (hd : P.degFr α = P.degFr β) {φ ψ : A' ⟶ B'}
    (h : P.toElem.map φ = P.toElem.map ψ) :
    P.toElem.map (frobTransport (F := F) α hα β hβ hd φ)
      = P.toElem.map (frobTransport (F := F) α hα β hβ hd ψ) := by
  have h1 := congrArg P.toElem.map (frobTransport_spec (F := F) α hα β hβ hd φ)
  have h2 := congrArg P.toElem.map (frobTransport_spec (F := F) α hα β hβ hd ψ)
  rw [P.toElem.map_comp, P.toElem.map_comp] at h1 h2
  rw [h] at h1
  refine elemFrob_cancel_of_baseIso (phi_isCancelAdd P) (P.toElem.map α) hα.2 ?_
  rw [← h1, ← h2]

/-! ## ★3. `Hom^pf` の帰納系の上の同値関係の族 -/

/-- ★★**`Hom^pf` の帰納系の上の unit-equivalence**。

★各段は `𝔽_Φ` での像の等しさ(`Proposition 3.3, (ii)` により unit-equivalence)。 -/
noncomputable def unTrQuotPf (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (A B : C) : HomColim.QuotData (homFunctorPf P F A B) where
  setoid Z :=
    { r := fun φ ψ => P.toElem.map φ.down = P.toElem.map ψ.down
      iseqv := ⟨fun _ => rfl, Eq.symm, Eq.trans⟩ }
  compat {Z W} u {φ ψ} h := by
    show P.toElem.map (idxTransport P F u φ.down) = P.toElem.map (idxTransport P F u ψ.down)
    exact frobTransport_toElem_congr _ u.right.property.1 _ u.right.property.2.1
      u.right.property.2.2 h

/-- ★各段の商は `Hom^un-tr`。 -/
noncomputable def stageEquiv (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    {A B : C} (Z : IdxPf P F A B) :
    Quotient ((unTrQuotPf P F A B).setoid Z)
      ≃ ULift.{u2} (HomUnTr P Z.right.obj.1 Z.right.obj.2) where
  toFun := Quotient.lift (fun φ => ULift.up (toHomUnTr P φ.down))
    (fun a b h => congrArg ULift.up ((toHomUnTr_eq_iff P a.down b.down).mpr h))
  invFun y := Quotient.liftOn y.down
    (fun φ => Quotient.mk ((unTrQuotPf P F A B).setoid Z) (ULift.up φ))
    (fun a b h => Quotient.sound ((toHomUnTr_eq_iff P a b).mp (Quotient.sound h)))
  left_inv z := by
    refine Quotient.inductionOn z (fun φ => ?_)
    rfl
  right_inv y := by
    obtain ⟨y⟩ := y
    refine Quotient.inductionOn y (fun φ => ?_)
    rfl

/-! ## ★4. 主定理 -/

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)**(un-tr の側、射の集合の層)——
`Hom^pf(A,B)` を unit-equivalence で割ったものは、
各段の `Hom^un-tr` の**帰納極限**に一致する。

★左辺が `(𝒞^pf)^un-tr` の射の集合、右辺が `(𝒞^un-tr)^pf` の射の集合である。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C -/
noncomputable def homPfUnTrEquiv (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (A B : C) :
    Quotient (HomColim.quotKer (unTrQuotPf P F A B))
      ≃ HomColim ((unTrQuotPf P F A B).functor) :=
  HomColim.quotEquiv (unTrQuotPf P F A B)

/-- ★★**核の完全な記述** —— `Hom^pf` の 2 元が `(𝒞^pf)^un-tr` で同じ射になるのは、
**ある段で `𝔽_Φ` の像が一致する**とき、ちょうどそのとき。 -/
theorem homPf_untr_eq_iff {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
    {A B : C} {Z W : IdxPf P F A B}
    (φ : ULift.{u2} (Z.right.obj.1 ⟶ Z.right.obj.2))
    (ψ : ULift.{u2} (W.right.obj.1 ⟶ W.right.obj.2)) :
    HomColim.toQuot (unTrQuotPf P F A B) (HomColim.mk _ Z φ)
        = HomColim.toQuot (unTrQuotPf P F A B) (HomColim.mk _ W ψ)
      ↔ ∃ (V : IdxPf P F A B) (u : Z ⟶ V) (u' : W ⟶ V),
          P.toElem.map (idxTransport P F u φ.down)
            = P.toElem.map (idxTransport P F u' ψ.down) :=
  HomColim.toQuot_eq_iff (unTrQuotPf P F A B)

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)`(un-tr の側)。 -/
def homPfUnTrEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — (𝒞^pf)^un-tr と (𝒞^un-tr)^pf の射の集合の自然な全単射",
    sectionId := "frdi-prop-5-5" }

end Transport

end ABC3.Found.FrdI
