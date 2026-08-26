/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Change

/-!
# [FrdI] model Frobenioid の 1-可換な四角形

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.103–104。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of

★★`Proposition 5.3` の 1-可換図式と `Corollary 5.4` の (4)(縦の矢印との 1-可換性)は
**どちらも同じ形**である ——

```
    E₁  --H₁-->  Obj M₁  --Fo-->  Obj M₂
     |                              ‖
     Ψ                              ‖
     v                              ‖
    E₂  --------- H₂ -------->  Obj M₂
```

★★model Frobenioid の射は **4 成分の明示的な組**(`base`, `div`, `deg`, `u`)なので、
自然同型を作るには**成分ごとの一致を与えればよい**。本ファイルはその一般形を置く。

| 宣言 | 中身 |
|---|---|
| `ModelData.objIsoOfBase` | ★底の同型＋類の一致 ⟹ 対象の同型(`u = 0` の場合) |
| `ModelData.objIsoOfBaseU` | ★★★同上、**単元 `u` を許す一般形**(`B` は group-like) |
| `ModelData.squareOfBase` | ★★★★1-可換な四角形(`u = 0` の場合) |
| `ModelData.squareOfBaseU` | ★★★★★同上、単元を許す一般形 |

## ★実務メモ

★`Theorem 5.2, (iv)` の圏同値がつくる関手(`modelType_equiv`)は
`PathCat` の同値の逆を経由するので、**成分が直接には計算できない**。
★したがって本ファイルの補題は「成分の一致」を**仮定として受け取る**形にしてある ——
その仮定を `Corollary 4.10` / `Corollary 4.11, (iii), (iv)` から導くのが残りの工程である。

★★対象の同型で `u ≠ 0` を許すのが要点である。model 型の Frobenioid では
対象の類 `cls` は**有理関数の単系 `B` の分だけ動く**ので、
`u = 0` に固定すると `Corollary 5.4` の図式には当たらない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u3 v3 u4 v4

namespace ModelData

/-! ## ★1. 対象の同型 -/

section ObjIso

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

/-- ★★**底の同型と類の一致から model Frobenioid の同型を作る**(`u = 0` の場合)。 -/
def objIsoOfBase {X Y : Obj M} (e : X.base ≅ Y.base)
    (hc : X.cls = M.phi.gpMapOn e.hom Y.cls) : X ≅ Y where
  hom :=
    { base := e.hom
      div := 0
      deg := 1
      u := 0
      cond := by
        show ((1 : ℕ+) : ℕ) • X.cls + toGpHom _ (0 : M.phi.val X.base)
          = M.phi.gpMapOn e.hom Y.cls + M.divB _ (0 : M.bmon.val X.base)
        rw [map_zero, map_zero, add_zero, add_zero]
        simpa using hc }
  inv :=
    { base := e.inv
      div := 0
      deg := 1
      u := 0
      cond := by
        show ((1 : ℕ+) : ℕ) • Y.cls + toGpHom _ (0 : M.phi.val Y.base)
          = M.phi.gpMapOn e.inv X.cls + M.divB _ (0 : M.bmon.val Y.base)
        rw [map_zero, map_zero, add_zero, add_zero, hc, ← M.phi.gpMapOn_comp,
          e.inv_hom_id, M.phi.gpMapOn_id]
        simp }
  hom_inv_id := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact e.hom_inv_id
    · show M.phi.map e.hom (0 : M.phi.val Y.base) + ((1 : ℕ+) : ℕ) • (0 : M.phi.val X.base) = 0
      simp
    · show (1 : ℕ+) * (1 : ℕ+) = 1
      simp
    · show M.bmon.map e.hom (0 : M.bmon.val Y.base)
        + ((1 : ℕ+) : ℕ) • (0 : M.bmon.val X.base) = 0
      simp
  inv_hom_id := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact e.inv_hom_id
    · show M.phi.map e.inv (0 : M.phi.val X.base) + ((1 : ℕ+) : ℕ) • (0 : M.phi.val Y.base) = 0
      simp
    · show (1 : ℕ+) * (1 : ℕ+) = 1
      simp
    · show M.bmon.map e.inv (0 : M.bmon.val X.base)
        + ((1 : ℕ+) : ℕ) • (0 : M.bmon.val Y.base) = 0
      simp

@[simp] theorem objIsoOfBase_hom_base {X Y : Obj M} (e : X.base ≅ Y.base)
    (hc : X.cls = M.phi.gpMapOn e.hom Y.cls) :
    (objIsoOfBase e hc).hom.base = e.hom := rfl

/-- ★★★**底の同型・単元・類の関係から model Frobenioid の同型を作る**(一般形)。

★`u` は `Theorem 5.2, (iv)` の `u_f` 成分に当たる。`B` は group-like なので
逆向きの単元は `neg` で作れる(`RatFnData` の `bneg` がそれを与える)。 -/
def objIsoOfBaseU (neg : ∀ d : D, M.bmon.val d → M.bmon.val d)
    (hneg : ∀ (d : D) (x : M.bmon.val d), neg d x + x = 0)
    {X Y : Obj M} (e : X.base ≅ Y.base) (u : M.bmon.val X.base)
    (hc : X.cls = M.phi.gpMapOn e.hom Y.cls + M.divB _ u) : X ≅ Y where
  hom :=
    { base := e.hom
      div := 0
      deg := 1
      u := u
      cond := by
        show ((1 : ℕ+) : ℕ) • X.cls + toGpHom _ (0 : M.phi.val X.base)
          = M.phi.gpMapOn e.hom Y.cls + M.divB _ u
        rw [map_zero, add_zero, PNat.one_coe, one_smul]
        exact hc }
  inv :=
    { base := e.inv
      div := 0
      deg := 1
      u := M.bmon.map e.inv (neg _ u)
      cond := by
        have h0 : M.divB X.base u + M.divB X.base (neg _ u) = 0 := by
          rw [← map_add, add_comm, hneg, map_zero]
        show ((1 : ℕ+) : ℕ) • Y.cls + toGpHom _ (0 : M.phi.val Y.base)
          = M.phi.gpMapOn e.inv X.cls + M.divB _ (M.bmon.map e.inv (neg _ u))
        rw [map_zero, add_zero, PNat.one_coe, one_smul, M.divB_nat, hc, map_add,
          ← M.phi.gpMapOn_comp, e.inv_hom_id, M.phi.gpMapOn_id, add_assoc, ← map_add, h0,
          map_zero, add_zero] }
  hom_inv_id := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact e.hom_inv_id
    · show M.phi.map e.hom (0 : M.phi.val Y.base) + ((1 : ℕ+) : ℕ) • (0 : M.phi.val X.base) = 0
      simp
    · show (1 : ℕ+) * (1 : ℕ+) = 1
      simp
    · show M.bmon.map e.hom (M.bmon.map e.inv (neg _ u)) + ((1 : ℕ+) : ℕ) • u = 0
      rw [← MonoidOn.map_comp, e.hom_inv_id, MonoidOn.map_id, PNat.one_coe, one_smul]
      exact hneg _ u
  inv_hom_id := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact e.inv_hom_id
    · show M.phi.map e.inv (0 : M.phi.val X.base) + ((1 : ℕ+) : ℕ) • (0 : M.phi.val Y.base) = 0
      simp
    · show (1 : ℕ+) * (1 : ℕ+) = 1
      simp
    · show M.bmon.map e.inv u + ((1 : ℕ+) : ℕ) • M.bmon.map e.inv (neg _ u) = 0
      rw [PNat.one_coe, one_smul, ← map_add, add_comm, hneg, map_zero]

@[simp] theorem objIsoOfBaseU_hom_base (neg : ∀ d : D, M.bmon.val d → M.bmon.val d)
    (hneg : ∀ (d : D) (x : M.bmon.val d), neg d x + x = 0)
    {X Y : Obj M} (e : X.base ≅ Y.base) (u : M.bmon.val X.base)
    (hc : X.cls = M.phi.gpMapOn e.hom Y.cls + M.divB _ u) :
    (objIsoOfBaseU neg hneg e u hc).hom.base = e.hom := rfl

/-- ★★locator —— model Frobenioid の対象の同型(`Theorem 5.2, (iv)` の射の形から)。 -/
def objIsoOfBaseU.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 100,
    item := "Theorem 5.2, (iv) — model Frobenioid の対象の同型(底の同型＋単元)",
    sectionId := "frdi-thm-5-2" }

end ObjIso

/-! ## ★2. 1-可換な四角形 -/

section Square

variable {D₁ : Type u} [Category.{v} D₁] {D₂ : Type u} [Category.{v} D₂]
  {M₁ : ModelData.{v, u, w} D₁} {M₂ : ModelData.{v, u, w} D₂} {ΨB : D₁ ⥤ D₂}

/-- ★★★★**1-可換性の一般形**(`u = 0` の場合)——
底の同型と 3 成分の一致から自然同型を作る。 -/
noncomputable def squareOfBase (Fo : ModelDataHomOver ΨB M₁ M₂)
    {E₁ : Type u3} [Category.{v3} E₁] {E₂ : Type u4} [Category.{v4} E₂]
    (H₁ : E₁ ⥤ Obj M₁) (Ψ : E₁ ⥤ E₂) (H₂ : E₂ ⥤ Obj M₂)
    (e : ∀ X : E₁, ((H₁ ⋙ Fo.functor).obj X).base ≅ ((Ψ ⋙ H₂).obj X).base)
    (hc : ∀ X : E₁, ((H₁ ⋙ Fo.functor).obj X).cls
      = M₂.phi.gpMapOn (e X).hom ((Ψ ⋙ H₂).obj X).cls)
    (hbase : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).base ≫ (e Y).hom = (e X).hom ≫ ((Ψ ⋙ H₂).map f).base)
    (hdiv : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).div = M₂.phi.map (e X).hom ((Ψ ⋙ H₂).map f).div)
    (hdeg : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).deg = ((Ψ ⋙ H₂).map f).deg)
    (hu : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).u = M₂.bmon.map (e X).hom ((Ψ ⋙ H₂).map f).u) :
    H₁ ⋙ Fo.functor ≅ Ψ ⋙ H₂ :=
  NatIso.ofComponents (fun X => objIsoOfBase (e X) (hc X)) (fun {X Y} f => by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact hbase f
    · show M₂.phi.map ((H₁ ⋙ Fo.functor).map f).base (0 : M₂.phi.val _)
          + ((1 : ℕ+) : ℕ) • ((H₁ ⋙ Fo.functor).map f).div
        = M₂.phi.map (e X).hom ((Ψ ⋙ H₂).map f).div
          + ((((Ψ ⋙ H₂).map f).deg : ℕ+) : ℕ) • (0 : M₂.phi.val _)
      rw [map_zero, zero_add, smul_zero, add_zero, PNat.one_coe, one_smul]
      exact hdiv f
    · show (1 : ℕ+) * ((H₁ ⋙ Fo.functor).map f).deg = ((Ψ ⋙ H₂).map f).deg * (1 : ℕ+)
      rw [one_mul, mul_one]
      exact hdeg f
    · show M₂.bmon.map ((H₁ ⋙ Fo.functor).map f).base (0 : M₂.bmon.val _)
          + ((1 : ℕ+) : ℕ) • ((H₁ ⋙ Fo.functor).map f).u
        = M₂.bmon.map (e X).hom ((Ψ ⋙ H₂).map f).u
          + ((((Ψ ⋙ H₂).map f).deg : ℕ+) : ℕ) • (0 : M₂.bmon.val _)
      rw [map_zero, zero_add, smul_zero, add_zero, PNat.one_coe, one_smul]
      exact hu f)

/-- ★★★★★**1-可換性の一般形**(単元 `u` を許す版)。

★★`Corollary 5.4` の (4) と `Proposition 5.3` の 1-可換図式は**この形**である ——
対象の類 `cls` は有理関数の単系 `B` の分だけ動くので、`u = 0` では当たらない。 -/
noncomputable def squareOfBaseU (Fo : ModelDataHomOver ΨB M₁ M₂)
    (neg : ∀ d : D₂, M₂.bmon.val d → M₂.bmon.val d)
    (hneg : ∀ (d : D₂) (x : M₂.bmon.val d), neg d x + x = 0)
    {E₁ : Type u3} [Category.{v3} E₁] {E₂ : Type u4} [Category.{v4} E₂]
    (H₁ : E₁ ⥤ Obj M₁) (Ψ : E₁ ⥤ E₂) (H₂ : E₂ ⥤ Obj M₂)
    (e : ∀ X : E₁, ((H₁ ⋙ Fo.functor).obj X).base ≅ ((Ψ ⋙ H₂).obj X).base)
    (uu : ∀ X : E₁, M₂.bmon.val ((H₁ ⋙ Fo.functor).obj X).base)
    (hc : ∀ X : E₁, ((H₁ ⋙ Fo.functor).obj X).cls
      = M₂.phi.gpMapOn (e X).hom ((Ψ ⋙ H₂).obj X).cls + M₂.divB _ (uu X))
    (hbase : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).base ≫ (e Y).hom = (e X).hom ≫ ((Ψ ⋙ H₂).map f).base)
    (hdiv : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).div = M₂.phi.map (e X).hom ((Ψ ⋙ H₂).map f).div)
    (hdeg : ∀ {X Y : E₁} (f : X ⟶ Y),
      ((H₁ ⋙ Fo.functor).map f).deg = ((Ψ ⋙ H₂).map f).deg)
    (hu : ∀ {X Y : E₁} (f : X ⟶ Y),
      M₂.bmon.map ((H₁ ⋙ Fo.functor).map f).base (uu Y) + ((H₁ ⋙ Fo.functor).map f).u
        = M₂.bmon.map (e X).hom ((Ψ ⋙ H₂).map f).u
          + ((((Ψ ⋙ H₂).map f).deg : ℕ+) : ℕ) • (uu X)) :
    H₁ ⋙ Fo.functor ≅ Ψ ⋙ H₂ :=
  NatIso.ofComponents (fun X => objIsoOfBaseU neg hneg (e X) (uu X) (hc X)) (fun {X Y} f => by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact hbase f
    · show M₂.phi.map ((H₁ ⋙ Fo.functor).map f).base (0 : M₂.phi.val _)
          + ((1 : ℕ+) : ℕ) • ((H₁ ⋙ Fo.functor).map f).div
        = M₂.phi.map (e X).hom ((Ψ ⋙ H₂).map f).div
          + ((((Ψ ⋙ H₂).map f).deg : ℕ+) : ℕ) • (0 : M₂.phi.val _)
      rw [map_zero, zero_add, smul_zero, add_zero, PNat.one_coe, one_smul]
      exact hdiv f
    · show (1 : ℕ+) * ((H₁ ⋙ Fo.functor).map f).deg = ((Ψ ⋙ H₂).map f).deg * (1 : ℕ+)
      rw [one_mul, mul_one]
      exact hdeg f
    · show M₂.bmon.map ((H₁ ⋙ Fo.functor).map f).base (uu Y)
          + ((1 : ℕ+) : ℕ) • ((H₁ ⋙ Fo.functor).map f).u
        = M₂.bmon.map (e X).hom ((Ψ ⋙ H₂).map f).u
          + ((((Ψ ⋙ H₂).map f).deg : ℕ+) : ℕ) • (uu X)
      rw [PNat.one_coe, one_smul]
      exact hu f)

/-- ★★★★★locator —— `Corollary 5.4` の (4)(縦の矢印との 1-可換性)の一般形
(★**条つき**: 成分の一致を仮定として受け取っている)。 -/
def squareOfBaseU.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印との 1-可換性(model Frobenioid の四角形の一般形)",
    sectionId := "frdi-cor-5-4" }

end Square

end ModelData

end ABC3.Found.FrdI
