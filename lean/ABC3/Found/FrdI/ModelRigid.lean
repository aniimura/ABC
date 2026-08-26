/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52
import ABC3.Found.FrdI.Prop113

/-!
# model Frobenioid への関手の rigidity —— 底が rigid なら十分

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.104。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of

## ★★`Corollary 5.4` の rigidity をどう出すか(測った結果)

`Corollary 4.11, (i)` の rigidity(`Cor411Rigid.lean`)は 3 段だった:
Div-identity → base-identity(Div-slim)→ 恒等(unit-trivial)。
★★**その第 1 段は `istrToUnTr` が充満であることに依っている** ——
`η` の成分が `𝒞^un-tr` の**すべての対象**で取れるからである。

★★★`Corollary 5.4` の縦の矢印 `𝒞 ⥤ 𝒞^rlf` は**充満でない**(係数拡大だから)。
そこで**別の道**を通る —— **model Frobenioid では射が 4 成分の明示的な組**なので、
自己同型の 4 成分を直接つぶせばよい:

| 成分 | 消える理由 |
|---|---|
| `base` | 底への合成 `F ⋙ (𝒞^rlf → 𝒟)` が rigid |
| `deg` | 同型なので `ℕ≥1` で `a·b = 1` |
| `div` | 同型なので `x + n·div = 0`、`Φ` が **sharp** |
| `u` | 上の 3 つを `cond` に入れると `Div_B(u) = 0`、`Div_B` が**単射** |

★★**`Div_B` の単射性は model Frobenioid `𝒞^rlf` では無料である** ——
`Div_B` は `ℝ·Φ^birat ⊆ (Φ^rlf)^gp` の**包含そのもの**だから。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ModelData.deg_eq_one_of_isIso` | 同型射は `deg = 1` |
| `ModelData.div_eq_zero_of_isIso` | 同型射は `div = 0`(`Φ` が sharp) |
| `ModelData.gpMapOn_id_apply` | `gpMapOn (𝟙 d)` は恒等 |
| `ModelData.hom_eq_id_of_base_id` | ★★4 成分をつぶす |
| `ModelData.isRigidFunctor_of_proj` | ★★★★**底が rigid なら関手も rigid** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u3 v3

namespace ModelData

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

/-! ## ★1. 同型射の 3 成分 -/

/-- ★**model Frobenioid の同型射は `deg = 1`** —— `ℕ≥1` は単元を持たない。 -/
theorem deg_eq_one_of_isIso {A B : Obj M} (φ : A ⟶ B) [IsIso φ] : φ.deg = 1 := by
  have h : ((inv φ).deg * φ.deg : ℕ+) = 1 :=
    congrArg (fun t : A ⟶ A => t.deg) (IsIso.hom_inv_id φ)
  have h' : ((inv φ).deg : ℕ) * (φ.deg : ℕ) = 1 := congrArg (fun t : ℕ+ => (t : ℕ)) h
  exact PNat.coe_injective (Nat.eq_one_of_mul_eq_one_left h')

/-- ★**model Frobenioid の同型射は `div = 0`** —— `Φ` が sharp だから。

★`x + n·div = 0` で `n ≥ 1` なので `div` は加法的単元。 -/
theorem div_eq_zero_of_isIso (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) [IsIso φ] :
    φ.div = 0 := by
  have h0 : M.phi.map φ.base (inv φ).div + (((inv φ).deg : ℕ+) : ℕ) • φ.div = 0 :=
    congrArg (fun t : A ⟶ A => t.div) (IsIso.hom_inv_id φ)
  set n : ℕ := (((inv φ).deg : ℕ+) : ℕ) with hn
  have hn1 : n - 1 + 1 = n := Nat.sub_add_cancel (inv φ).deg.2
  have hsplit : n • φ.div = φ.div + (n - 1) • φ.div := by
    conv_lhs => rw [← hn1]
    rw [succ_nsmul']
  rw [hsplit] at h0
  refine (h.divisorial _).2 φ.div (isAddUnit_iff_exists_neg.mpr
    ⟨(n - 1) • φ.div + M.phi.map φ.base (inv φ).div, ?_⟩)
  rw [← add_assoc]
  rw [add_comm (M.phi.map φ.base (inv φ).div)] at h0
  exact h0

/-! ## ★2. `gpMapOn (𝟙 d)` は恒等 -/

@[simp] theorem gpMapOn_id_apply (Φ : MonoidOn.{v, u, w} D) (d : D) (x : Gp (Φ.val d)) :
    Φ.gpMapOn (𝟙 d) x = x := by
  have h : Φ.map (𝟙 d) = AddMonoidHom.id (Φ.val d) := by
    refine AddMonoidHom.ext fun a => ?_
    exact Φ.map_id d a
  show gpMap _ (Φ.map (𝟙 d)) x = x
  rw [h, gpMap_id]
  rfl

/-! ## ★3. 4 成分をつぶす -/

/-- ★★★**底が恒等な自己同型は恒等** —— `deg`・`div` は同型性から、
`u` は `cond` と `Div_B` の単射性から。 -/
theorem hom_eq_id_of_base_id (h : Hyp M) (hinj : ∀ d : D, Function.Injective (M.divB d))
    {A : Obj M} (φ : A ⟶ A) [IsIso φ] (hb : φ.base = 𝟙 A.base) : φ = 𝟙 A := by
  have hdeg : φ.deg = 1 := deg_eq_one_of_isIso φ
  have hdiv : φ.div = 0 := div_eq_zero_of_isIso h φ
  have hu : φ.u = 0 := by
    have hc := φ.cond
    rw [hdeg, hdiv, hb] at hc
    have hc' : A.cls = A.cls + M.divB _ φ.u := by
      simpa using hc
    refine hinj A.base ?_
    rw [map_zero]
    have hz : A.cls + 0 = A.cls + M.divB _ φ.u := by rw [add_zero]; exact hc'
    exact (add_left_cancel hz).symm
  refine Hom.ext ?_ ?_ ?_ ?_
  · exact hb
  · exact hdiv
  · exact hdeg
  · exact hu

/-! ## ★4. 底が rigid なら関手も rigid -/

/-- ★★★★★★**model Frobenioid への関手は、`𝒟` への合成が rigid なら rigid**。

★★これが `Corollary 5.4` の rigidity の実体である ——
`Corollary 4.11, (i)` の議論(充満性に依る)は使えないが、
**model Frobenioid の射が 4 成分の明示的な組であること**で代わりになる。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of -/
theorem isRigidFunctor_of_proj {X : Type u3} [Category.{v3} X] (h : Hyp M)
    (hinj : ∀ d : D, Function.Injective (M.divB d)) (F : X ⥤ Obj M)
    (hrig : IsRigidFunctor (F ⋙ (modelPre h).proj)) : IsRigidFunctor F := by
  intro η
  haveI hiso : ∀ V : X, IsIso (η.hom.app V) := fun V =>
    ⟨η.inv.app V, η.hom_inv_id_app V, η.inv_hom_id_app V⟩
  have hb : ∀ V : X, (η.hom.app V).base = 𝟙 ((F.obj V).base) := by
    intro V
    have hr := hrig (NatIso.ofComponents
      (fun W => ((modelPre h).proj).mapIso (asIso (η.hom.app W)))
      (fun {W W'} f => by
        show ((modelPre h).proj).map (F.map f) ≫ ((modelPre h).proj).map (η.hom.app W')
          = ((modelPre h).proj).map (η.hom.app W) ≫ ((modelPre h).proj).map (F.map f)
        rw [← Functor.map_comp, ← Functor.map_comp, η.hom.naturality f]))
    have hV := congrArg (fun t : (F ⋙ (modelPre h).proj) ≅ (F ⋙ (modelPre h).proj) =>
      t.hom.app V) hr
    exact hV
  refine Iso.ext (NatTrans.ext (funext fun V => ?_))
  exact hom_eq_id_of_base_id h hinj (η.hom.app V) (hb V)

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Corollary 5.4` の rigidity の実体
(★**条つき**: `Corollary 5.4` 全体はまだ揃っていない)。 -/
def isRigidFunctor_of_proj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 合成関手の rigidity(model Frobenioid の 4 成分でつぶす)",
    sectionId := "frdi-cor-5-4" }

end ModelData

end ABC3.Found.FrdI
