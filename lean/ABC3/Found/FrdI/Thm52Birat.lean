/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Frob
import ABC3.Found.FrdI.Prop44Otri

/-!
# model Frobenioid の `𝒞^birat` —— `𝒪^▷(A^birat)` は可換

★`Proposition 4.4, (ii)` に残る 1 条 `otriBase` は
**`𝒪^▷(A^birat)` が可換**であることと同値だった(`Gap/FrdI/Prop44.lean`)。

★★**model Frobenioid では Ore の四角形を明示的に作れる**ので、可換性が直接出る:

`x = [a]⁻¹ ≫ [f]`・`y = [a]⁻¹ ≫ [g]`(共通の添字、`f`・`g`・`a` は次数 1・同じ底)に取り、
**同じ `p := preStepDown E (Div a)`** で 2 つの Ore の四角形

  `p ≫ f = m ≫ a`(`Div m = Div f`)、`p ≫ g = m' ≫ a`(`Div m' = Div g`)

を作ると、`𝒞` の中で **`m ≫ g = m' ≫ f`**(4 成分すべて一致)が成り立つ。
あとは `[a]`・`[p]` が同型なので割り算するだけである。

★★したがって **model Frobenioid については `Proposition 4.4, (ii)` が閉じる**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

namespace ModelData

/-- ★★**Ore の 2 本の脚は交換する** —— 同じ `p` で `f`・`g` に四角形を当てたときの
もう一方の脚 `m`・`m'` について `m ≫ g = m' ≫ f`。
★4 成分すべてが `Φ`・`B` の可換性で潰れる。 -/
theorem ore_legs_comm (h : Hyp M) {E A : Obj M} {K : M.phi.val E.base}
    {a f g : E ⟶ A} {m m' : objDown E K ⟶ E}
    (hfd : f.deg = 1) (hgd : g.deg = 1) (hfb : f.base = a.base) (hgb : g.base = a.base)
    (hmb : m.base = 𝟙 E.base) (hmd : m.div = f.div) (hmn : m.deg = 1)
    (hmu : m.u = f.u + bneg h a.u)
    (hm'b : m'.base = 𝟙 E.base) (hm'd : m'.div = g.div) (hm'n : m'.deg = 1)
    (hm'u : m'.u = g.u + bneg h a.u) :
    m ≫ g = m' ≫ f := by
  refine Hom.ext ?_ ?_ ?_ ?_
  · show m.base ≫ g.base = m'.base ≫ f.base
    rw [hmb, hm'b, hgb, hfb]
  · have h1 : M.phi.map m.base g.div = g.div := by
      rw [hmb]; exact MonoidOn.map_id M.phi E.base g.div
    have h2 : M.phi.map m'.base f.div = f.div := by
      rw [hm'b]; exact MonoidOn.map_id M.phi E.base f.div
    show M.phi.map m.base g.div + (g.deg : ℕ) • m.div
      = M.phi.map m'.base f.div + (f.deg : ℕ) • m'.div
    rw [h1, h2, hgd, hfd, hmd, hm'd]
    simp only [PNat.one_coe, one_smul]
    exact add_comm _ _
  · show g.deg * m.deg = f.deg * m'.deg
    rw [hgd, hfd, hmn, hm'n]
  · have h3 : M.bmon.map m.base g.u = g.u := by
      rw [hmb]; exact MonoidOn.map_id M.bmon E.base g.u
    have h4 : M.bmon.map m'.base f.u = f.u := by
      rw [hm'b]; exact MonoidOn.map_id M.bmon E.base f.u
    show M.bmon.map m.base g.u + (g.deg : ℕ) • m.u
      = M.bmon.map m'.base f.u + (f.deg : ℕ) • m'.u
    rw [h3, h4, hgd, hfd, hmu, hm'u]
    simp only [PNat.one_coe, one_smul]
    abel

/-! ## ★2. `𝒪^▷(A^birat)` の可換性 -/

/-- ★★★**model Frobenioid では `𝒪^▷(A^birat)` は可換**。

★★`x = [a]⁻¹ ≫ [f]`・`y = [a]⁻¹ ≫ [g]` に取り、**同じ `p`** で 2 つの Ore の四角形を作ると
`m ≫ g = m' ≫ f`(`ore_legs_comm`)。あとは `[a]`・`[p]` で割るだけ。 -/
theorem model_birat_otri_comm (h : Hyp M)
    (X : BiratCat (modelPre h) (model_frobenioid h))
    (x y : OTri (biratPre (modelPre h) (model_frobenioid h)) X) : x * y = y * x := by
  refine Subtype.ext ?_
  show ((y : End X) : X ⟶ X) ≫ ((x : End X) : X ⟶ X)
    = ((x : End X) : X ⟶ X) ≫ ((y : End X) : X ⟶ X)
  -- ★共通の代表元
  obtain ⟨W, f, g, hfx, hgy⟩ :=
    HomBirat.exists_rep_pair (P := modelPre h) (G := model_frobenioid h)
      (x : End X) (y : End X)
  have hac : IsCoAngular (modelPre h) W.unop.hom.hom := W.unop.hom.property.1
  have has : IsPreStep (modelPre h) W.unop.hom.hom := W.unop.hom.property.2
  have had : W.unop.hom.hom.deg = 1 := has.1
  haveI hab : IsIso W.unop.hom.hom.base := has.2
  -- ★次数
  have hfd : f.deg = 1 := by
    have h1 : biratDeg (HomBirat.mk W f) = 1 := by rw [hfx]; exact x.2.2
    rwa [biratDeg_mk] at h1
  have hgd : g.deg = 1 := by
    have h1 : biratDeg (HomBirat.mk W g) = 1 := by rw [hgy]; exact y.2.2
    rwa [biratDeg_mk] at h1
  -- ★底
  have hbase : ∀ (z : End X) (hz : IsBaseIdentity (biratPre (modelPre h) (model_frobenioid h)) z)
      (φ : W.unop.left.obj ⟶ _), HomBirat.mk W φ = z → φ.base = W.unop.hom.hom.base := by
    intro z hz φ hφ
    have h1 : biratBase (HomBirat.mk W φ) = 𝟙 _ := by
      rw [hφ]
      have := hz
      show (biratPre (modelPre h) (model_frobenioid h)).Base z = 𝟙 _
      rw [this, (biratPre (modelPre h) (model_frobenioid h)).Base_id]
    rw [biratBase_mk, sliceBaseOf_eq] at h1
    have h2 := congrArg (fun t => (modelPre h).Base W.unop.hom.hom ≫ t) h1
    simp only [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at h2
    exact h2.trans (Category.comp_id _)
  have hfb : f.base = W.unop.hom.hom.base := hbase _ x.2.1 f hfx
  have hgb : g.base = W.unop.hom.hom.base := hbase _ y.2.1 g hgy
  -- ★2 つの Ore の四角形(同じ `p`)
  obtain ⟨m, hmb, hmd, hmn, hmu, horef⟩ :=
    exists_ore_square h W.unop.hom.hom f had hfb W.unop.hom.hom.div f.div (by
      rw [hfd]; simp only [PNat.one_coe, one_smul]; exact add_comm _ _)
  obtain ⟨m', hm'b, hm'd, hm'n, hm'u, horeg⟩ :=
    exists_ore_square h W.unop.hom.hom g had hgb W.unop.hom.hom.div g.div (by
      rw [hgd]; simp only [PNat.one_coe, one_smul]; exact add_comm _ _)
  have hlegs : m ≫ g = m' ≫ f :=
    ore_legs_comm h hfd hgd hfb hgb hmb hmd (hmn.trans hfd) hmu
      hm'b hm'd (hm'n.trans hgd) hm'u
  -- ★`[a]`・`[p]` は同型
  haveI hia : IsIso ((toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom) :=
    birat_isIso_of_coaPre W.unop.hom.hom hac has
  haveI hip : IsIso ((toBiratCat (modelPre h) (model_frobenioid h)).map
      (preStepDown W.unop.left.obj W.unop.hom.hom.div)) :=
    birat_isIso_of_coaPre _ (model_coAngular h _) (preStepDown_isPreStep h _ _)
  -- ★`[a] ≫ x = [f]`、`[a] ≫ y = [g]`
  have hax : (toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
      ≫ ((x : End X) : X ⟶ X)
      = (toBiratCat (modelPre h) (model_frobenioid h)).map f := by
    rw [← hfx]
    exact birat_toHom_comp_mk W.unop.hom.hom hac has f
  have hay : (toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
      ≫ ((y : End X) : X ⟶ X)
      = (toBiratCat (modelPre h) (model_frobenioid h)).map g := by
    rw [← hgy]
    exact birat_toHom_comp_mk W.unop.hom.hom hac has g
  -- ★割り算
  refine (cancel_epi ((toBiratCat (modelPre h) (model_frobenioid h)).map
    W.unop.hom.hom)).mp ?_
  have e1 : (toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
      ≫ (((y : End X) : X ⟶ X) ≫ ((x : End X) : X ⟶ X))
      = (toBiratCat (modelPre h) (model_frobenioid h)).map g ≫ ((x : End X) : X ⟶ X) :=
    ((Category.assoc _ _ _).symm).trans (congrArg (fun t => t ≫ ((x : End X) : X ⟶ X)) hay)
  have e2 : (toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
      ≫ (((x : End X) : X ⟶ X) ≫ ((y : End X) : X ⟶ X))
      = (toBiratCat (modelPre h) (model_frobenioid h)).map f ≫ ((y : End X) : X ⟶ X) :=
    ((Category.assoc _ _ _).symm).trans (congrArg (fun t => t ≫ ((y : End X) : X ⟶ X)) hax)
  refine e1.trans (Eq.trans ?_ e2.symm)
  refine (cancel_epi ((toBiratCat (modelPre h) (model_frobenioid h)).map
    (preStepDown W.unop.left.obj W.unop.hom.hom.div))).mp ?_
  have e3 : (toBiratCat (modelPre h) (model_frobenioid h)).map
        (preStepDown W.unop.left.obj W.unop.hom.hom.div)
      ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map g ≫ ((x : End X) : X ⟶ X))
      = (toBiratCat (modelPre h) (model_frobenioid h)).map (m' ≫ f) :=
    calc (toBiratCat (modelPre h) (model_frobenioid h)).map
          (preStepDown W.unop.left.obj W.unop.hom.hom.div)
        ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map g ≫ ((x : End X) : X ⟶ X))
        = ((toBiratCat (modelPre h) (model_frobenioid h)).map
            (preStepDown W.unop.left.obj W.unop.hom.hom.div)
          ≫ (toBiratCat (modelPre h) (model_frobenioid h)).map g)
            ≫ ((x : End X) : X ⟶ X) := (Category.assoc _ _ _).symm
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map
            (preStepDown W.unop.left.obj W.unop.hom.hom.div ≫ g)
            ≫ ((x : End X) : X ⟶ X) := by rw [Functor.map_comp]
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map (m' ≫ W.unop.hom.hom)
            ≫ ((x : End X) : X ⟶ X) :=
          congrArg (fun t => (toBiratCat (modelPre h) (model_frobenioid h)).map t
            ≫ ((x : End X) : X ⟶ X)) horeg
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map m'
            ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
              ≫ ((x : End X) : X ⟶ X)) := by
          rw [Functor.map_comp, Category.assoc]
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map m'
            ≫ (toBiratCat (modelPre h) (model_frobenioid h)).map f :=
          congrArg (fun t => (toBiratCat (modelPre h) (model_frobenioid h)).map m' ≫ t) hax
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map (m' ≫ f) :=
          (Functor.map_comp _ _ _).symm
  have e4 : (toBiratCat (modelPre h) (model_frobenioid h)).map
        (preStepDown W.unop.left.obj W.unop.hom.hom.div)
      ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map f ≫ ((y : End X) : X ⟶ X))
      = (toBiratCat (modelPre h) (model_frobenioid h)).map (m ≫ g) :=
    calc (toBiratCat (modelPre h) (model_frobenioid h)).map
          (preStepDown W.unop.left.obj W.unop.hom.hom.div)
        ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map f ≫ ((y : End X) : X ⟶ X))
        = ((toBiratCat (modelPre h) (model_frobenioid h)).map
            (preStepDown W.unop.left.obj W.unop.hom.hom.div)
          ≫ (toBiratCat (modelPre h) (model_frobenioid h)).map f)
            ≫ ((y : End X) : X ⟶ X) := (Category.assoc _ _ _).symm
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map
            (preStepDown W.unop.left.obj W.unop.hom.hom.div ≫ f)
            ≫ ((y : End X) : X ⟶ X) := by rw [Functor.map_comp]
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map (m ≫ W.unop.hom.hom)
            ≫ ((y : End X) : X ⟶ X) :=
          congrArg (fun t => (toBiratCat (modelPre h) (model_frobenioid h)).map t
            ≫ ((y : End X) : X ⟶ X)) horef
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map m
            ≫ ((toBiratCat (modelPre h) (model_frobenioid h)).map W.unop.hom.hom
              ≫ ((y : End X) : X ⟶ X)) := by
          rw [Functor.map_comp, Category.assoc]
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map m
            ≫ (toBiratCat (modelPre h) (model_frobenioid h)).map g :=
          congrArg (fun t => (toBiratCat (modelPre h) (model_frobenioid h)).map m ≫ t) hay
      _ = (toBiratCat (modelPre h) (model_frobenioid h)).map (m ≫ g) :=
          (Functor.map_comp _ _ _).symm
  exact e3.trans (Eq.trans
    (congrArg (fun t => (toBiratCat (modelPre h) (model_frobenioid h)).map t) hlegs.symm)
    e4.symm)

/-! ## ★3. したがって model については `Proposition 4.4, (ii)` が閉じる -/

/-- ★★★★★**model Frobenioid の `𝒞^birat` は `Definition 1.3` の 21 条をすべて満たす**。

★★`Proposition 4.4, (ii)` に残っていた 1 条 `otriBase` が、
model については `𝒪^▷(A^birat)` の可換性(`model_birat_otri_comm`)から出る。 -/
theorem model_birat_frobenioidCore (h : Hyp M) :
    FrobenioidCore (biratPre (modelPre h) (model_frobenioid h)) :=
  birat_frobenioidCore_of_comm (modelPre h) (model_frobenioid h) (model_birat_otri_comm h)

/-- ★**locator** —— `Proposition 4.4, (ii)`(★**条つき**: model Frobenioid について)。 -/
def model_birat_frobenioidCore.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85,
    item := "Proposition 4.4, (ii) — Definition 1.3 の 21 条(model Frobenioid について)",
    sectionId := "frdi-prop-4-4" }

end ModelData

end ABC3.Found.FrdI
