/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop48Cpt
import ABC3.Found.FrdI.Prop114
import ABC3.Found.FrdI.Prop44Gp

/-!
# [FrdI] Proposition 4.8, (iii), (b) —— `Frobenius-compact` の第 3 条

## ★★★★★★以前の見立ての誤り

台帳は「`σ` が恒等 ⟹ `Div^gp` の核に入る ⟹ **核が torsion であること**が要る」と
書いていたが、これは**遠回り**であった(`frobeniusCompact_cond3_of_eq` を見よ)。

★仮定は `∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^(d*k) = u^(c*k)` であり、
結論は `∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^k = u^k` である。
★★**`c = d` なら仮定の式そのものが結論**である。
★★★したがって `𝒪^×(A)` の torsion 性は**まったく関係が無い**。

## ★残る問い —— `c = d`

`Div^gp` を当てると `(d*k) • T(x) = (c*k) • x`(`T := Φ^gp(Base θ.inv)`、`x := Div^gp u`)。
★`T = id` が言えれば `(d*k) • x = (c*k) • x` となり、
`x` が torsion でなければ `c = d` が出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. `non-dilating` から「恒等」を出す —— 一般のモノイド版 -/

section MonoidGen

variable {M : Type w} [AddCommMonoid M]

/-- ★★★★**`isDivIdentity_of_forall_mprec` の一般形** ——
`Frobenioid` を経由しない、素の加法モノイド準同型についての主張。

★元の証明と同じ 2 段: `charMap` を `non-dilating` で恒等にし、
`toChar` の単射性(sharp)で `M` へ戻す。 -/
theorem addMonoidHom_eq_id_of_forall_mprec (hsharp : IsSharp M) (f : M →+ M)
    (hnd : IsNonDilating f) (h : ∀ x : M, x ≠ 0 → MPrec (f x) x) :
    f = AddMonoidHom.id M := by
  have hchar : charMap f = AddMonoidHom.id (MChar M) := by
    refine hnd (fun a ha => ?_)
    obtain ⟨x, rfl⟩ := toChar_surjective _ a
    have hx : x ≠ 0 := fun hz => ha.1 (by rw [hz, map_zero])
    rw [charMap_toChar]
    exact MPrec.map toChar (h x hx)
  refine AddMonoidHom.ext (fun x => ?_)
  refine toChar_injective_of_isSharp hsharp ?_
  have := congrArg (fun g : MChar M →+ MChar M => g (toChar x)) hchar
  rw [charMap_toChar] at this
  exact this

def addMonoidHom_eq_id_of_forall_mprec.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 41,
    item := "Proposition 1.14, (v) — non-dilating から恒等（一般形）",
    sectionId := "frdi-prop-1-14" }

end MonoidGen

/-! ## ★2. `Div^gp` と共役 -/

section BiratConj

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★`𝒞^birat` では同型の Frobenius 次数は 1。 -/
theorem biratDeg_of_isIso {A B : BiratCat P G} (f : A ⟶ B) [IsIso f] :
    ((biratDeg f : ℕ+) : ℕ) = 1 := by
  have h : (biratPre P G).degFr f = 1 := degFr_of_isIso (biratPre P G) f
  rw [biratPre_degFr] at h
  rw [h]
  rfl

/-- ★★★★★**共役の `Div^gp`** —— `Div^gp (endConj θ u) = Φ^gp(Base θ.inv) (Div^gp u)`。

★`θ.inv ≫ θ.hom = 𝟙` から `Div^gp θ.inv` が消える。 -/
theorem biratDivGp_endConj {A : BiratCat P G} (θ : A ≅ A) {u : End A}
    (hu : u ∈ OTimes (biratPre P G) A) :
    biratDivGp ((endConj θ u : End A) : A ⟶ A)
      = gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A))) := by
  -- ★`u ≫ θ.hom`
  have h1 : biratDivGp ((u : A ⟶ A) ≫ θ.hom)
      = biratDivGp θ.hom + biratDivGp ((u : A ⟶ A)) := by
    have h := biratDivGp_comp' ((u : A ⟶ A)) θ.hom
    rw [gpMap_biratBase_of_baseIdentity hu.1.1, biratDeg_of_isIso θ.hom, one_smul] at h
    exact h
  -- ★`θ.inv ≫ (u ≫ θ.hom)`
  have h2 : biratDivGp (θ.inv ≫ ((u : A ⟶ A) ≫ θ.hom))
      = gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A) ≫ θ.hom))
        + biratDivGp θ.inv := by
    have h := biratDivGp_comp' θ.inv ((u : A ⟶ A) ≫ θ.hom)
    have hdeg : ((biratDeg ((u : A ⟶ A) ≫ θ.hom) : ℕ+) : ℕ) = 1 := by
      have hiso : IsIso ((u : A ⟶ A) ≫ θ.hom) := by
        haveI : IsIso (u : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso u).mp hu.2
        infer_instance
      exact biratDeg_of_isIso _
    rw [hdeg, one_smul] at h
    exact h
  -- ★`𝟙 = θ.inv ≫ θ.hom` から `Div^gp θ.inv` を消す
  have h3 : gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp θ.hom) + biratDivGp θ.inv = 0 := by
    have h := biratDivGp_comp' θ.inv θ.hom
    rw [biratDeg_of_isIso θ.hom, one_smul, θ.inv_hom_id, biratDivGp_id] at h
    exact h.symm
  show biratDivGp (θ.inv ≫ ((u : A ⟶ A) ≫ θ.hom)) = _
  rw [h2, h1, map_add]
  rw [add_assoc, add_comm (gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A)))),
    ← add_assoc, h3, zero_add]

def biratDivGp_endConj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 共役の Div^gp",
    sectionId := "frdi-prop-4-4" }

end BiratConj

end ABC3.Found.FrdI
