/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Frob

/-!
# model Frobenioid の `𝒪^×(A)` は `{u ∈ B(A_𝒟) | Div_B(u) = 0}`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.100。

原文 (FrdI p.100):
> refer to as the projection to D to

## ★★何のためか

`Example 6.1` の `𝒪^×(A) = 𝒪^▷(A) = k_L^×`(鎖 `normalize` の `ex61-units`)を
出すには、まず **model Frobenioid 一般**で `𝒪^×(A)` が何かを確定させておく必要がある。

★model Frobenioid の射は `(Base, Div, deg_Fr, u)` の 4 つ組なので、
`𝒪^×(A)` の条件(`Base = 𝟙`、`deg_Fr = 1`、可逆)を入れると

| 成分 | 値 |
|---|---|
| `Base` | `𝟙` (`IsBaseIdentity`) |
| `deg_Fr` | `1` (`IsLinear`) |
| `Div` | `0` (可逆性 —— `model_isIso_iff`) |
| `u` | **自由**、ただし条件式から `Div_B(u) = 0` |

★★すなわち `𝒪^×(A) ≃ {u ∈ B(A_𝒟) | Div_B(u) = 0}` である。

## ★★★合成が加法になること

`End A` の積は `φ * ψ = ψ ≫ φ` で、`Hom.comp ψ φ` の `u` 成分は
`B(ψ.base)^*(φ.u) + φ.deg • ψ.u`。`ψ.base = 𝟙`、`φ.deg = 1` なので
**`φ.u + ψ.u`** になる。★だから乗法群としての同型が出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory ABC3.Meta

universe v u w

/-! ## ★0. 恒等射での `gpMapOn` -/

/-- ★**`gpMapOn` は恒等射で恒等**。 -/
theorem gpMapOn_id {D : Type u} [Category.{v} D] (Φ : MonoidOn.{v, u, w} D) (A : D)
    (x : Gp (Φ.val A)) : Φ.gpMapOn (𝟙 A) x = x := by
  have hmap : Φ.map (𝟙 A) = AddMonoidHom.id _ := AddMonoidHom.ext (MonoidOn.map_id Φ A)
  show (gpMap (Φ.val A) (Φ.map (𝟙 A))) x = x
  rw [hmap, gpMap_id]
  rfl

namespace ModelData

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

/-! ## ★1. `Div_B` の核 -/

variable (M) in
/-- ★★**`𝒪^×(A)` に対応する `B(A_𝒟)` の部分単系** —— `Div_B` の核。 -/
def unitSub (A : Obj M) : AddSubmonoid (M.bmon.val A.base) where
  carrier := {u | M.divB A.base u = 0}
  add_mem' := by
    intro x y hx hy
    show M.divB A.base (x + y) = 0
    rw [map_add, hx, hy, add_zero]
  zero_mem' := by show M.divB A.base 0 = 0; simp

variable (M) in
@[simp] theorem mem_unitSub {A : Obj M} {x : M.bmon.val A.base} :
    x ∈ unitSub M A ↔ M.divB A.base x = 0 := Iff.rfl

/-! ## ★2. 両向きの構成 -/

/-- ★`Div_B(u) = 0` から作る自己射。 -/
def endOfU {A : Obj M} (u : unitSub M A) : End A where
  base := 𝟙 A.base
  div := 0
  deg := 1
  u := (u : M.bmon.val A.base)
  cond := by
    have hu : M.divB A.base (u : M.bmon.val A.base) = 0 := u.2
    rw [hu, gpMapOn_id]
    simp

/-- ★`endOfU` は `𝒪^×(A)` に入る。 -/
theorem endOfU_mem (h : M.Hyp) {A : Obj M} (u : unitSub M A) :
    endOfU u ∈ OTimes (modelPre h) A := by
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · show (endOfU u).base = (𝟙 A : End A).base
    rfl
  · show (modelPre h).degFr (endOfU u) = 1
    rfl
  · refine (isUnit_iff_isIso _).mpr ((model_isIso_iff h _).mpr ⟨?_, rfl, rfl⟩)
    show IsIso (𝟙 A.base)
    infer_instance

/-- ★★`𝒪^×(A)` の元は `Div_B(u) = 0` を満たす。

★可逆性から `Div = 0`(`model_isIso_iff`)、`Base = 𝟙`・`deg_Fr = 1` と合わせて
射の条件式が `A.cls = A.cls + Div_B(u)` になり、`Gp` は群なので消せる。 -/
theorem divB_eq_zero_of_mem_otimes (h : M.Hyp) {A : Obj M}
    (φ : OTimes (modelPre h) A) : M.divB A.base (φ : End A).u = 0 := by
  obtain ⟨⟨hb, hlin⟩, hu⟩ := φ.2
  have hbase : (φ : End A).base = 𝟙 A.base := hb
  have hdeg : (φ : End A).deg = 1 := hlin
  have hiso : IsIso ((φ : End A) : A ⟶ A) := (isUnit_iff_isIso _).mp hu
  have hdiv : (φ : End A).div = 0 := ((model_isIso_iff h _).mp hiso).2.1
  have hc := (φ : End A).cond
  rw [hbase, hdiv, hdeg] at hc
  simp only [PNat.one_coe, one_smul, map_zero, add_zero] at hc
  rw [gpMapOn_id] at hc
  exact add_eq_left.mp hc.symm

/-! ## ★3. 同型 -/

/-- ★★★★★★**model Frobenioid の `𝒪^×(A)`** ——
`𝒪^×(A) ≃ {u ∈ B(A_𝒟) | Div_B(u) = 0}`(乗法群としての同型)。

★これが `Example 6.1` の `𝒪^×(A) = 𝒪^▷(A) = k_L^×` の出発点である。 -/
noncomputable def otimesModelEquiv (h : M.Hyp) (A : Obj M) :
    OTimes (modelPre h) A ≃* Multiplicative (unitSub M A) where
  toFun φ := Multiplicative.ofAdd ⟨(φ : End A).u, divB_eq_zero_of_mem_otimes h φ⟩
  invFun u := ⟨endOfU (Multiplicative.toAdd u), endOfU_mem h _⟩
  left_inv := by
    intro φ
    obtain ⟨⟨hb, hlin⟩, hu⟩ := φ.2
    have hbase : (φ : End A).base = 𝟙 A.base := hb
    have hdeg : (φ : End A).deg = 1 := hlin
    have hiso : IsIso ((φ : End A) : A ⟶ A) := (isUnit_iff_isIso _).mp hu
    have hdiv : (φ : End A).div = 0 := ((model_isIso_iff h _).mp hiso).2.1
    exact Subtype.ext (Hom.ext hbase.symm hdiv.symm hdeg.symm rfl)
  right_inv := by intro u; rfl
  map_mul' := by
    intro φ ψ
    obtain ⟨⟨hb, hlin⟩, hu⟩ := φ.2
    obtain ⟨⟨hb2, hlin2⟩, hu2⟩ := ψ.2
    have hbase2 : (ψ : End A).base = 𝟙 A.base := hb2
    have hdeg : (φ : End A).deg = 1 := hlin
    refine congrArg Multiplicative.ofAdd (Subtype.ext ?_)
    show ((φ : End A) * (ψ : End A)).u = (φ : End A).u + (ψ : End A).u
    show M.bmon.map (ψ : End A).base (φ : End A).u
        + ((φ : End A).deg : ℕ) • (ψ : End A).u = _
    rw [hbase2, hdeg]
    simp [MonoidOn.map_id]

@[simp] theorem otimesModelEquiv_apply (h : M.Hyp) (A : Obj M) (φ : OTimes (modelPre h) A) :
    ((Multiplicative.toAdd (otimesModelEquiv h A φ) : unitSub M A) : M.bmon.val A.base)
      = (φ : End A).u := rfl

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def otimesModelEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 100,
    item := "Theorem 5.2 — model Frobenioid の 𝒪^×(A) は Div_B の核",
    sectionId := "frdi-thm-5-2" }

def otimesModelEquiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "model_isIso_iff(model Frobenioid の同型の判定)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.model_isIso_iff") 100,
    .citation "[ABC3]" "OTimes(𝒪^× の定義)"
      (.inProject "ABC3" "ABC3.Found.FrdI.OTimes") 100,
    .derivation
      "Base = 𝟙・deg_Fr = 1・Div = 0 を射の条件式に入れると Div_B(u) = 0 が残る。合成の u 成分は φ.u + ψ.u" 100 ]

def unitSub.src : Source :=
  { paper := "FrdI", pdfPage := 100,
    item := "Theorem 5.2 — Div_B の核",
    sectionId := "frdi-thm-5-2" }

end ModelData

end ABC3.Found.FrdI
