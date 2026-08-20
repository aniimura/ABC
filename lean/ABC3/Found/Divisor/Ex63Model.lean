/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Bmon

/-!
# `Example 6.3` の model Frobenioid `C_{F̄/F}`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★最後の 1 本 —— `div : B(L) → Φ(L)^gp`

`Ex63Datum.lean` が因子の側(`Φ`)、`Ex63Bmon.lean` が有理関数の側(`B = L^×`)を
与えたので、残るのは 2 つを繋ぐ `div` とその関手性である。

★`arithDiv : L^× → arithDivGroup L`(`ArithDivisor.lean`)を
`effRGpEquiv`(`ArithGp.lean`)の逆で `Φ(L)^gp` へ運ぶ。

★★関手性の中身は **`arithExtend_arithDiv`**(`ArithFunctor.lean`)——
「体の射で送ってから因子を取る」と「因子を取ってから引き戻す」が一致する。
★それを `Gp` の層へ持ち上げるのが `effRGpHom_gpMapOn` である
(`gp_eq_sub` で差に分解し、`gpLift_toGp` で降ろすだけ)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `arithDivGroupHom` | `Additive (L^×) →+ arithDivGroup L` |
| `effRGpHom_gpMapOn` | ★`Φ^gp ≃ Γ` は引き戻しと両立する |
| `divBGalois` / `divBGalois_nat` | ★★`div : B(L) → Φ(L)^gp` とその関手性 |
| `ex63Frobenioid` | ★★★★★★**`Example 6.3` の model Frobenioid `C_{F̄/F}`** |
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI

open scoped NumberField

/-! ## ★1. `div : L^× → Γ` -/

/-- ★**`arithDiv` を加法的な準同型として書いたもの**。 -/
noncomputable def arithDivGroupHom (L : Type*) [Field L] [NumberField L] :
    Additive Lˣ →+ arithDivGroup L :=
  AddMonoidHom.mk' (fun x => ⟨arithDiv (Additive.toMul x), arithDiv_mem_arithDivGroup _⟩)
    (fun _ _ => Subtype.ext (arithDiv_mul _ _))

@[simp] theorem arithDivGroupHom_coe (L : Type*) [Field L] [NumberField L] (x : Additive Lˣ) :
    ((arithDivGroupHom L x : arithDivGroup L) : ArithPlace L →₀ ℝ)
      = arithDiv (Additive.toMul x) := rfl

/-! ## ★2. `Φ^gp ≃ Γ` は引き戻しと両立する -/

universe v u w

variable {D : Type u} [Category.{v} D] (Δ : ArithDatum.{v, u, w} D) (hD : IsOfFSMType D)

/-- ★`Φ^gp → Γ`(束ねた型の上で)。 -/
noncomputable def phiGpHom {A : D} : Gp ((Δ.phi hD).val A) →+ Δ.grp A :=
  effRGpHom (Δ.grp A)

@[simp] theorem phiGpHom_toGp {A : D} (a : effR (Δ.grp A)) :
    ((phiGpHom Δ hD (toGp _ a) : Δ.grp A) : Δ.primes A →₀ ℝ)
      = ((a : Δ.primes A →₀ ℝ)) := by
  show ((effRGpHom (Δ.grp A) (toGp _ a) : Δ.grp A) : Δ.primes A →₀ ℝ) = _
  rw [effRGpHom, gpLift_toGp]
  rfl

/-- ★★★**`Φ^gp ≃ Γ` は引き戻しと両立する**。

★`gp_eq_sub` で差に分解し、`gpLift_toGp` で降ろすだけ。 -/
theorem phiGpHom_gpMapOn {A B : D} (α : B ⟶ A) (z : Gp ((Δ.phi hD).val A)) :
    ((phiGpHom Δ hD ((Δ.phi hD).gpMapOn α z) : Δ.grp B) : Δ.primes B →₀ ℝ)
      = Δ.pull α ((phiGpHom Δ hD z : Δ.grp A) : Δ.primes A →₀ ℝ) := by
  obtain ⟨a, b, rfl⟩ := gp_eq_sub z
  rw [map_sub, map_sub]
  have h1 : ((Δ.phi hD).gpMapOn α) (toGp ((Δ.phi hD).val A) a)
      = toGp ((Δ.phi hD).val B) ((Δ.phi hD).map α a) := by
    exact MonoidOn.gpMapOn_toGpHom (Δ.phi hD) α a
  have h2 : ((Δ.phi hD).gpMapOn α) (toGp ((Δ.phi hD).val A) b)
      = toGp ((Δ.phi hD).val B) ((Δ.phi hD).map α b) := by
    exact MonoidOn.gpMapOn_toGpHom (Δ.phi hD) α b
  rw [h1, h2, map_sub, AddSubgroup.coe_sub, AddSubgroup.coe_sub,
    phiGpHom_toGp, phiGpHom_toGp, phiGpHom_toGp, phiGpHom_toGp, map_sub]
  rfl

theorem phiGpHom_injective {A : D} :
    Function.Injective (phiGpHom Δ hD (A := A)) :=
  effRGpHom_injective (Δ.grp A)

/-! ## ★3. `Example 6.3` の `div : B(L) → Φ(L)^gp` -/

section Ex63

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-- ★★**`div : B(L) → Φ(L)^gp`** —— `arithDiv` を `effRGpEquiv` の逆で運んだもの。 -/
noncomputable def divBGalois (A : (FinSub F Kbar)ᵒᵖ) :
    (bmonGalois F Kbar).val A →+
      Gp (((arithDatumGalois F Kbar).phi finSubOp_isOfFSMType).val A) :=
  ((effRGpEquiv (arithDivGroup A.unop.toIF)
      isGenSubgroupR_arithDivGroup).symm.toAddMonoidHom).comp
    (arithDivGroupHom A.unop.toIF)

theorem phiGpHom_divBGalois (A : (FinSub F Kbar)ᵒᵖ)
    (x : (bmonGalois F Kbar).val A) :
    Subtype.val (phiGpHom (arithDatumGalois F Kbar) finSubOp_isOfFSMType (divBGalois A x))
      = arithDiv (Additive.toMul x) := by
  show Subtype.val ((effRGpEquiv (arithDivGroup A.unop.toIF) isGenSubgroupR_arithDivGroup)
      ((effRGpEquiv (arithDivGroup A.unop.toIF) isGenSubgroupR_arithDivGroup).symm
        (arithDivGroupHom A.unop.toIF x))) = _
  rw [AddEquiv.apply_symm_apply]
  rfl

/-- ★★★**`div` は関手的** —— 中身は `arithExtend_arithDiv` 1 本。 -/
theorem divBGalois_nat {A B : (FinSub F Kbar)ᵒᵖ} (f : A ⟶ B)
    (x : (bmonGalois F Kbar).val B) :
    divBGalois A ((bmonGalois F Kbar).map f x)
      = ((arithDatumGalois F Kbar).phi finSubOp_isOfFSMType).gpMapOn f (divBGalois B x) := by
  refine phiGpHom_injective (arithDatumGalois F Kbar) finSubOp_isOfFSMType (Subtype.ext ?_)
  rw [phiGpHom_divBGalois, phiGpHom_gpMapOn, phiGpHom_divBGalois]
  letI := algOfHom f.unop
  show arithDiv (Units.map
      (algebraMap (Opposite.unop B).toIF (Opposite.unop A).toIF).toMonoidHom (Additive.toMul x))
    = arithExtend (L := (Opposite.unop B).toIF) (arithDiv (Additive.toMul x))
  exact arithExtend_arithDiv (Additive.toMul x)

/-! ## ★4. `Example 6.3` の model Frobenioid -/

/-- ★★★★★★**[FrdI] `Example 6.3`** —— `C_{F̄/F}` は **Frobenioid** である。

★`Theorem 5.2, (ii)` を、`𝒟 = B(G)⁰` の上で実際に組んだデータに当てたもの。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model] -/
noncomputable def ex63Frobenioid :
    Frobenioid (ModelData.modelPre
      ((arithDatumGalois F Kbar).modelHyp finSubOp_isOfFSMType (bmonGalois F Kbar)
        (fun A => divBGalois A) (fun f x => divBGalois_nat f x)
        (bmonGalois_isGroupLike F Kbar) finSubOp_totallyEpimorphic finSubOp_isConnected)) :=
  (arithDatumGalois F Kbar).arithFrobenioid finSubOp_isOfFSMType (bmonGalois F Kbar)
    (fun A => divBGalois A) (fun f x => divBGalois_nat f x)
    (bmonGalois_isGroupLike F Kbar) finSubOp_totallyEpimorphic finSubOp_isConnected

end Ex63

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Example 6.3` の model Frobenioid `C_{F̄/F}`。 -/
def ex63Frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
