/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRep
import ABC3.Meta.Claim

/-!
# 基底の取り替えは共役である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か（2026-08-31、第 824 の測定）

`galRep W l e` は Tate 加群の基底 `e` に依存する。★基底を取り替えると

    `galRep W l e σ = C · (galRep W l e₀ σ) · C⁻¹`

——**共役**になる（`C` は基底変換の行列）。

## ★★★なぜ要るか

葉 3（`alpha_in_modl_image`）は「mod `l` 像が `α = (1 1 / 0 1)` を含む」であるが、
★これは **Tate 一意化に適合した基底 `(ζ_l, q^{1/l})` で**成り立つ主張である。
一方 `Found/GenEll/EllModuliGalois.lean` の `imageContainsSL2J_of_alpha` は
仮説を `∀ e` の形で持っていた——★★**それでは循環する**（`α` がすべての共役に
入ることは `SL₂ ⊆ 像` を既に知っていないと言えない）。

★★★本ファイルは共役の関係を与え、`∃ e₀` から `∀ e` を出せるようにする。
`SL₂` は `GL₂` の**正規部分群**なので、共役で保たれるからである。

☆同じ形の測定: `Check/GenEll/EllModuliDegInfPos.lean`（第 745）、
`Check/GenEll/LcyclicExcTooStrong.lean`（第 754-755）、
`Check/GenEll/ImageSL2NeedsL5.lean`（第 776）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep

universe u

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★**基底変換の行列**——`e₀` から `e` へ。 -/
noncomputable def basisChange (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :
    Matrix (Fin 2) (Fin 2) ℤ_[l] :=
  Classical.choose (exists_matrix_of_addHom (e₀.symm.trans e).toAddMonoidHom)

theorem basisChange_apply (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (y : tateModule (W.baseChange L) l) :
    e y = Matrix.mulVec (basisChange W l e₀ e) (e₀ y) := by
  have h := Classical.choose_spec
    (exists_matrix_of_addHom (e₀.symm.trans e).toAddMonoidHom) (e₀ y)
  simpa only [basisChange, AddEquiv.coe_toAddMonoidHom, AddEquiv.coe_trans,
    Function.comp_apply, AddEquiv.symm_apply_apply] using h

/-- ★★★★**基底変換は可逆**——逆向きの行列との積は単位行列。 -/
theorem basisChange_mul (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :
    basisChange W l e₀ e * basisChange W l e e₀ = 1 := by
  refine matrix_ext_of_mulVec (fun x => ?_)
  have h1 := basisChange_apply W l e e₀ (e.symm x)
  have h2 := basisChange_apply W l e₀ e (e.symm x)
  simp only [AddEquiv.apply_symm_apply] at h1 h2
  rw [← Matrix.mulVec_mulVec, ← h1, ← h2, Matrix.one_mulVec]

/-- ★★★★★★★★★★★★★★★★
**基底を取り替えると `galMat` は共役になる**。 -/
theorem galMat_basisChange (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L) :
    galMat W l e σ
      = basisChange W l e₀ e * galMat W l e₀ σ * basisChange W l e e₀ := by
  refine matrix_ext_of_mulVec (fun x => ?_)
  have hy := basisChange_apply W l e e₀ (e.symm x)
  have h1 := galMat_apply W l e σ (e.symm x)
  have h2 := galMat_apply W l e₀ σ (e.symm x)
  have h3 := basisChange_apply W l e₀ e (galTate W l σ (e.symm x))
  simp only [AddEquiv.apply_symm_apply] at hy h1 h2 h3
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, ← hy, ← h2, ← h3, ← h1]

/-! ## ★★★★★★★★★★★★★★★★`GL` の水準と `SL₂` の移植 -/

open scoped MatrixGroups

/-- ★★★★★★**基底変換の `GL` 元**。 -/
noncomputable def basisChangeGL (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) : GL (Fin 2) ℤ_[l] :=
  ⟨basisChange W l e₀ e, basisChange W l e e₀,
    basisChange_mul W l e₀ e, basisChange_mul W l e e₀⟩

/-- ★★★★★★★★★★★★★★★★★★
**`galRep` は基底の取り替えで共役になる**。 -/
theorem galRep_basisChange (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L) :
    galRep W l e σ
      = basisChangeGL W l e₀ e * galRep W l e₀ σ * (basisChangeGL W l e₀ e)⁻¹ := by
  refine Units.ext ?_
  show galMat W l e σ = _
  rw [galMat_basisChange W l e₀ e σ]
  rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`SL₂ ⊆ 像` は基底に依らない**——1 つの基底で言えればすべての基底で言える。

★★★`SL₂` は `GL₂` の**正規部分群**なので、共役で保たれる。
★★これで葉 3 の仮説を `∀ e` から **`∃ e₀`** へ弱められる。 -/
theorem sl2_range_basis_transfer (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (h : ∀ g : SL(2, ℤ_[l]), (Matrix.SpecialLinearGroup.toGL g : GL (Fin 2) ℤ_[l])
      ∈ (galRep W l e₀).range)
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (g : SL(2, ℤ_[l])) :
    (Matrix.SpecialLinearGroup.toGL g : GL (Fin 2) ℤ_[l]) ∈ (galRep W l e).range := by
  classical
  set C : GL (Fin 2) ℤ_[l] := basisChangeGL W l e₀ e with hC
  set M : GL (Fin 2) ℤ_[l] := C⁻¹ * (Matrix.SpecialLinearGroup.toGL g) * C with hM
  have hdet : (M : Matrix (Fin 2) (Fin 2) ℤ_[l]).det = 1 := by
    have h1 : Matrix.GeneralLinearGroup.det M
        = Matrix.GeneralLinearGroup.det (Matrix.SpecialLinearGroup.toGL g) := by
      rw [hM, map_mul, map_mul, map_inv]
      simp [mul_comm, mul_left_comm, mul_assoc]
    have h3 : Matrix.GeneralLinearGroup.det (Matrix.SpecialLinearGroup.toGL g)
        = (1 : ℤ_[l]ˣ) := Units.ext (by simpa using g.2)
    rw [h3] at h1
    have h4 := congrArg Units.val h1
    simpa using h4
  obtain ⟨σ, hσ⟩ := h ⟨(M : Matrix (Fin 2) (Fin 2) ℤ_[l]), hdet⟩
  refine ⟨σ, ?_⟩
  have hMe : (Matrix.SpecialLinearGroup.toGL
      (⟨(M : Matrix (Fin 2) (Fin 2) ℤ_[l]), hdet⟩ : SL(2, ℤ_[l]))) = M := by
    ext i j
    rfl
  rw [galRep_basisChange W l e₀ e σ, hσ, hMe, hM, ← hC]
  group

/-! ## ★出典の紐付け(`.src`) -/

def basisChangeGL.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基底変換の GL 元)",
    sectionId := "genell-thm-3-8" }

def galRep_basisChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galRep は基底の取り替えで共役になる。★無条件)",
    sectionId := "genell-thm-3-8" }

def sl2_range_basis_transfer.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(SL₂ ⊆ 像 は基底に依らない。★無条件)",
    sectionId := "genell-thm-3-8" }

def basisChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 加群の基底変換の行列)",
    sectionId := "genell-thm-3-8" }

def basisChange_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基底変換の行列の作用。★無条件)",
    sectionId := "genell-thm-3-8" }

def basisChange_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基底変換は可逆。★無条件)",
    sectionId := "genell-thm-3-8" }

def galMat_basisChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基底を取り替えると galMat は共役になる。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
