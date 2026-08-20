import ABC3.Found.GaloisRep.PadicLinear
import ABC3.Interface.GaloisRep.Representation

/-!
# Galois (G3) 第 79 ブロック —— **★★★★★★★★`ρ_{E,l} : Gal(L/K) → GL₂(ℤ_l)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★原文が名指しする表現が書けた

    Gal(L/K) → GL₂(ℤ_l),   σ ↦ (T_l E への σ の作用の行列)

★作用そのものは mathlib の `Affine.Point.map`(関手性つき)で書ける。
★★行列になるのは第 78 ブロック(加法写像は自動的に `ℤ_l` 線型)による。
★★★可逆性は `σ` と `σ⁻¹` の行列の積が `1` になることから出る。

## ★★★残るのは `det = 円分指標`

G3 の obligation はもう一つ **`det ρ(σ)` が 1 の `l` 冪根への作用の指数である**
ことを要求する(§9-403 で訂正した形)。★これは **Weil 対**の内容であり、
mathlib に在庫が無い(2026-08-17 実測)。★★したがって G3 はまだ閉じない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `galPoint_mul` / `galPoint_one` | ★作用が群作用であること |
| `galTate_mul` / `galTate_one` | ★★`T_l E` の上でも同じ |
| `galMat` ほか | ★★★行列の割り当てと乗法性 |
| `galRep` | ★★★★**`Gal(L/K) →* GL₂(ℤ_l)`** |
| `exists_galRep` | ★★★★★★★★**作用の行列表示の存在**(G3 の `rep_action`) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-! ## ★点への作用は群作用である -/

theorem galPoint_mul (W : WeierstrassCurve K) (σ τ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point) :
    galPoint W (σ * τ) P = galPoint W σ (galPoint W τ P) := by
  rw [galPoint, galPoint, galPoint, Point.map_map]
  rfl

theorem galPoint_one (W : WeierstrassCurve K) (P : (W.baseChange L).toAffine.Point) :
    galPoint W (1 : L ≃ₐ[K] L) P = P := by
  cases P with
  | zero => exact Point.map_zero _
  | some h =>
    rw [galPoint, Point.map_some]
    rfl

theorem galTate_mul (W : WeierstrassCurve K) (l : ℕ) (σ τ : L ≃ₐ[K] L) :
    galTate W l (σ * τ) = (galTate W l σ).comp (galTate W l τ) := by
  refine AddMonoidHom.ext (fun x => Subtype.ext (funext fun n => Subtype.ext ?_))
  exact galPoint_mul W σ τ _

theorem galTate_one (W : WeierstrassCurve K) (l : ℕ) :
    galTate W l (1 : L ≃ₐ[K] L) = AddMonoidHom.id (tateModule (W.baseChange L) l) := by
  refine AddMonoidHom.ext (fun x => Subtype.ext (funext fun n => Subtype.ext ?_))
  exact galPoint_one W _

/-! ## ★★★行列表示 -/

/-- ★★`σ` の `T_l E` への作用の行列。 -/
noncomputable def galMat (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L) :
    Matrix (Fin 2) (Fin 2) ℤ_[l] := matOf e (galTate W l σ)

theorem galMat_apply (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (x : tateModule (W.baseChange L) l) :
    e (galTate W l σ x) = Matrix.mulVec (galMat W l e σ) (e x) :=
  matOf_apply e (galTate W l σ) x

theorem galMat_mul (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ τ : L ≃ₐ[K] L) :
    galMat W l e (σ * τ) = galMat W l e σ * galMat W l e τ := by
  rw [galMat, galMat, galMat, galTate_mul]
  exact matOf_comp e _ _

theorem galMat_one (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :
    galMat W l e (1 : L ≃ₐ[K] L) = 1 := by
  rw [galMat, galTate_one]
  exact matOf_id e

/-- ★★★★**Galois 表現** `ρ : Gal(L/K) → GL₂(ℤ_l)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
noncomputable def galRep (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :
    (L ≃ₐ[K] L) →* Matrix.GeneralLinearGroup (Fin 2) ℤ_[l] where
  toFun := fun σ => ⟨galMat W l e σ, galMat W l e σ⁻¹,
    by rw [← galMat_mul, mul_inv_cancel, galMat_one],
    by rw [← galMat_mul, inv_mul_cancel, galMat_one]⟩
  map_one' := Units.ext (galMat_one W l e)
  map_mul' := fun σ τ => Units.ext (galMat_mul W l e σ τ)

theorem galRep_coe (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L) :
    ((galRep W l e σ : Matrix.GeneralLinearGroup (Fin 2) ℤ_[l]) :
      Matrix (Fin 2) (Fin 2) ℤ_[l]) = galMat W l e σ := rfl

/-- ★★★★★★★★**Galois 作用は `GL₂(ℤ_l)` の行列で書ける**——G3 の `rep_action`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_galRep [IsAlgClosed L] [CharZero K] (W : WeierstrassCurve K) (hell : W.IsElliptic)
    (l : ℕ) [Fact l.Prime] :
    ∃ (ρ : (L ≃ₐ[K] L) →* Matrix.GeneralLinearGroup (Fin 2) ℤ_[l])
      (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])),
      ∀ (σ : L ≃ₐ[K] L) (x : tateModule (W.baseChange L) l),
        e (galTate W l σ x)
          = Matrix.mulVec ((ρ σ : Matrix.GeneralLinearGroup (Fin 2) ℤ_[l]) :
              Matrix (Fin 2) (Fin 2) ℤ_[l]) (e x) := by
  haveI : CharZero L := charZero_of_injective_algebraMap (algebraMap K L).injective
  haveI hE : (W.baseChange L).IsElliptic := by
    haveI := hell
    unfold WeierstrassCurve.baseChange
    infer_instance
  have hΔ : (W.baseChange L).Δ ≠ 0 := (W.baseChange L).isUnit_Δ.ne_zero
  obtain ⟨e0⟩ := torsion_padic_addEquiv (W.baseChange L) hΔ l
  refine ⟨galRep W l (e0.trans (prodEquivPiTwo l)), e0.trans (prodEquivPiTwo l), ?_⟩
  intro σ x
  rw [galRep_coe]
  exact galMat_apply W l _ σ x

/-! ## ★出典の紐付け(`.src`) -/

def exists_galRep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進 Galois 表現 rho_{E,l} の構成)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
