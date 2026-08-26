import ABC3.Found.GaloisRep.WeilNondegFull
import ABC3.Found.GaloisRep.BasisFree
import ABC3.Found.GaloisRep.ModLWitness

/-!
# Galois (G5) 第 328 ブロック —— **★★★★★★★★★★G5 達成**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★到達点

> **`FullImageData` の witness が組めた**(`FullImageData.nonvacuous`)——**(G5) 達成**

★★★これで Galois 表現論の 8 件のうち **7 件**(G1-G7)が埋まる。

## ★★★★★★配線の中身

    第 327 `weilPairing_nondegenerate`(非退化性)
      → 第 210 `det_cyclotomic_of_nondeg`(`BasisFree.lean`)
      → `det_cyclotomic_full`(本ブロック)
      → `FullImageData.det_cyclotomic`

★`IsDedekindDomain F[W]` は第 137 の `isDedekindDomain_coordinateRing` で供給する
(`IsAlgClosed L` と `IsUnit (2:L)` から)。
★★`Infinite L` は `CharZero L` から、`CharZero L` は `CharZero K` と
`algebraMap K L` の単射性から出る。

## ★★★★`ImageContainsSL2` の型が `Fact l.Prime` を持たない件

界面の `ImageContainsSL2` は `... → WeierstrassCurve K → ℕ → Prop` であって
`[Fact l.Prime]` を取らない。★しかし `ℤ_[l]` は `Fact l.Prime` が無いと環ですらない。
★★そこで述語の本体を **`∀ hp : Nat.Prime l, ...`** の形にし、
中身は `Fact` 付きの補助述語 `ImageSL2Aux` に流した。
★★★`imageContainsSL2_iff` は `Fact` が場にあるときの同値で、
`⟨Fact.out⟩ ≡ inst`(構造の η と証明の非関与)で両向きとも通る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ImageSL2Aux` | ★像が `SL₂` を含むという述語の中身 |
| `isElliptic_baseChange_affine` | ★係数拡大は楕円性を保つ |
| `det_cyclotomic_full` | ★★★★★★★★★★**`det ρ = 円分指標`(仮定なし)** |
| `fullImageDataWitness` | ★★★★★★★★★★**`FullImageData` の実装** |
| `FullImageData.nonvacuous` | ★★★★★★★★★★**G5 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Interface.GaloisRep

/-! ## ★述語の中身 -/

/-- ★像が `SL₂` を含むという述語の中身(`Fact l.Prime` つき)。 -/
def ImageSL2Aux {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] : Prop :=
  ∀ A : Matrix.GeneralLinearGroup (Fin 2) ℤ_[l],
    (A : Matrix (Fin 2) (Fin 2) ℤ_[l]).det = 1 →
    ∃ σ : L ≃ₐ[K] L, modLRepDataWitness.toGaloisRepData.rep W l σ = A

/-- ★係数拡大は楕円性を保つ。 -/
theorem isElliptic_baseChange_affine {K L : Type} [Field K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (hell : W.IsElliptic) :
    WeierstrassCurve.IsElliptic ((W.baseChange L).toAffine) := by
  haveI := hell
  refine ⟨?_⟩
  show IsUnit (W.baseChange L).Δ
  rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_Δ]
  exact hell.isUnit.map (algebraMap K L)

/-! ## ★★★★★★★★★★`det ρ = 円分指標` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`det ρ(σ)` は円分指標である**——仮定なしの形。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 210 の `det_cyclotomic_of_nondeg` に第 327 の非退化性を渡すだけ。 -/
theorem det_cyclotomic_full {K L : Type} [Field K] [DecidableEq K] [CharZero K]
    [Field L] [DecidableEq L] [Algebra K L] [IsAlgClosed L]
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic (W.baseChange L).toAffine]
    (l : ℕ) [Fact l.Prime]
    (e : (ABC3.Interface.GaloisRep.tateModule (W.baseChange L) l) ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (n : ℕ) (ζ : L) (hζ : ζ ^ (l ^ n) = 1) :
    σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val) := by
  haveI : CharZero L := charZero_of_injective_algebraMap (algebraMap K L).injective
  haveI : Infinite L := Infinite.of_injective (Nat.cast : ℕ → L) Nat.cast_injective
  have h2 : IsUnit (2 : L) := isUnit_iff_ne_zero.2 two_ne_zero
  have h4 : (4 : L) ≠ 0 := by norm_num
  haveI : IsDedekindDomain ((W.baseChange L).toAffine).CoordinateRing :=
    isDedekindDomain_coordinateRing h2 _
  refine det_cyclotomic_of_nondeg W l e σ n ζ hζ ?_
  intro S hS hall
  exact weilPairing_nondegenerate _ h2 h4 (l ^ n)
    (Nat.one_le_pow _ _ (Fact.out : l.Prime).pos) S hS hall

/-! ## ★★★★★★★★★★G5 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`FullImageData` の実装**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
noncomputable def fullImageDataWitness : FullImageData where
  toModLRepData := modLRepDataWitness
  ImageContainsSL2 := fun {K L} _ _ _ _ _ W l =>
    ∀ hp : Nat.Prime l, @ImageSL2Aux K L _ _ _ _ _ W l ⟨hp⟩
  imageContainsSL2_iff := by
    intro K L _ _ _ _ _ W l hp
    haveI := hp
    exact ⟨fun h A hA => h Fact.out A hA, fun h _ A hA => h A hA⟩
  det_cyclotomic := by
    intro K L _ _ _ _ _ _ _ W hell l hp σ n ζ hζ
    haveI := hp
    haveI : WeierstrassCurve.IsElliptic ((W.baseChange L).toAffine) :=
      isElliptic_baseChange_affine W hell
    have hne : Nonempty (tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) :=
      nonempty_tate_basis W hell l
    have hrep : repChoice W l = galRep W l (Classical.choice hne) := repChoice_eq W l hne
    show σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((repChoice W l σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val)
    rw [hrep]
    exact det_cyclotomic_full W l _ σ n ζ hζ

/-- ★★★★★★★★★★**`FullImageData` は非空虚である**——G5 達成。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが Galois 表現論の 8 件のうち **G5** である。 -/
theorem FullImageData.nonvacuous : Nonempty FullImageData :=
  ⟨fullImageDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def det_cyclotomic_full.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式が円分指標であること)",
    sectionId := "genell-thm-3-8" }

def fullImageDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(像が SL₂ を含むという結論の定式化のみ)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
