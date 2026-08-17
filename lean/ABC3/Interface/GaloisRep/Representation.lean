import ABC3.Interface.GaloisRep.Torsion
import Mathlib.LinearAlgebra.Matrix.GeneralLinearGroup.Defs
import Mathlib.FieldTheory.Galois.Basic
import Mathlib.NumberTheory.Padics.RingHoms

/-!
# Galois 表現のスケルトン(2/3)—— **表現そのものと像**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★本ファイルが受ける 3 件

| # | 受けるもの | mathlib の在庫(2026-08-17 実測) |
|---|---|---|
| G3 | `ρ_{E,l} : Gal(K̄/K) → GL₂(ℤ_l)` | ★**無い**((G2) に従属) |
| G4 | `mod l` 表現 `Gal(K̄/K) → GL₂(𝔽_l)` | ★**無い**((G1) に従属) |
| G5 | 像が `SL₂` を含むこと | ★**無い**。原文 `Theorem 3.8` の結論そのもの |

★★★**G5 が `[GenEll] Theorem 3.8` の主張である**——
`Lemma 3.1`(`SL₂` の群論、我々は**実装済み**)は
その `mod l` 版を出すための道具である。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve

/-! ## ★★G3 —— `l` 進表現 -/

/-- **(G3)** `l` 進 Galois 表現 `ρ_{E,l} : Gal(K̄/K) → GL₂(ℤ_l)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★**原文が名指ししているのはこれである。** -/
structure GaloisRepData where
  /-- 台となる Tate 加群。 -/
  toTateModuleData : TateModuleData
  /-- `ρ_{E,l}`。 -/
  rep : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    (W : WeierstrassCurve K) → (l : ℕ) → [Fact l.Prime] → (L ≃ₐ[K] L) →
    Matrix.GeneralLinearGroup (Fin 2) ℤ_[l]
  /-- ★群準同型であること。 -/
  rep_mul : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ τ : L ≃ₐ[K] L),
    rep W l (σ * τ) = rep W l σ * rep W l τ
  /-- ★★**行列式は円分指標**(Weil 対から出る)。 -/
  det_eq_cyclotomic : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    (W : WeierstrassCurve K) → (l : ℕ) → [Fact l.Prime] → ((L ≃ₐ[K] L) → ℤ_[l]ˣ)

def GaloisRepData.waiting : WaitingFor :=
  { what := "(G3) l 進 Galois 表現 rho_{E,l} : Gal(K̄/K) → GL_2(Z_l) と、その行列式が円分指標であること"
    trackB := "Found/GaloisRep — ★(G2) の Tate 加群に従属する。★★mathlib は `Matrix.GeneralLinearGroup` と `AlgEquiv`(Galois 群)を持つので、**行き先と定義域は書ける**——無いのは Tate 加群と、その上の Galois 作用である。★★★行列式が円分指標であることには **Weil 対**が要る(mathlib に 0 件、2026-08-17 実測)" }

/-! ## ★★G4 —— `mod l` 表現 -/

/-- **(G4)** `mod l` 表現 `Gal(K̄/K) → GL₂(𝔽_l)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★`Lemma 3.1`(`SL₂` の群論)が効くのはこちらである。 -/
structure ModLRepData where
  /-- 台となる `l` 進表現。 -/
  toGaloisRepData : GaloisRepData
  /-- `mod l` 表現。 -/
  repMod : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    (W : WeierstrassCurve K) → (l : ℕ) → [Fact l.Prime] → (L ≃ₐ[K] L) →
    Matrix.GeneralLinearGroup (Fin 2) (ZMod l)
  /-- ★群準同型であること。 -/
  repMod_mul : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ τ : L ≃ₐ[K] L),
    repMod W l (σ * τ) = repMod W l σ * repMod W l τ
  /-- ★★`l` 進表現の還元であること。 -/
  repMod_eq_reduction : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ : L ≃ₐ[K] L),
    (repMod W l σ : Matrix (Fin 2) (Fin 2) (ZMod l))
      = (toGaloisRepData.rep W l σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).map
          (fun a => (PadicInt.toZMod a : ZMod l))

def ModLRepData.waiting : WaitingFor :=
  { what := "(G4) mod l 表現 Gal(K̄/K) → GL_2(F_l) と、それが l 進表現の還元であること"
    trackB := "Found/GaloisRep — ★(G1)(G3) に従属する。★★`PadicInt.toZMod` は mathlib にあるので還元写像は書ける。★★★`Lemma 3.1`(SL_2 の群論)は **我々が Found/GenEll/Lemma31.lean に実装済み**(sorry 0)——効くのはこの mod l の側である" }

/-! ## ★★★G5 —— 像が `SL₂` を含むこと(`Theorem 3.8` の結論) -/

/-- **(G5)** Galois 表現の像が `SL₂(ℤ_l)` を含むこと。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★**これが `[GenEll] Theorem 3.8` の結論そのものである。** -/
structure FullImageData where
  /-- 台となる `mod l` 表現。 -/
  toModLRepData : ModLRepData
  /-- 像が `SL₂` を含むという述語。 -/
  ImageContainsSL2 : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    WeierstrassCurve K → ℕ → Prop
  /-- ★★**述語の内容**——`SL₂(ℤ_l)` のどの元も像に入ること。 -/
  imageContainsSL2_iff : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime],
    @ImageContainsSL2 K L _ _ _ _ W l ↔
      ∀ A : Matrix.GeneralLinearGroup (Fin 2) ℤ_[l],
        (A : Matrix (Fin 2) (Fin 2) ℤ_[l]).det = 1 →
        ∃ σ : L ≃ₐ[K] L, toModLRepData.toGaloisRepData.rep W l σ = A

def FullImageData.waiting : WaitingFor :=
  { what := "(G5) Galois 表現の像が SL_2(Z_l) を含むこと —— [GenEll] Theorem 3.8 の結論そのもの"
    trackB := "Found/GaloisRep — ★(G3)(G4) に従属する。★★原文の証明は `Lemma 3.1`(**我々は実装済み**)・`Lemma 3.2`(Tate 曲線、= G6)・`Lemma 3.5` / `Lemma 3.7`(Faltings 高さ、= G8)を使う。★★★したがって S3 の最後の段であり、G1-G4 と G6-G8 がすべて揃ってから着手する" }

/-! ## ★出典の紐付け(`.src`) -/

def GaloisRepData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(l 進表現の構成のみ——像の主張は含まない)",
    sectionId := "genell-thm-3-8" }

def ModLRepData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(mod l 表現のみ——像の主張は含まない)",
    sectionId := "genell-thm-3-8" }

def FullImageData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(像が SL₂ を含むという結論の定式化のみ)",
    sectionId := "genell-thm-3-8" }

end ABC3.Interface.GaloisRep
