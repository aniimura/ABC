import ABC3.Found.GaloisRep.GalRepWitness

/-!
# Galois (G4) 第 93 ブロック —— **★★★★★★★★G4 達成 `mod l` 表現**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★`mod l` 表現は還元そのものである

    GL₂(ℤ_l) → GL₂(𝔽_l)

★mathlib の `PadicInt.toZMod` を行列に広げ、単元群へ持ち上げるだけである。
★★`Lemma 3.1`(`SL₂` の群論、**我々は実装済み**)が効くのはこちら側である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `redGL` | ★`GL₂(ℤ_l) → GL₂(𝔽_l)` |
| `modLRepDataWitness` | ★`ModLRepData` の実装 |
| `ModLRepData.nonvacuous` | ★★★★★★★★**G4 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Interface.GaloisRep

/-- ★`GL₂(ℤ_l) → GL₂(𝔽_l)` —— 係数の還元。 -/
noncomputable def redGL (l : ℕ) [Fact l.Prime] :
    Matrix.GeneralLinearGroup (Fin 2) ℤ_[l] →* Matrix.GeneralLinearGroup (Fin 2) (ZMod l) :=
  Units.map ((PadicInt.toZMod (p := l)).mapMatrix (m := Fin 2)).toMonoidHom

/-- ★★★★★★★★**`ModLRepData` の実装**。 -/
noncomputable def modLRepDataWitness : ModLRepData where
  toGaloisRepData := galoisRepDataWitness
  repMod := by
    intro K L _ _ _ _ _ W l hp σ
    haveI := hp
    exact redGL l (galoisRepDataWitness.rep W l σ)
  repMod_mul := by
    intro K L _ _ _ _ _ W l hp σ τ
    haveI := hp
    show redGL l (galoisRepDataWitness.rep W l (σ * τ))
      = redGL l (galoisRepDataWitness.rep W l σ) * redGL l (galoisRepDataWitness.rep W l τ)
    rw [galoisRepDataWitness.rep_mul W l σ τ, map_mul]
  repMod_eq_reduction := by
    intro K L _ _ _ _ _ W l hp σ
    haveI := hp
    rfl

/-- ★★★★★★★★**`ModLRepData` は非空虚である**——G4 達成。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが Galois 表現論の 8 件のうち **G4** である。 -/
theorem ModLRepData.nonvacuous : Nonempty ModLRepData :=
  ⟨modLRepDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def modLRepDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 表現の構成)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
