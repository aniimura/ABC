import ABC3.Found.GaloisRep.OmegaNum

/-!
# Galois (G1) 第 2 ブロック —— **`ω_n` の分子は底変換と可換**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★普遍環から降ろすための鍵

`ω_n := omegaNum / 2` を定義するには `2 ∣ omegaNum` が要る。
mathlib の docstring が示す筋は:

    標数 0 の普遍環 ℤ[A₁,…,A₆][X,Y] で 2 ∣ omegaNum を示し、
    普遍射 ℤ[A₁,…,A₆] → R(Aᵢ ↦ aᵢ)で降ろす

★★この「降ろす」段に要るのが**本ブロック**である:

    omegaNum (W.map f) n = (omegaNum W n).map (mapRingHom f)

## ★★部品はすべて mathlib に在った

| 部品 | 名前 |
|---|---|
| `complEDS₂` の写像可換性 | ★`map_complEDS₂` |
| `ψ` の写像可換性 | ★`WeierstrassCurve.map_ψ` |
| `φ` の写像可換性 | ★`WeierstrassCurve.map_φ` |
| `ψ₂`/`Ψ₃`/`preΨ₄` の写像可換性 | ★`map_ψ₂` / `map_Ψ₃` / `map_preΨ₄` |

★★★**`ω_n` 以外はすべて写像可換性が整備されている**——
mathlib は `ω_n` を TODO にしたが、**その周りは完成している**。

## ★★摩擦は「二重多項式環の綴り」だった

EDS の係数環は `R[X][Y] = (R[X])[X]` で、写像は `mapRingHom (mapRingHom f)` である。
★`map_complEDS₂` の暗黙引数 `b c d` は section 変数なので、
**位置引数で渡す**必要があった(`map_complEDS₂ b c d f n` の順)。
★★`Polynomial.map g` と `(mapRingHom g)` が定義的に等しいことを `show` で橋渡しする。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `map_psiComp` | ★★★`psiComp` は底変換と可換 |
| `map_omegaNum` | ★★★★**`ω_n` の分子は底変換と可換** |
-/

universe u v

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} {S : Type v} [CommRing R] [CommRing S] (W : WeierstrassCurve R) (f : R →+* S)

/-- ★★★**`psiComp` は底変換と可換である**。 -/
theorem map_psiComp (n : ℤ) :
    ABC3.Found.GaloisRep.psiComp (W.map f) n
      = (ABC3.Found.GaloisRep.psiComp W n).map (mapRingHom f) := by
  have h : (mapRingHom (mapRingHom f)) (complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) n)
      = complEDS₂ ((mapRingHom (mapRingHom f)) W.ψ₂)
          ((mapRingHom (mapRingHom f)) (C W.Ψ₃))
          ((mapRingHom (mapRingHom f)) (C W.preΨ₄)) n :=
    map_complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) (mapRingHom (mapRingHom f)) n
  rw [ABC3.Found.GaloisRep.psiComp, ABC3.Found.GaloisRep.psiComp]
  show complEDS₂ _ _ _ n = (mapRingHom (mapRingHom f)) (complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) n)
  rw [h]
  congr 1
  · exact WeierstrassCurve.map_ψ₂ W f
  · simp [WeierstrassCurve.map_Ψ₃]
  · simp [WeierstrassCurve.map_preΨ₄]


/-- ★★★★**`ω_n` の分子は底変換と可換である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが「普遍環で示して降ろす」の**降ろす**段である。 -/
theorem map_omegaNum (n : ℤ) :
    omegaNum (W.map f) n = (omegaNum W n).map (mapRingHom f) := by
  rw [omegaNum, omegaNum, map_psiComp W f n, WeierstrassCurve.map_ψ W f n,
    WeierstrassCurve.map_φ W f n]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def map_omegaNum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——ω_n の分子は底変換と可換)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
