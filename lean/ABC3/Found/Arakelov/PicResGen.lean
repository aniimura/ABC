import ABC3.Found.Arakelov.PicLiftGen
import ABC3.Found.Arakelov.PicSectionEquiv

/-!
# Arakelov (B1) 第 131 ブロック —— **切断側での「制限した生成元」の全単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 130 を切断の言葉へ運ぶ

第 130 ブロックは **`M_{t·g}` の言葉**で

    Bijective (c : R_{t·g} ↦ c • liftAwayMapA (e.symm 1))

を与えた。★★本ブロックはこれを **`Γ(tilde M, D(t·g))` の言葉**へ運ぶ。

## ★★運び方は 3 段の合成である

    𝒪(D(t·g)) --awayRingEquiv--> R_{t·g} --(· • y)--> M_{t·g}
      --tildeAwayEquivScalar--> Γ(tilde M, D(t·g))

★★★真ん中が第 130、両端は**全単射**(環同型・線型同型)だから合成も全単射である。

★可換性は `tildeAwayEquivScalar` が `𝒪(D(t·g))` 線型であること(第 125)から
`map_smul` **一発**で出る——`𝒪` の作用は `awayRingEquiv` 経由で定義されているので `rfl`。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `bijective_smul_restrictGen` | ★★★★★★**切断側での制限した生成元の全単射性** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (g t : (R : Type u))

/-- ★★★★★★**切断側での「制限した生成元」による乗法は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが §9-145 の残る `sorry` の**切断側の形**である。 -/
theorem bijective_smul_restrictGen
    (e : LocalizedModule (Submonoid.powers g) M
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    Function.Bijective (fun c : (Γ(Spec R, PrimeSpectrum.basicOpen (t * g)) : Type u) =>
      c • (tildeAwayEquivScalar R M (t * g)
        (liftAwayMapA (R : Type u) g t (M : Type u) (e.symm 1)))) := by
  have h1 := bijective_smul_liftGen (R : Type u) g t (M : Type u) e
  have hcomp : (fun c : (Γ(Spec R, PrimeSpectrum.basicOpen (t * g)) : Type u) =>
      c • (tildeAwayEquivScalar R M (t * g)
        (liftAwayMapA (R : Type u) g t (M : Type u) (e.symm 1))))
      = (tildeAwayEquivScalar R M (t * g))
        ∘ (fun d : Localization (Submonoid.powers (t * g)) =>
            d • liftAwayMapA (R : Type u) g t (M : Type u) (e.symm 1))
        ∘ (awayRingEquiv R (t * g)) := by
    funext c
    show c • (tildeAwayEquivScalar R M (t * g)) _ = _
    rw [← map_smul]
    rfl
  rw [hcomp]
  exact ((tildeAwayEquivScalar R M (t * g)).bijective.comp h1).comp
    (awayRingEquiv R (t * g)).bijective

/-! ## ★出典の紐付け(`.src`) -/

def bijective_smul_restrictGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断側での制限した生成元の全単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
