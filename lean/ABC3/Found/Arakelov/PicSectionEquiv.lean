import ABC3.Found.Arakelov.PicAwayTransport

/-!
# Arakelov (B1) 第 127 ブロック —— **切断の同型と生成元の乗法**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★局所自明性の全単射が出た

    Γ(tilde M, D f)  ≃ₗ[𝒪(D f)]  𝒪(D f)

★2 つの同型の合成である:

| 段 | 出典 |
|---|---|
| `Γ(tilde M, D f) ≅ M_f` | ★第 125(`𝒪(D f)` 線型) |
| `M_f ≅ 𝒪(D f)` | ★第 126(`𝒪(D f)` 線型) |

★★★これに第 103(生成元の乗法は全単射)を当てると、
**`c ↦ c • (生成元の切断)` が全単射**である——これが第 115 の仮定そのものである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `sectionEquivScalar` | ★★★★★**`Γ(tilde M, D f) ≃ₗ[𝒪(D f)] 𝒪(D f)`** |
| `bijective_smul_genSection` | ★★★★★★**生成元の切断による乗法は全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★★★★★**`Γ(tilde M, D f) ≃ₗ[𝒪(D f)] 𝒪(D f)`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 125(`Γ ≅ M_f`)と第 126(`M_f ≅ 𝒪(D f)`)の合成である。 -/
noncomputable def sectionEquivScalar
    (e : LocalizedModule (Submonoid.powers f) M
      ≃ₗ[Localization (Submonoid.powers f)] Localization (Submonoid.powers f)) :
    (Γ(tilde M, PrimeSpectrum.basicOpen f) : Type u)
      ≃ₗ[(Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)]
        (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u) :=
  (tildeAwayEquivScalar R M f).symm ≪≫ₗ (awayEquivScalar R M f e)

/-- ★★★★★★**生成元の切断による乗法は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが第 115 ブロック(切断から `𝟙_ ≅ P`)の仮定そのものである。 -/
theorem bijective_smul_genSection
    (e : LocalizedModule (Submonoid.powers f) M
      ≃ₗ[Localization (Submonoid.powers f)] Localization (Submonoid.powers f)) :
    Function.Bijective (fun c : (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u) =>
      c • generatorOf (sectionEquivScalar R M f e)) :=
  bijective_smul_generator _

/-! ## ★出典の紐付け(`.src`) -/

def bijective_smul_genSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の切断による乗法は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
