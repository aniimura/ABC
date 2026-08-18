import ABC3.Found.Arakelov.PicSecUnit

/-!
# Arakelov (B1) 第 120 ブロック —— **制限と局所化の可換図式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★「生成元の制限は生成元」の中身

    M_g  --tildeAwayEquiv-->  Γ(tilde M, D(g))
     |                              |
     | 局所化(第 98)                | 制限
     v                              v
    M_{h·g} --tildeAwayEquiv--> Γ(tilde M, D(h·g))

★★この図式が**可換**であることが、局所自明性の最後の一点である。

## ★★機構 —— `IsLocalizedModule.ext`

★両辺を `mk` と合成すると、どちらも `toOpen M (D(h·g))` になる:

| 辺 | 経路 |
|---|---|
| 上→右 | `iso_mk_one` で `toOpen M (D g)`、`toOpen_res'`(第 118)で制限 |
| 左→下 | `liftAwayMap_comp`(第 98)、`iso_mk_one` |

★★`powers g` が終域に可逆に作用する(第 119)ので `IsLocalizedModule.ext` が使える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `basicOpenMul_le` | ★`D(h·g) ≤ D(g)` |
| `tildeAwayEquiv_res` | ★★★★★★**可換図式** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (g h : (R : Type u))

/-- ★**`D(h·g) ≤ D(g)`**。 -/
theorem basicOpenMul_le : PrimeSpectrum.basicOpen (h * g) ≤ PrimeSpectrum.basicOpen g := by
  rw [PrimeSpectrum.basicOpen_mul]; exact inf_le_right

/-- ★★★★★★**制限と局所化の可換図式**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「生成元の制限は生成元」の中身である。 -/
theorem tildeAwayEquiv_res :
    (((modulesSpecToSheaf.obj (tilde M)).presheaf.map
          (homOfLE (basicOpenMul_le R g h)).op).hom)
        ∘ₗ (tildeAwayEquiv R M g).toLinearMap
      = (tildeAwayEquiv R M (h * g)).toLinearMap ∘ₗ (liftAwayMap R g h M) := by
  refine IsLocalizedModule.ext (Submonoid.powers g)
    (LocalizedModule.mkLinearMap (Submonoid.powers g) M)
    (isUnit_end_powers_section R M g h) ?_
  refine LinearMap.ext fun m => ?_
  show ((modulesSpecToSheaf.obj (tilde M)).presheaf.map _).hom
      ((tildeAwayEquiv R M g) (LocalizedModule.mk m 1)) = _
  rw [tildeAwayEquiv, IsLocalizedModule.iso_mk_one]
  show _ = (tildeAwayEquiv R M (h * g)) (liftAwayMap R g h M (LocalizedModule.mk m 1))
  rw [show liftAwayMap R g h M (LocalizedModule.mk m 1) = LocalizedModule.mk m 1 from
    DFunLike.congr_fun (liftAwayMap_comp R g h M) m,
    tildeAwayEquiv, IsLocalizedModule.iso_mk_one]
  exact DFunLike.congr_fun (congrArg ModuleCat.Hom.hom
    (toOpen_res' R M _ _ (homOfLE (basicOpenMul_le R g h)))) m

/-! ## ★出典の紐付け(`.src`) -/

def tildeAwayEquiv_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限と局所化の可換図式)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
