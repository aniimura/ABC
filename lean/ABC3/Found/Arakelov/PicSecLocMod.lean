import ABC3.Found.Arakelov.PicPointwise

/-!
# Arakelov (B1) 第 129 ブロック —— **切断に `Localization` 加群構造を入れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★最後の `sorry` へ向けた橋

§9-145 で組み立て全体が型検査を通り、残るのは
「**制限した生成元が生成元である**」ことの同定だけになった。

★★それを第 112(局所化は全単射を保つ)で言うには、
**切断が `Localization (powers f)` 加群でなければならない**。

★★★`Γ(tilde M, D f)` は `𝒪(D f)` 加群であり、
`𝒪(D f) ≅ Localization (powers f)`(第 124 の `awayRingEquiv`)なので
`Module.compHom` で入る。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `locModOnSection` | ★★切断は `Localization (powers f)` 加群 |
| `locTowerOnSection` | ★★`R → Localization → 切断` は足場 |

## ★★★これは第 124 と同じ手の**逆向き**である

第 124 は `M_f` に `𝒪(D f)` 構造を入れた。本ブロックは切断に `Localization` 構造を入れる。
★どちらも `awayRingEquiv` を `Module.compHom` に通すだけである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★★**切断は `Localization (powers f)` 加群である**。 -/
noncomputable instance locModOnSection : Module (Localization (Submonoid.powers f))
    ((Γ(tilde M, PrimeSpectrum.basicOpen f) : Type u)) :=
  Module.compHom _ (awayRingEquiv R f).symm.toRingHom

/-- ★★**`R → Localization (powers f) → 切断` は足場をなす**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 112 ブロック(局所化は全単射を保つ)が切断の側で使える。 -/
instance locTowerOnSection : IsScalarTower (R : Type u) (Localization (Submonoid.powers f))
    ((Γ(tilde M, PrimeSpectrum.basicOpen f) : Type u)) := by
  refine IsScalarTower.of_algebraMap_smul fun r x => ?_
  show ((awayRingEquiv R f).symm
    (algebraMap (R : Type u) (Localization (Submonoid.powers f)) r)) • x = r • x
  rw [AlgEquiv.commutes, IsScalarTower.algebraMap_smul]

/-! ## ★出典の紐付け(`.src`) -/

def locTowerOnSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断に Localization 加群構造を入れる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
