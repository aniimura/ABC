import ABC3.Found.Arakelov.PicLocSurj
import ABC3.Found.Arakelov.PicFiniteCover

/-!
# Arakelov (B1) 第 143 ブロック —— ★★★★★★★★**局所自明 ⟹ `F ≅ (Γ F)~`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★mathlib に無かった段が埋まった

`Tilde.lean` 547 行の TODO——「準連接 ⟹ `fromTildeΓ` 同型」——のうち、
★**局所自明(局所自由)な場合**を自前で埋めた。

    IsLocallyTrivial (Spec R) F.val  ⟹  IsLocalizing (modulesSpecToSheaf.obj F)
                                     ⟹  IsIso F.fromTildeΓ
                                     ⟹  tilde (Γ F) ≅ F

## ★★組み上げ(第 134–143、10 ブロック)

| # | 内容 |
|---|---|
| 134 | 自明化を基本開集合へ落とす |
| 135 | 有限被覆(`span = ⊤`) |
| 136 | 構造層の制限は局所化(mathlib の `isLocalization_of_eq_basicOpen`) |
| 137 | 自明化の四角形で `F` 側へ(`naturality` 一発) |
| 138 | `f` 倍は単射 |
| 139 | `f` 倍は全射(層で貼る) |
| 140 | `map_units` と `exists_of_eq` |
| 141 | 制限写像の短縮記法 |
| 142 | ★`surj'`(2 段の最大値) |
| 143 | ★★組み立て |

## ★★★これで `equivPicRing` の**逆向き**の土台が出来た

残るは「`Γ F` が可逆加群であること」と群同型の組み立てである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)

/-- ★★★★★★★**局所自明なら `IsLocalizing`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★mathlib の `Tilde.lean` が TODO に置いた段である。 -/
theorem isLocalizing_of_isLocallyTrivial (h : IsLocallyTrivial (Spec R) F.val) :
    AlgebraicGeometry.IsLocalizing (modulesSpecToSheaf.obj F) := by
  obtain ⟨n, g, hspan, htriv⟩ := exists_finite_basicOpen_trivial R F.val h
  intro f
  exact
    { map_units := fun s => isUnit_pow_smul_specD R F g hspan htriv f s
      surj := fun y => by
        obtain ⟨a, N, hN⟩ := exists_pow_smul_eq_res R F g hspan htriv f y
        exact ⟨⟨a, ⟨f ^ N, N, rfl⟩⟩, hN⟩
      exists_of_eq := fun {x₁ x₂} hx => by
        obtain ⟨N, hN⟩ := exists_pow_smul_eq R F g hspan htriv f x₁ x₂ hx
        exact ⟨⟨f ^ N, N, rfl⟩, hN⟩ }

theorem isIso_fromTildeΓ_of_isLocallyTrivial (h : IsLocallyTrivial (Spec R) F.val) :
    IsIso F.fromTildeΓ :=
  (AlgebraicGeometry.isIso_fromTildeΓ_iff_isLocalizing F).2
    (isLocalizing_of_isLocallyTrivial R F h)

noncomputable def tildeGammaIsoOfTrivial (h : IsLocallyTrivial (Spec R) F.val) :
    tilde (AlgebraicGeometry.moduleSpecΓFunctor.obj F) ≅ F :=
  @asIso _ _ _ _ F.fromTildeΓ (isIso_fromTildeΓ_of_isLocallyTrivial R F h)


/-! ## ★出典の紐付け(`.src`) -/

def isLocalizing_of_isLocallyTrivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所自明 ⟹ F ≅ (Γ F)~)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
