import ABC3.Found.Arakelov.PicPointwise

/-!
# Arakelov (B1) 第 134 ブロック —— **自明化を基本開集合へ落とす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★逆向き(可逆層 ⟹ 可逆加群)の第 1 歩

第 133 ブロックで「可逆加群 ⟹ 可逆層」が通った。★逆向きには

    IsLocallyTrivial (Spec R) F.val  ⟹  IsLocalizing (modulesSpecToSheaf.obj F)

が要る(`isIso_fromTildeΓ_iff_isLocalizing`)。★★これは mathlib に**無い**
——`Tilde.lean` の TODO で「準連接 ⟹ `fromTildeΓ` 同型」が未着手と明記されている。

## ★★古典的な筋(Hartshorne II.5.1 / Stacks 01I3)

    Spec R = ⋃ D(gᵢ)(有限)、各 D(gᵢ) 上で F ≅ 𝒪
      ⟹ Γ(F,⊤) → Γ(F,D f) は powers f の局所化

★本ブロックはその**第 1 歩**——被覆を**基本開集合**で取り直す。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `trivialOfLe` | ★自明化は小さい開集合へ制限できる(第 57 の言い換え) |
| `exists_basicOpen_trivial` | ★★★各点に**基本開集合**の自明化近傍が取れる |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

/-- ★**自明化は小さい開集合へ制限できる**——第 57 ブロック(制限の推移律)の言い換え。 -/
noncomputable def trivialOfLe {X : Scheme.{u}} (P : X.PresheafOfModules) {V W : X.Opens}
    (h : W ≤ V) (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    (restrictPresheafFunctor X W).obj P ≅ 𝟙_ (PresheafModulesOn X W) := by
  have hiso := (restrictOnFunctor h).mapIso e
  rw [restrict_trans h P, restrictOnUnit h] at hiso
  exact hiso

/-- ★★★**各点に基本開集合の自明化近傍が取れる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Opens` の位相は点ごとなので篩から近傍が取れ、
基本開集合は基底をなすのでその中に基本開集合が入る。 -/
theorem exists_basicOpen_trivial (R : CommRingCat.{u}) (P : (Spec R).PresheafOfModules)
    (h : IsLocallyTrivial (Spec R) P) (x : (Spec R)) :
    ∃ g : (R : Type u), x ∈ PrimeSpectrum.basicOpen g ∧
      Nonempty ((restrictPresheafFunctor (Spec R) (PrimeSpectrum.basicOpen g)).obj P
        ≅ 𝟙_ (PresheafModulesOn (Spec R) (PrimeSpectrum.basicOpen g))) := by
  obtain ⟨S, hS, htriv⟩ := h ⊤
  obtain ⟨V, i, hi, hxV⟩ := hS x trivial
  obtain ⟨U', ⟨g, hg⟩, hxU', hU'V⟩ :=
    (Opens.isBasis_iff_nbhd.1 (PrimeSpectrum.isBasis_basic_opens (R := (R : Type u)))) hxV
  subst hg
  exact ⟨g, hxU', ⟨trivialOfLe P hU'V (htriv i hi).some⟩⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_basicOpen_trivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明化を基本開集合へ落とす)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
