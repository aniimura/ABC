import ABC3.Found.Arakelov.PicPredSieve

/-!
# Arakelov (B2) 第 237 ブロック —— **自明化は制限で降り、点ごとに取れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★第 236 の述語を作るための 2 つの道具

引き戻しの比較射が全単射になる開集合 `V` を**点ごとに**作るには

| 要る性質 | 道具 |
|---|---|
| 自明化開の中でさらに小さく取れる | ★`restrict_trivial_of_le` |
| 局所自明性から点の近傍を取り出す | ★`exists_trivial_nbhd` |

★★`restrict_trivial_of_le` は第 57 ブロック(制限の推移律が `rfl`)の**帰結**である
——`eqToIso` を 2 つ噛ませるだけで済む。

★★★`exists_trivial_nbhd` は `Opens` の Grothendieck 位相が**点ごと**であること
(第 58 ブロックと同じ事実)をそのまま使う。

| 定義・定理 | 内容 |
|---|---|
| `restrict_trivial_of_le` | ★★自明化は小さい開集合へ降りる |
| `exists_trivial_nbhd` | ★★局所自明性から点の近傍を取り出す |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★★**自明化は小さい開集合へ降りる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 57 ブロックで `(M|_A)|_B = M|_B` と `𝟙_|_B = 𝟙_` がどちらも `rfl` と
測ってあるので、`eqToIso` を 2 つ噛ませるだけである。 -/
noncomputable def restrict_trivial_of_le {M : X.PresheafOfModules} {A B : X.Opens} (h : B ≤ A)
    (e : (restrictPresheafFunctor X A).obj M ≅ 𝟙_ (PresheafModulesOn X A)) :
    (restrictPresheafFunctor X B).obj M ≅ 𝟙_ (PresheafModulesOn X B) :=
  eqToIso (restrict_trans h M).symm ≪≫ (restrictOnFunctor h).mapIso e
    ≪≫ eqToIso (restrictOnUnit h)

/-- ★★**局所自明性から点の近傍を取り出す**——位相が点ごとだから。 -/
theorem exists_trivial_nbhd {M : X.PresheafOfModules} (hM : IsLocallyTrivial X M)
    (U : X.Opens) (x : X) (hx : x ∈ U) :
    ∃ V : X.Opens, x ∈ V ∧ V ≤ U ∧
      Nonempty ((restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) := by
  obtain ⟨S, hS, hSt⟩ := hM U
  obtain ⟨V, i, hi, hxV⟩ := hS x hx
  exact ⟨V, hxV, leOfHom i, hSt i hi⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_trivial_nbhd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明化は制限で降り、点ごとに取れる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
