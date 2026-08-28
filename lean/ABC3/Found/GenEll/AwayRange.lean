/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayHom
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★★段 E3c の橋 —— 「`g·s_i = s_j` なら像に入る」（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 段 E3c の最後の橋

`§9-833` は段 E3c を「**有限個**の生成元が `A⁰_{x_i} → Γ(X, X_{s_i})` の像に入れば全射」に
落とした。★本ファイルは**その「像に入る」を確認する道具**である:

> `g · s_i = s_j`（座標の言葉で）なら、`g` は像に入る——それは `x_j/x_i` の像である。

## ★★★★★機構 —— 比の一意性だけ

`homogRatio_unique`（`§9-811`）が「`r·s_i^n = a(s)` を満たす `r` は `a(s)/s_i^n` に限る」と言う。
★`a = x_j`、`n = 1` と取れば本ファイルの主張になる。

★★`awayToSectionHom … (projCoord N R i j) = homogRatio … 1 (x_j)` は **`rfl`** である
（2026-08-28 実測）——`Away.mk` の `deg` が `1·1 = 1`、`num` が `x_j` だからである。

## ★★★段 E3a-3 との噛み合い

段 E3a-3（`§9-831`）が出すのは「`g · f_i^n = a_g`」（`n` 乗）である。
★一方本ファイルが要求するのは**1 乗**である。
★★これは `M` を `M^{⊗n}` に、`s_i` を `s_i^{⊗n}` に取り替えれば一致する
——`trivValue_secPow`（`§9-825`）が `trivValue (s^{⊗n}) = (trivValue s)^n` を与えるからである。
★★★**冪を取る段（§9-825）がここで効く**。

## ★残っている段（明示）

★★**座標の選び方**は本ファイルに無い——`s : Fin (N+1) → M(⊤)` に
「各チャートの生成元を延ばして得た切断」を**含める**ように取る段である。
★それは有限個（チャート有限×生成元有限）なので `N` を大きく取れば済むが、
**その台帳付けは書いていない**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov
open MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★正規化座標の像 -/

/-- ★**`x_j/x_i` の像は比 `s_j/s_i` である**（`rfl`）。

★`Away.mk` の `deg` が `1·1 = 1`、`num` が `x_j` なので定義的に等しい（2026-08-28 実測）。 -/
theorem awayToSectionHom_projCoord (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) :
    awayToSectionHom M V e φ s i (projCoord N R i j)
      = homogRatio M V e φ s i 1 (MvPolynomial.X j) := rfl

/-- ★★★**`x_j/x_i` の像に `s_i` を掛けると `s_j` になる**——比の特徴づけ。 -/
theorem awayToSectionHom_projCoord_mul (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) :
    awayToSectionHom M V e φ s i (projCoord N R i j)
        * X.presheaf.map (homOfLE (inf_le_right :
            nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s i))
      = X.presheaf.map (homOfLE (inf_le_right :
          nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s j)) := by
  rw [awayToSectionHom_projCoord]
  have h := homogRatio_mul M V e φ s i 1 (MvPolynomial.X j)
  rw [pow_one] at h
  rw [h, homogValue_X]

/-! ## ★★★★★★★★段 E3c の橋 -/

/-- ★★★★★★★★**`g·s_i = s_j` なら `g` は像に入る** —— 段 E3c の最後の橋。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`§9-833` が段 E3c を「有限個の生成元が像に入れば全射」に落としたので、
★本補題があれば**確認は「その生成元に `s_i` を掛けたものが座標の 1 つになるか」だけ**になる。

★★機構は `homogRatio_unique`（`§9-811`——比は一意）である。

★★★段 E3a-3（`§9-831`）が出すのは `n` 乗の形だが、
`M` を `M^{⊗n}`・`s_i` を `s_i^{⊗n}` に取り替えれば 1 乗になる
（`trivValue_secPow`、`§9-825`）——**冪を取る段がここで効く**。 -/
theorem mem_range_awayToSectionHom_of_mul_eq (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    (g : (Γ(X, nonVanishing M (s i) ⊓ V) : Type))
    (hg : g * X.presheaf.map (homOfLE (inf_le_right :
            nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s i))
      = X.presheaf.map (homOfLE (inf_le_right :
          nonVanishing M (s i) ⊓ V ≤ V)).op (trivValue M V e (s j))) :
    g ∈ Set.range (awayToSectionHom M V e φ s i) := by
  refine ⟨projCoord N R i j, ?_⟩
  rw [awayToSectionHom_projCoord]
  refine (homogRatio_unique M V e φ s i 1 (MvPolynomial.X j) g ?_).symm
  rw [pow_one, hg, homogValue_X]

/-! ## ★出典の紐付け(`.src`) -/

def awayToSectionHom_projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(x_j/x_i の像は比 s_j/s_i である)",
    sectionId := "genell-prop-1-4" }

def awayToSectionHom_projCoord_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(x_j/x_i の像に s_i を掛けると s_j になる)",
    sectionId := "genell-prop-1-4" }

def mem_range_awayToSectionHom_of_mul_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c の橋——g·s_i = s_j なら像に入る)",
    sectionId := "genell-prop-1-4" }

def mem_range_awayToSectionHom_of_mul_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "homogRatio_unique(比は一意、§9-811)"
      (.inProject "ABC3" "ABC3.Found.GenEll.homogRatio_unique") 2,
    .citation "[ABC3]" "exists_finset_surjective_criterion(段 E3c の還元、§9-833)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_surjective_criterion") 2,
    .citation "[ABC3]" "trivValue_secPow(冪と座標は可換、§9-825)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_secPow") 2,
    .implicitStep
      ("★段 E3a-3(§9-831)が出すのは n 乗の形だが、M を M^{⊗n}・s_i を s_i^{⊗n} に" ++
       "取り替えれば 1 乗になる(trivValue_secPow、§9-825)。**冪を取る段がここで効く**") 6,
    .implicitStep
      ("★★**座標の選び方**は本ファイルに無い——s : Fin (N+1) → M(⊤) に" ++
       "「各チャートの生成元を延ばして得た切断」を**含める**ように取る段である。" ++
       "★有限個(チャート有限×生成元有限)なので N を大きく取れば済むが、その台帳付けは書いていない") 7 ]

end ABC3.Found.GenEll
