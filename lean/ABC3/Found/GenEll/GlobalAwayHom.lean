/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalRatio
import ABC3.Found.GenEll.AwayNumDen
import ABC3.Found.GenEll.AwayToSection
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E3c-2 —— 大域のチャート写像 `A⁰_{x_i} →+* Γ(X, X_{s_i})`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 「`⊓ V`」が消えた

`§9-841` で見つけた欠陥——`§9-815` の `awayToSectionHom` は
`A⁰_{x_i} →+* Γ(X, X_{s_i} **⊓ V**)` という形で、段 E3d が要る
`Γ(X, X_{s_i})`（`⊓ V` なし）と噛み合っていなかった——を本ファイルが埋める。

★**`trivValue` の代わりに `globalRatio`（`§9-841`）を使う**だけで `⊓ V` が消える。

## ★★★★★★★機構 —— むしろ簡単になった

★★`Ψ ≝ eval₂Hom φ (fun j ↦ s_j/s_i)` は**多項式環全体の環準同型**であり、
`Ψ(x_i) = s_i/s_i = 1` である。★★★したがって

* **分母は自動で消える**（`den` は `x_i` の冪、`§9-813`）
* **well-defined 性は `Ψ(x_i) = 1` だけで出る**——
  `x_i^{k}·num c = x_i^{k'}·num c'`（`§9-814`）に `Ψ` を当てるだけ

★`§9-814`〜`§9-815` は `trivValue` の言葉で比を作っていたので
`homogRatio_congr_of_mul_eq'` のような補題が要ったが、
**大域の比を使うと環準同型 1 本で済む**——`n` 次斉次式の非同次化そのものだからである。

## ★これが段 E3c-2 の中身である

★★`globalAwayHom … (x_j/x_i) = s_j/s_i`（`globalAwayHom_projCoord`）なので、
`§9-838` の座標条件（`g·s_i = s_j`）は「`g = s_j/s_i`」と読み替えられる。

## ★残っている段（明示）

★★★段 E3c の全射性を**この大域版で言い直す**段が残る
——`§9-833`（有限生成 ⟹ 全射判定）はそのまま使えるが、
`§9-837`（比の判定）を `globalRatio_unique` で置き換える必要がある。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★自分自身との比は `1` -/

/-- ★**自分自身との比は `1` である**。 -/
theorem globalRatio_self (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (t : (M.obj (op ⊤) : Type)) : globalRatio M hM t t = 1 := by
  refine (globalRatio_unique M hM t t 1 (fun i => ?_)).symm
  rw [map_one, sectionRatio_self]

/-! ## ★★★★★★大域の同次評価 -/

/-- ★★★★★★**大域の同次評価** `a ↦ a(s_0/s_i, …, s_N/s_i)`。

★これは**多項式環全体の環準同型**である——`trivValue` 版と違い `⊓ V` が付かない。
★★`n` 次斉次式に当てれば `a(s)/s_i^n` の非同次化そのものである。 -/
noncomputable def globalHomogHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) :
    MvPolynomial σ R →+* (Γ(X, nonVanishing M (s i)) : Type) :=
  MvPolynomial.eval₂Hom
    (((X.presheaf.map (homOfLE (le_top : nonVanishing M (s i) ≤ ⊤)).op).hom).comp φ)
    (fun j => globalRatio M hM (s j) (s i))

theorem globalHomogHom_X (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i j : σ) :
    globalHomogHom M hM φ s i (MvPolynomial.X j) = globalRatio M hM (s j) (s i) := by
  rw [globalHomogHom, MvPolynomial.eval₂Hom_X']

/-- ★★★**割る変数は `1` に行く** —— これが分母を消す。 -/
theorem globalHomogHom_X_self (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) :
    globalHomogHom M hM φ s i (MvPolynomial.X i) = 1 := by
  rw [globalHomogHom_X, globalRatio_self]

/-- ★★★★**`A⁰_{x_i}` の分母はすべて `1` に行く**（`den` は `x_i` の冪、`§9-813`）。 -/
theorem globalHomogHom_den (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (c : NumDenSameDeg (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
    globalHomogHom M hM φ s i (c.den : MvPolynomial (Fin (N + 1)) R) = 1 := by
  obtain ⟨k, hk, -⟩ := exists_pow_of_numDenSameDeg i c
  rw [hk, map_pow, globalHomogHom_X_self, one_pow]

/-! ## ★★★★★★★★★★大域のチャート写像 -/

/-- ★★★★★**`A⁰_{x_i}` の元を大域の切断に送る**（写像そのもの）。

★well-defined 性は **`Ψ(x_i) = 1` だけ**で出る
——`§9-814` の `x_i^k·num c = x_i^{k'}·num c'` に `Ψ` を当てるだけである。 -/
noncomputable def globalAwayToSection (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (z : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) :
    (Γ(X, nonVanishing M (s i)) : Type) :=
  Quotient.liftOn' z
    (fun c => globalHomogHom M hM φ s i (c.num : MvPolynomial (Fin (N + 1)) R))
    (fun c c' h => by
      obtain ⟨j, hj⟩ := exists_mul_eq_of_embedding_eq i c c' h
      have h2 := congrArg (globalHomogHom M hM φ s i) hj
      simpa [map_mul, map_pow, globalHomogHom_X_self] using h2)

/-- ★★★★★★★★★★**大域のチャート写像** `A⁰_{x_i} →+* Γ(X, X_{s_i})` —— 段 E3c-2。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★**`⊓ V` が付かない**——これが段 E3d（`§9-834`）の消費する形である。

★★機構は `globalHomogHom`（多項式環全体の環準同型）で、
分母は `Ψ(x_i) = 1` により自動で消える。
★★★`§9-815` の `trivValue` 版より**むしろ簡単になった**
——大域の比を使うと環準同型 1 本で済むからである。 -/
noncomputable def globalAwayHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)
      →+* (Γ(X, nonVanishing M (s i)) : Type) where
  toFun := globalAwayToSection M hM φ s i
  map_one' := by
    show globalHomogHom M hM φ s i
      ((NumDenSameDeg.num (1 : NumDenSameDeg
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
          MvPolynomial (Fin (N + 1)) R)) = 1
    rw [NumDenSameDeg.num_one]
    exact map_one _
  map_zero' := by
    show globalHomogHom M hM φ s i
      ((NumDenSameDeg.num (0 : NumDenSameDeg
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
          MvPolynomial (Fin (N + 1)) R)) = 0
    rw [NumDenSameDeg.num_zero, ZeroMemClass.coe_zero]
    exact map_zero _
  map_mul' := by
    rintro ⟨c⟩ ⟨c'⟩
    show globalHomogHom M hM φ s i (((c * c').num : MvPolynomial (Fin (N + 1)) R)) = _
    rw [NumDenSameDeg.num_mul, map_mul]
    rfl
  map_add' := by
    rintro ⟨c⟩ ⟨c'⟩
    show globalHomogHom M hM φ s i (((c + c').num : MvPolynomial (Fin (N + 1)) R)) = _
    rw [NumDenSameDeg.num_add, map_add, map_mul, map_mul,
      globalHomogHom_den M hM φ s i c, globalHomogHom_den M hM φ s i c', one_mul, one_mul]
    show _ = globalHomogHom M hM φ s i (c.num : MvPolynomial (Fin (N + 1)) R)
      + globalHomogHom M hM φ s i (c'.num : MvPolynomial (Fin (N + 1)) R)
    ring

/-- ★★★★★★**正規化座標 `x_j/x_i` は大域の比 `s_j/s_i` に行く**。 -/
theorem globalAwayHom_projCoord (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) :
    globalAwayHom M hM φ s i (projCoord N R i j) = globalRatio M hM (s j) (s i) := by
  show globalHomogHom M hM φ s i (MvPolynomial.X j) = _
  rw [globalHomogHom_X]

/-- ★★★★★★★**大域の比は像に入る** —— 段 E3c を `⊓ V` なしで言うための橋。 -/
theorem mem_range_globalAwayHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    (g : (Γ(X, nonVanishing M (s i)) : Type)) (hg : g = globalRatio M hM (s j) (s i)) :
    g ∈ Set.range (globalAwayHom M hM φ s i) :=
  ⟨projCoord N R i j, by rw [globalAwayHom_projCoord, hg]⟩

/-! ## ★出典の紐付け(`.src`) -/

def globalHomogHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(大域の同次評価——⊓ V なし)",
    sectionId := "genell-prop-1-4" }

def globalAwayHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c-2——大域のチャート写像 A⁰_{x_i} → Γ(X, X_{s_i}))",
    sectionId := "genell-prop-1-4" }

def globalAwayHom_projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(正規化座標は大域の比に行く)",
    sectionId := "genell-prop-1-4" }

def mem_range_globalAwayHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(大域の比は像に入る)",
    sectionId := "genell-prop-1-4" }

def globalAwayHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "globalRatio(大域の比、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio") 2,
    .citation "[ABC3]" "exists_pow_of_numDenSameDeg(分母は x_i の冪、§9-813)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_of_numDenSameDeg") 2,
    .citation "[ABC3]" "exists_mul_eq_of_embedding_eq(§9-814)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_mul_eq_of_embedding_eq") 2,
    .implicitStep
      ("★大域の比を使うと**むしろ簡単になった**——Ψ ≝ eval₂Hom φ (fun j ↦ s_j/s_i) は" ++
       "多項式環全体の環準同型で Ψ(x_i) = 1 なので、分母は自動で消え、" ++
       "well-defined 性も Ψ(x_i) = 1 だけで出る。" ++
       "★§9-814〜§9-815 は trivValue の言葉で比を作っていたので homogRatio_congr_of_mul_eq' の" ++
       "ような補題が要ったが、ここでは環準同型 1 本で済む——n 次斉次式の**非同次化そのもの**だからである") 5,
    .implicitStep
      ("★★★段 E3c の全射性を**この大域版で言い直す**段が残る" ++
       "——§9-833(有限生成 ⟹ 全射判定)はそのまま使えるが、" ++
       "§9-837(比の判定)を globalRatio_unique で置き換える必要がある") 6 ]

end ABC3.Found.GenEll
