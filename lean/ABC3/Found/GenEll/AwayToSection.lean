/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayNumDen
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★`A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)` の写像 —— 段 E1c（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か

大域切断 `s_0, …, s_N` が定める射 `X ⟶ ℙᴺ_R`（段 E1）は、
チャートごとに **`X_{s_i} ⟶ D₊(x_i) ≅ Spec A⁰_{x_i}`** を作って貼り合わせる。
★その各チャートの射は、環準同型

    `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)`,   `a / x_i^n ↦ a(s) / s_i^n`

に対応する。★★本ファイルはその**写像**を作る（環準同型の欄は次段）。

## ★★★★★★機構 —— `Quotient.liftOn'` の 2 つの入力

mathlib の `HomogeneousLocalization 𝒜 x` は
`Quotient (Setoid.ker (NumDenSameDeg.embedding 𝒜 x))` である。したがって要るのは:

| 入力 | 場所 |
|---|---|
| `NumDenSameDeg` から値を作る | `homogRatio ... c.deg c.num`（`§9-811`） |
| 関係を潰す | ★本ファイル（`exists_mul_eq_of_embedding_eq` ＋ `homogRatio_congr_of_mul_eq'`） |

★`embedding c = embedding c'` は `Localization.r_iff_exists` で
`∃ u ∈ powers(x_i), u·(den'·num) = u·(den·num')` になり、
`den = x_i^{deg}`（`§9-813`）を入れると

    `x_i^{deg' + j} · num = x_i^{deg + j} · num'`

になる。★★そこに `homogRatio_congr_of_mul_eq'`（本ファイル）を当てる。

## ★★★★★局所化の余分な `x_i^j` を吸収する

`§9-812` の `homogRatio_congr_of_mul_eq` は `x_i^m·a = x_i^n·b` の形だったが、
局所化の関係には**余分な `x_i^j`** が付く。★本ファイルの `'` 版はそれを込みで扱う
——機構は「両辺に `t_i^{n+m+j}` を掛けて `linear_combination` で潰す」だけである。

## ★残っている段（明示）

★★**環準同型の欄（`map_add` / `map_mul` / `map_zero` / `map_one`）は本ファイルに無い**。
`NumDenSameDeg` の `0`・`1`・`+`・`*` の `deg`/`num` を追って
`homogValue` の環準同型性（`homogValue_mul` 等）に落とす段である。
★★★そこが済めば段 E1c は閉じ、次は E1d（貼り合わせ）になる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization

variable {X : Scheme.{0}}

/-! ## ★★★★★★★(1) 局所化の余分な `x_i^j` を吸収する -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★**`x_i^{m+j}·a = x_i^{n+j}·b` なら比は等しい**。

★`§9-812` の版に**局所化の余分な `x_i^j`** を込めたものである。
★★機構は「両辺に `t_i^{n+m+j}` を掛けて `linear_combination` で潰す」だけである。 -/
theorem homogRatio_congr_of_mul_eq' (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) {n m j : ℕ}
    {a b : MvPolynomial σ R}
    (h : (MvPolynomial.X i) ^ (m + j) * a = (MvPolynomial.X i) ^ (n + j) * b) :
    homogRatio M V e φ s i n a = homogRatio M V e φ s i m b := by
  set ρ := X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op with hρ
  have hu : IsUnit (ρ.hom (trivValue M V e (s i))) := isUnit_trivValue_res M V e (s i)
  have hv : homogValue M V e φ s ((MvPolynomial.X i) ^ (m + j) * a)
      = homogValue M V e φ s ((MvPolynomial.X i) ^ (n + j) * b) := by rw [h]
  rw [homogValue_mul, homogValue_mul, homogValue_pow_X, homogValue_pow_X] at hv
  have hvr : ρ.hom (trivValue M V e (s i)) ^ (m + j) * ρ.hom (homogValue M V e φ s a)
      = ρ.hom (trivValue M V e (s i)) ^ (n + j) * ρ.hom (homogValue M V e φ s b) := by
    have h0 := congrArg (fun z => ρ.hom z) hv
    simpa [map_mul, map_pow] using h0
  have h1 : homogRatio M V e φ s i n a * ρ.hom (trivValue M V e (s i)) ^ n
      = ρ.hom (homogValue M V e φ s a) := homogRatio_mul M V e φ s i n a
  have h2 : homogRatio M V e φ s i m b * ρ.hom (trivValue M V e (s i)) ^ m
      = ρ.hom (homogValue M V e φ s b) := homogRatio_mul M V e φ s i m b
  refine IsUnit.mul_right_cancel (hu.pow (n + m + j)) ?_
  linear_combination (ρ.hom (trivValue M V e (s i)) ^ (m + j)) * h1
    - (ρ.hom (trivValue M V e (s i)) ^ (n + j)) * h2 + hvr

/-! ## ★★★★★★(2) 局所化の関係を多項式の等式に直す -/

/-- ★★★★★★**`embedding c = embedding c'` を多項式の等式に直す**。

★`Localization.r_iff_exists` ＋ `exists_pow_of_numDenSameDeg`（`§9-813`）である。 -/
theorem exists_mul_eq_of_embedding_eq {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (i : Fin (N + 1))
    (c c' : NumDenSameDeg (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R)))
    (h : NumDenSameDeg.embedding _ _ c = NumDenSameDeg.embedding _ _ c') :
    ∃ j : ℕ, (MvPolynomial.X i) ^ (c'.deg + j) * (c.num : MvPolynomial (Fin (N + 1)) R)
      = (MvPolynomial.X i) ^ (c.deg + j) * (c'.num : MvPolynomial (Fin (N + 1)) R) := by
  rw [NumDenSameDeg.embedding, NumDenSameDeg.embedding, Localization.mk_eq_mk_iff,
    Localization.r_iff_exists] at h
  obtain ⟨u, hu⟩ := h
  obtain ⟨j, hj⟩ := u.2
  obtain ⟨k, hk, hdk⟩ := exists_pow_of_numDenSameDeg i c
  obtain ⟨k', hk', hdk'⟩ := exists_pow_of_numDenSameDeg i c'
  refine ⟨j, ?_⟩
  simp only at hu
  rw [hk, hk', ← hj, ← hdk, ← hdk'] at hu
  calc (MvPolynomial.X i) ^ (c'.deg + j) * (c.num : MvPolynomial (Fin (N + 1)) R)
      = (MvPolynomial.X i) ^ j
        * ((MvPolynomial.X i) ^ c'.deg * (c.num : MvPolynomial (Fin (N + 1)) R)) := by
        rw [pow_add]; ring
    _ = (MvPolynomial.X i) ^ j
        * ((MvPolynomial.X i) ^ c.deg * (c'.num : MvPolynomial (Fin (N + 1)) R)) := hu
    _ = (MvPolynomial.X i) ^ (c.deg + j) * (c'.num : MvPolynomial (Fin (N + 1)) R) := by
        rw [pow_add]; ring

/-! ## ★★★★★★★★★(3) 写像そのもの -/

/-- ★★★★★★★★★**`A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)`** —— `a/x_i^n ↦ a(s)/s_i^n`。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これが「大域切断が定める `ℙᴺ` への射」の**チャートごとの中身**である。
★★★環準同型の欄（`map_add` 等）は本ファイルには無い——次段である。 -/
noncomputable def awayToSection (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (z : HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.X i)) : Γ(X, nonVanishing M (s i) ⊓ V) :=
  Quotient.liftOn' z
    (fun c => homogRatio M V e φ s i c.deg (c.num : MvPolynomial (Fin (N + 1)) R))
    (fun c c' h => by
      obtain ⟨j, hj⟩ := exists_mul_eq_of_embedding_eq i c c' h
      exact homogRatio_congr_of_mul_eq' M V e φ s i hj)

/-! ## ★出典の紐付け(`.src`) -/

def homogRatio_congr_of_mul_eq'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(局所化の余分な x_i^j を吸収した形)",
    sectionId := "genell-prop-1-4" }

def exists_mul_eq_of_embedding_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(局所化の関係を多項式の等式に直す)",
    sectionId := "genell-prop-1-4" }

def awayToSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V) の写像。環準同型の欄は含まない)",
    sectionId := "genell-prop-1-4" }

def awayToSection.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "homogRatio(斉次式の比、§9-811)"
      (.inProject "ABC3" "ABC3.Found.GenEll.homogRatio") 7,
    .citation "[ABC3]" "exists_pow_of_numDenSameDeg(分母は x_i^k で次数も k、§9-813)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_of_numDenSameDeg") 7,
    .citation "[mathlib]" "Localization.r_iff_exists / HomogeneousLocalization の Quotient 表示"
      (.inMathlib "Localization.r_iff_exists") 7,
    .implicitStep
      ("★**環準同型の欄(map_add / map_mul / map_zero / map_one)は本ファイルに無い**。" ++
       "NumDenSameDeg の 0・1・+・* の deg/num を追って homogValue の環準同型性" ++
       "(homogValue_mul 等)に落とす段である。" ++
       "★★そこが済めば段 E1c は閉じ、次は E1d(貼り合わせ)になる") 7 ]

end ABC3.Found.GenEll
