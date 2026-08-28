/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayToSection
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★環準同型 `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)` —— **段 E1c が閉じた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これで段 E1c が閉じた

`§9-814` で写像 `awayToSection` は作った。本ファイルはその**環準同型の欄**を埋める:

| 欄 | 中身 |
|---|---|
| `map_one` | `deg 1 = 0`・`num 1 = 1` なので `homogRatio 0 1 = 1` |
| `map_zero` | `num 0 = 0` なので比も `0` |
| `map_mul` | ★`homogRatio (n+m) (a·b) = homogRatio n a · homogRatio m b` |
| `map_add` | ★★`homogRatio (n+m) (x_i^m·a + x_i^n·b) = homogRatio n a + homogRatio m b` |

★`map_add` だけは分母 `den = x_i^{deg}`（`§9-813`）を使う
——`NumDenSameDeg` の和の分子が `den₁·num₂ + den₂·num₁` だからである。

## ★★★★★★機構 —— どれも `homogRatio_unique` に落ちる

比は分母が単元なので**一意**である（`homogRatio_unique`、`§9-811`）。
★したがって「候補に `s_i^{n+m}` を掛けると `a(s)` になる」ことだけ確かめればよく、
それは `homogValue` が環準同型であること（`eval₂`）から出る。

## ★測定の記録

`NumDenSameDeg.num_add` / `deg_add` は **`x`（分母の submonoid）が明示引数**である
（`{𝒜}` は暗黙、`x` は明示）。★`num_add c1 c2` と書くと `c1` が `x` に食われる（2026-08-28 実測）。

★★`(c1 + c2).deg` は `num` の**型に現れる**ので、先に `num` を書き換えてから
`deg` を書き換える必要がある（順序を逆にすると motive が通らない）。

## ★残っている段（明示）

★★★次は **E1d（貼り合わせ）**である——各チャートの環準同型から
`X_{s_i} ⟶ Spec A⁰_{x_i} ≅ D₊(x_i)` を作り、`X` の上へ貼る。
材料は `Scheme.homSpecIso`（環準同型 ⟷ `Spec` への射）と
`Proj.awayι`（`Spec (A⁰_f) ⟶ Proj 𝒜` の開埋め込み、mathlib）である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★(1) 比は積と和を保つ -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**比は積を保つ**: `a·b / s_i^{n+m} = (a/s_i^n)·(b/s_i^m)`。 -/
theorem homogRatio_mul_mul (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n m : ℕ)
    (a b : MvPolynomial σ R) :
    homogRatio M V e φ s i (n + m) (a * b)
      = homogRatio M V e φ s i n a * homogRatio M V e φ s i m b := by
  set ρ := X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op with hρ
  refine (homogRatio_unique M V e φ s i (n + m) (a * b) _ ?_).symm
  have h1 : homogRatio M V e φ s i n a * ρ.hom (trivValue M V e (s i)) ^ n
      = ρ.hom (homogValue M V e φ s a) := homogRatio_mul M V e φ s i n a
  have h2 : homogRatio M V e φ s i m b * ρ.hom (trivValue M V e (s i)) ^ m
      = ρ.hom (homogValue M V e φ s b) := homogRatio_mul M V e φ s i m b
  have hv : ρ.hom (homogValue M V e φ s (a * b))
      = ρ.hom (homogValue M V e φ s a) * ρ.hom (homogValue M V e φ s b) := by
    rw [homogValue_mul, map_mul]
  rw [hv, ← h1, ← h2]
  ring

theorem homogValue_add (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (p q : MvPolynomial σ R) :
    homogValue M V e φ s (p + q) = homogValue M V e φ s p + homogValue M V e φ s q :=
  MvPolynomial.eval₂_add _ _

set_option maxHeartbeats 1000000 in
/-- ★★★★★**比は和を保つ**: `(x_i^m·a + x_i^n·b)/s_i^{n+m} = a/s_i^n + b/s_i^m`。 -/
theorem homogRatio_add_add (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n m : ℕ)
    (a b : MvPolynomial σ R) :
    homogRatio M V e φ s i (n + m)
        ((MvPolynomial.X i) ^ m * a + (MvPolynomial.X i) ^ n * b)
      = homogRatio M V e φ s i n a + homogRatio M V e φ s i m b := by
  set ρ := X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op with hρ
  refine (homogRatio_unique M V e φ s i (n + m) _ _ ?_).symm
  have h1 : homogRatio M V e φ s i n a * ρ.hom (trivValue M V e (s i)) ^ n
      = ρ.hom (homogValue M V e φ s a) := homogRatio_mul M V e φ s i n a
  have h2 : homogRatio M V e φ s i m b * ρ.hom (trivValue M V e (s i)) ^ m
      = ρ.hom (homogValue M V e φ s b) := homogRatio_mul M V e φ s i m b
  have hv : ρ.hom (homogValue M V e φ s
        ((MvPolynomial.X i) ^ m * a + (MvPolynomial.X i) ^ n * b))
      = ρ.hom (trivValue M V e (s i)) ^ m * ρ.hom (homogValue M V e φ s a)
        + ρ.hom (trivValue M V e (s i)) ^ n * ρ.hom (homogValue M V e φ s b) := by
    rw [homogValue_add, homogValue_mul, homogValue_mul, homogValue_pow_X, homogValue_pow_X,
      map_add, map_mul, map_mul, map_pow, map_pow]
  rw [hv, ← h1, ← h2]
  ring

/-! ## ★★★★★★★★(2) 環準同型の欄 -/

theorem awayToSection_one (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    awayToSection M V e φ s i 1 = 1 := by
  show homogRatio M V e φ s i 0 (1 : MvPolynomial (Fin (N + 1)) R) = 1
  show _ * _ ^ (0 : ℕ) = 1
  rw [pow_zero, mul_one]
  have h : homogValue M V e φ s (1 : MvPolynomial (Fin (N + 1)) R) = 1 :=
    MvPolynomial.eval₂_one _ _
  rw [h, map_one]

theorem awayToSection_zero (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    awayToSection M V e φ s i 0 = 0 := by
  show homogRatio M V e φ s i 0 (0 : MvPolynomial (Fin (N + 1)) R) = 0
  show _ * _ = 0
  have h : homogValue M V e φ s (0 : MvPolynomial (Fin (N + 1)) R) = 0 :=
    MvPolynomial.eval₂_zero _ _
  rw [h, map_zero, zero_mul]

set_option maxHeartbeats 1000000 in
theorem awayToSection_mul (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (z w : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) :
    awayToSection M V e φ s i (z * w)
      = awayToSection M V e φ s i z * awayToSection M V e φ s i w := by
  induction z using Quotient.ind with
  | _ c1 =>
    induction w using Quotient.ind with
    | _ c2 =>
      show homogRatio M V e φ s i (c1.deg + c2.deg)
          ((c1.num : MvPolynomial (Fin (N + 1)) R) * (c2.num : MvPolynomial (Fin (N + 1)) R))
        = homogRatio M V e φ s i c1.deg (c1.num : MvPolynomial (Fin (N + 1)) R)
          * homogRatio M V e φ s i c2.deg (c2.num : MvPolynomial (Fin (N + 1)) R)
      exact homogRatio_mul_mul M V e φ s i c1.deg c2.deg _ _

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**和の欄**——分子が `den₁·num₂ + den₂·num₁` なので
分母の同定（`§9-813`）が要る。 -/
theorem awayToSection_add (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (z w : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) :
    awayToSection M V e φ s i (z + w)
      = awayToSection M V e φ s i z + awayToSection M V e φ s i w := by
  induction z using Quotient.ind with
  | _ c1 =>
    induction w using Quotient.ind with
    | _ c2 =>
      obtain ⟨k1, hk1, hd1⟩ := exists_pow_of_numDenSameDeg i c1
      obtain ⟨k2, hk2, hd2⟩ := exists_pow_of_numDenSameDeg i c2
      have hk1' : (c1.den : MvPolynomial (Fin (N + 1)) R)
          = MvPolynomial.X i ^ c1.deg := by rw [hk1, hd1]
      have hk2' : (c2.den : MvPolynomial (Fin (N + 1)) R)
          = MvPolynomial.X i ^ c2.deg := by rw [hk2, hd2]
      have hnum := NumDenSameDeg.num_add
        (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R)) c1 c2
      have hdeg := NumDenSameDeg.deg_add
        (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R)) c1 c2
      show homogRatio M V e φ s i (c1 + c2).deg
          ((c1 + c2).num : MvPolynomial (Fin (N + 1)) R) = _
      rw [hnum, hdeg, hk1', hk2', Nat.add_comm c1.deg c2.deg,
        homogRatio_add_add M V e φ s i c2.deg c1.deg]
      exact add_comm _ _

/-! ## ★★★★★★★★★★(3) 環準同型 -/

/-- ★★★★★★★★★★**環準同型 `A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V)`**——**段 E1c が閉じた**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これが「大域切断が定める `ℙᴺ` への射」の**チャートごとの中身**である。
★★★次は `Scheme.homSpecIso` で `X_{s_i} ⟶ Spec A⁰_{x_i}` に直し、
`Proj.awayι` で `ℙᴺ` へ入れて貼り合わせる（段 E1d）。 -/
noncomputable def awayToSectionHom (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] [Nontrivial R] (φ : R →+* (Γ(X, V) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.X i) →+* (Γ(X, nonVanishing M (s i) ⊓ V) : Type) where
  toFun := awayToSection M V e φ s i
  map_one' := awayToSection_one M V e φ s i
  map_mul' := awayToSection_mul M V e φ s i
  map_zero' := awayToSection_zero M V e φ s i
  map_add' := awayToSection_add M V e φ s i

/-! ## ★出典の紐付け(`.src`) -/

def homogRatio_mul_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比は積を保つ)",
    sectionId := "genell-prop-1-4" }

def homogRatio_add_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比は和を保つ)",
    sectionId := "genell-prop-1-4" }

def awayToSectionHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(環準同型 A⁰_{x_i} → Γ(X, X_{s_i} ⊓ V))",
    sectionId := "genell-prop-1-4" }

def awayToSectionHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "homogRatio_unique(比は一意、§9-811)"
      (.inProject "ABC3" "ABC3.Found.GenEll.homogRatio_unique") 7,
    .citation "[ABC3]" "awayToSection(写像そのもの、§9-814)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayToSection") 7,
    .citation "[ABC3]" "exists_pow_of_numDenSameDeg(分母は x_i^{deg}、§9-813)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_of_numDenSameDeg") 7,
    .implicitStep
      ("★NumDenSameDeg.num_add / deg_add は **x(分母の submonoid)が明示引数**である" ++
       "(2026-08-28 実測——num_add c1 c2 と書くと c1 が x に食われる)。" ++
       "★★(c1 + c2).deg は num の**型に現れる**ので、先に num を書き換えてから" ++
       "deg を書き換える必要がある") 7,
    .implicitStep
      ("★★★次は段 E1d(貼り合わせ)——Scheme.homSpecIso で X_{s_i} ⟶ Spec A⁰_{x_i} に直し、" ++
       "Proj.awayι(mathlib)で ℙᴺ へ入れて貼る") 7 ]

end ABC3.Found.GenEll
