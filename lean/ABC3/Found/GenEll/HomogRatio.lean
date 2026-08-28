/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HomogValue
import ABC3.Meta.Claim

/-!
# ★★★★★★★★**斉次式の比 `a(s)/s_i^n`** —— 段 E1b（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは台帳の段 E1b である

| 段 | 内容 | 状態 |
|---|---|---|
| E1a | 斉次多項式を切断で評価する（`homogValue`） | ★済（`§9-810`） |
| **E1b** | **`a(s)/s_i^n` が `Γ(X, X_{s_i} ⊓ V)` の元を定める** | ★★本ファイル |
| E1c | 環準同型 `A⁰_{x_i} → Γ(X, X_{s_i})` の well-defined 性 | 残る |
| E1d | 貼り合わせて `X ⟶ ℙᴺ_R` | 残る |

★段 D3 の `sectionRatio` は**次数 1** の比（`s/t`、どちらも `Γ(X,M)`）だけを扱っていた。
★★本ファイルはそれを **`n` 次斉次式 `a` に対する `a(s)/s_i^n`** へ広げる。

## ★★★★★★機構 —— `sectionRatio` と同じ形

`X_{s_i} ⊓ V` の上では `t_i = trivValue M V e (s i)` の制限が**単元**である
（`isUnit_trivValue_res`、段 D3）。★そこで

    `homogRatio ≝ (a(t) の制限) · ((t_i の制限)⁻¹)^n`

と置く。★★**型が `e` を含まない**のが要点である
——`nonVanishing_inf`（段 D2）で `X.basicOpen (t_i)` を `X_{s_i} ⊓ V` と書けるからで、
これで「自明化に依らない」が**転送なしの等式**として言える（`homogRatio_congr`）。

## ★★★★★★★自明化に依らないこと

`§9-810` の `homogValue_congr` が「分子は `u^n` 倍」を、
`trivValue_congr'` が「分母 `t_i` は `u` 倍」を与える。★分母の `n` 乗も `u^n` 倍なので、
比では `u^n` が消える。★★機構は `sectionRatio_congr`（段 D3）と**同じ**である。

## ★★★★★これが `northcott_of_projModel` の同次座標である

原文の「射影埋め込みの座標」は `s_i/s_j` であり、その `n` 次版が `a(s)/s_i^n` である。
★`homogRatio_self`（割る成分そのものの比は `1`）は `projCoord_self`（段 C2a）の類似である。

## ★残っている段（明示）

★★環準同型 `A⁰_{x_i} → Γ(X, X_{s_i})`（段 E1c）は本ファイルに無い。
`HomogeneousLocalization.Away.mk` の同値関係（`mk n a = mk m b`）を潰す段が要る
——`MvPolynomial` は整域なので `x_i` は非零因子であり、消去できるはずである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-- ★★★★★★★★**斉次式の比 `a(s)/s_i^n`** —— `X_{s_i} ⊓ V` の上の正則関数。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★型が `e` を含まないのが要点である（`sectionRatio` と同じ設計）。 -/
noncomputable def homogRatio (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n : ℕ)
    (a : MvPolynomial σ R) : Γ(X, nonVanishing M (s i) ⊓ V) :=
  X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op
      (homogValue M V e φ s a)
    * (↑(isUnit_trivValue_res M V e (s i)).unit⁻¹ :
        (Γ(X, nonVanishing M (s i) ⊓ V) : Type)) ^ n

/-- ★★★**比の特徴づけ**: `(a(s)/s_i^n) · s_i^n = a(s)`。 -/
theorem homogRatio_mul (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n : ℕ)
    (a : MvPolynomial σ R) :
    homogRatio M V e φ s i n a
        * (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op
            (trivValue M V e (s i))) ^ n
      = X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op
          (homogValue M V e φ s a) := by
  have h := (isUnit_trivValue_res M V e (s i)).val_inv_mul
  simp only [homogRatio, mul_assoc, ← mul_pow, h, one_pow, mul_one]

/-- ★★★**比は一意である**（分母が単元だから）。 -/
theorem homogRatio_unique (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n : ℕ)
    (a : MvPolynomial σ R) (r : Γ(X, nonVanishing M (s i) ⊓ V))
    (hr : r * (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op
            (trivValue M V e (s i))) ^ n
      = X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op
          (homogValue M V e φ s a)) :
    r = homogRatio M V e φ s i n a := by
  refine IsUnit.mul_right_cancel ((isUnit_trivValue_res M V e (s i)).pow n) ?_
  rw [hr, homogRatio_mul]

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**比は自明化の取り方に依らない**（`n` 次斉次式の場合）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★分子は `u^n` 倍（`homogValue_congr`、`§9-810`）、分母 `t_i^n` も `u^n` 倍なので消える。
★★機構は `sectionRatio_congr`（段 D3）と**同じ**である。 -/
theorem homogRatio_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) {n : ℕ}
    {a : MvPolynomial σ R} (ha : a.IsHomogeneous n) :
    homogRatio M V e' φ s i n a = homogRatio M V e φ s i n a := by
  obtain ⟨u, hu, hall⟩ := trivValue_congr' M V e e'
  set ρ := X.presheaf.map (homOfLE (inf_le_right : nonVanishing M (s i) ⊓ V ≤ V)).op with hρ
  have huρ : IsUnit (ρ.hom u) := hu.map ρ.hom
  have hval : homogValue M V e' φ s a = homogValue M V e φ s a * u ^ n := by
    show MvPolynomial.eval₂ φ (fun j => trivValue M V e' (s j)) a = _
    have hx : (fun j => trivValue M V e' (s j)) = (fun j => trivValue M V e (s j) * u) :=
      funext (fun j => hall (s j))
    rw [hx]
    exact eval₂_mul_const_of_isHomogeneous φ _ u ha
  refine homogRatio_unique M V e φ s i n a (homogRatio M V e' φ s i n a) ?_
  have h1 := homogRatio_mul M V e' φ s i n a
  rw [hall (s i), hval, map_mul, map_mul, map_pow, mul_pow] at h1
  refine IsUnit.mul_right_cancel (huρ.pow n) ?_
  rw [mul_assoc]
  exact h1

/-- ★★★**割る成分そのものの比は `1`**（`projCoord_self`（段 C2a）の類似）。 -/
theorem homogRatio_self (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (i : σ) (n : ℕ) :
    homogRatio M V e φ s i n (MvPolynomial.X i ^ n) = 1 := by
  refine (homogRatio_unique M V e φ s i n _ 1 ?_).symm
  rw [one_mul]
  have hv : homogValue M V e φ s (MvPolynomial.X i ^ n) = trivValue M V e (s i) ^ n := by
    show MvPolynomial.eval₂ φ (fun j => trivValue M V e (s j)) (MvPolynomial.X i ^ n) = _
    rw [MvPolynomial.eval₂_pow, MvPolynomial.eval₂_X]
  rw [hv, map_pow]

/-! ## ★出典の紐付け(`.src`) -/

def homogRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(斉次式の比 a(s)/s_i^n——射影埋め込みの座標の n 次版)",
    sectionId := "genell-prop-1-4" }

def homogRatio_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比は自明化の取り方に依らない——n 次斉次式の場合)",
    sectionId := "genell-prop-1-4" }

def homogRatio_self.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(割る成分そのものの比は 1)",
    sectionId := "genell-prop-1-4" }

def homogRatio_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "homogValue_congr(自明化を替えると n 次斉次式の値は u^n 倍、§9-810)"
      (.inProject "ABC3" "ABC3.Found.GenEll.homogValue_congr") 7,
    .citation "[ABC3]" "isUnit_trivValue_res / nonVanishing_inf(段 D2・D3)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isUnit_trivValue_res") 6,
    .implicitStep
      ("★環準同型 A⁰_{x_i} → Γ(X, X_{s_i})(段 E1c)は本ファイルに無い。" ++
       "HomogeneousLocalization.Away.mk の同値関係を潰す段が要る" ++
       "——MvPolynomial は整域なので x_i は非零因子であり、消去できるはずである") 7 ]

end ABC3.Found.GenEll
