/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DivisorOfSection
import Mathlib.RingTheory.MvPolynomial.Homogeneous
import ABC3.Meta.Claim

/-!
# ★★★★★★**斉次多項式を切断で評価する** —— 段 E1a（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは台帳の段 E1a である

`ResearchPaper/mathlib-gap.json` の `ample-and-projective-embedding` は
段 E1（大域切断が定める `ℙᴺ` への射）を 4 つに割った（2026-08-28）:

| 段 | 内容 |
|---|---|
| **E1a** | **斉次多項式を切断で評価する** ← ★本ファイル |
| E1b | `a(s)/s_i^n` が `Γ(X, X_{s_i})` の元を定める（`sectionRatio` の次数 `n` 版） |
| E1c | 環準同型 `A⁰_{x_i} → Γ(X, X_{s_i})` の well-defined 性 |
| E1d | 貼り合わせて `X ⟶ ℙᴺ_R` |

★mathlib の `HomogeneousLocalization.Away.mk_surjective` が `A⁰_{x_i}` の元を
**`a / x_i^n`（`a` は `n` 次斉次）**に落とすので、環準同型を作るには
まず「`a` を切断 `s` で評価する」段が要る。

## ★★★★★★機構 —— 自明化を替えると `u^n` 倍になる

自明化 `e` のもとで切断 `s_j` は関数 `t_j = trivValue M V e s_j` になる。
★自明化を `e'` に替えると**すべての `t_j` が同じ単元 `u` 倍**される
（`trivValue_congr'`——`u` は切断に依らない）。

★★したがって `n` 次斉次式 `a` の値は

    `a(t·u) = a(t) · u^n`

となる（`eval₂_mul_const_of_isHomogeneous`）。★★★これが
「比 `a(s)/s_i^n` が自明化に依らない」ことの根拠である
——分子が `u^n` 倍、分母も `u^n` 倍だからである。

## ★測定の記録

mathlib に「斉次式の評価はスカラー倍で `u^n` 倍になる」補題は**無い**
（2026-08-28 実測——`Mathlib/RingTheory/MvPolynomial/Homogeneous.lean` の
`IsHomogeneous.eval₂` は**次数の合成**であって値のスカラー則ではない）。
★本ファイルで入れた（`eval₂_mul_const_of_isHomogeneous`）。

## ★残っている段（明示）

★★**比そのもの（E1b）は本ファイルに無い**。`sectionRatio`（段 D3）と同じ形で

    `homogRatio ≝ (a(t) の制限) * ((t_i の制限)⁻¹)^n : Γ(X, X_{s_i} ⊓ V)`

と置けばよく、自明化に依らないことは本ファイルの `homogValue_congr` から出る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

/-! ## ★★★★★(1) 斉次式の評価のスカラー則 -/

open MvPolynomial in
/-- ★★★★★**`n` 次斉次式の評価は変数を `u` 倍すると `u^n` 倍になる**。

    `a(x·u) = a(x) · u^n`

★mathlib に無い（2026-08-28 実測）。★★機構は
「`d` 次単項式は `∏ (x_i u)^{d_i} = (∏ x_i^{d_i}) · u^{deg d}`」＋
「斉次なら台の各 `d` で `deg d = n`」だけである。 -/
theorem eval₂_mul_const_of_isHomogeneous {σ R S : Type} [CommRing R] [CommRing S]
    (φ : R →+* S) (x : σ → S) (u : S) {n : ℕ} {a : MvPolynomial σ R}
    (ha : a.IsHomogeneous n) :
    eval₂ φ (fun j => x j * u) a = (eval₂ φ x a) * u ^ n := by
  rw [eval₂_eq, eval₂_eq, Finset.sum_mul]
  refine Finset.sum_congr rfl (fun d hd => ?_)
  have hdeg : Finsupp.degree d = n := by
    by_contra h
    exact (MvPolynomial.mem_support_iff.1 hd) (ha.coeff_eq_zero h)
  have hprod : (∏ i ∈ d.support, (x i * u) ^ d i)
      = (∏ i ∈ d.support, x i ^ d i) * u ^ n := by
    rw [Finset.prod_congr rfl (fun i _ => mul_pow (x i) u (d i)), Finset.prod_mul_distrib,
      Finset.prod_pow_eq_pow_sum, ← hdeg]
    rfl
  rw [hprod]
  ring

/-! ## ★★★★★★(2) 切断での評価 -/

variable {X : Scheme.{0}}

/-- ★★★★★★**斉次多項式を切断で評価する**（自明化 `e` のもとで）。

★係数環 `R` から `Γ(X, V)` への環準同型 `φ` を受ける
——`X` が `Spec R` 上のとき、構造射がそれを与える。 -/
noncomputable def homogValue (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (a : MvPolynomial σ R) : Γ(X, V) :=
  MvPolynomial.eval₂ φ (fun j => trivValue M V e (s j)) a

/-- ★★★★★★★**自明化を替えると `n` 次斉次式の値は `u^n` 倍になる**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★`u` は**切断に依らない**（`trivValue_congr'`）のが要点である
——だから比 `a(s)/s_i^n` を取ると `u^n` が分子と分母で消える。
★★★これが「射影埋め込みの座標が自明化に依らない」ことの中身である。 -/
theorem homogValue_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) :
    ∃ u : Γ(X, V), IsUnit u ∧ ∀ {n : ℕ} {a : MvPolynomial σ R}, a.IsHomogeneous n →
      homogValue M V e' φ s a = homogValue M V e φ s a * u ^ n := by
  obtain ⟨u, hu, hval⟩ := trivValue_congr' M V e e'
  refine ⟨u, hu, fun {n} {a} ha => ?_⟩
  show MvPolynomial.eval₂ φ (fun j => trivValue M V e' (s j)) a = _
  have hx : (fun j => trivValue M V e' (s j))
      = (fun j => trivValue M V e (s j) * u) := funext (fun j => hval (s j))
  rw [hx]
  exact eval₂_mul_const_of_isHomogeneous φ _ u ha

/-- ★★次数 `1` の場合は `trivValue` そのものである（健全性）。 -/
theorem homogValue_X (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    {R : Type} [CommRing R] (φ : R →+* (Γ(X, V) : Type))
    {σ : Type} (s : σ → (M.obj (op ⊤) : Type)) (j : σ) :
    homogValue M V e φ s (MvPolynomial.X j) = trivValue M V e (s j) := by
  show MvPolynomial.eval₂ φ (fun j => trivValue M V e (s j)) (MvPolynomial.X j) = _
  rw [MvPolynomial.eval₂_X]

/-! ## ★出典の紐付け(`.src`) -/

def eval₂_mul_const_of_isHomogeneous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(n 次斉次式の評価は変数を u 倍すると u^n 倍になる)",
    sectionId := "genell-prop-1-4" }

def homogValue.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——斉次多項式を切断で評価する)",
    sectionId := "genell-prop-1-4" }

def homogValue_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(自明化を替えると n 次斉次式の値は u^n 倍になる)",
    sectionId := "genell-prop-1-4" }

def homogValue_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_congr'(単元 u は切断に依らない、段 D1)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_congr'") 6,
    .citation "[mathlib]" "MvPolynomial.eval₂_eq / IsHomogeneous.coeff_eq_zero"
      (.inMathlib "MvPolynomial.eval₂_eq") 7,
    .implicitStep
      ("★mathlib に「斉次式の評価はスカラー倍で u^n 倍になる」補題は**無い**" ++
       "(2026-08-28 実測——IsHomogeneous.eval₂ は次数の合成であって値のスカラー則ではない)。" ++
       "本ファイルで入れた") 7,
    .implicitStep
      ("★★比そのもの(段 E1b)は本ファイルに無い。sectionRatio(段 D3)と同じ形で " ++
       "homogRatio ≝ (a(t) の制限) * ((t_i の制限)⁻¹)^n と置けばよく、" ++
       "自明化に依らないことは homogValue_congr から出る") 7 ]

end ABC3.Found.GenEll
