/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonDegree
import ABC3.Meta.Claim

/-!
# ★★★★★アフィン開の上で分母を払う —— 段 E3a の局所版（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★これは段 E3a の局所版である

段 E3（`ψ : X ⟶ ℙᴺ_R` が immersion であること）の葉は

> **[Stacks] Lemma 01PW** —— `g ∈ Γ(X, X_s)` なら、ある `n` で `g·s^n` が
> 大域切断の制限になる

である（2026-08-28 の測定）。★その**局所版**——「アフィン開 `U` の上では
`g ∈ Γ(X, X_s ⊓ U)` の分母を `trivValue s` の冪で払える」——が本ファイルである。

## ★★★★機構 —— mathlib の局所化がそのまま効く

`X_s ⊓ U = X.basicOpen (trivValue M U e s)`（`nonVanishing_inf`、段 D2）であり、
mathlib の

    `IsAffineOpen.isLocalization_basicOpen : IsLocalization.Away f Γ(X, X.basicOpen f)`

が `Γ(X, X_s ⊓ U)` を `Γ(X, U)` の `f = trivValue s` による局所化だと言う。
★あとは `IsLocalization.Away.surj` を呼ぶだけである。

## ★残っている段（明示）

★★**大域化は本ファイルに無い**。要るのは:

1. `X` を有限個のアフィン開 `U_j`（`M` が自明化する）で覆う（段 E2 の材料が使える）
2. 各 `j` で本ファイルの補題を当て、指数 `n_j` を**最大値で揃える**
3. 得られた `a_j ∈ Γ(X, U_j)` を貼り合わせて `Γ(X, M^{⊗?})` の大域切断にする

★★★3 が本体である——重なりで一致することを見る必要があり、
そこでは `trivValue` の遷移単元（段 D1）が効く。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-- ★★★★★**アフィン開の上では分母を `f` の冪で払える**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★mathlib の `IsAffineOpen.isLocalization_basicOpen`（`Γ(X, D(f))` は `Γ(X,U)` の
`f` による局所化）＋ `IsLocalization.Away.surj` そのものである。 -/
theorem exists_pow_mul_eq_res (U : X.Opens) (hU : IsAffineOpen U) (f : Γ(X, U))
    (g : Γ(X, X.basicOpen f)) :
    ∃ (n : ℕ) (a : Γ(X, U)),
      g * (algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ n
        = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a := by
  haveI := hU.isLocalization_basicOpen f
  exact IsLocalization.Away.surj f g

/-- ★★★★★★**切断の非消失軌跡の側で読んだ形**。

`X_s ⊓ U = X.basicOpen (trivValue M U e s)`（`nonVanishing_inf`、段 D2）なので、
★`X_s ⊓ U` の上の関数は `trivValue s` の冪を掛ければ `U` へ延びる。 -/
theorem exists_pow_mul_eq_res_nonVanishing (M : X.PresheafOfModules) (U : X.Opens)
    (hU : IsAffineOpen U)
    (e : (restrictPresheafFunctor X U).obj M ≅ 𝟙_ (PresheafModulesOn X U))
    (s : M.obj (op ⊤)) (g : Γ(X, X.basicOpen (trivValue M U e s))) :
    ∃ (n : ℕ) (a : Γ(X, U)),
      g * (algebraMap (Γ(X, U) : Type)
            (Γ(X, X.basicOpen (trivValue M U e s)) : Type) (trivValue M U e s)) ^ n
        = algebraMap (Γ(X, U) : Type)
            (Γ(X, X.basicOpen (trivValue M U e s)) : Type) a :=
  exists_pow_mul_eq_res U hU (trivValue M U e s) g

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_mul_eq_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(アフィン開の上では分母を f の冪で払える)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_res_nonVanishing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X_s ⊓ U の上の関数は trivValue s の冪で U へ延びる)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_res.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsAffineOpen.isLocalization_basicOpen / IsLocalization.Away.surj"
      (.inMathlib "AlgebraicGeometry.IsAffineOpen.isLocalization_basicOpen") 7,
    .citation "[ABC3]" "nonVanishing_inf(X_s ⊓ U = basicOpen (trivValue)、段 D2)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_inf") 6,
    .citation "[Stacks]" "Lemma 01PW(ample な可逆層の切断の延長——本ファイルはその局所版)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)。段 E3a の大域版は未着手") 7,
    .implicitStep
      ("★**大域化は本ファイルに無い**。要るのは (1) X を有限個の自明化するアフィン開で覆う" ++
       "(段 E2 の材料が使える)、(2) 各々で本補題を当て指数を最大値で揃える、" ++
       "(3) 得られた a_j を貼り合わせて大域切断にする、の 3 段である。" ++
       "★★3 が本体で、重なりでの一致に trivValue の遷移単元(段 D1)が効く") 7 ]

end ABC3.Found.GenEll
