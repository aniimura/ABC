/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalChartSurjective
import ABC3.Found.GenEll.GluedTrivValue
import ABC3.Meta.Claim

/-!
# ★★★★★★★大域の比を**被覆で**同定する —— `§9-839` を大域版へ渡す道具（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★これは何か

`§9-843`（段 E3c 大域版）が要求するのは「試験元 `g` が `s_j/s_i` の形である」ことである。
★一方 `§9-839` が出すのは「**チャートごとの**等式 `g·trivValue(s^{⊗n}) = trivValue(t)`」である。

★★両者を繋ぐには「**チャートごとの等式から大域の比を同定する**」道具が要る。
それが本ファイルの `globalRatio_unique_of_cover` である
——`§9-841` の `globalRatio_unique` は**すべての**自明化つき開について要求していたが、
**被覆になっている部分族だけ**で足りる（構造層が層だから）。

## ★★★機構 —— 3 本

| 補題 | 役割 |
|---|---|
| `globalRatio_unique_of_cover` | ★被覆で一致すれば大域の比である（`eq_of_locally_eq'`） |
| `nonVanishing_unit_secPow` | ★★`X_{unit(s^{⊗n})} = X_s`（`§9-824` ＋ `§9-818`） |
| `basicOpen_pow_succ` | `D(f^{n+1}) = D(f)`（`n ≥ 1` のとき冪は基本開集合を変えない） |

## ★残っている段（明示）

★★★**開集合の同一視に沿った切断の輸送**が残る。
`§9-839` は `g ∈ Γ(X, X_s)` の言葉で条件を出し、`§9-843` は
`Γ(X, X_{unit(s^{⊗n})})` の言葉で `g` を要求する。★2 つの開は**等しい**
（`nonVanishing_unit_secPow`）が、**型は違う**ので `X.presheaf.map (eqToHom …)` で運ぶ段が要る。
★★そこは純粋な配管であり、数学は無い。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★冪は基本開集合を変えない（`n ≥ 1`） -/

/-- ★**`D(f^{n+1}) = D(f)`** —— 冪は基本開集合を変えない。

★`X.basicOpen_mul` と `inf_idem` の再帰だけである。 -/
theorem basicOpen_pow_succ (V : X.Opens) (f : (Γ(X, V) : Type)) :
    ∀ n : ℕ, X.basicOpen (f ^ (n + 1)) = X.basicOpen f
  | 0 => by rw [zero_add, pow_one]
  | n + 1 => by
      rw [pow_succ, X.basicOpen_mul, basicOpen_pow_succ V f n, inf_idem]

/-! ## ★★`X_{unit(s^{⊗n})} = X_s` -/

/-- ★★**層化した冪の非消失軌跡はもとの切断のそれである**。

★機構は `nonVanishing_sheafify`（`§9-824`）＋ `nonVanishing_secPow`（`§9-818`）である。
★★`n ≥ 1` が要る——`n = 0` では `X_{1} = ⊤` になってしまう。 -/
theorem nonVanishing_unit_secPow (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s : (M.obj (op ⊤) : Type)) (k : ℕ) :
    nonVanishing ((sheafifyFunctor X).obj (presheafTensorPow M (k + 1))).val
        (((sheafifyUnit X (presheafTensorPow M (k + 1))).app (op ⊤)).hom
          (secPow M s (k + 1)))
      = nonVanishing M s := by
  rw [nonVanishing_sheafify _ (isLocallyTrivial_presheafTensorPow hM (k + 1)),
    nonVanishing_secPow M hM s k]

/-! ## ★★★★★★★被覆で大域の比を同定する -/

/-- ★★★★★★★**被覆で一致すれば大域の比である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`§9-841` の `globalRatio_unique` は**すべての**自明化つき開について
`sectionRatio` に一致することを要求していた。
★本補題は**被覆になっている部分族だけ**で足りることを言う——構造層が層だからである。

★★これで `§9-839`（チャートごとの等式）を `§9-843`（大域の比の形）へ渡せる。
★★★機構は mathlib の `TopCat.Sheaf.eq_of_locally_eq'`（分離性）だけである。 -/
theorem globalRatio_unique_of_cover (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) {ι : Type} (W : ι → X.Opens)
    (eW : ∀ i, (restrictPresheafFunctor X (W i)).obj M ≅ 𝟙_ (PresheafModulesOn X (W i)))
    (hcov : nonVanishing M t ≤ ⨆ i, (nonVanishing M t ⊓ W i))
    (r : (Γ(X, nonVanishing M t) : Type))
    (hr : ∀ i, X.presheaf.map (homOfLE (inf_le_left :
        nonVanishing M t ⊓ W i ≤ nonVanishing M t)).op r
      = sectionRatio M (W i) (eW i) s t) :
    r = globalRatio M hM s t := by
  refine TopCat.Sheaf.eq_of_locally_eq' X.sheaf (fun i => nonVanishing M t ⊓ W i)
    (nonVanishing M t) (fun i => homOfLE inf_le_left) hcov r _ ?_
  intro i
  show X.presheaf.map (homOfLE _).op r = X.presheaf.map (homOfLE _).op (globalRatio M hM s t)
  rw [hr i, ← globalRatio_res M hM s t ⟨W i, eW i⟩]

/-! ## ★出典の紐付け(`.src`) -/

def basicOpen_pow_succ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(冪は基本開集合を変えない)",
    sectionId := "genell-prop-1-4" }

def nonVanishing_unit_secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(層化した冪の非消失軌跡はもとの切断のそれ)",
    sectionId := "genell-prop-1-4" }

def globalRatio_unique_of_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(被覆で一致すれば大域の比である)",
    sectionId := "genell-prop-1-4" }

def globalRatio_unique_of_cover.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "TopCat.Sheaf.eq_of_locally_eq(分離性)"
      (.inMathlib "TopCat.Sheaf.eq_of_locally_eq'") 2,
    .citation "[ABC3]" "globalRatio_res(大域の比はチャートで sectionRatio に戻る、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_res") 2,
    .citation "[ABC3]" "nonVanishing_sheafify(§9-824) / nonVanishing_secPow(§9-818)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.nonVanishing_sheafify") 2,
    .implicitStep
      ("★★★**開集合の同一視に沿った切断の輸送**が残る。§9-839 は g ∈ Γ(X, X_s) の言葉で条件を出し、" ++
       "§9-843 は Γ(X, X_{unit(s^{⊗n})}) の言葉で g を要求する。" ++
       "★2 つの開は**等しい**(nonVanishing_unit_secPow)が**型は違う**ので " ++
       "X.presheaf.map (eqToHom …) で運ぶ段が要る。★★そこは純粋な配管であり数学は無い") 5 ]

end ABC3.Found.GenEll
