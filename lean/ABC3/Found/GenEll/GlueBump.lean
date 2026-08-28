/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GluedTrivValue
import ABC3.Meta.Claim

/-!
# ★★★★★★★★指数揃え —— 貼った大域切断の冪を上げる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 段 E3c の最後の帳簿

`§9-839` は「試験元 `g` ごとに `n` と大域切断 `t` を作る」まで来た。
★しかし**`g` ごとに `n` が違う**——座標族は 1 つの `n` で揃っていなければならない。

★★本ファイルがその**指数揃え**である:

> 指数 `n` で貼り合う族 `a_i` は、`a_i · f_i^k` として**指数 `n+k` でも貼り合う**。
> しかも新しい大域切断の座標は `a_i · f_i^k` である。

## ★★★機構 —— `§9-831` の `overlap_bump` がそのまま効く

一致条件は `res(a_i)·u^n = res(a_j)` である。`a_i → a_i·f_i^k` と替えると、
`f_j| = f_i|·u`（`§9-829`）なので

    `res(a_i)·f_i^k·u^{n+k} = res(a_j)·f_i^k·u^k = res(a_j)·f_j^k`

★すなわち `u^k` が両辺に現れて**ちょうど相殺する**——`overlap_bump`（`§9-831`）そのものである。

## ★★これで段 E3c の帳簿が閉じる

1. 試験元は有限個（`§9-832`）、チャートも有限個（`§9-817`）
2. 各試験元に `§9-839` を当てて `(n_g, t_g)` を得る
3. ★**本ファイルで `n := max n_g` に揃える**
4. `t_g` たちを座標族に並べ、`§9-838` に渡す

★★★残るのは 1・2・4 の**並べ替えの帳簿だけ**である（数学は無い）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★★指数を上げても貼り合わせは保たれる -/

/-- ★★★★★★★**指数を上げても貼り合わせは保たれる**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

指数 `n` で貼り合う族 `a_i` は、`a_i · f_i^k` として**指数 `n+k` でも貼り合う**。

★機構は `overlap_bump`（`§9-831`）そのものである
——`f_j| = f_i|·u` なので `u^k` が両辺に現れてちょうど相殺する。 -/
theorem exists_glue_bump {ι : Type} (M : X.PresheafOfModules)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (n k : ℕ) (a : ∀ i, (Γ(X, U i) : Type))
    (hagree : ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) :
    ∃! t : (((sheafifyFunctor X).obj (presheafTensorPow M (n + k))).val.obj
        (op (⊤ : X.Opens)) : Type),
      ∀ i, ((sheafifyFunctor X).obj (presheafTensorPow M (n + k))).val.map
          (homOfLE (le_top : U i ≤ ⊤)).op t
        = ((sheafifyUnit X (presheafTensorPow M (n + k))).app (op (U i))).hom
            (secOfFun M (U i) (e i) (n + k)
              (a i * (trivValue M (U i) (e i) s) ^ k)) := by
  refine exists_glue_of_agree M U hcov (n + k) e
    (fun i => a i * (trivValue M (U i) (e i) s) ^ k) ?_
  intro i j
  simp only [map_mul, map_pow]
  rw [trivValue_res_transUnit M (U i) (U j) (e i) (e j) s]
  refine overlap_bump _ _ _ _ n k ?_
  rw [hagree i j]

/-- ★★★★★★★★**指数を上げた大域切断の座標は `a_i · f_i^k` である**。

★これで「`n` を揃えたあとも座標が読める」——`§9-838` の座標条件が扱える形のまま残る。 -/
theorem exists_glue_bump_trivValue {ι : Type} (M : X.PresheafOfModules)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (n k : ℕ) (a : ∀ i, (Γ(X, U i) : Type))
    (hagree : ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) :
    ∃ t : (((sheafifyFunctor X).obj (presheafTensorPow M (n + k))).val.obj
        (op (⊤ : X.Opens)) : Type), ∀ i,
      trivValue ((sheafifyFunctor X).obj (presheafTensorPow M (n + k))).val (U i)
          (sheafifyTriv (presheafTensorPow M (n + k)) (tensorPowTriv (e i) (n + k))) t
        = a i * (trivValue M (U i) (e i) s) ^ k := by
  obtain ⟨t, ht, -⟩ := exists_glue_bump M U hcov e s n k a hagree
  exact ⟨t, fun i => trivValue_glued M (U i) (n + k) (e i) _ t (ht i)⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_glue_bump.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(指数を上げても貼り合わせは保たれる)",
    sectionId := "genell-prop-1-4" }

def exists_glue_bump_trivValue.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(指数を上げた大域切断の座標は a_i · f_i^k)",
    sectionId := "genell-prop-1-4" }

def exists_glue_bump.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "overlap_bump(単元倍のずれの吸収、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.overlap_bump") 2,
    .citation "[ABC3]" "trivValue_res_transUnit(隣のチャートの座標は遷移単元で移る、§9-829)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_res_transUnit") 2,
    .citation "[ABC3]" "exists_glue_of_agree(層化の側で貼る、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_glue_of_agree") 2,
    .implicitStep
      ("★★これで段 E3c の帳簿が閉じる: (1) 試験元は有限個(§9-832)、チャートも有限個(§9-817)、" ++
       "(2) 各試験元に §9-839 を当てて (n_g, t_g) を得る、" ++
       "(3) **本ファイルで n := max n_g に揃える**、(4) t_g たちを座標族に並べて §9-838 に渡す。" ++
       "★★★残るのは 1・2・4 の**並べ替えの帳簿だけ**である(数学は無い)") 6 ]

end ABC3.Found.GenEll
