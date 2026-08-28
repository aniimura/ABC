/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalizeStep
import ABC3.Meta.Claim

/-!
# ★★★★★★★★重なりの一致を実際に作る —— 段 E3a-3 の配管（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 抽象形を具体の状況へ落とす

`§9-829` は「単元倍のずれは冪で吸収できる」を**抽象形**で持っている。
★本ファイルはそれを**チャートの言葉**（`trivValue`・遷移単元）へ落としたものである:

> `X_s` の上で一致する 2 つの関数 `g_i, g_j` から分母を払って得た `a_i, a_j` は、
> **ある `m` で** `f_i|^m · (a_i| · u^N) = f_i|^m · a_j|` を満たす。

★★これが `§9-828` の判定基準（`a| · u^N = a'|`）を**実際に満たす形へ持っていく**段である。

## ★★★★機構 —— 3 つの同定

| 何を | どう |
|---|---|
| `D := X.basicOpen (f_i\|_W)` は `D(f_i)` にも `D(f_j)` にも入る | `basicOpen_res` ＋ `basicOpen_res_trivValue_eq`（§9-829） |
| `D` の上で `g_i` と `g_j` は同じ関数 `g` になる | 仮定 `hg` を `D ≤ D(f_i) ⊓ D(f_j)` へ制限 |
| `D` の上で `a_j` は `g·(f_i\|·u)^N` に見える | `trivValue_res_transUnit`（§9-829） |

★あとは `§9-829` の抽象形に渡すだけである。

## ★測定の記録

★`W = U_i ⊓ U_j` が**アフィンであること**を仮定に置いた。
★★これは `X` が分離的なら自動である（アフィン開の交わりはアフィン）が、
**そこは本ファイルの担当ではない**——消費側が `IsSeparated X` から供給する。

## ★残っている段（明示）

★★**最終の組み上げは本ファイルに無い**。残るのは:

1. 有限被覆の全ペアで `m` を揃える（`exists_common_pow`、§9-826）
2. 指数を `N → N + m` に上げ、`a_i → a_i · f_i^m` と取り替える
3. `§9-828` の判定基準に通す
4. `sheafifyUnit` で層化の側へ送り、`§9-827` の `exists_unique_glue_top` で貼る
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-- ★★★★★★★★**重なりの一致を実際に作る** —— 段 E3a-3 の配管。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

仮定は次のとおり（`f_i = trivValue M V_i e_i s` と書く）:

* `hai` —— `D(f_i)` の上で `g_i · f_i^N = a_i`（分母を払った、§9-822）
* `haj` —— `D(f_j)` の上で `g_j · f_j^N = a_j`
* `hg`  —— `g_i` と `g_j` は `D(f_i) ⊓ D(f_j)` で一致する（同じ `X_s` 上の関数から来た）

結論は「**ある `m` で** `f_i|^m · (a_i| · u^N) = f_i|^m · a_j|`」である。

★★これで `§9-828` の判定基準が（指数を `N+m` に上げれば）満たせる形になる。 -/
theorem exists_pow_mul_eq_overlap (M : X.PresheafOfModules) (Vi Vj : X.Opens)
    (hW : IsAffineOpen (Vi ⊓ Vj))
    (ei : (restrictPresheafFunctor X Vi).obj M ≅ 𝟙_ (PresheafModulesOn X Vi))
    (ej : (restrictPresheafFunctor X Vj).obj M ≅ 𝟙_ (PresheafModulesOn X Vj))
    (s : (M.obj (op ⊤) : Type)) (N : ℕ)
    (gi : (Γ(X, X.basicOpen (trivValue M Vi ei s)) : Type))
    (gj : (Γ(X, X.basicOpen (trivValue M Vj ej s)) : Type))
    (ai : (Γ(X, Vi) : Type)) (aj : (Γ(X, Vj) : Type))
    (hai : gi * (X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M Vi ei s))).op
              (trivValue M Vi ei s)) ^ N
        = X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M Vi ei s))).op ai)
    (haj : gj * (X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M Vj ej s))).op
              (trivValue M Vj ej s)) ^ N
        = X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M Vj ej s))).op aj)
    (hg : X.presheaf.map (homOfLE (inf_le_left :
            X.basicOpen (trivValue M Vi ei s) ⊓ X.basicOpen (trivValue M Vj ej s)
              ≤ X.basicOpen (trivValue M Vi ei s))).op gi
        = X.presheaf.map (homOfLE (inf_le_right :
            X.basicOpen (trivValue M Vi ei s) ⊓ X.basicOpen (trivValue M Vj ej s)
              ≤ X.basicOpen (trivValue M Vj ej s))).op gj) :
    ∃ m : ℕ,
      (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op (trivValue M Vi ei s)) ^ m
          * (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op ai
            * (transUnit M (Vi ⊓ Vj)
                (trivialOfLe M (inf_le_left : Vi ⊓ Vj ≤ Vi) ei)
                (trivialOfLe M (inf_le_right : Vi ⊓ Vj ≤ Vj) ej)) ^ N)
        = (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op (trivValue M Vi ei s)) ^ m
          * X.presheaf.map (homOfLE (inf_le_right : Vi ⊓ Vj ≤ Vj)).op aj := by
  have hD := basicOpen_res_trivValue_eq M Vi Vj ei ej s
  have hDi : X.basicOpen (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op
      (trivValue M Vi ei s)) ≤ X.basicOpen (trivValue M Vi ei s) := by
    rw [Scheme.basicOpen_res X (trivValue M Vi ei s)
      (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op]
    exact inf_le_right
  have hDj : X.basicOpen (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op
      (trivValue M Vi ei s)) ≤ X.basicOpen (trivValue M Vj ej s) := by
    rw [hD, Scheme.basicOpen_res X (trivValue M Vj ej s)
      (homOfLE (inf_le_right : Vi ⊓ Vj ≤ Vj)).op]
    exact inf_le_right
  refine exists_pow_mul_eq_of_unit_scaling (Vi ⊓ Vj) hW
    (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op (trivValue M Vi ei s))
    (transUnit M (Vi ⊓ Vj)
      (trivialOfLe M (inf_le_left : Vi ⊓ Vj ≤ Vi) ei)
      (trivialOfLe M (inf_le_right : Vi ⊓ Vj ≤ Vj) ej)) N
    (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op ai)
    (X.presheaf.map (homOfLE (inf_le_right : Vi ⊓ Vj ≤ Vj)).op aj)
    (X.presheaf.map (homOfLE hDi).op gi) ?_ ?_
  · -- ★`D` の上で `a_i` は `g · f_i^N` に見える
    have h := congrArg (fun z => X.presheaf.map (homOfLE hDi).op z) hai
    simp only [map_mul, map_pow] at h
    show X.presheaf.map (homOfLE (X.basicOpen_le _)).op
      (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op ai) = _
    rw [res_trans]
    rw [res_trans, res_trans] at h
    rw [← h, algebraMap_basicOpen_eq_res, res_trans]
  · -- ★★`D` の上で `a_j` は `g · (f_i·u)^N` に見える
    have h := congrArg (fun z => X.presheaf.map (homOfLE hDj).op z) haj
    simp only [map_mul, map_pow] at h
    rw [res_trans, res_trans] at h
    have hgg : X.presheaf.map (homOfLE hDj).op gj = X.presheaf.map (homOfLE hDi).op gi := by
      have h2 := congrArg (fun z => X.presheaf.map
        (homOfLE (le_inf hDi hDj :
          X.basicOpen (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op
            (trivValue M Vi ei s))
          ≤ X.basicOpen (trivValue M Vi ei s) ⊓ X.basicOpen (trivValue M Vj ej s))).op z) hg
      simp only at h2
      rw [res_trans, res_trans] at h2
      exact h2.symm
    rw [hgg] at h
    show X.presheaf.map (homOfLE (X.basicOpen_le _)).op
      (X.presheaf.map (homOfLE (inf_le_right : Vi ⊓ Vj ≤ Vj)).op aj) = _
    rw [res_trans, ← h, algebraMap_basicOpen_eq_res,
      ← trivValue_res_transUnit M Vi Vj ei ej s, res_trans]

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_mul_eq_overlap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの一致を実際に作る——チャートの言葉)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_overlap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pow_mul_eq_of_unit_scaling(単元倍のずれは冪で吸収、§9-829)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_mul_eq_of_unit_scaling") 3,
    .citation "[ABC3]" "basicOpen_res_trivValue_eq / trivValue_res_transUnit(§9-829)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_res_transUnit") 3,
    .citation "[mathlib]" "Scheme.basicOpen_res(制限した切断の非消失軌跡)"
      (.inMathlib "AlgebraicGeometry.Scheme.basicOpen_res") 3,
    .citation "[Stacks]" "Lemma 01PW(の大域化の段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★W = U_i ⊓ U_j が**アフィンであること**を仮定に置いた。" ++
       "★★X が分離的なら自動である(アフィン開の交わりはアフィン)が、" ++
       "そこは本ファイルの担当ではない——消費側が IsSeparated X から供給する") 6,
    .implicitStep
      ("★★★**最終の組み上げは本ファイルに無い**。残るのは " ++
       "(1) 全ペアで m を揃える(§9-826)、(2) 指数を N → N+m に上げ a_i → a_i·f_i^m と取り替える、" ++
       "(3) §9-828 の判定基準に通す、(4) sheafifyUnit で層化の側へ送り §9-827 で貼る") 8 ]

end ABC3.Found.GenEll
