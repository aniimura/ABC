/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalizeOverlap
import ABC3.Found.GenEll.SheafifyGlue
import ABC3.Found.GenEll.OverlapCriterion
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 E3a-3 —— 大域化が閉じた（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— [Stacks] Lemma 01PW の大域版

> `X_s` の上の関数（チャートごとに与えられ、重なりで一致するもの）は、
> **`M` を十分高いテンソル冪へ上げれば** `X` 全体の切断に延びる。

★これが原典の「**[some positive tensor power of]**」という括弧の中身である。

## ★★★★★★★★組み立て —— 7 本の合流

| 段 | 出典 | 役割 |
|---|---|---|
| 分母を払う | `§9-822` `exists_pow_mul_eq_res` | アフィン開の上で `g·f^n` を延ばす |
| 指数を揃える | `§9-826` `exists_common_pow` | `Finset.sup` で有限個まとめる |
| 冪と座標 | `§9-825` `trivValue_secPow` | `(s^{⊗n})` の座標は `(trivValue s)^n` |
| 遷移単元の `N` 乗 | `§9-828` `transUnit_tensorPowTriv` | `M^{⊗N}` の遷移函数は `N` 乗 |
| 一致判定 | `§9-828` `res_secOfFun_eq_iff` | 幾何の一致 ⟺ `a\|·u^N = a'\|` |
| 単元倍の吸収 | `§9-829`・`§9-830` `exists_pow_mul_eq_overlap` | ずれを冪で潰す |
| 貼り合わせ | `§9-827` `exists_unique_glue_top` | 層化の側で貼る（前層では貼れない） |

★★本ファイルはこの 7 本を **3 つの定理**に組む:

1. `exists_common_exponent` —— 分母と重なりの**両方**に効く単一の指数 `n` を取る
2. `exists_glue_of_agree` —— 一致する族を層化の側で貼る
3. `exists_global_section_of_localData` —— ★★★**段 E3a-3 の本体**

## ★★★測定 —— 原典の括弧には冪が 3 回畳まれている

`[some positive tensor power of]` の中身を型で開くと、**冪を取る理由が 3 つ**あった:

* (a) **分母を払う冪**（`§9-822`）—— `X_s` 上の関数を `U` へ延ばすため
* (b) **重なりを合わせる冪**（`§9-826`・`§9-830`）—— `X_s` の外で一致させるため
* (c) **遷移函数の `N` 乗**（`§9-828`）—— `M^{⊗N}` の遷移が `N` 乗になること

★合図 1 つを節点 1 つに対応させる規律で数えると、この括弧は**節点 3 個**であった。

## ★逸脱の記録

★`U i ⊓ U j` が**アフィンであること**を仮定に置いた（`hUij`）。
★★`X` が分離的なら自動だが、本ファイルはそれを導かない——**消費側が供給する**。
★★★また被覆は `Fintype ι` で表した（`ample` から来る有限被覆はこの形である、`§9-817`）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★底上げ —— 指数を上げても形は保たれる -/

/-- ★**分母側の底上げ**——指数を `N → N+m` に上げると分子が `a·f^m` になる。 -/
theorem res_pow_bump (U : X.Opens) (f : (Γ(X, U) : Type))
    (g : (Γ(X, X.basicOpen f) : Type)) (a : (Γ(X, U) : Type)) (N m : ℕ)
    (h : g * (algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ N
      = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a) :
    g * (algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ (N + m)
      = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) (a * f ^ m) := by
  rw [pow_add, ← mul_assoc, h, map_mul, map_pow]

/-- ★★**重なり側の底上げ**——`f^m` を掛けた形が、そのまま指数 `N+m` の一致条件になる。

★★★ここが「1 回冪を上げれば重なりも揃う」の心臓である
——`f_j = f_i·u` なので `a_j·f_j^m = a_j·f_i^m·u^m` となり、`u^m` が両辺に現れて相殺する。 -/
theorem overlap_bump {W : X.Opens} (fi u ai aj : (Γ(X, W) : Type)) (N m : ℕ)
    (hm : fi ^ m * (ai * u ^ N) = fi ^ m * aj) :
    (ai * fi ^ m) * u ^ (N + m) = aj * (fi * u) ^ m := by
  have h2 : u ^ m * (fi ^ m * (ai * u ^ N)) = u ^ m * (fi ^ m * aj) := by rw [hm]
  rw [pow_add, mul_pow]
  linear_combination h2

/-! ## ★★★★★★★★分母と重なりの両方に効く単一の指数 -/

open scoped Classical in
/-- ★★★★★★★★**分母と重なりの両方に効く単一の指数を取る**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

チャートごとの関数 `g i`（`X_s ⊓ U i` の上）が重なりで一致するなら、
★**ある `n` と `a i ∈ Γ(X, U i)` があって**

* `g i · f_i^n = a_i`（`X_s ⊓ U i` の上）——分母が払えている
* `a_i\|_{U_i⊓U_j} · u_{ij}^n = a_j\|_{U_i⊓U_j}`——重なりで一致する（`§9-828` の判定基準）

★★機構は「分母を揃える（`§9-829`）→ 重なりのずれを潰す `m` を全ペアで揃える（`§9-826`）→
`n := N + m` に底上げする（`res_pow_bump` / `overlap_bump`）」の 3 手である。 -/
theorem exists_common_exponent {ι : Type} (I : Finset ι)
    (M : X.PresheafOfModules) (U : ι → X.Opens)
    (hU : ∀ i ∈ I, IsAffineOpen (U i))
    (hUij : ∀ i ∈ I, ∀ j ∈ I, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type))
    (g : ∀ i, (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type))
    (hg : ∀ i ∈ I, ∀ j ∈ I,
      X.presheaf.map (homOfLE (inf_le_left :
        X.basicOpen (trivValue M (U i) (e i) s) ⊓ X.basicOpen (trivValue M (U j) (e j) s)
          ≤ X.basicOpen (trivValue M (U i) (e i) s))).op (g i)
        = X.presheaf.map (homOfLE (inf_le_right :
            X.basicOpen (trivValue M (U i) (e i) s) ⊓ X.basicOpen (trivValue M (U j) (e j) s)
              ≤ X.basicOpen (trivValue M (U j) (e j) s))).op (g j)) :
    ∃ (n : ℕ) (a : ∀ i, (Γ(X, U i) : Type)),
      (∀ i ∈ I, g i * (algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
            (trivValue M (U i) (e i) s)) ^ n
          = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a i)) ∧
      (∀ i ∈ I, ∀ j ∈ I,
        X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
            * (transUnit M (U i ⊓ U j)
                (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
                (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
          = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) := by
  obtain ⟨N, hN⟩ := exists_common_denominator I U hU
    (fun i => trivValue M (U i) (e i) s) g
  choose! a ha using hN
  have hpair : ∀ p ∈ I ×ˢ I, ∃ n : ℕ,
      (X.presheaf.map (homOfLE (inf_le_left : U p.1 ⊓ U p.2 ≤ U p.1)).op
        (trivValue M (U p.1) (e p.1) s)) ^ n
        * (X.presheaf.map (homOfLE (inf_le_left : U p.1 ⊓ U p.2 ≤ U p.1)).op (a p.1)
          * (transUnit M (U p.1 ⊓ U p.2)
              (trivialOfLe M (inf_le_left : U p.1 ⊓ U p.2 ≤ U p.1) (e p.1))
              (trivialOfLe M (inf_le_right : U p.1 ⊓ U p.2 ≤ U p.2) (e p.2))) ^ N)
      = (X.presheaf.map (homOfLE (inf_le_left : U p.1 ⊓ U p.2 ≤ U p.1)).op
          (trivValue M (U p.1) (e p.1) s)) ^ n
        * X.presheaf.map (homOfLE (inf_le_right : U p.1 ⊓ U p.2 ≤ U p.2)).op (a p.2) := by
    intro p hp
    rw [Finset.mem_product] at hp
    exact exists_pow_mul_eq_overlap M (U p.1) (U p.2) (hUij p.1 hp.1 p.2 hp.2)
      (e p.1) (e p.2) s N (g p.1) (g p.2) (a p.1) (a p.2)
      (ha p.1 hp.1) (ha p.2 hp.2) (hg p.1 hp.1 p.2 hp.2)
  obtain ⟨m, hm⟩ := exists_common_pow (I ×ˢ I) _
    (fun p n₁ n₂ hle h => pow_mul_eq_of_le (U p.1 ⊓ U p.2) _ _ _ h hle) hpair
  refine ⟨N + m, fun i => a i * (trivValue M (U i) (e i) s) ^ m, ?_, ?_⟩
  · intro i hi
    exact res_pow_bump (U i) _ (g i) (a i) N m (ha i hi)
  · intro i hi j hj
    have h := hm (i, j) (Finset.mem_product.2 ⟨hi, hj⟩)
    simp only [map_mul, map_pow]
    rw [trivValue_res_transUnit M (U i) (U j) (e i) (e j) s]
    exact overlap_bump _ _ _ _ N m h

/-! ## ★★★★★★層化の側へ送って貼る -/

/-- ★**層化の単位は制限と可換である**（単位の自然性）。

★★これで「前層で一致していれば層化の側でも一致する」が言える。 -/
theorem sheafifyUnit_res (P : X.PresheafOfModules) {V W : X.Opens} (h : W ≤ V)
    (x : (P.obj (op V) : Type)) :
    ((sheafifyUnit X P).app (op W)).hom (P.map (homOfLE h).op x)
      = ((sheafifyFunctor X).obj P).val.map (homOfLE h).op
          (((sheafifyUnit X P).app (op V)).hom x) :=
  PresheafOfModules.naturality_apply (sheafifyUnit X P) (homOfLE h).op x

/-- ★★★★★★**一致する族は層化の側で貼り合う**。

★機構は `§9-828` の判定基準（幾何の一致 ⟺ `a|·u^n = a'|`）＋
単位の自然性＋`§9-827` の `exists_unique_glue_top` だけである。 -/
theorem exists_glue_of_agree {ι : Type} (M : X.PresheafOfModules)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤) (n : ℕ)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (a : ∀ i, (Γ(X, U i) : Type))
    (hagree : ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) :
    ∃! t : (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.obj
        (op (⊤ : X.Opens)) : Type),
      ∀ i, ((sheafifyFunctor X).obj (presheafTensorPow M n)).val.map
          (homOfLE (le_top : U i ≤ ⊤)).op t
        = ((sheafifyUnit X (presheafTensorPow M n)).app (op (U i))).hom
            (secOfFun M (U i) (e i) n (a i)) := by
  refine exists_unique_glue_top (presheafTensorPow M n) U hcov
    (fun i => ((sheafifyUnit X (presheafTensorPow M n)).app (op (U i))).hom
      (secOfFun M (U i) (e i) n (a i))) ?_
  intro i j
  rw [← sheafifyUnit_res (presheafTensorPow M n) (inf_le_left : U i ⊓ U j ≤ U i),
    ← sheafifyUnit_res (presheafTensorPow M n) (inf_le_right : U i ⊓ U j ≤ U j)]
  congr 1
  exact (res_secOfFun_eq_iff (e i) (e j) n (a i) (a j)).2 (hagree i j)

/-! ## ★★★★★★★★★★段 E3a-3 の本体 -/

/-- ★★★★★★★★★★**段 E3a-3 —— 大域化**（[Stacks] Lemma 01PW の大域版）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X` が有限個のアフィン開 `U i`（重なりもアフィン）で覆われ、そこで `M` が自明化し、
`X_s ⊓ U i` の上の関数 `g i` が重なりで一致するとき、
★**ある `n` と `a i ∈ Γ(X, U i)` と `sheafify (M^{⊗n})` の大域切断 `t` があって**

* `g i · f_i^n = a_i` —— 分母が払えている
* `t|_{U i}` は `a_i` が定める切断（の層化）である

★★★これが原典の「**[some positive tensor power of]**」の中身である。

## ★逸脱の記録

`U i ⊓ U j` がアフィンであること（`hUij`）を仮定に置いた。
★`X` が分離的なら自動だが、本定理はそれを導かない——**消費側が供給する**。 -/
theorem exists_global_section_of_localData {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type))
    (g : ∀ i, (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type))
    (hg : ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left :
        X.basicOpen (trivValue M (U i) (e i) s) ⊓ X.basicOpen (trivValue M (U j) (e j) s)
          ≤ X.basicOpen (trivValue M (U i) (e i) s))).op (g i)
        = X.presheaf.map (homOfLE (inf_le_right :
            X.basicOpen (trivValue M (U i) (e i) s) ⊓ X.basicOpen (trivValue M (U j) (e j) s)
              ≤ X.basicOpen (trivValue M (U j) (e j) s))).op (g j)) :
    ∃ (n : ℕ) (a : ∀ i, (Γ(X, U i) : Type))
      (t : (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.obj
        (op (⊤ : X.Opens)) : Type)),
      (∀ i, g i * (algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
            (trivValue M (U i) (e i) s)) ^ n
          = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a i)) ∧
      (∀ i, ((sheafifyFunctor X).obj (presheafTensorPow M n)).val.map
          (homOfLE (le_top : U i ≤ ⊤)).op t
        = ((sheafifyUnit X (presheafTensorPow M n)).app (op (U i))).hom
            (secOfFun M (U i) (e i) n (a i))) := by
  obtain ⟨n, a, hden, hagree⟩ := exists_common_exponent Finset.univ M U
    (fun i _ => hU i) (fun i _ j _ => hUij i j) e s g
    (fun i _ j _ => hg i j)
  obtain ⟨t, ht, -⟩ := exists_glue_of_agree M U hcov n e a
    (fun i j => hagree i (Finset.mem_univ i) j (Finset.mem_univ j))
  exact ⟨n, a, t, fun i => hden i (Finset.mem_univ i), ht⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_common_exponent.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母と重なりの両方に効く単一の指数)",
    sectionId := "genell-prop-1-4" }

def exists_glue_of_agree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(一致する族は層化の側で貼り合う)",
    sectionId := "genell-prop-1-4" }

def exists_global_section_of_localData.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3a-3——大域化、[Stacks] 01PW の大域版)",
    sectionId := "genell-prop-1-4" }

def exists_global_section_of_localData.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[Stacks]" "Lemma 01PW(ample な可逆層の切断の延長——本定理がその形式化)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)。本定理は前層加群の言葉で独自に組んだ") 7,
    .citation "[ABC3]" "exists_pow_mul_eq_res(分母を払う、§9-822)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_mul_eq_res") 3,
    .citation "[ABC3]" "res_secOfFun_eq_iff(重なりの一致判定、§9-828)"
      (.inProject "ABC3" "ABC3.Found.GenEll.res_secOfFun_eq_iff") 3,
    .citation "[ABC3]" "exists_unique_glue_top(層化の側で貼る、§9-827)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_unique_glue_top") 3,
    .implicitStep
      ("★★★原典の「[some positive tensor power of]」には**冪を取る理由が 3 つ**畳まれていた: " ++
       "(a) 分母を払う冪(§9-822)、(b) 重なりを合わせる冪(§9-826・§9-830)、" ++
       "(c) 遷移函数の N 乗(§9-828)。★合図 1 つを節点 1 つに対応させる規律で数えると節点 3 個である") 8,
    .implicitStep
      ("★U i ⊓ U j が**アフィンであること**を仮定に置いた(hUij)。" ++
       "★★X が分離的なら自動だが、本定理はそれを導かない——消費側が供給する。" ++
       "★★★また被覆は Fintype ι で表した(ample から来る有限被覆はこの形である、§9-817)") 7,
    .implicitStep
      ("★★大域切断は **sheafify (M^{⊗n}) の側**にある(前層 M^{⊗n} の側ではない)。" ++
       "§9-824 のとおり「局所自明な前層加群は層である」は偽だからである") 8 ]

end ABC3.Found.GenEll
