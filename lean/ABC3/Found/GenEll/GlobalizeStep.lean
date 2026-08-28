/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.OverlapCriterion
import ABC3.Found.GenEll.ClearDenomGlue
import ABC3.Meta.Claim

/-!
# ★★★★★★★大域化の 3 段 —— 指数を揃え、単元のずれを吸収する（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 段 E3a-3 の残り 3 段

`§9-828`（`OverlapCriterion.lean`）は「重なりで一致するか」を

    `a|_{V⊓V'} · u^N = a'|_{V⊓V'}`

という**関数の等式**に落とした。★本ファイルはその等式を**実際に作る**ための 3 段である:

| 段 | 補題 | 内容 |
|---|---|---|
| 1 | `trivValue_res_transUnit` | ★隣り合うチャートの座標は遷移単元で移る: `f_j\| = f_i\| · u` |
| 2 | `exists_common_denominator` | ★分母を払う指数は**有限個まとめて揃えられる** |
| 3 | `exists_pow_mul_eq_of_unit_scaling` | ★★**単元倍のずれは冪で吸収できる** |

★★★段 3 が心臓である。段 2 で揃えた `a_i` は、重なりでは
「同じ `g` に `f_i^N` を掛けたもの」対「同じ `g` に `f_j^N = (f_i·u)^N` を掛けたもの」なので、
**`u^N` だけずれている**。★そのずれは `X_s` の上では等式だが `V ⊓ V'` 全体では等式でない
——そこで `§9-826`（`exists_pow_mul_eq_of_res_eq`）をもう一度使って冪で潰す。

## ★★測定の記録

★`basicOpen_res_trivValue_eq`——**2 つのチャートは重なりで同じ非消失軌跡を見る**。
`f_j| = f_i|·u` で `u` が単元だからである。★これがないと「どちらの `D(f)` で
一致を見るのか」が定まらない。

## ★残っている段（明示）

★★**組み上げは本ファイルに無い**。要るのは、具体の状況（有限被覆・`trivValue`）を
段 3 の**抽象形**へ落とす配管である:

1. `a_i ∈ Γ(X, U_i)` を `W = U_i ⊓ U_j` へ制限して `b_i` にする
2. `D = X.basicOpen (f_i|_W)` の上で `b_i| = g·f_i|^N`、`b_j| = g·(f_i|·u)^N` を出す
3. 段 3 を当て、`§9-826` の `exists_common_pow` で `m` を有限個まとめて揃える
4. `§9-828` の判定基準に通し、`§9-827` の `exists_unique_glue_top` で貼る
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★制限の合成 -/

/-- ★**制限を 2 回するのは 1 回するのと同じ**（構造層、元に当てた形）。

★`Opens` の射は `Subsingleton` なので、合成の証明項を気にしなくてよい。 -/
theorem res_trans {V W Z : X.Opens} (h1 : Z ≤ W) (h2 : W ≤ V) (x : (Γ(X, V) : Type)) :
    X.presheaf.map (homOfLE h1).op (X.presheaf.map (homOfLE h2).op x)
      = X.presheaf.map (homOfLE (h1.trans h2)).op x := by
  rw [show (homOfLE (h1.trans h2) : Z ⟶ V) = homOfLE h1 ≫ homOfLE h2 from
    Subsingleton.elim _ _, op_comp, X.presheaf.map_comp]
  rfl

/-! ## ★★★★★段 1 —— 隣り合うチャートの座標は遷移単元で移る -/

/-- ★★★★★**隣り合うチャートの座標は遷移単元で移る**。

    `f_j|_{V⊓V'} = f_i|_{V⊓V'} · u`   （`f_i = trivValue M V e s` 等）

★★これが「重なりで `u^N` だけずれる」の出どころである。
★機構は `trivValue_restrict`（座標は制限と可換）＋ `trivValue_transUnit`（在庫）だけである。 -/
theorem trivValue_res_transUnit (M : X.PresheafOfModules) (V V' : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (e' : (restrictPresheafFunctor X V').obj M ≅ 𝟙_ (PresheafModulesOn X V'))
    (s : (M.obj (op ⊤) : Type)) :
    X.presheaf.map (homOfLE (inf_le_right : V ⊓ V' ≤ V')).op (trivValue M V' e' s)
      = X.presheaf.map (homOfLE (inf_le_left : V ⊓ V' ≤ V)).op (trivValue M V e s)
        * transUnit M (V ⊓ V')
            (trivialOfLe M (inf_le_left : V ⊓ V' ≤ V) e)
            (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e') := by
  rw [← trivValue_restrict M (inf_le_right : V ⊓ V' ≤ V') e' s,
    ← trivValue_restrict M (inf_le_left : V ⊓ V' ≤ V) e s,
    trivValue_transUnit]

/-- ★★★★★★**2 つのチャートは重なりで同じ非消失軌跡を見る**。

★これがないと「どちらの `D(f)` の上で一致を見るのか」が定まらない。
★★機構は「`f_j| = f_i|·u` で `u` は単元」だけである。 -/
theorem basicOpen_res_trivValue_eq (M : X.PresheafOfModules) (Vi Vj : X.Opens)
    (ei : (restrictPresheafFunctor X Vi).obj M ≅ 𝟙_ (PresheafModulesOn X Vi))
    (ej : (restrictPresheafFunctor X Vj).obj M ≅ 𝟙_ (PresheafModulesOn X Vj))
    (s : (M.obj (op ⊤) : Type)) :
    X.basicOpen (X.presheaf.map (homOfLE (inf_le_left : Vi ⊓ Vj ≤ Vi)).op
        (trivValue M Vi ei s))
      = X.basicOpen (X.presheaf.map (homOfLE (inf_le_right : Vi ⊓ Vj ≤ Vj)).op
        (trivValue M Vj ej s)) := by
  rw [trivValue_res_transUnit M Vi Vj ei ej s, X.basicOpen_mul,
    X.basicOpen_of_isUnit (isUnit_transUnit M (Vi ⊓ Vj) _ _)]
  exact (inf_eq_left.2 (X.basicOpen_le _)).symm

/-! ## ★★★★★段 2 —— 分母を払う指数は揃えられる -/

/-- ★★**分母を払う指数は上げてよい**——余分の `f^k` は分子に押し込めばよい。 -/
theorem exists_res_pow_mono (U : X.Opens) (f : (Γ(X, U) : Type))
    (g : (Γ(X, X.basicOpen f) : Type)) {n : ℕ}
    (h : ∃ a : (Γ(X, U) : Type),
      g * (algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ n
        = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a)
    {m : ℕ} (hnm : n ≤ m) :
    ∃ a : (Γ(X, U) : Type),
      g * (algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ m
        = algebraMap (Γ(X, U) : Type) (Γ(X, X.basicOpen f) : Type) a := by
  obtain ⟨a, ha⟩ := h
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hnm
  refine ⟨a * f ^ k, ?_⟩
  rw [pow_add, ← mul_assoc, ha, map_mul, map_pow]

/-- ★★★★★**分母を払う指数は有限個まとめて揃えられる** —— 段 E3a-3 の第 2 段。

★機構は `exists_pow_mul_eq_res`（§9-822）＋ `exists_common_pow`（§9-826、`Finset.sup`）
＋ 単調性（`exists_res_pow_mono`）だけである。 -/
theorem exists_common_denominator {ι : Type} (I : Finset ι) (U : ι → X.Opens)
    (hU : ∀ i ∈ I, IsAffineOpen (U i)) (f : ∀ i, (Γ(X, U i) : Type))
    (g : ∀ i, (Γ(X, X.basicOpen (f i)) : Type)) :
    ∃ N : ℕ, ∀ i ∈ I, ∃ a : (Γ(X, U i) : Type),
      g i * (algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) (f i)) ^ N
        = algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) a := by
  refine exists_common_pow I
    (fun i N => ∃ a : (Γ(X, U i) : Type),
      g i * (algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) (f i)) ^ N
        = algebraMap (Γ(X, U i) : Type) (Γ(X, X.basicOpen (f i)) : Type) a)
    (fun i n m hnm h => exists_res_pow_mono (U i) (f i) (g i) h hnm) ?_
  intro i hi
  obtain ⟨n, a, ha⟩ := exists_pow_mul_eq_res (U i) (hU i hi) (f i) (g i)
  exact ⟨n, a, ha⟩

/-! ## ★★★★★★★★段 3 —— 単元倍のずれは冪で吸収できる -/

/-- ★★★★★★★★**単元倍のずれは冪で吸収できる** —— 段 E3a-3 の心臓（抽象形）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`D(f)` の上で `b_i` が `g·f^N`、`b_j` が `g·(f·u)^N` に見えるなら、
★**ある `m` で `f^m·(b_i·u^N) = f^m·b_j`**。

★★これが `§9-828` の判定基準（`a|·u^N = a'|`）を**実際に作る**段である。
★★★機構は「`D(f)` の上では等式である」（単元 `u` の `N` 乗が両辺に現れて相殺する）
＋ `§9-826`（`exists_pow_mul_eq_of_res_eq`）だけである。

★抽象形にしてあるのは、具体の `trivValue` の言葉を持ち込むと式が読めなくなるからである
——`f = f_i|_W`、`u` = 遷移単元、`b_i = a_i|_W` と読む。 -/
theorem exists_pow_mul_eq_of_unit_scaling (W : X.Opens) (hW : IsAffineOpen W)
    (f u : (Γ(X, W) : Type)) (N : ℕ) (bi bj : (Γ(X, W) : Type))
    (g : (Γ(X, X.basicOpen f) : Type))
    (hbi : algebraMap (Γ(X, W) : Type) (Γ(X, X.basicOpen f) : Type) bi
        = g * (algebraMap (Γ(X, W) : Type) (Γ(X, X.basicOpen f) : Type) f) ^ N)
    (hbj : algebraMap (Γ(X, W) : Type) (Γ(X, X.basicOpen f) : Type) bj
        = g * (algebraMap (Γ(X, W) : Type) (Γ(X, X.basicOpen f) : Type) (f * u)) ^ N) :
    ∃ m : ℕ, f ^ m * (bi * u ^ N) = f ^ m * bj := by
  refine exists_pow_mul_eq_of_res_eq W hW f (bi * u ^ N) bj ?_
  rw [map_mul, map_pow, hbi, hbj, map_mul, mul_pow]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def trivValue_res_transUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(隣り合うチャートの座標は遷移単元で移る)",
    sectionId := "genell-prop-1-4" }

def basicOpen_res_trivValue_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(2 つのチャートは重なりで同じ非消失軌跡を見る)",
    sectionId := "genell-prop-1-4" }

def exists_common_denominator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母を払う指数は有限個まとめて揃えられる)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_of_unit_scaling.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(単元倍のずれは冪で吸収できる)",
    sectionId := "genell-prop-1-4" }

def res_trans.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(制限の合成——構造層、元に当てた形)",
    sectionId := "genell-prop-1-4" }

def exists_pow_mul_eq_of_unit_scaling.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pow_mul_eq_of_res_eq(重なりの一致は冪で強制できる、§9-826)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_mul_eq_of_res_eq") 3,
    .citation "[ABC3]" "trivValue_transUnit / transUnit_restrict(在庫、LocalMetric.lean)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_transUnit") 3,
    .citation "[Stacks]" "Lemma 01PW(の大域化の段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★重なりでは「同じ g に f_i^N を掛けたもの」対「同じ g に f_j^N = (f_i·u)^N を掛けたもの」" ++
       "なので **u^N だけずれている**。★★そのずれは D(f) の上では等式だが " ++
       "V ⊓ V' 全体では等式でない——そこで §9-826 をもう一度使って冪で潰す") 8,
    .implicitStep
      ("★★★**組み上げは本ファイルに無い**。具体の状況(有限被覆・trivValue)を段 3 の" ++
       "抽象形へ落とす配管が要る: (1) a_i を W = U_i ⊓ U_j へ制限、" ++
       "(2) D = basicOpen (f_i|_W) の上で b_i| = g·f_i|^N と b_j| = g·(f_i|·u)^N を出す、" ++
       "(3) 段 3 と exists_common_pow で m を揃える、" ++
       "(4) §9-828 の判定基準に通し §9-827 で貼る") 8 ]

end ABC3.Found.GenEll
