/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim

/-!
# 第 1206 ブロック —— **座標が部分体に入る点は、その部分体の点から来る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`L(H)` の道の最後の原始関数

第 1205 で `HasLCyclicJ` から `L̄` 上の位数 `l` の点 `Q` を作った。
第 1199 が要るのは**有限拡大 `L''` の上で有理な** `Q` である。

☆第 1195（`exists_finite_subextension`）が「有限個の座標は有限次拡大に入る」
を与えるので、あとは**座標が部分体に入る点はその部分体の点から来る**が要る。
★それが本ブロックの `exists_rhPoint_eq` である。

| 定理 | 内容 |
|---|---|
| `exists_rhPoint_eq` | ★★★座標が `f` の像に入る点は `rhPoint f` の像 |
| `baseChange_map_intermediate` | ★塔 `L → M → L̄` で `(W ⊗ M) ⊗ L̄ = W ⊗ L̄` |

★★`rhPoint` は加法を保ち単射なので、位数は `addOrderOf_rhPoint`（在庫）で移る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

/-! ## ★★★★★★★★★★★★座標から点を戻す -/

/-- ★★★★★★★★★★★★
**座標が `f` の像に入る点は `rhPoint f` の像である**——★**無条件**（第 1206）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`O` は `0` から来る。`(x, y)` の場合は `x = f a`・`y = f b` を取り、
非特異性は `map_nonsingular` の逆向きで降ろす。

★★★これが `L̄` の点を有限拡大へ降ろす原始関数である
——第 1195（座標は有限次拡大に入る）と対で使う。 -/
theorem exists_rhPoint_eq {F K : Type*} [Field F] [Field K] (f : F →+* K)
    (W : WeierstrassCurve F) (P : (W.map f).toAffine.Point)
    (hx : (pointCoords P).1 ∈ Set.range f) (hy : (pointCoords P).2 ∈ Set.range f) :
    ∃ P' : W.toAffine.Point, rhPoint f W P' = P := by
  cases P with
  | zero => exact ⟨0, rfl⟩
  | some x y h =>
    obtain ⟨a, ha⟩ := hx
    obtain ⟨b, hb⟩ := hy
    simp only [pointCoords_some] at ha hb
    subst ha
    subst hb
    exact ⟨.some a b ((W.toAffine.map_nonsingular (f := f) f.injective a b).1 h), rfl⟩

/-! ## ★塔での底変換 -/

/-- ★**塔 `L → M → L̄` で `(W ⊗ M) ⊗ L̄ = W ⊗ L̄`**——★**無条件**（第 1206）。 -/
theorem baseChange_map_intermediate {L Lbar : Type*} [Field L] [Field Lbar] [Algebra L Lbar]
    (W : WeierstrassCurve L) (M : Type*) [Field M] [Algebra L M] [Algebra M Lbar]
    [IsScalarTower L M Lbar] :
    (W.baseChange M).map (algebraMap M Lbar) = W.baseChange Lbar := by
  show (W.map (algebraMap L M)).map (algebraMap M Lbar) = W.map (algebraMap L Lbar)
  rw [WeierstrassCurve.map_map]
  congr 1
  exact (IsScalarTower.algebraMap_eq L M Lbar).symm

/-! ## ★出典の紐付け(`.src`) -/

def exists_rhPoint_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標が f の像に入る点は rhPoint f の像。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_rhPoint_eq.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1206）**——`L̄` の点を有限拡大へ降ろす原始関数である。" ++
       "☆第 1195（`exists_finite_subextension`——有限個の座標は有限次拡大に入る）と" ++
       "対で使うと、第 1205 が作った `L̄` 上の位数 `l` の点が" ++
       "**有限拡大 `L''` の上で有理**になり、第 1199 に入る。" ++
       "★位数は `addOrderOf_rhPoint`（在庫）で移る。") 2 ]

def baseChange_map_intermediate.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(塔 L → M → L̄ で (W ⊗ M) ⊗ L̄ = W ⊗ L̄。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
