/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.Artinian.Module
import Mathlib.RingTheory.DiscreteValuationRing.TFAE
import Mathlib.RingTheory.DedekindDomain.Basic
import ABC3.Meta.Claim

/-!
# 第 1363 ブロック —— **道の残り 2 段（Artin 環と DVR）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1357 の道の第 5 段と第 7 段

「完備 DVR `R` の有限整閉包 `R′` は完備 DVR」の 7 段のうち、

| 段 | 内容 | 第 |
|---|---|---|
| 1 | 完備なイデアルに沿う冪等元の持ち上げ | 1357 |
| 2 | 完備な整域の商に非自明な冪等元なし | 1358 |
| 3 | `I`-進完備性は線型同型で移る | 1359 |
| 4 | 有限直積・係数環の移送 | 1360-1361 |
| 6 | Artin ＋ 冪等元なし ⟹ 局所 | 1362 |

はすでに済んでいる。★本ブロックは残る

* **段 5**——`R′/m_R R′` は Artin 環（剰余体上有限次元）
* **段 7**——Dedekind ＋ 局所 ⟹ DVR

を作る。☆どちらも mathlib の在庫の組み合わせ 1 行である。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★
**Artin 環の上の加群有限代数は Artin 環**——★**無条件**（第 1363、段 5）。

☆剰余体 `k` は Artin 環なので、`k` 上有限次元の環は Artin 環である。 -/
theorem isArtinianRing_of_module_finite (k : Type*) {B : Type*} [CommRing k] [IsArtinianRing k]
    [Ring B] [Algebra k B] [Module.Finite k B] : IsArtinianRing B :=
  isArtinian_of_tower k inferInstance

set_option linter.overlappingInstances false in
/-- ★★★★★★★★★★★★
**Dedekind 整域が局所で体でなければ DVR**——★**無条件**（第 1363、段 7）。

☆`IsDiscreteValuationRing.TFAE`（在庫）の 3 番目と 1 番目を繋ぐだけである。 -/
theorem isDiscreteValuationRing_of_isDedekindDomain (R : Type*) [CommRing R] [IsDomain R]
    [IsLocalRing R] [IsDedekindDomain R] (h : ¬ IsField R) : IsDiscreteValuationRing R :=
  ((IsDiscreteValuationRing.TFAE R h).out 2 0).mp ‹IsDedekindDomain R›

/-! ## ★出典の紐付け(`.src`) -/

def isArtinianRing_of_module_finite.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Artin 環の上の加群有限代数は Artin 環。★無条件)",
    sectionId := "genell-thm-3-8" }

def isDiscreteValuationRing_of_isDedekindDomain.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Dedekind 整域が局所で体でなければ DVR。★無条件)",
    sectionId := "genell-thm-3-8" }

def isDiscreteValuationRing_of_isDedekindDomain.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1363）**——第 1357 の道の**第 5 段と第 7 段**である。" ++
       "☆これで 7 段すべてが揃った——" ++
       "残るのは `IsIntegralClosure` の言葉で組み立てる配管だけである。") 19 ]

end ABC3.Found.GenEll
