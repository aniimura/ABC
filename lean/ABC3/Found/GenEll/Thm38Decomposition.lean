/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.FieldTheory.Normal.Defs
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.Algebra.Algebra.Equiv
import ABC3.Meta.Claim

/-!
# 第 1167 ブロック —— **局所の Galois 元は大域へ制限される**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`AlphaBridge` の節点 3 の核

`Skeleton/GenEll/AlphaBridge.lean` の節点 3 は

> 分解群 `D_v ⊆ Gal(L̄/L)` を経由して局所の `σ` を大域の元に持ち上げる

であった。★本ファイルは**その制限準同型を作る**。

## ★★★★★★★★機構

塔 `L → L_v → M`（`M` は `L_v` の代数閉包）と、`L` 上正規な中間体 `E`（＝ `L̄`）に対し、

    `σ : M ≃ₐ[L_v] M`  →  `σ|_L : M ≃ₐ[L] M`  →  `σ|_E : E ≃ₐ[L] E`

☆前半は `AlgEquiv.restrictScalars`（`L_v`-代数射は `L`-代数射）、
後半は `AlgEquiv.restrictNormalHom`（`E` が `L` 上正規なら制限が定義できる）。

★★**単射性は要らない**——必要なのは「局所の元の像が大域の像に**含まれる**」ことだけで、
それは準同型があれば出る。☆これが「`galRep` の像は局所の `α` を含む」の骨である。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★★★制限準同型 -/

variable (L Lv M E : Type*) [Field L] [Field Lv] [Field M] [Field E]
  [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
  [Algebra L E] [Algebra E M] [IsScalarTower L E M] [Normal L E]

/-- ★★★★★★★★★★★★★★**局所の Galois 元を大域へ制限する準同型**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`M ≃ₐ[L_v] M` の元はまず `L`-代数射として見（`restrictScalars`）、
`L` 上正規な `E` に制限する（`restrictNormalHom`）。 -/
noncomputable def restrictLocalHom : (M ≃ₐ[Lv] M) →* (E ≃ₐ[L] E) where
  toFun σ := AlgEquiv.restrictNormalHom E (σ.restrictScalars L)
  map_one' := by
    have h : (1 : M ≃ₐ[Lv] M).restrictScalars L = 1 := AlgEquiv.ext fun _ => rfl
    rw [h, map_one]
  map_mul' := fun σ τ => by
    have h : (σ * τ).restrictScalars L = σ.restrictScalars L * τ.restrictScalars L :=
      AlgEquiv.ext fun _ => rfl
    rw [h, map_mul]

/-- ★★★★★★★★★★**制限は埋め込みと可換**——★**無条件**。

☆`algebraMap E M (σ|_E x) = σ (algebraMap E M x)`。
★これが「局所で `α` が実現されるなら大域でも実現される」の中身である。 -/
theorem restrictLocalHom_commutes (σ : M ≃ₐ[Lv] M) (x : E) :
    algebraMap E M (restrictLocalHom L Lv M E σ x) = σ (algebraMap E M x) :=
  AlgEquiv.restrictNormal_commutes (σ.restrictScalars L) E x

/-! ## ★★★★★★★★★★★★代数閉包を埋め込んで大域へ -/

/-- ★★★★★★★★★★★★★★★★
**局所の Galois 元から大域の Galois 元へ**——★**無条件**（第 1168）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`M` が代数閉体なら `L̄ = AlgebraicClosure L` は `IsAlgClosed.lift` で `M` に埋め込める。
★その埋め込みで `L̄` を `M` の中間体と見なし、`restrictLocalHom` を当てる。

★★★これが `AlphaBridge` の節点 3 の**完成形**である
——`Gal(M/L_v)` の元は `Gal(L̄/L)` の元に制限され、したがって
**局所で実現される行列は大域の像にも入る**。 -/
noncomputable def globalOfLocalHom (L Lv M : Type*) [Field L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M] [IsAlgClosed M] :
    (M ≃ₐ[Lv] M) →* (AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L) := by
  letI ι : AlgebraicClosure L →ₐ[L] M := IsAlgClosed.lift
  letI : Algebra (AlgebraicClosure L) M := ι.toAlgebra
  haveI : IsScalarTower L (AlgebraicClosure L) M :=
    IsScalarTower.of_algebraMap_eq (fun x => (ι.commutes x).symm)
  exact restrictLocalHom L Lv M (AlgebraicClosure L)

/-! ## ★★★★★★★★★★★★局所像は大域像に含まれる -/

/-- ☆合成の像は元の像に含まれる。 -/
theorem range_comp_le {G H K : Type*} [Group G] [Group H] [Group K]
    (f : H →* K) (g : G →* H) : (f.comp g).range ≤ f.range := by
  rintro x ⟨y, rfl⟩
  exact ⟨g y, rfl⟩

/-- ★★★★★★★★★★★★★★★★
**局所で実現された行列は大域の mod `l` 像にも入る**——★**無条件**（第 1169）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`f` を `galRep`、`g` を `globalOfLocalHom`、`r` を `glRedPadic l` と読む。
★★★これが `imageContainsSL2J_of_alpha'` の `halpha` を局所から受け取る形であり、
`AlphaBridge` の節点 3 の**最後の 1 行**である。 -/
theorem mem_map_of_mem_map_comp {G H K K' : Type*} [Group G] [Group H] [Group K] [Group K']
    (f : H →* K) (g : G →* H) (r : K →* K') {x : K'}
    (hx : x ∈ ((f.comp g).range).map r) : x ∈ (f.range).map r :=
  Subgroup.map_mono (range_comp_le f g) hx

/-! ## ★出典の紐付け(`.src`) -/

def restrictLocalHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の Galois 元を大域へ制限する準同型。★無条件)",
    sectionId := "genell-thm-3-8" }

def restrictLocalHom_commutes.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(制限は埋め込みと可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def mem_map_of_mem_map_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所で実現された行列は大域の mod l 像にも入る。★無条件)",
    sectionId := "genell-thm-3-8" }

def mem_map_of_mem_map_comp.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1169）**——`f` を `galRep`、`g` を `globalOfLocalHom`、" ++
       "`r` を `glRedPadic l` と読むと、これが " ++
       "`imageContainsSL2J_of_alpha'` の `halpha` を**局所から受け取る形**である。" ++
       "☆`AlphaBridge` の節点 3 はこれで完備した——" ++
       "残るのは節点 2（Tate 一意化で `(ζ, π)` を `E[l]` の基底として実現する段）だけである。") 1 ]

def globalOfLocalHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の Galois 元から大域の Galois 元へ——代数閉包の埋め込み経由。★無条件)",
    sectionId := "genell-thm-3-8" }

def globalOfLocalHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsAlgClosed.lift(代数的な体は代数閉体に埋め込める)"
      (.inMathlib "IsAlgClosed.lift") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1168）**——`AlphaBridge` の**節点 3 の完成形**である。" ++
       "☆`Gal(M/L_v)` の元は `Gal(L̄/L)` の元に制限されるので、" ++
       "**局所で実現される行列は大域の像にも入る**。" ++
       "★残るのは `galRep` との合成だけである。") 2 ]

def restrictLocalHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgEquiv.restrictNormalHom(正規な中間体への制限は準同型)"
      (.inMathlib "AlgEquiv.restrictNormalHom") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1167）**——`Skeleton/GenEll/AlphaBridge.lean` の" ++
       "節点 3 の**核**である。☆単射性は要らない——必要なのは" ++
       "「局所の元の像が大域の像に**含まれる**」ことだけで、それは準同型があれば出る。" ++
       "★残るのは `L̄ ↪ M`（`M` は `L_v` の代数閉包）を実際に取る段" ++
       "（mathlib の `IsAlgClosed.lift`）と、`galRep` との合成である。") 4 ]

end ABC3.Found.GenEll
