/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRep
import ABC3.Meta.Claim

/-!
# 第 1201 ブロック —— **`T_l E` から `E[l]` への射影**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`HasLCyclicJ` を点の言葉へ運ぶ第 1 段

`HasLCyclicJ`（`Found/GenEll/EllModuliGalois.lean`）は
**`T_l E` の `mod l` 還元の中の `Gal`-安定な直線**として述べられている。
一方 `Lemma 3.5`・`Lemma 3.7` が要るのは**位数 `l` の点 `Q`** である。

☆`tateModule` は「両立する系 `f : ∀ n, E[l^n]`、`l • f(n+1) = f(n)`」として
**具体的に定義されている**（`Interface/GaloisRep/Torsion.lean`）ので、
第 1 層を取る射影 `f ↦ f 1` がそのまま書ける。

| 定理 | 内容 |
|---|---|
| `tateProj` | ★`T_l E →+ E(L)`、`f ↦ f 1` |
| `nsmul_tateProj` | ★`l • (tateProj t) = 0`——像は `E[l]` に入る |
| `tateProj_galTate` | ★★**`Gal` と可換**（定義から `rfl`） |
| `addOrderOf_tateProj` | ★★★`≠ 0` なら**位数はちょうど `l`** |

★★★これで「`T_l E` の側の言明」から「位数 `l` の点」への橋の
**`Gal` 同変性の段**が済む。
☆残るのは `tateProj` の核が `l • T_l E` であること（`T_l E / l ≅ E[l]`）で、
これは `l •` が `E(L̄)` の上で全射であることを要る——次のブロックで取る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {K : Type} [Field K] [DecidableEq K]

/-! ## ★★★★★★★★第 1 層への射影 -/

/-- ★★★★★★★★**`T_l E` から `E(L)` への第 1 層の射影**——★**無条件**（第 1201）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`tateModule` の元は両立する系 `f : ∀ n, E[l^n]` なので、`f 1` を取るだけである。 -/
def tateProj (L : Type) [Field L] [DecidableEq L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) :
    tateModule (W.baseChange L) l →+ (W.baseChange L).toAffine.Point where
  toFun := fun f =>
    ((f : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) 1 :
      (W.baseChange L).toAffine.Point)
  map_zero' := rfl
  map_add' := fun _ _ => rfl

variable {L : Type} [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★**像は `l`-捩れに入る**——★**無条件**（第 1201）。 -/
theorem nsmul_tateProj (W : WeierstrassCurve K) (l : ℕ)
    (t : tateModule (W.baseChange L) l) :
    l • tateProj L W l t = 0 := by
  have h : (l ^ 1) • (tateProj L W l t) = 0 :=
    ((t : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) 1).2
  simpa using h

/-- ★★★★★★★★**射影は `Gal` と可換**——★**無条件**（第 1201）。

☆`galTate` は各層で `galPoint` を当てるだけなので、定義から `rfl` である。
★★★これが `HasLCyclicJ`（`T_l E` の側の安定直線）を
`E[l]` の中の安定直線へ運ぶ**要**である。 -/
theorem tateProj_galTate (W : WeierstrassCurve K) (l : ℕ) (σ : L ≃ₐ[K] L)
    (t : tateModule (W.baseChange L) l) :
    tateProj L W l (galTate W l σ t) = galPoint W σ (tateProj L W l t) := rfl

/-- ★★★★★★★★★★**`≠ 0` なら位数はちょうど `l`**——★**無条件**（第 1201）。

☆`l` が素数で `l • P = 0` なら位数は `1` か `l`、`P ≠ 0` なら `l` である。
★★★これが `Lemma 3.5`・`Lemma 3.7` が要る `addOrderOf Q = l` の形である。 -/
theorem addOrderOf_tateProj (W : WeierstrassCurve K) {l : ℕ} (hl : l.Prime)
    (t : tateModule (W.baseChange L) l) (h : tateProj L W l t ≠ 0) :
    addOrderOf (tateProj L W l t) = l := by
  have hdvd : addOrderOf (tateProj L W l t) ∣ l :=
    addOrderOf_dvd_of_nsmul_eq_zero (nsmul_tateProj W l t)
  rcases (Nat.Prime.eq_one_or_self_of_dvd hl _ hdvd) with h1 | h2
  · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h1) h
  · exact h2

/-! ## ★出典の紐付け(`.src`) -/

def tateProj.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(T_l E から E(L) への第 1 層の射影。★無条件)",
    sectionId := "genell-thm-3-8" }

def nsmul_tateProj.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(第 1 層の射影の像は l-捩れに入る。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateProj_galTate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(第 1 層の射影は Gal と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def addOrderOf_tateProj.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(第 1 層の射影が 0 でなければ位数はちょうど l。★無条件)",
    sectionId := "genell-thm-3-8" }

def addOrderOf_tateProj.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1201）**——`HasLCyclicJ`（`T_l E` の `mod l` 還元の" ++
       "中の `Gal`-安定な直線）を `Lemma 3.5`・`Lemma 3.7` が要る" ++
       "**位数 `l` の点**へ運ぶ第 1 段である。" ++
       "☆残るのは `tateProj` の核が `l • T_l E` であること" ++
       "（`T_l E / l ≅ E[l]`）で、これは `l •` が `E(L̄)` の上で全射であることを要る。") 3 ]

end ABC3.Found.GaloisRep
