/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Meta.Claim

/-!
# 第 1204 ブロック —— **ベクトルの `mod l` 還元**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★これは何か——`HasLCyclicJ` の直線を `ℤ_l²` へ持ち上げる配管

`HasLCyclicJ`（`Found/GenEll/EllModuliGalois.lean`）は
`v : Fin 2 → ZMod l`（`≠ 0`）と、`glRedPadic` で落とした像の行列 `M` について
`M.mulVec v = c • v` を要求する。

一方 `T_l E` の側の言明は `ℤ_l²` の上にある。両者を繋ぐのが本ブロックの `redVec`:

| 定理 | 内容 |
|---|---|
| `redVec` | ★`(Fin 2 → ℤ_[l]) → (Fin 2 → ZMod l)`、成分ごとの `PadicInt.toZMod` |
| `redVec_mulVec` | ★★行列の作用と可換 |
| `coe_glRedPadic` | ★`glRedPadic` の行列は成分ごとの還元（定義から `rfl`） |
| `redVec_mulVec_glRed` | ★★★**`HasLCyclicJ` が使う形** |
| `redVec_eq_zero_iff` | ★★★**核はちょうど `l · ℤ_l²`** |
| `exists_redVec_eq` | ★★還元は**全射**——`v` は持ち上がる |

★★★`redVec_eq_zero_iff` が第 1203（`tateProj` の核）と対になって、
「`mod l` で `≠ 0`」と「第 1 層の点が `≠ 0`」を結ぶ。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable (l : ℕ) [Fact (Nat.Prime l)]

/-! ## ★★★★★★成分ごとの還元 -/

/-- ★★★★★★**ベクトルの `mod l` 還元**（第 1204）。 -/
noncomputable def redVec (w : Fin 2 → ℤ_[l]) : Fin 2 → ZMod l := fun i => PadicInt.toZMod (w i)

@[simp] theorem redVec_apply (w : Fin 2 → ℤ_[l]) (i : Fin 2) :
    redVec l w i = PadicInt.toZMod (w i) := rfl

/-- ★★★★**還元は行列の作用と可換**——★**無条件**（第 1204）。 -/
theorem redVec_mulVec (M : Matrix (Fin 2) (Fin 2) ℤ_[l]) (w : Fin 2 → ℤ_[l]) :
    redVec l (M.mulVec w) = (M.map PadicInt.toZMod).mulVec (redVec l w) := by
  funext i
  show PadicInt.toZMod (∑ j, M i j * w j)
    = ∑ j, PadicInt.toZMod (M i j) * PadicInt.toZMod (w j)
  rw [map_sum]
  exact Finset.sum_congr rfl (fun j _ => map_mul _ _ _)

/-- ★★**`glRedPadic` の行列は成分ごとの還元**（定義から `rfl`）。 -/
theorem coe_glRedPadic (g : GL (Fin 2) ℤ_[l]) :
    ((glRedPadic l g : GL (Fin 2) (ZMod l)) : Matrix (Fin 2) (Fin 2) (ZMod l))
      = (g : Matrix (Fin 2) (Fin 2) ℤ_[l]).map PadicInt.toZMod := rfl

/-- ★★★★★★**`HasLCyclicJ` が使う形**——★**無条件**（第 1204）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`ℤ_l` の上で行列を当ててから還元しても、還元してから当てても同じ。 -/
theorem redVec_mulVec_glRed (g : GL (Fin 2) ℤ_[l]) (w : Fin 2 → ℤ_[l]) :
    redVec l ((g : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec w)
      = ((glRedPadic l g : GL (Fin 2) (ZMod l)) :
          Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec (redVec l w) := by
  rw [redVec_mulVec, coe_glRedPadic]

/-- ★★**還元は和と可換**（第 1233）。 -/
theorem redVec_add (a b : Fin 2 → ℤ_[l]) :
    redVec l (a + b) = redVec l a + redVec l b := by
  funext i
  show PadicInt.toZMod (a i + b i) = PadicInt.toZMod (a i) + PadicInt.toZMod (b i)
  exact map_add _ _ _

/-- ★★**還元は差と可換**（第 1205）。 -/
theorem redVec_sub (a b : Fin 2 → ℤ_[l]) :
    redVec l (a - b) = redVec l a - redVec l b := by
  funext i
  show PadicInt.toZMod (a i - b i) = PadicInt.toZMod (a i) - PadicInt.toZMod (b i)
  exact map_sub _ _ _

/-- ★★**還元は自然数倍と可換**（第 1205）。 -/
theorem redVec_nsmul (k : ℕ) (w : Fin 2 → ℤ_[l]) :
    redVec l (k • w) = k • redVec l w := by
  funext i
  show PadicInt.toZMod (k • w i) = k • PadicInt.toZMod (w i)
  exact map_nsmul _ _ _

/-! ## ★★★★★★★★★★核と全射性 -/

/-- ★★★★★★★★★★**還元の核はちょうど `l · ℤ_l²`**——★**無条件**（第 1204）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`PadicInt.ker_toZMod` と `maximalIdeal_eq_span_p` で成分ごとに落とす。

★★★これが第 1203（`tateProj` の核は `l · T_l E`）と**対**になって、
「`mod l` で `≠ 0`」と「第 1 層の点が `≠ 0`」を結ぶ。 -/
theorem redVec_eq_zero_iff (w : Fin 2 → ℤ_[l]) :
    redVec l w = 0 ↔ ∃ u : Fin 2 → ℤ_[l], w = l • u := by
  constructor
  · intro h
    have hdvd : ∀ i, (l : ℤ_[l]) ∣ w i := by
      intro i
      have hi : w i ∈ RingHom.ker (PadicInt.toZMod (p := l)) := congrFun h i
      rwa [PadicInt.ker_toZMod, PadicInt.maximalIdeal_eq_span_p,
        Ideal.mem_span_singleton] at hi
    choose u hu using hdvd
    refine ⟨u, funext fun i => ?_⟩
    show w i = (l : ℕ) • u i
    rw [nsmul_eq_mul]
    exact hu i
  · rintro ⟨u, rfl⟩
    funext i
    show PadicInt.toZMod ((l : ℕ) • u i) = (0 : Fin 2 → ZMod l) i
    rw [nsmul_eq_mul, map_mul, map_natCast]
    simp

/-- ★★★★**還元は全射**——★**無条件**（第 1204）。

☆`v i` の代表 `(v i).val` を `ℤ_l` の自然数として取ればよい。 -/
theorem exists_redVec_eq (v : Fin 2 → ZMod l) : ∃ w : Fin 2 → ℤ_[l], redVec l w = v := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  refine ⟨fun i => (((v i).val : ℕ) : ℤ_[l]), funext fun i => ?_⟩
  show PadicInt.toZMod (((v i).val : ℕ) : ℤ_[l]) = v i
  rw [map_natCast]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def redVec.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ベクトルの mod l 還元の定義のみ)",
    sectionId := "genell-thm-3-8" }

def redVec_mulVec.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(還元は行列の作用と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_add.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 還元は和と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_sub.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 還元は差と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_nsmul.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 還元は自然数倍と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_mulVec_glRed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(glRedPadic の像への作用は還元と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_eq_zero_iff.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 還元の核はちょうど l · Z_l²。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_eq_zero_iff.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1204）**——第 1203（`tateProj` の核は `l · T_l E`）と" ++
       "**対**になって、「`mod l` で `≠ 0`」と「第 1 層の点が `≠ 0`」を結ぶ。" ++
       "☆これが `HasLCyclicJ` の直線（`mod l`）を" ++
       "**点の直線**へ運ぶための最後の配管である。") 2 ]

def exists_redVec_eq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(mod l 還元は全射。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
