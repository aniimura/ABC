/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.Section2Converse
import ABC3.Meta.Claim

/-!
# ★★★節点 1（`exists_belyi_noncritical_general`）の**現在の主張は偽**である

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★何を見つけたか（2026-09-02、第 1389）

`Skeleton/GenEll/Section2Converse.lean` の**節点 1**

    exists_belyi_noncritical_general
      (X : Curve) (V) (Xi : V → Set (Pt X)) (hfin : ∀ v, (Xi v).Finite) (UP : Set P1) :
      ∃ φ : Pt X → P1, ∀ v, ∀ x ∈ Xi v, φ x ∈ UP

は **`UP` を任意の集合にしている**ので、`UP = ∅` かつ `Ξ_v ≠ ∅` で**偽**である。
☆さらに `φ` は**単なる函数**であって射ではないので、
`[NCBelyi] Theorem 2.5` の内容（`ℙ¹ ∖ {0,1,∞}` の上で不分岐な射）は
まったく写っていない。

★★★**この節点は `sorry` を埋めようとしても永久に埋まらない**——
主張が偽だからである。★先に**忠実な主張へ書き直す**必要がある。

## ★★★★何に書き直すべきか

`[NCBelyi] Theorem 2.5` の内容は
「滑らかな固有連結曲線 `X` と有限個の点 `S` に対し、
`ℙ¹ ∖ {0,1,∞}` の上で不分岐な射 `φ : X → ℙ¹` で `φ(S) ⊆ {0,1,∞}` なるものが取れる」
である。☆少なくとも

* `φ` を**射**（あるいは `Interface` の射の語彙）にすること
* `UP` を `ℙ¹ ∖ {0,1,∞}` の形に固定すること
* 不分岐性を主張に含めること

が要る。★`Found/NCBelyi/` の `ℙ¹` 版（21 ファイル、`sorry` 0）が
そのための語彙を既に持っている。
-/

namespace ABC3.Check.GenEll

open ABC3.Meta

/-- ★★★★★★★★**節点 1 の現在の主張は偽である**——機械検証（第 1389）。

☆`Curve = Pt X = P1 = Unit`、`Ξ_v = univ`、`UP = ∅` で反例になる。 -/
theorem exists_belyi_noncritical_general_false :
    ¬ (∀ (Curve : Type) (Pt : Curve → Type) (P1 : Type) (X : Curve)
        (V : Type) (Xi : V → Set (Pt X)), (∀ v, (Xi v).Finite) →
        ∀ UP : Set P1, ∃ φ : Pt X → P1, ∀ v, ∀ x ∈ Xi v, φ x ∈ UP) := by
  intro h
  obtain ⟨φ, hφ⟩ :=
    h Unit (fun _ => Unit) Unit () Unit (fun _ => (Set.univ : Set Unit))
      (fun _ => Set.finite_univ) (∅ : Set Unit)
  exact (hφ () () (Set.mem_univ ()))

/-! ## ★出典の紐付け(`.src`) -/

def exists_belyi_noncritical_general_false.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 1 の現在の主張は偽である——機械検証)",
    sectionId := "genell-thm-2-1" }

def exists_belyi_noncritical_general_false.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1389）**——`UP` が任意の集合なので " ++
       "`UP = ∅` かつ `Ξ_v ≠ ∅` で偽である。" ++
       "☆`φ` も単なる函数で射ではないので `[NCBelyi] Thm 2.5` の内容が写っていない。" ++
       "★先に忠実な主張へ書き直すこと——`Found/NCBelyi/` の `ℙ¹` 版が語彙を持っている。") 11 ]

end ABC3.Check.GenEll
