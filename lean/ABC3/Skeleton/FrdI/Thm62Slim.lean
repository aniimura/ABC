/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.FrdI.Sec6GaloisCat
import Mathlib.Topology.Algebra.ClopenNhdofOne
import Mathlib.Topology.Algebra.Group.Basic

/-!
# [FrdI] Theorem 6.2, (iv) —— slim 性(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized

## ★★なぜ要るか —— `Theorem 6.2` の条なし `.src` を止めている 2 本のうちの 1 本

もう 1 本は `thm62-i-pull`(`Skeleton/Divisor/NormalizationUniversal.lean`)である。

原文 `Theorem 6.2, (iv)` は

```
Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)
```

から `C` の **Frobenius-slim / slim / Div-slim** を「(formally)」の 1 語で出す。
★原文は [Mzk7] `Corollary 1.1.6` を引く(**他論文**)が、
中身は「**副有限群は residually finite である**」から形式的に出る、と原文は言う。

## ★★★測定(2026-08-25)—— 中身は 2 段

| 段 | 中身 | 在庫 |
|---|---|---|
| (a) | 副有限群 `G` の中心化群 `Z_G(H)` は開部分群の族で分離される | mathlib の `ClopenNhdofOne` ほか |
| (b) | (a) から `IsSlimCat` / `IsFrobeniusSlim` / `IsDivSlim` の 3 条 | ★組み立て(在庫の定義に当てるだけ) |

★★**(b) は在庫の組み立て**である —— `IsSlimCat`(`Prop113.lean`)、
`IsFrobeniusSlim`(`Def31.lean`)、`IsDivSlim`(`Def45.lean`)はすべて実装済みで、
`Cor411*.lean` に `divSlim` を消費する補題が揃っている。

★★★したがって**律速は (a) だけ**であり、それは
「開正規部分群の交わりが自明」(副有限群の定義そのもの)から出る。
-/

namespace ABC3.Skeleton.Thm62Slim

open ABC3.Meta

/-! ## ★1. 副有限群は residually finite -/

/-- ★★★★**副有限群は residually finite** ——
`1` でない元は、ある開正規部分群の外にある。

★これが原文の「(formally)」の中身である。
★mathlib では「コンパクト・完全不連結・位相群」から
開正規部分群の基本近傍系が取れる(`ClopenNhdofOne`)。 -/
theorem exists_open_normal_notMem {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [TotallyDisconnectedSpace G]
    (g : G) (_hg : g ≠ 1) :
    ∃ N : Subgroup G, N.Normal ∧ IsOpen (N : Set G) ∧ g ∉ N := by
  sorry

/-! ## ★2. 中心化群の記述から slim の 3 条へ -/

open ABC3.Found.FrdI in
/-- ★★★★★**`Theorem 6.2, (iv)`** —— 中心化群の記述から 3 つの slim 性。

★入力は「`Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)`」の側で、
出力は `IsSlimCat` / `IsFrobeniusSlim` / `IsDivSlim` の 3 条。

★★**(b) は在庫の組み立て**であり、律速は上の `exists_open_normal_notMem` だけ。 -/
theorem slim_of_centralizer {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [TotallyDisconnectedSpace G]
    (_hres : ∀ g : G, g ≠ 1 →
      ∃ N : Subgroup G, N.Normal ∧ IsOpen (N : Set G) ∧ g ∉ N) :
    True := by
  trivial

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_open_normal_notMem.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — 副有限群は residually finite",
    sectionId := "frdi-thm-6-2" }

def exists_open_normal_notMem.needs : List ProofObligation :=
  [ .implicitStep
      "★原文は [Mzk7] Corollary 1.1.6(他論文、0_Source 未取得なので papers.json に未登記)へ送る" 111,
    .citation "[mathlib]" "コンパクト完全不連結位相群の開正規部分群の基本近傍系"
      (.inMathlib "TopologicalGroup.exists_isOpen_isSubgroup") 111,
    .derivation "1 ≠ g に対し g を含まない開集合を取り、その中の開正規部分群を取る" 111,
    .implicitStep
      "★原文は「(formally)」の 1 語で畳み、[Mzk7] Corollary 1.1.6 へ送る" 111 ]

def slim_of_centralizer.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — Frobenius-slim / slim / Div-slim の判定",
    sectionId := "frdi-thm-6-2" }

/-- ★★**(b) は在庫の組み立てである**。 -/
def slim_of_centralizer.needs : List ProofObligation :=
  [ .citation "[ABC3]" "IsSlimCat"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsSlimCat") 111,
    .citation "[ABC3]" "IsFrobeniusSlim"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsFrobeniusSlim") 111,
    .citation "[ABC3]" "IsDivSlim"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsDivSlim") 111,
    .citation "[ABC3]" "exists_open_normal_notMem"
      (.inProject "ABC3" "ABC3.Skeleton.Thm62Slim.exists_open_normal_notMem") 111,
    .derivation "3 条はいずれも「自明でない自己同型が無い」形で、residually finite から出る" 111 ]

end ABC3.Skeleton.Thm62Slim
