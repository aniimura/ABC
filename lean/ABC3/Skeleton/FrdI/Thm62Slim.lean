/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.FrdI.Sec6GaloisCat
import ABC3.Found.FrdI.Profinite
import ABC3.Found.FrdI.Thm62Slim
import Mathlib.Topology.Algebra.ClopenNhdofOne
import Mathlib.Topology.Algebra.Group.Basic

/-!
# [FrdI] Theorem 6.2, (iv) —— slim 性(★★★★★★**閉じた**、2026-08-25)

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
★原文は [Mzk7] `Corollary 1.1.6`(**他論文**、`papers.json` に登記済み)を引くが、
中身は「**副有限群は residually finite である**」から形式的に出る、と原文は言う。

## ★★★測定(2026-08-25)—— 中身は 2 段、どちらも閉じた

| 段 | 中身 | 状態 |
|---|---|---|
| (a) | 副有限群は residually finite | ★★**済**(`Found/FrdI/Profinite.lean`、`sorry` 無し) |
| (b) | (a) から `IsFrobeniusSlim` / `IsSlimCat` / `IsDivSlim` の 3 条 | ★★**済**(`Found/FrdI/Thm62Slim.lean`、`sorry` 無し) |

★★律速だった (a) の中身は「**開部分群の正規核は開**」で、これは
**mathlib に無かった**(2026-08-25 実測)。
`H.normalCore = ⋂_{q : G ⧸ H} (a ↦ (out q)⁻¹ a (out q))⁻¹ (H)` と書き直し、
`G ⧸ H` が有限であることから有限交叉に落として閉じた。

★★★(b) は在庫の組み立てだった —— `IsSlimCat`(`Prop113.lean`)、
`IsFrobeniusSlim`(`Def31.lean`)、`IsDivSlim`(`Def45.lean`)、
`Remark 3.1.2`(`Remark312.lean`)がすべて実装済みだった。

★★★★**残る他論文の入力は `Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)` だけ**であり、
これは [Mzk7] `Corollary 1.1.6`(物理 p.14)なので**仮定として受ける**。
-/

namespace ABC3.Skeleton.Thm62Slim

open ABC3.Meta CategoryTheory ABC3.Found.FrdI

universe v u w

/-! ## ★1. 副有限群は residually finite -/

/-- ★★★★**副有限群は residually finite** ——
`1` でない元は、ある開正規部分群の外にある。

★これが原文の「(formally)」の中身である。
★★★★★**閉じた**(2026-08-25)—— 中身は `Found/FrdI/Profinite.lean` にある。 -/
theorem exists_open_normal_notMem {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [T2Space G] [TotallyDisconnectedSpace G]
    (g : G) (hg : g ≠ 1) :
    ∃ N : Subgroup G, N.Normal ∧ IsOpen (N : Set G) ∧ g ∉ N :=
  ABC3.Found.FrdI.Profinite.exists_open_normal_notMem g hg

/-! ## ★2. 中心化群の記述から slim の 3 条へ

原文 (FrdI p.112):
> Finally, we consider assertion (iv). First, we observe that if L ⊆ K is a finite
-/

/-- ★★★★★**`Theorem 6.2, (iv)`** —— 中心化群の記述から slim 性。

★入力は [Mzk7] `Corollary 1.1.6` の「`Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)`」の側
(仮定 `Hsub`、`e`)で、出力は

* `𝒟` は **Frobenius-slim**(無条件)
* `𝒟` が **slim** ⟺ すべての `A` で `Z_G(H_A) = {1}`(原文の「`Z = {1}`」)

★★★★★**閉じた**(2026-08-25)—— 中身は `Found/FrdI/Thm62Slim.lean` にある。
★`Div-slim` の同値(`isDivSlim_iff_of_mulEquiv`)も同ファイルにある。 -/
theorem slim_of_centralizer {E : Type u} [Category.{v} E]
    {G : Type w} [Group G] [TopologicalSpace G] [IsTopologicalGroup G] [CompactSpace G]
    [T2Space G] [TotallyDisconnectedSpace G]
    (Hsub : E → Subgroup G) (e : ∀ A : E, Aut (Over.forget A) ≃* (Hsub A)) :
    IsFrobeniusSlim E ∧ (IsSlimCat E ↔ ∀ A : E, Hsub A = ⊥) :=
  ⟨isFrobeniusSlim_of_mulEquiv_subgroup Hsub e, isSlimCat_iff_of_mulEquiv Hsub e⟩

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_open_normal_notMem.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — 副有限群は residually finite",
    sectionId := "frdi-thm-6-2" }

def exists_open_normal_notMem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Found 側の本体(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.FrdI.Profinite.exists_open_normal_notMem") 111,
    .citation "[mathlib]" "コンパクト完全不連結位相群の clopen 近傍の中の開部分群"
      (.inMathlib "IsTopologicalGroup.exist_openSubgroup_sub_clopen_nhds_of_one") 111,
    .derivation "1 ≠ g に対し g を含まない clopen を取り、その中の開部分群の正規核を取る" 111 ]

def slim_of_centralizer.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — Frobenius-slim / slim / Div-slim の判定",
    sectionId := "frdi-thm-6-2" }

/-- ★★**(b) は在庫の組み立てだった**。 -/
def slim_of_centralizer.needs : List ProofObligation :=
  [ .otherPaper "[Mzk7]" "Corollary 1.1.6 — Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)。★仮定 e として受ける" 14,
    .citation "[ABC3]" "Found 側の本体(Frobenius-slim、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.FrdI.isFrobeniusSlim_of_mulEquiv_subgroup") 111,
    .citation "[ABC3]" "Found 側の本体(slim ⟺ Z = {1}、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.FrdI.isSlimCat_iff_of_mulEquiv") 111,
    .citation "[ABC3]" "Found 側の本体(Div-slim の同値、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.FrdI.isDivSlim_iff_of_mulEquiv") 111,
    .citation "[ABC3]" "exists_open_normal_notMem"
      (.inProject "ABC3" "ABC3.Skeleton.Thm62Slim.exists_open_normal_notMem") 111 ]

end ABC3.Skeleton.Thm62Slim
