import Mathlib.Data.Set.Basic
import Mathlib.Order.Monotone.Basic
import ABC3.Meta.Claim

/-!
# [GenEll] §2 —— `Theorem 2.1` の **(ii) ⟹ (i)** が要る幾何を受ける `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.11–p.13。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★これは何か——**段 B の幾何だけ**を欄にする

`Theorem 2.1` は同値 `(i) ⟺ (ii)` である。

* `(i) ⟹ (ii)` は原文が「immediate from the definitions」と言うとおり、
  **BD-類の不等式は部分集合へ制限できる**という 1 行である。
* `(ii) ⟹ (i)` は原文 p.11–p.13 の 3 ページで、背理法である:
  `(i)` が偽なら「どんどん悪くなる点の列」が取れ、
  コンパクト性でその**部分列**が収束し、
  極限のまわりに noncritical Belyi 写像で compactly bounded subset `K_V` を作ると、
  部分列は（有限個を除いて）`K_V` に入ってしまう。★そこで `(ii)` と矛盾する。

☆★★**その最後の 1 文だけ**が幾何である。本 `Interface` はそれを `cover` 欄で受ける。

## ★★受けているものと、受けていないもの

| 段 | 場所 | 状態 |
|---|---|---|
| 数値の鎖(段 A・段 B) | `Found/GenEll/Thm21Chain.lean` | ✅ `sorry` 0 |
| コンパクト性で部分列を取る | `Found/GenEll/Thm21Extract.lean` | ✅ `sorry` 0 |
| 次数の段(Riemann–Hurwitz の帳簿) | `Found/GenEll/Thm21DegRatio.lean` | ✅ `sorry` 0 |
| noncritical Belyi(一般曲線版) | `Skeleton/GenEll/Section2Converse.lean` の `NoncriticalBelyiData` | ❌ 欄 |
| **列を呑み込む `K_V` が取れる** | ★**本ファイルの `cover` 欄** | ❌ 欄 |
| 組み立て | `Found/GenEll/Thm21Equiv.lean` | ✅ `sorry` 0 |

★★**`cover` は「両方の欄を合わせた出力」である**——
`Found/GenEll/Thm21Extract.lean` の `exists_bad_seq_tendsto` が部分列を出し、
`NoncriticalBelyiData` がその極限のまわりに `K_V` を作る。

## ★★★★posit の向きに注意

`cover` は**仮定**であるから、強く取れば定理は弱くなる。
☆そこで `cover` は「どんな列にも」ではなく
**「`degLe`（＝ `X(ℚ̄)^{≤d}`）の中の列に対して」**だけ要求する
——原文のコンパクト性は「次数が `d` 以下」から来るので、これが原文の強さである。
★また結論は「列そのもの」ではなく**部分列**であり、これも原文どおりである
（`Found/GenEll/Thm21Extract.lean` が出すのは部分列である）。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **`Theorem 2.1` の (ii) ⟹ (i) が要る幾何**を受ける `Interface`（第 1446）。

☆`P` は「点」の型（原文の `X(ℚ̄)`）である。 -/
structure Thm21Data (P : Type) where
  /-- 原文「`X(ℚ̄)^{≤d}`」——次数が `d` 以下の点の集合。 -/
  degLe : Set P
  /-- 原文「a **compactly bounded subset** … whose **support contains** `Σ`」。
  ☆`Interface/GenEll/AbcSetup.lean` の `SupportContains` を合わせたものを 1 つの述語で持つ。 -/
  CB : Set P → Prop
  /-- ★★★★**段 B の幾何**——`degLe` の中の列は、部分列を取れば
  ある compactly bounded subset `K` に**丸ごと入る**。

  ☆原文 p.12 の「コンパクト性で `Ξ` と `Ξ_v` を取り、そのまわりに `K_V` を作る」段である。
  ★材料は `Found/GenEll/Thm21Extract.lean` の `exists_bad_seq_tendsto`（`sorry` 0）と
  `Skeleton/GenEll/Section2Converse.lean` の `NoncriticalBelyiData` である。 -/
  cover : ∀ u : ℕ → P, (∀ n, u n ∈ degLe) →
    ∃ (K : Set P) (φ : ℕ → ℕ), StrictMono φ ∧ CB K ∧ ∀ n, u (φ n) ∈ K

/-- ★Track B は何を作らねばならないか。 -/
def Thm21Data.waiting : WaitingFor :=
  { what := "★段 B の幾何 1 本だけ——「X(ℚ̄)^{≤d} の中の列は、部分列を取ればある compactly bounded subset に丸ごと入る」。☆コンパクト性の側(Found/GenEll/Thm21Extract.lean の exists_bad_seq_tendsto)は既に sorry 0 で手元にあり、残るのは noncritical Belyi 写像で極限のまわりに K_V を作る段である"
    trackB := "Found/GenEll — ★mathlib に `Belyi` は 0 件、`LineBundle` は 0 件(2026-08-16 実測)。☆一般曲線への帰着に Riemann–Roch が要る。★本定理は [IUTchIV] を一切使わないので IUT 側の進捗と独立に進められる" }

/-! ## ★出典の紐付け(`.src`) -/

def Thm21Data.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(段 B の幾何——列を呑み込む compactly bounded subset)",
    sectionId := "genell-thm-2-1" }

end ABC3.Interface.GenEll
