/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm21Chain
import Mathlib.Topology.Sequences
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`Theorem 2.1` のコンパクト性の段（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.12。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

原文 p.12 は、`ht_{ω_X} ≲ (1+ϵ)·log-diff_X` が `X(ℚ̄)^{=d}` の上で**偽**だとしたとき、

> it follows immediately from the compactness of the set of rational points of X over
> any finite extension of Qv for v ∈V that there exists a subset Ξ ⊆X(Q)=d, together
> with a(n) [unordered] d-tuple of points Ξv of X(Qv) for each v ∈V , such that ...

と言う。★★★**本ファイルはこの段を取る**——`§9-976` が残した幾何の 4 点のうちの
**第 3 点（(c) 局所コンパクト性から `Ξ`・`Ξ_v` を取る段）**である。

## ★機構 —— 2 つに分かれる

1. **`¬ BDge α β` から「悪い点の列」を作る**（`exists_seq_of_not_bdge`）。
   `BDge α β ≔ ∃ C, ∀ x, α x − β x ≤ C` だから、その否定は
   「任意の `C` に**超える点がある**」——選択公理で `α(x_n) − β(x_n) > n` なる列が取れる。
2. **有限個のコンパクト空間の中で同時に収束する部分列を取る**
   （`exists_subseq_tendsto_finite_family`）。
   ★有限積はコンパクト（Tychonoff）かつ第一可算なので、`IsCompact.tendsto_subseq` が使える。
   ★★`tendsto_pi_nhds` で成分ごとの収束に戻す。

★★★この 2 つを合わせたのが `exists_bad_seq_tendsto` である。

## ★原文の `Ξ_v` は「順序なし `d`-組」である

原文の `Ξ_v` は `X(ℚ_v)` の**順序なし `d`-組**であり、`Ξ` の点の `ℚ̄`-共役の組が
そこへ収束する。★本ファイルの `M i` にその「順序なし `d`-組の空間」を、
`emb i` に「点をその共役の組へ送る写像」を入れれば、原文そのものになる。
★★**コンパクト性（`CompactSpace (M i)`）と第一可算性が仮定である**
——`X(ℚ_v)` が固有なら成り立つ（`v` がアルキメデスでも非アルキメデスでも）。

## ★★これで `Theorem 2.1` に残るのは 3 点

| 入力 | 状態 |
|---|---|
| (a) 分岐指数がちょうど `e` の連結有限エタール Galois 被覆 | ☆残る |
| (b) noncritical Belyi 写像（[NCBelyi] `Theorem 2.5`） | ☆残る |
| ★(c) 局所コンパクト性から `Ξ`・`Ξ_v` を取る段 | ★★**本ファイル** |
| (d) `deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)` | ☆残る |
-/

namespace ABC3.Found.GenEll

open Filter Topology

/-! ## ★★★★★★★★★★★悪い点の列 -/

/-- ★★★★★★★★★★★**`≲` が成り立たないなら「どんどん悪くなる点の列」がある**。

`BDge α β ≔ ∃ C, ∀ x, α x − β x ≤ C` の否定は
「任意の `C` に対して `α x − β x > C` なる `x` がある」である。
★選択公理で `C = n` と取れば列になる。 -/
theorem exists_seq_of_not_bdge {Pt : Type} (α β : Pt → ℝ) (h : ¬ BDge α β) :
    ∃ x : ℕ → Pt, ∀ n : ℕ, (n : ℝ) < α (x n) - β (x n) := by
  have hall : ∀ C : ℝ, ∃ p : Pt, C < α p - β p := by
    intro C
    by_contra hc
    exact h ⟨C, by intro p; by_contra hp; exact hc ⟨p, lt_of_not_ge hp⟩⟩
  choose g hg using hall
  exact ⟨fun n => g (n : ℝ), fun n => hg _⟩

/-! ## ★★★★★★★★★★★★★★有限個のコンパクト空間で同時に収束させる -/

/-- ★★★★★★★★★★★★★★**有限個のコンパクト第一可算空間の中で同時に収束する
部分列が取れる**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★原文 p.12 の「`the compactness of the set of rational points of X over any finite
extension of Q_v for v ∈ V`」がこれである。
★★有限積はコンパクト（Tychonoff）かつ第一可算だから、
`IsCompact.tendsto_subseq` を積の上で 1 回使えばよい。 -/
theorem exists_subseq_tendsto_finite_family
    {ι : Type} [Fintype ι] {M : ι → Type}
    [∀ i, TopologicalSpace (M i)] [∀ i, CompactSpace (M i)]
    [∀ i, FirstCountableTopology (M i)]
    (f : ∀ i, ℕ → M i) :
    ∃ φ : ℕ → ℕ, StrictMono φ ∧ ∃ a : ∀ i, M i,
      ∀ i, Tendsto (fun n => f i (φ n)) atTop (𝓝 (a i)) := by
  have hcpt : IsCompact (Set.univ : Set (∀ i, M i)) := isCompact_univ
  obtain ⟨a, -, φ, hφ, htend⟩ := hcpt.tendsto_subseq (x := fun n i => f i n)
    (fun n => Set.mem_univ _)
  exact ⟨φ, hφ, a, fun i => (tendsto_pi_nhds.mp htend) i⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★原文 p.12 の段そのもの -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**原文 p.12 のコンパクト性の段**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

`ht_{ω_X} ≲ (1+ϵ)·log-diff_X` が偽（`¬ BDge α β`）なら、

* 不等式が**破れ続ける**点の列 `x : ℕ → Pt`（`α(x_n) − β(x_n) > n`）と
* 各 `v ∈ V`（有限）での**極限** `a_v`

が取れて、`x_n` の `v` での像が `a_v` に収束する。

★これが原文の `Ξ`（点の集合）と `Ξ_v`（順序なし `d`-組の極限）である。
★★`M i` に「`X(ℚ_v)` の順序なし `d`-組の空間」を、`emb i` に
「点をその `ℚ̄`-共役の組へ送る写像」を入れれば原文そのものになる。 -/
theorem exists_bad_seq_tendsto {Pt : Type} (α β : Pt → ℝ) (h : ¬ BDge α β)
    {ι : Type} [Fintype ι] {M : ι → Type}
    [∀ i, TopologicalSpace (M i)] [∀ i, CompactSpace (M i)]
    [∀ i, FirstCountableTopology (M i)]
    (emb : ∀ i, Pt → M i) :
    ∃ (x : ℕ → Pt) (a : ∀ i, M i),
      (∀ n : ℕ, (n : ℝ) < α (x n) - β (x n))
      ∧ ∀ i, Tendsto (fun n => emb i (x n)) atTop (𝓝 (a i)) := by
  obtain ⟨y, hy⟩ := exists_seq_of_not_bdge α β h
  obtain ⟨φ, hφ, a, ha⟩ :=
    exists_subseq_tendsto_finite_family (M := M) (fun i n => emb i (y n))
  refine ⟨fun n => y (φ n), a, fun n => ?_, ha⟩
  have h1 : (n : ℝ) ≤ (φ n : ℝ) := by exact_mod_cast hφ.le_apply
  exact lt_of_le_of_lt h1 (hy (φ n))

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def exists_seq_of_not_bdge.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(≲ が偽なら、どんどん悪くなる点の列がある)",
    sectionId := "genell-thm-2-1" }

def exists_subseq_tendsto_finite_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(有限個のコンパクト空間で同時に収束する部分列)",
    sectionId := "genell-thm-2-1" }

def exists_bad_seq_tendsto.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(段 B のコンパクト性の段——Ξ と Ξ_v を取る)",
    sectionId := "genell-thm-2-1" }

def exists_bad_seq_tendsto.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsCompact.tendsto_subseq(コンパクト第一可算なら収束部分列)"
      (.inMathlib "IsCompact.tendsto_subseq") 2,
    .citation "[mathlib]" "tendsto_pi_nhds(積の収束は成分ごと)"
      (.inMathlib "tendsto_pi_nhds") 1,
    .implicitStep
      ("★★★★★測定(2026-08-29): 原文 p.12 の『the compactness of the set of rational " ++
       "points of X over any finite extension of Q_v』の段は、" ++
       "**有限個のコンパクト第一可算空間の積で収束部分列を取ること**に尽きる。" ++
       "★Tychonoff と IsCompact.tendsto_subseq で 1 行である") 5,
    .implicitStep
      ("★★原文の Ξ_v は X(ℚ_v) の**順序なし d-組**であり、Ξ の点の ℚ̄-共役の組が" ++
       "そこへ収束する。★本ファイルの M i にその空間を、emb i に共役の組へ送る写像を" ++
       "入れれば原文そのものになる——**その 2 つは受けている**" ++
       "(コンパクト性と第一可算性が仮定)") 4,
    .implicitStep
      ("★★★これで Theorem 2.1 に残るのは 3 点である: " ++
       "(a) 分岐指数がちょうど e の連結有限エタール Galois 被覆、" ++
       "(b) noncritical Belyi 写像([NCBelyi] Theorem 2.5)、" ++
       "(d) deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)(deg が語彙の外)") 6 ]

end ABC3.Found.GenEll
